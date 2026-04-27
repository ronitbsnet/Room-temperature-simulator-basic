

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>


typedef struct {
    double room_volume_m3;
    double air_density;
    double air_cp;
    double thermal_R;
    double heater_power;
    double setpoint;
    double amb_temp;
    double dt;
    double sim_time_s;
    double hysteresis;
    double Kp;
    double Ki;
    double Kd;
    double pid_window_s;
    double sensor_noise_std;
    double sensor_sample_interval;
} SimParams;


void SimParams_init(SimParams *p) {
    p->room_volume_m3 = 50.0;
    p->air_density = 1.2;
    p->air_cp = 1005.0;
    p->thermal_R = 1.5;
    p->heater_power = 2000.0;
    p->setpoint = 22.0;
    p->amb_temp = 18.0;
    p->dt = 1.0;
    p->sim_time_s = 3 * 3600.0;
    p->hysteresis = 1.0;
    p->Kp = 60.0;
    p->Ki = 0.02;
    p->Kd = 0.0;
    p->pid_window_s = 30.0;
    p->sensor_noise_std = 0.05;
    p->sensor_sample_interval = 1.0;
}


typedef struct {
    SimParams *params;
    double C_air;
    double T;
} ThermalModel;

void ThermalModel_init(ThermalModel *m, SimParams *p) {
    m->params = p;
    m->C_air = p->room_volume_m3 * p->air_density * p->air_cp;
    m->T = p->amb_temp; 
}

void ThermalModel_step(ThermalModel *m, double dt, double Q_heater, double Tamb) {
    double Q_loss = (Tamb - m->T) / m->params->thermal_R;
    double dT = (Q_loss + Q_heater) * dt / m->C_air;
    m->T += dT;
}

double ThermalModel_temperature(ThermalModel *m) {
    return m->T;
}


static double randn_boxmuller() {
    static bool have_spare = false;
    static double spare;
    if (have_spare) {
        have_spare = false;
        return spare;
    }
    double u, v, s;
    do {
        u = (rand() + 1.0) / (RAND_MAX + 2.0); 
        v = (rand() + 1.0) / (RAND_MAX + 2.0);
        u = 2.0 * u - 1.0; 
        v = 2.0 * v - 1.0;
        s = u*u + v*v;
    } while (s == 0.0 || s >= 1.0);
    double mul = sqrt(-2.0 * log(s) / s);
    spare = v * mul;
    have_spare = true;
    return u * mul;
}


static double randn(double mean, double std) {
    return mean + randn_boxmuller() * std;
}


typedef struct {
    SimParams *params;
    double last_sample_time;
    double last_meas;
} Sensor;

void Sensor_init(Sensor *s, SimParams *p) {
    s->params = p;
    s->last_sample_time = -p->sensor_sample_interval; 
    s->last_meas = p->amb_temp;
}

double Sensor_read(Sensor *s, double true_temp, double time_s) {
    if (time_s - s->last_sample_time >= s->params->sensor_sample_interval - 1e-9) {
        double noise = randn(0.0, s->params->sensor_noise_std);
        s->last_meas = true_temp + noise;
        s->last_sample_time = time_s;
    }
    return s->last_meas;
}


typedef struct {
    SimParams *params;
    bool state;
} BangBangController;

void BangBang_init(BangBangController *c, SimParams *p) {
    c->params = p;
    c->state = false;
}

bool BangBang_update(BangBangController *c, double measured_temp) {
    double sp = c->params->setpoint;
    double h = c->params->hysteresis;
    if (measured_temp <= sp - h/2.0) c->state = true;
    else if (measured_temp >= sp + h/2.0) c->state = false;
    return c->state;
}


typedef struct {
    SimParams *params;
    double integral;
    double prev_error;
} PIDController;

void PID_init(PIDController *pc, SimParams *p) {
    pc->params = p;
    pc->integral = 0.0;
    pc->prev_error = 0.0;
}

double PID_update(PIDController *pc, double measured_temp, double dt) {
    double error = pc->params->setpoint - measured_temp;
    pc->integral += error * dt;
    double derivative = (error - pc->prev_error) / dt;
    pc->prev_error = error;
    double u = pc->params->Kp * error + pc->params->Ki * pc->integral + pc->params->Kd * derivative;
    double frac = u / pc->params->heater_power;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;
    
    if (frac <= 0.0 && error < 0) pc->integral -= error * dt;
    if (frac >= 1.0 && error > 0) pc->integral -= error * dt;
    return frac;
}


typedef struct {
    SimParams *params;
    int window_steps;
    int step_counter;
    int on_steps;
} TimeProportionalRelay;

void TPR_init(TimeProportionalRelay *r, SimParams *p) {
    r->params = p;
    r->window_steps = (int)fmax(1.0, p->pid_window_s / p->dt);
    r->step_counter = 0;
    r->on_steps = 0;
}

bool TPR_update(TimeProportionalRelay *r, double frac) {
    if (r->step_counter == 0) {
        int os = (int)lround(frac * r->window_steps);
        if (os < 0) os = 0;
        if (os > r->window_steps) os = r->window_steps;
        r->on_steps = os;
    }
    bool out = (r->step_counter < r->on_steps);
    r->step_counter = (r->step_counter + 1) % r->window_steps;
    return out;
}


typedef struct {
    double time_s;
    double T_actual;
    double T_measured;
    int u_bin;
    double u_frac;
} LogRow;


int main(int argc, char **argv) {
    SimParams params;
    SimParams_init(&params);

    printf("Enter setpoint temperature (°C) [default %.1f]: ", params.setpoint);
    if (scanf("%lf", &params.setpoint) != 1) { params.setpoint = 22.0; }
    
    int c;
    while ((c = getchar()) != '\n' && c != EOF) { }

    printf("Enter ambient temperature (°C) [default %.1f]: ", params.amb_temp);
    if (scanf("%lf", &params.amb_temp) != 1) { params.amb_temp = 18.0; }
    while ((c = getchar()) != '\n' && c != EOF) { }

    printf("Enter heater power (W) [default %.0f]: ", params.heater_power);
    if (scanf("%lf", &params.heater_power) != 1) { params.heater_power = 2000.0; }
    while ((c = getchar()) != '\n' && c != EOF) { }

    printf("Enter simulation time (s) [default %.0f]: ", params.sim_time_s);
    if (scanf("%lf", &params.sim_time_s) != 1) { params.sim_time_s = 3*3600.0; }
    while ((c = getchar()) != '\n' && c != EOF) { }

    char ctrl[16] = "pid";
    if (argc >= 2) {
        strncpy(ctrl, argv[1], sizeof(ctrl)-1);
        ctrl[sizeof(ctrl)-1] = '\0';
    }

    
    srand((unsigned int)time(NULL));

    ThermalModel room;
    ThermalModel_init(&room, &params);

    Sensor sensor;
    Sensor_init(&sensor, &params);

    BangBangController bb;
    BangBang_init(&bb, &params);

    PIDController pid;
    PID_init(&pid, &params);

    TimeProportionalRelay tpr;
    TPR_init(&tpr, &params);

    double dt = params.dt;
    int steps = (int)ceil(params.sim_time_s / dt);

    LogRow *log = (LogRow*)malloc((size_t)(steps + 5) * sizeof(LogRow));
    if (!log) {
        fprintf(stderr, "Failed to allocate log buffer\n");
        return 1;
    }
    int log_idx = 0;

    double time_s = 0.0;
    
    double T_true = ThermalModel_temperature(&room);
    double T_meas = Sensor_read(&sensor, T_true, time_s);

    for (int i = 0; i < steps; ++i) {
        T_true = ThermalModel_temperature(&room);
        T_meas = Sensor_read(&sensor, T_true, time_s);

        double frac = 0.0;
        int u_bin = 0;

        if (strcmp(ctrl, "bang") == 0) {
            bool on = BangBang_update(&bb, T_meas);
            u_bin = on ? 1 : 0;
            frac = (double)u_bin;
        } else {
            frac = PID_update(&pid, T_meas, dt);
            bool on = TPR_update(&tpr, frac);
            u_bin = on ? 1 : 0;
        }

        double Q_heater = (u_bin ? params.heater_power : 0.0);
        ThermalModel_step(&room, dt, Q_heater, params.amb_temp);

        double T_actual = ThermalModel_temperature(&room);
        log[log_idx].time_s = time_s;
        log[log_idx].T_actual = T_actual;
        log[log_idx].T_measured = T_meas;
        log[log_idx].u_bin = u_bin;
        log[log_idx].u_frac = frac;
        log_idx++;

        time_s += dt;
    }

    
    const char *fname = "sim_output.csv";
    FILE *ofs = fopen(fname, "w");
    if (!ofs) {
        fprintf(stderr, "Failed to open %s for writing\n", fname);
        free(log);
        return 1;
    }
    fprintf(ofs, "time_s,T_actual,T_measured,u_bin,u_frac\n");
    for (int i = 0; i < log_idx; ++i) {
        fprintf(ofs, "%.6f,%.6f,%.6f,%d,%.6f\n",
                log[i].time_s,
                log[i].T_actual,
                log[i].T_measured,
                log[i].u_bin,
                log[i].u_frac);
    }
    fclose(ofs);

    if (log_idx > 0) {
        printf("Simulation complete. Controller=%s\n", ctrl);
        printf("Final temp: %.3f °C (setpoint %.3f)\n", log[log_idx-1].T_actual, params.setpoint);
        printf("CSV written to: %s\n", fname);
    } else {
        printf("No data logged.\n");
    }

    free(log);
    return 0;
}
