#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 512

/* Initialization variables */
const int    mc_steps      = 10000;
const int    output_steps  = 100;
const double density       = 0.8;
const double delta         = 0.1;
const double r_cut         = 2.5;
const double beta          = 0.5;
const char*  init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double e_cut;
double r[N][NDIM];
double box[NDIM];

double energy = 0.0;
double virial = 0.0;

typedef struct{
    double energy;
    double virial;
}particle_info_t;

typedef struct{
    double average_pressure;
    double mu_excess;
}measurement_t;

particle_info_t particle_energy_and_virial(int);

/* Functions */
measurement_t measure(void){
    measurement_t result;
    /*--------- Your code goes here -----------*/

    return result;
}

particle_info_t particle_energy_and_virial(int pid){
    particle_info_t info;
    info.energy = 0.0;
    info.virial = 0.0;
    int n, d;
    for(n = 0; n < n_particles; ++n){
        if(n == pid) continue;
        double dist2 = 0.0;
        for(d = 0; d < NDIM; ++d){
            double min_d = r[pid][d] - r[n][d];
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            dist2 += min_d * min_d;
        }

        if(dist2 <= r_cut * r_cut){
            double temp = 1.0 / (dist2 * dist2 * dist2);
            info.energy += 4.0 * temp * (temp - 1.0) - e_cut;
            info.virial += 24.0 * temp * (2.0 * temp - 1.0);
        }
    }

    return info;
}

void read_data(void){
    FILE* fp = fopen(init_filename, "r");
    int n, d;
    double dmin,dmax;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
        double diameter;
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
}

int move_particle(void){
    int rpid = n_particles * dsfmt_genrand();


    particle_info_t info = particle_energy_and_virial(rpid);

    double old_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    particle_info_t new_info = particle_energy_and_virial(rpid);

    double dE = new_info.energy - info.energy;
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        virial += new_info.virial - info.virial;
        return 1;
    }

    for(d = 0; d < NDIM; ++d) r[rpid][d] = old_pos[d];

    return 0;
}

void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", 1.0);
    }
    fclose(fp);
}

void set_density(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

int main(int argc, char* argv[]){

    assert(delta > 0.0);

    e_cut = 4.0 * (pow(1.0 / r_cut, 12.0) - pow(1.0 / r_cut, 6.0));

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_density();

    int d;
    for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box[d]);

    int step, n;
    for(n = 0; n < n_particles; ++n){
        particle_info_t info = particle_energy_and_virial(n);
        energy += info.energy;
        virial += info.virial;
    }
    energy *= 0.5;
    virial *= 0.5;

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    double volume = 1.0;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    printf("Starting volume: %f\n", volume);
    printf("Starting energy: %f\n", energy);
    printf("Starting virial: %f\n", virial);
    printf("Starting seed: %lu\n", seed);

    FILE* fp = fopen("measurements.dat", "w");

    int accepted = 0;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }

        measurement_t ms = measure();

        fprintf(fp, "%d\t%f\t%f\n", step, ms.average_pressure, ms.mu_excess);

        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %f.\n",
                step, (double)accepted / (n_particles * output_steps)
            );
            accepted = 0;
            write_data(step);
        }
    }

    fclose(fp);

    return 0;
}
