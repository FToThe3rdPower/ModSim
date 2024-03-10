#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 1000;
const int output_steps = 100;
const double packing_fraction = 0.65;
const double diameter = 0.1;
const double delta = 0.1;
const char* init_filename = "fcc5.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];


/* Functions */

/* Read the initial configuration provided */
void read_data(void){
    FILE* fp = fopen(init_filename, "r");
    int n, d;
    double dmin, dmax;//, dia;
    fscanf(fp, "%d\n", &n_particles);
    assert(n_particles<10000);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
        //fscanf(fp, "%lf\n", &dia);
        //assert(dia==diameter);
    }
    fclose(fp);
}

/* Move a particle randomly */
int move_particle(void){
    //Choose a random particle
    int rpid = n_particles * dsfmt_genrand();

    double new_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        //Displace by a random value between -delta and delta
        new_pos[d] = r[rpid][d] + delta * (2.0 * dsfmt_genrand() - 1.0);
        //Apply periodic boundaries
        //new_pos[d] -= floor(new_pos[d] / box[d]) * box[d];
        if(new_pos[d]>box[d]) {
            new_pos[d]-=box[d];
        }
        if(new_pos[d]<0.0) {
            new_pos[d]+=box[d];
        }
        assert(new_pos[d] < box[d]);
        assert(new_pos[d] >= 0.0);
    }

    int n;
    //Check for overlaps
    for(n = 0; n < n_particles; ++n){
        if(n == rpid) continue;
        double dist = 0.0;
        for(d = 0; d < NDIM; ++d){
            double min_d = new_pos[d] - r[n][d];
            // Find the distance with the Nearest Image Convension
            //min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            if(min_d>0.5*box[d]) {
                min_d -= box[d];
            }
            if(min_d<-0.5*box[d]) {
                min_d += box[d];
            }
            assert(min_d <= 0.5 * box[d]);
            assert(min_d >= -0.5 * box[d]);
            dist += min_d * min_d;
        }
        if(dist <= diameter * diameter){
            //reject the move
            return 0;
        }
    }

    //Accept the move if reached here
    for(d = 0; d < NDIM; ++d) r[rpid][d] = new_pos[d];

    return 1;
}

/* Write the configuration files */
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
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

/* Scales the box volume and coordinates read from 'init_fileaname'
 * according to the value of 'packing fraction'
 * */
void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

int main(int argc, char* argv[]){

    assert(packing_fraction > 0.0 && packing_fraction < 1.0);
    assert(diameter > 0.0);
    assert(delta > 0.0);

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }
    
    printf("Particle Volume: %lf\n",particle_volume);

    //Read the input configuration
    read_data();
    printf("Box: %lf %lf %lf\n",box[0], box[1], box[2]);

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    //Set the packing fraction according to the set variable
    set_packing_fraction();
    printf("Box: %lf %lf %lf\n",box[0], box[1], box[2]);

    dsfmt_seed(time(NULL));

    int accepted = 0;
    int step, n;
    //Perform MC moves
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }

        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %f.\n", step, (double)accepted / (n_particles * output_steps));
            accepted = 0;
            write_data(step);
        }
    }

    return 0;
}