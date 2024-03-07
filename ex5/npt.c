#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;
const double packing_fraction = 0.6;
double diameter = 1.0;
double delta  = 0.1;
/* Volume change -deltaV, delta V */
double deltaV = 2.0;
/* Reduced pressure \beta P */
const double betaP = 3.0;
const char* init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];


/* Functions */
int change_volume(void){
    /*--------- Your code goes here -----------*/
    return 0;
}

void read_data(void)
{
    //open the file
    FILE *file;
    file = fopen(init_filename, "r");

    //some vars for getting the box size
    float len0;
    float len1;

    //scan the first line for the number of particles
    fscanf(file, "%i\n", &n_particles);

    //set n_particles right
    if (n_particles > N) n_particles = N;

    //sanity check
    printf("num of particles: %i", n_particles);

    //get the length of the box's side
    for(int b=0; b<NDIM; b++)
    {
        //where am I
        printf("\nb: %i\t", b);

        //scannin the second line for the box start and stop
        fscanf(file, "%f\t%f\n", &len0, &len1);

        //storin the box len
        box[b] = len1 - len0;

        //check
        printf("\tbox len %lf", box[b]);
    }

    //reading the coords
    for (int a = 0 ; a < n_particles ; a++) //reading coordinates of particles
    {
        //print so we know what's goin on
        //printf("%i\n", a);

        //scan the line for the coords
        fscanf(file ,"%lf\t%lf\t%lf\n", &(r[a][0]), &(r[a][1]), &(r[a][2]));//had another %lf\n in the string, &diameter);
    }

    //close the file, we're done here boys
    fclose (file);
}

int move_particle(void){
    //pick a random particle
    int index = (int) (dsfmt_genrand() * n_particles);

    //arrays for the new position and move
    double new_pos[NDIM];
    double delta_r[NDIM];

    //loop over dims and gen new coords
    for(int i = 0; i < NDIM; i++)
    {
        //generate the new coord
        delta_r[i] = (dsfmt_genrand() - 0.5) * delta;
        new_pos[i] = r[index][i] + delta_r[i];

        //make sure it's a position in the box
        if(new_pos[i] < 0.0) new_pos[i] += box[i];
        if(new_pos[i] >= box[i]) new_pos[i] -= box[i];
    }

    //overlap checker
    for(int q = 0; q < n_particles; q++)
    {
        if(q == index) continue;

        //some dubs for storage
        double dr[NDIM];
        double dr2 = 0.0;

        //dim loop for new positions
        for(int d = 0; d < NDIM; d++)
        {
            dr[d] = new_pos[d] - r[q][d];
            if(dr[d] > 0.5 * box[d]) dr[d] -= box[d];
            else if(dr[d] < -0.5 * box[d]) dr[d] += box[d];
            dr2 += dr[d]*dr[d];
        }

        //if (overlap){don't}
        if (sqrt(dr2) < diameter)
        {
        // Overlap detected, return 0 (not accepted)
        return 0;
        }
    }

    //loop through dimensions to save the new coord
    for(int l = 0; l < NDIM; l++)
    {
        //process the move
        r[index][l] = new_pos[l];
    }
    //apparently returning 1 is good in this context...
    return 1;
}

int adjust_delta(int move_accepted, int vol_accepted){
    //Adjusting delta and deltaV by 5% each cycle of output_steps to keep the respective acceptance ratios between 45% and 55%
    if(((double)move_accepted / (n_particles * output_steps)) < 0.45) delta *= 0.95;
    else if(((double)move_accepted / (n_particles * output_steps)) > 0.55) delta *= 1.05;
    if(((double)vol_accepted / output_steps) < 0.45) deltaV *= 0.95;
    else if(((double)vol_accepted / output_steps) > 0.55) deltaV *= 1.05;
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
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

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

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("\nNumber of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("\nError: Number of particles, n_particles = 0.");
        return 0;
    }

    set_packing_fraction();

    dsfmt_seed(time(NULL));
            
    printf("\n\n#Step \t Volume \t Move-acceptance\t Volume-acceptance");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            move_accepted += move_particle();
        }
        vol_accepted += change_volume();

        if(step % output_steps == 0){
            printf("\n%d \t %lf \t %lf \t %lf", 
                    step, box[0] * box[1] * box[2], 
                    (double)move_accepted / (n_particles * output_steps), 
                    (double)vol_accepted /  output_steps);
            move_accepted = 0;
            vol_accepted = 0;
            write_data(step);
        }
    }

    return 0;
}