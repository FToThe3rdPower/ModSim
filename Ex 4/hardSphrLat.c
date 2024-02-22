//
//  hardSphrLat.c
//  
//
//
//Goal: save a file that has the total number of hard spheres and where they are on a lattice.

#include <stdio.h>
#include <math.h>

/* First we're going to define the main function,
 before any constants or vars, so they're all local
 and will be cleared from the memory when the funtion
 is done.*/
int main() /* The "int" is saying this func will return integers*/
{
    /* Let's define our consts */
    int N=20; /*number of spheres per axis*/
    int n=0; //sphere id?
    int numSphr=4*N*N*N; /*the total number of sphere's we'll be worrying about*/
    float a=10.0; /*the distance between points???*/
    float d=5.0; /*the DIAMETER of the sphere, the unit length for this exercise*/
    float b=1.1*d; // ??? the upper bound for "l", the lattice spacing in the innermost loop

    /* Now lets define the axis arrays, they need to be big enough to store a coord for every sphere*/
    double x[numSphr]; //C doesn't have an exponential operation built in unfortunately. The one in math.c requires double floats, and N is an int
    double y[numSphr];
    double z[numSphr];
    
    /* This needs to be output as a file "fcc_lat.dat" we're gonna call coords*/
    FILE *coords;
    coords = fopen("lat.dat", "w+"); //w+ mode lets us read & write to & from the file

    /*put the stuff into the file*/
   fprintf(coords, "%i\n", numSphr); //this makes the first line in the file the total number of spheres
   //fprintf(coords, "Box dimensions x:(%lf, %lf) \t y:(%lf, %lf) \t z:(%lf, %lf)", 0.0, N*a, 0.0, N*a, 0.0, N*a); //and this is line 2, the min, max of the box in each dim
   fprintf(coords, "%lf\t%lf\n", -a, N*a);
   fprintf(coords, "%lf\t%lf\n", -a, N*a);
   fprintf(coords, "%lf\t%lf\n", -a, N*a);

   for(int i=0; i<N; i++) /*this is the x coord loop,
   it keeps track of its count with "i", it will run from 0 to N, so if N=3 it will do i=0, i=1, and finish after i=2*/
   {
        for(int j=0; j<N; j++) //this is the y loop
        {
            for(int k=0; k<N; k++) //this is the z loop
            {
                for(int l=0; l<b; l++)
                {
                    switch(l) //this is like the whole mess of if statements in the example, but cleaner and more efficient
                    {
                        case(0) : printf("l = 0   "); //it will create a ¿sphere? at the point (i,j,k)*a
                        x[n+1] = i*a;
                        y[n+1] = j*a;
                        z[n+1] = k*a;
                        printf("\n");
                        break;

                        case(1) : printf("l = 1   "); //this one scoots over .5a in y & z to make a new
                        x[n+1] = i*a;
                        y[n+1] = j*a + 0.5*a;
                        z[n+1] = k*a + 0.5*a;
                        printf("\n");
                        break;

                        case(2) : printf("l = 2   ");
                        x[n+1] = i*a + 0.5*a;
                        y[n+1] = j*a;
                        z[n+1] = k*a + 0.5*a;
                        printf("\n");
                        break;

                        case(3) : printf("l = 3   ");
                        x[n+1] = i*a + 0.5*a;
                        y[n+1] = j*a + 0.5*a;
                        z[n+1] = k*a;
                        printf("\n");
                        break;
                    }
                    printf("x = %lf, y = %lf, z = %lf\n", x[n], y[n], z[n]); //printing to the console to spot errors
                    fprintf(coords, "%lf, %lf, %lf\n", x[n], y[n], z[n]);    //printing to the file
                }
                n++;
            }
        }
   }

   //printf("X array: %lf, %lf \n", x[30], x[1]);
   fclose(coords);
}
