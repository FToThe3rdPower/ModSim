// Include standard input and output library for IO operations
#include <stdio.h>
// Include the standard library for general purpose functions like memory allocation, process control, conversions, etc.
#include <stdlib.h>
// import standard bool for storing only the unique data before writin it to the file
#include <stdbool.h>

//check for duplicates to avoid russian nesting doll particles that violate the Pauli exc. prin.
//They will appear to multiply with the movement code
/*
bool is_duplicate(Coordinate* arr, size_t size, Coordinate coord)
{
  for (size_t i = 0; i < size; ++i)
  {
    if (arr[i].x == coord.x && arr[i].y == coord.y && arr[i].z == coord.z)
    {
      return true;
    }
  }
  return false;
}
*/


// Entry point of the program with command line arguments
// argc: argument count, argv: argument vector (list of arguments)
int main(int argc, char *argv[]){
    // Convert the first command line argument to an int
    //used determine the number of points along one edge of the cube
    //now it's the number of unit cells
    int N = atoi(argv[1]);

    //determining the total number of particle coords. 1 whole particle per fcc face, 6 faces, N^3 for the dimensions
    int totalParticles = 14*N*N*N;//N*N*N;


    // Set the spacing between points in the lattice
    float spacing = 2.0;
    //float spacing = atoi(argv[2]);
    float halfSpacing = spacing / 2;


    // Declare loop counter vars
    int i, j, k;

    // Declare arrays to store the x, y, and z coordinates of the points
    //in the lattice
    float x[totalParticles];
    float y[totalParticles];
    float z[totalParticles]; // These arrays can store N^3 points as it's a 3D cubic lattice

    //An empty array to store only unique coordinates
 /*   size_t max_unique_coords = totalParticles;
    Coordinate unique_coords[max_unique_coords];
    size_t unique_count = 0;
*/

    // Initialize a counter for indexing the arrays
    int m = 0;

    // Declare a file pointer for writing the coordinates to a file
    // Open the file "simple_cubic.xyz" for writing.
    //w modes means if the file doesn't exist, it will be created.
    //an err will be caught by the if
    FILE *dataFilePntr = fopen("fcc.dat", "w");
    if (dataFilePntr == NULL)
    {
    printf("Error opening file!\n");
    return 1;
    }

    //The first 4 lines are info for the WebGL visualizer to know # of 'particles' and box size
    // Write the total number of points to the file as the first line
    fprintf(dataFilePntr, "%i\n", totalParticles);
    // Write the lattice dimensions to the file.
    //These lines are placeholders and may not represent actual physical dimensions.
    fprintf(dataFilePntr, "%f %f\n", 0.0, N*spacing);
    fprintf(dataFilePntr, "%f %f\n", 0.0, N*spacing);
    fprintf(dataFilePntr, "%f %f\n", 0.0, N*spacing);


    // Generate the coordinates for each point in the cell
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            for(k = 0; k < N; k++){
                // Calculate the coordinates for the current point and store them in the arrays
                //1st middle
                x[m] = (float)i * spacing;
                y[m] = (float)j * spacing + halfSpacing;
                z[m] = (float)k * spacing + halfSpacing;

                //second middle
                x[m+1] = (float)i * spacing + halfSpacing;
                y[m+1] = (float)j * spacing + halfSpacing;
                z[m+1] = (float)k * spacing ;//+ spacing;

                //third middle
                x[m+2] = (float)i * spacing + spacing;
                y[m+2] = (float)j * spacing + halfSpacing;
                z[m+2] = (float)k * spacing + halfSpacing;

                //fourth middle
                x[m+3] = (float)i * spacing + halfSpacing;
                y[m+3] = (float)j * spacing + halfSpacing;
                z[m+3] = (float)k * spacing + spacing;

                //top, 5th
                x[m+4] = (float)i * spacing + halfSpacing;
                y[m+4] = (float)j * spacing + spacing;
                z[m+4] = (float)k * spacing + halfSpacing;

                //bottom, 6th
                x[m+5] = (float)i * spacing + halfSpacing;
                y[m+5] = (float)j * spacing;
                z[m+5] = (float)k * spacing + halfSpacing;
                

                //corner 1 
                x[m+6] = (float)i * spacing;
                y[m+6] = (float)j * spacing;
                z[m+6] = (float)k * spacing;

                //corner 2
                x[m+7] = (float)i * spacing + spacing;
                y[m+7] = (float)j * spacing + spacing;
                z[m+7] = (float)k * spacing + spacing;

                //corner 3
                x[m+8] = (float)i * spacing;
                y[m+8] = (float)j * spacing + spacing;
                z[m+8] = (float)k * spacing + spacing;

                //corner 4
                x[m+9] = (float)i * spacing + spacing;
                y[m+9] = (float)j * spacing;
                z[m+9] = (float)k * spacing + spacing;

                //corner 5
                x[m+10] = (float)i * spacing + spacing;
                y[m+10] = (float)j * spacing + spacing;
                z[m+10] = (float)k * spacing;

                //corner 6
                x[m+11] = (float)i * spacing;
                y[m+11] = (float)j * spacing;
                z[m+11] = (float)k * spacing + spacing;

                //corner 7
                x[m+12] = (float)i * spacing + spacing;
                y[m+12] = (float)j * spacing;
                z[m+12] = (float)k * spacing;

                //corner 8
                x[m+13] = (float)i * spacing;
                y[m+13] = (float)j * spacing + spacing;
                z[m+13] = (float)k * spacing;
                
                
                // Write the coordinates to the file
                fprintf(dataFilePntr, "%f %f %f\n", x[m], y[m], z[m]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+1], y[m+1], z[m+1]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+2], y[m+2], z[m+2]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+3], y[m+3], z[m+3]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+4], y[m+4], z[m+4]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+5], y[m+5], z[m+5]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+6], y[m+6], z[m+6]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+7], y[m+7], z[m+7]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+8], y[m+8], z[m+8]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+9], y[m+9], z[m+9]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+10], y[m+10], z[m+10]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+11], y[m+11], z[m+11]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+12], y[m+12], z[m+12]);
                fprintf(dataFilePntr, "%f %f %f\n", x[m+13], y[m+13], z[m+13]);



                // Print the coordinates to the console for verification
                printf("\n%f %f %f\n", x[m], y[m], z[m]);
                printf("%f %f %f\n", x[m+1], y[m+1], z[m+1]);
                printf("%f %f %f\n", x[m+2], y[m+2], z[m+2]);
                printf("%f %f %f\n", x[m+3], y[m+3], z[m+3]);
                printf("%f %f %f\n", x[m+4], y[m+4], z[m+4]);
                printf("%f %f %f\n", x[m+5], y[m+5], z[m+5]);
                printf("%f %f %f\n", x[m+6], y[m+6], z[m+6]);
                printf("%f %f %f\n", x[m+7], y[m+7], z[m+7]);
                printf("%f %f %f\n", x[m+8], y[m+8], z[m+8]);
                printf("%f %f %f\n", x[m+9], y[m+9], z[m+9]);
                printf("%f %f %f\n", x[m+10], y[m+10], z[m+10]);
                printf("%f %f %f\n", x[m+11], y[m+11], z[m+11]);
                printf("%f %f %f\n", x[m+12], y[m+12], z[m+12]);
                printf("%f %f %f\n", x[m+13], y[m+13], z[m+13]);


                // Increment the counter for the next cell
                m+=14;
                
            }
        }
    }

    //output some relevant info
    printf("total coords: %i, space: %f, hlfSpace: %f\n", totalParticles, spacing, halfSpacing);

    // Close the file after writing all the coordinates
    fclose(dataFilePntr);

    // Return 0 to indicate successful execution of the program
    return 0;
}