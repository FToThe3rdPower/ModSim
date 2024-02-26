// Include standard input and output library for IO operations
#include <stdio.h>
// Include the standard library for general purpose functions like memory allocation, process control, conversions, etc.
#include <stdlib.h>

// Entry point of the program with command line arguments
// argc: argument count, argv: argument vector (list of arguments)
int main(int argc, char *argv[]){
    // Convert the first command line argument to an integer to
    //determine the number of points along one edge of the cube
    int N = atoi(argv[1]);

    //determining the total number of particles
    int totalParticles = N*N*N;

    // Set the spacing between points in the lattice
    float spacing = 2.0;
    //float spacing = atoi(argv[2]);
    float halfSpacing = spacing / 2;

    //number of full particles per unit cell, 4 for FCC
    int particlePerCell = 4;

    //determine the number of unit cells based the total particles
    int numOfCells = totalParticles / particlePerCell; 


    // Declare loop counter vars
    int c, i, j, k;

    // Declare arrays to store the x, y, and z coordinates of the points
    //in the lattice
    float x[totalParticles];
    float y[totalParticles];
    float z[totalParticles]; // These arrays can store N^3 points as it's a 3D cubic lattice

    // Initialize a counter for indexing the arrays, counts the corners
    int m = 0;

    // Declare a file pointer for writing the coordinates to a file
    FILE *print_coordinates;
    // Open the file "simple_cubic.xyz" for writing.
    //w modes means if the file doesn't exist, it will be created.
    print_coordinates = fopen("fcc.dat", "w");

    // Write the total number of points to the file as the first line
    fprintf(print_coordinates, "totalParticles: %i, space: %f, hlfSpace: %f\n", totalParticles, spacing, halfSpacing);
    // Write the lattice dimensions to the file.
    //These lines are placeholders and may not represent actual physical dimensions.
    fprintf(print_coordinates, "%f %f\n", 0.0, N*spacing);
    fprintf(print_coordinates, "%f %f\n", 0.0, N*spacing);
    fprintf(print_coordinates, "%f %f\n", 0.0, N*spacing);


    // Generate the coordinates for each point in the cell
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            for(k = 0; k < N; k++){
                // Calculate the coordinates for the current point and store them in the arrays
                //1st corner
                x[m] = (float)i * spacing;
                y[m] = (float)j * spacing + halfSpacing;
                z[m] = (float)k * spacing + halfSpacing;

                //second corner
                x[m+1] = (float)i * spacing + halfSpacing;
                y[m+1] = (float)j * spacing +halfSpacing;;
                z[m+1] = (float)k * spacing ;//+ spacing;

                //third corner
                x[m+2] = (float)i * spacing + spacing;//+ spacing;
                y[m+2] = (float)j * spacing + halfSpacing ;//- halfSpacing;
                z[m+2] = (float)k * spacing + halfSpacing;// + halfSpacing;

                //fourth corner
                x[m+3] = (float)i * spacing + halfSpacing;
                y[m+3] = (float)j * spacing + halfSpacing;
                z[m+3] = (float)k * spacing + spacing;
                
                
                // Write the coordinates to the file
                fprintf(print_coordinates, "%f %f %f\n", x[m], y[m], z[m]);
                fprintf(print_coordinates, "%f %f %f\n", x[m+1], y[m+1], z[m+1]);
                fprintf(print_coordinates, "%f %f %f\n", x[m+2], y[m+2], z[m+2]);
                fprintf(print_coordinates, "%f %f %f\n", x[m+3], y[m+3], z[m+3]);
                // Print the coordinates to the console for verification
                printf("\n%f %f %f\n", x[m], y[m], z[m]);
                printf("%f %f %f\n", x[m+1], y[m+1], z[m+1]);
                printf("%f %f %f\n", x[m+2], y[m+2], z[m+2]);
                printf("%f %f %f\n", x[m+3], y[m+3], z[m+3]);
                // Increment the counter for the next cell
                m++;
                
            }
        }
    }

    // Close the file after writing all the coordinates
    fclose(print_coordinates);

    // Return 0 to indicate successful execution of the program
    return 0;
}