//standard input and output
#include <stdio.h>

//main func def, doesn't need any args so ()
int main()
{
	//initializing a variable to store the number in and count the loop
	int numba;

	//da loop
	for (numba = 1; numba <=100; numba++)
	{


		//conditional practice
		(numba % 2 == 0 ) ? printf("even %d\n", numba) : printf("oddball\n");

		//if/else practice
		//only print evens
		/*if (numba % 2 != 0  &&  numba % 3 != 0
			&&  numba % 5 != 0  && numba % 7 != 0)
		{
			//printing to the console
			printf("primetime %d \n", numba);
		}
*/	}
	//returning 0 if all went well
	return 0;
}