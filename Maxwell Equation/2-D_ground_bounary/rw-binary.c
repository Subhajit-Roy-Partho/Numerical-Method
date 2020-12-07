//test program for binary read-write and file positioning

#include <stdio.h>

void main (void)
{
     //array should be large enough
	int i, values[20];
	FILE *fp;
	
     printf ("Assigned values:\n");
	for (i = 0; i < 20; i++){
		values[i] = i * i;
          printf ("%i ", values[i]);
	}
	printf ("\n");
     
	//write out the values to data file
	fp = fopen ("test-rw.out", "w");
	fwrite (&values, sizeof (int), 20, fp);
	fclose (fp);

	//clear the values
	for (i = 0; i < 20; i++){
		values[i] = 0.0;
	}
	
	//reopen the file for reading
	fp = fopen ("test-rw.out", "r");
     //position the file pointer at the second half from origin
	fseek (fp, 10*sizeof(int), SEEK_SET);
	//and read in 10 values
	fread (&values, sizeof(int), 10, fp);
	//and print them out
     printf ("Second half of data file:\n");
	for (i = 0; i < 10; i++){
		printf ("%i ", values[i]);
	}
	printf ("\n");
	fclose (fp);

     //clear the values
     for (i = 0; i < 20; i++){
          values[i] = 0.0;
     }
     
	//reopen the file for reading
	fp = fopen ("test-rw.out", "r");
	//read in all values from the beginning
	fread (&values, sizeof(int), 20, fp);
	//and print them out
     printf ("All values:\n");
	for (i = 0; i < 20; i++){
		printf ("%i ", values[i]);
	}
	printf ("\n");
	fclose (fp);
}