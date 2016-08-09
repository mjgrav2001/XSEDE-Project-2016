/* 
 * Classic trapezoid integration in an MPI contect for the SCS technical topics 
 * seminar. We calculate the integral of a function using a composite 
 * trapezium rule and each node calculates it's subpart of the entire domain.
 * 
 * Paul van der Mark
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

/* 
 * Integrate function f over the domain a to b with n intervals and 
 * trapezium width h. 
 * 
 * The function returns an aproximation of int f_a^b
 *
 */

double calcTrapezoid(double (*f)(double),double a, double b, int n, double  h) {

    double sum;
    double x;
    int    i;
        
    sum = (f(a) + f(b))/2.0;
    x = a;
    for (i = 1; i < n; i++) {
        x = x + h;
        sum = sum + f(x);
    }
    sum = sum*h;
    return sum;
}


int main(int argc, char** argv) {

    int         rank;      /* My process rank */
    int         numproc;   /* The number of processes */
    double      a = 0.0;   /* Startpoint */
    double      b = 1.0;   /* Endpoint */
    int         n = 1024;  /* Trapezoids steps */
    double      h;         /* Trapezoid base length*/
    double      local_a;   /* Local startpoint */
    double      local_b;   /* Local endpoint */
    int         local_n;   /* Local trapezoids steps */
    double      integral;  /* Integral over my interval */
    double      total;     /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;

    /* Initialize */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    if(rank == 0)  
        printf("Number of processors %d\n",numproc);

    if(argc==3) {
	float na = (float) atoi(argv[1]);
	float nb = (float) atoi(argv[2]);
	if(nb>na) {
		a = na;
		b = nb;
	}	
    }

    /* h and n are the same for each processor */
    h = (b-a)/n;    
    local_n = n/numproc; 

    /* Calculate the subdomain */
    local_a = a + rank*local_n*h;
    local_b = local_a + local_n*h;

    /* In this example we use the exponential function to integrate */
    integral = calcTrapezoid(exp, local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    if (rank == 0) {
        total = integral;
        for (source = 1; source < numproc; source++) {
            MPI_Recv(&integral, 1, MPI_DOUBLE, source, tag,
                MPI_COMM_WORLD, &status);
            total = total + integral;
        }
    } else {  
        MPI_Send(&integral, 1, MPI_DOUBLE, dest,
            tag, MPI_COMM_WORLD);
    }

    /* Print the result */
    if (rank == 0) {
        printf("With n = %d trapezoids, our estimate\n"
               "of the integral from %f to %f = %1.8f\n",
                n, a, b, total);
    }
    /* Shut down MPI */
    MPI_Finalize();

} /*  main  */






