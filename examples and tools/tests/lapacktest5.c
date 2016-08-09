#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <complex.h>

/* Parameters */
#define N 10

typedef complex dcomp;
/* Main program */
int main() {
        /* Locals */
        int n = N, info,i,j;
		int iter = 1;
        /* Local arrays */
	
		static complex A2[N][N];
		static complex B2[N][N];
		static complex A[N][N];
		static complex U[N][N];
		static complex b[N];
		
		
			
	    srand (1337);
        static int ipiv[N];

		int max = 0;
		int min = 0;
		int mu = 3;
		int ml = 3;
		int md = ml+mu+1;
		complex sumA = 0;
		complex sumU = 0;
		complex ratio = 0;
		
//Initialize random band matrix		
		for (i=0;i<n;i++){
			for (j=0; j<n; j++){
				A2[i][j]=0;
				B2[i][j]=0;
				U[i][j]=0;				
			}
		}

		for (j=0;j<n;j++){
			max = 0 > j-mu ? 0    : j-mu ;
			min = n < j+ml+1 ? n    : j+ml+1 ;
			for(i=max;i<min;i++){
				A2[i][j] =  ((1 - ( 20.0 * rand() / ( RAND_MAX + 1.0 ) )),(1 - ( 20.0 * rand() / ( RAND_MAX + 1.0 ))));	
			                                       }
						}


		for (i=0;i<n;i++){
				for (j=0;j<n;j++){
					A[i][j]=A2[i][j];
								}			
						}


		for (i=0;i<n;i++){
			U[i][i]=1;	
						}
/// B2(ML+MU+1+i-j,j) = A(i,j) for max(1,j-MU)<=i<=min(N,j+ML)						
		for (j=0;j<n;j++){
			max = (0 < j-mu) ? j-mu : 0    ;
			min = (n < j+ml+1) ? n    : j+ml+1 ;
			for(i=max;i<min;i++){
				B2[j][md+i-j-1]= A[j][i];			
			                     }
						}
						
				 printf( "Matrix A\n" );	
				 						for( i = 0; i < n; i++ ) {
				                for( j = 0; j < n; j++ ) printf( " (%6.2f ,%6.2f) ", A[j][i] );
				                printf( "\n" );
				        }
				 						printf( "Compressed A\n" );
				        for( i = 0; i < n; i++ ) {
				                for( j = 0; j < n; j++ ) printf( " (%6.2f, %6.2f)", B2[j][i] );
				                printf( "\n" );
				        }				
		printf( "Enter General Factor\n" );				
						for (i=0;i<iter;i++){	
//		printf( "Enter General Factor\n" );				
		zgetrf_(&n, &n, A2, &n, ipiv, &info);
//		printf( "%d",info );
//		printf( "Enter Inversion\n" );
		zgetri_( &n, A2, &n, ipiv, b, &n, &info );
//		printf( "%d",info );
//		printf( "Exit inversion\n" );
			}
		printf( "\nExit inversion\n" );				
						
		 printf( "Entering Banded Solve\n" );
		for (i=0;i<iter;i++){			
		zgbsv_( &n, &ml, &mu, &n, B2, &n, ipiv, U, &n, &info );
//		printf( "%d",info );
	     }
		 printf( "Exited Banded Solve\n" );	
		
		for (i=0;i<N;i++){
			for (j=0;j<N;j++){
				sumA = A2[j][i]+sumA;
				sumU = U[j][i]+sumU;

			}
		}
		ratio = sumA/sumU;
		printf( "Accuracy - %f\n",ratio );
		printf( "%d\n",info );

		printf( "Inversed A\n" );
		        for( i = 0; i < n; i++ ) {
		                for( j = 0; j < n; j++ ) printf( "   (%6.8f, %6.2f)", A2[j][i] );
		                printf( "\n" );
		        }
		printf( "Inversed A (solve)\n" );
		        for( i = 0; i < n; i++ ) {
		                for( j = 0; j < n; j++ ) printf( "   (%6.8f, %6.2f)", U[j][i] );
		                printf( "\n" );
		        }

        exit( 0 );
}
