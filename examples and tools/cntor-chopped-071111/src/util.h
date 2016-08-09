/*
 * util.h
 *
 *  Created on: Jul 13, 2010
 *      Author: leond
 *
 *   	Tests whether or not the matrix was set up correctly by copying it to a 3600 by 3600 matrix and inverting that.
 *   	Results can then be compared to the old blas code.
 *
 */
#ifndef INVERT_H_
#define INVERT_H_
#include <iostream>
#include <complex>
using namespace std;

extern "C" void zgetrf_( int*, int* , complex<double>* , int*, int* , int* );
extern "C" void zgetri_( int*, complex<double>* , int*, int* , complex<double>*, int* , int* );

void test_invert (block hmod[300][3], complex <double> ener, int isl[], int isr[],
		double bf, complex <double> gl[4][4],complex <double> gr[4][4])
{
	//	cout <<"Entered Inversion!" << endl;
	int n = 3600, lda = 3600, info, n2=10,lda2=10;
	int ipiv[n];
	static complex <double> exmat[3600][3600];
	static complex <double> exmatc[3600][3600];
	static std::complex <double> exmat2[10][10];
	complex <double> b[n];
	complex <double> gaml[4][4],gamrr[4][4];
	complex <double> gtl[4][4],gtr[4][4];
	complex <double> prodl[4][4],prodr[4][4],zres[4][4];
	double dos=0;
	complex <double> doscomplex = 0;
	complex <double> doscomplexover = 0;
	double half=.5e0;
	double therm=.03e0;
	complex <double> doscomplexdown =0;
	complex <double> fermi1,fermi2,fermif, tofe;
	double imagpart, emu1, emu2, vsdmax;
	double thop;
	thop = -0.25;


	vsdmax = .1;
    emu1 =-vsdmax/2e0;
    emu2 = vsdmax/2e0;


		for (int j=0;j<300;j++){
			for (int k=0;k<12;k++){
				for (int l=0;l<12;l++){
				exmat[k+12*j][l+12*j]=(hmod[j][1].getData(k,l)) ;
				}
			}
		}
;
		for (int j=0;j<299;j++){
			for (int k=0;k<12;k++){
				for (int l=0;l<12;l++){
				exmat[k+12*j][l+12*(j-1)]=hmod[j][2].getData(k,l);

				if (j==298){
				exmat[n-12][n-23] = hmod[j][2].getData(0,1);
				exmat[n-9][n-22] = hmod[j][2].getData(3,2);
				exmat[n-8][n-19] = hmod[j][2].getData(4,5);
				exmat[n-5][n-18] = hmod[j][2].getData(7,6);
				exmat[n-4][n-15] = hmod[j][2].getData(8,9);
				exmat[n-1][n-14] = hmod[j][2].getData(11,10);;
				}

				}
			}
		}

		for (int j=1;j<300;j++){
			for (int k=0;k<12;k++){
				for (int l=0;l<12;l++){
				exmat[k+12*(j-1)][l+12*j]=hmod[j][0].getData(k,l);
				}
			}
		}

		for (int k=0;k<12;k++){
			for (int l=0;l<12;l++){
			exmat[k][l+12*299]=hmod[0][0].getData(k,l);
			exmat[k+12*299][l]=hmod[299][2].getData(k,l);
			}
		}

		for (int i=0;i<n-11;i++){
			doscomplex= exmatc[i][i]+doscomplex;
			doscomplexover= exmatc[i][i+11]+doscomplexover;
			doscomplexdown= exmatc[i+11][i]+doscomplexdown;
		}
/////////////////////////////////////////////////////
//		long k = 0;
//		FILE* file = fopen("hmodc.out", "w");
//		for(int i = 0; i < n; i++) {
//			for(int j = 0; j < n; j++) {
//				double re = exmat[i][j].real();
//				double im = exmat[i][j].imag();
//
//				const double A = 1e-10;
//				if(fabs(re)<A) re=0e0;
//				if(fabs(im)<A) im=0e0;
//
//				if(re==0 && im==0) continue;
//				fprintf(file, "%4d,%4d: %15.5E , %15.5E\n", i+1, j+1, re, im );
//				//printf( "%15.8E , %15.8E\n", exmatc[i][j].real(), exmatc[i][j].imag() );
//				//if(++k>10) return;
//			}
//		}
//		fclose(file);
//		return;
/////////////////////////////////////////////////////

		dos = 0;

		zgetrf_(&n, &n, &(exmat[0][0]), &n, ipiv, &info);
		if (info == 0 ){
//			printf("Factorization successful!\n");
			} else {cout <<" Factorization didn't work- Info-" << info <<""<< endl;}

		zgetri_( &n, &(exmat[0][0]), &lda, ipiv, b, &n, &info );

		if (info == 0 ){
//			printf("Inverse successful!\n");
			} else {cout <<" Inversion didn't work- Info-" << info <<""<< endl;}

		for (int i=0;i<n;i++){
			imagpart = exmat[i][i].imag();
			dos = -1.0e0*imagpart/3.141592635e0+dos;
		}

		for (int i=0;i<4;i++){
			for (int j=0;j<4;j++){
				gaml[i][j]  = exmat[isl[i]-1][isr[j]-1];
				gamrr[i][j] = conj(exmat[isr[j]-1][isl[i]-1]);
//        		cout <<"gaml-"<<gaml[i][j]<<endl;
// 	        	cout <<"gamrr-"<<gamrr[i][j]<<endl;
				gtr[i][j] = -gr[i][j].imag();
				gtl[i][j] = -gl[i][j].imag();
			}
		}

		//NOTE: We need to multiply the previous four 4x4 matrices together to get tofe
				for (int i=0;i<4;i++){
					for (int j=0;j<4;j++){
						prodl[i][j]=std::complex <double>(0,0);
						for (int k=0;k<4;k++){
							prodl[i][j]=gtl[i][k]*gaml[k][j]+prodl[i][j];
						}
				}
				}
				for (int i=0;i<4;i++){
					for (int j=0;j<4;j++){
						prodr[i][j]=std::complex <double>(0,0);
						for (int k=0;k<4;k++){
							prodr[i][j]=gtr[i][k]*gamrr[k][j]+prodr[i][j];
						}
				}
				}

				for (int i=0;i<4;i++){
					for (int j=0;j<4;j++){
						zres[i][j]=std::complex <double>(0,0);
						for (int k=0;k<4;k++){
							zres[i][j]=prodl[i][k]*prodr[k][j]+zres[i][j];
						}
				}
				}

				for (int i=0;i<4;i++){
						     tofe=tofe+zres[i][i];
				}

		//
				fermi1 = 1.0e0/(exp((ener-emu1)/therm)+1.0e0);
				fermi2 = 1.0e0/(exp((ener-emu2)/therm)+1.0e0);
				fermif=2e0*abs(tofe)*(fermi1-fermi2);


		cout <<ener <<"  " << dos <<" "<< 4*thop*thop*thop*thop*tofe<<  fermif<<endl;
//		cout <<"Second="<< 4*thop*thop*thop*thop*tofe<<endl;
//		cout <<"fermif="<< fermif<<endl;
} //end test_invert

#endif /* INVERT_H_ */
