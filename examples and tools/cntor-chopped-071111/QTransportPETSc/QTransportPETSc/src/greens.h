/*
 * greens.h
 *
 *  Created on: Jul 8, 2010
 *      Author: leond
 *
 *      This header file should contain all of the functions that are essential to creating the greens functions.
 *
 */
#include <iostream>
#include <complex>
#include <complex.h>
#include <math.h>

#ifndef GREENS_H_
#define GREENS_H_


#endif /* GREENS_H_ */



	void wiregf(float a, double e,float mtop,float ntop,int isl180[],int isr180[],int isl[], int isr[], complex <double> gl[4][4],complex <double> gr[4][4],double coord[3600][3]){

	int ig=4;
	int i,j,m,n;
	complex <double> zfunc;
	double yl180[ig], yr180[ig],zl180[ig],zr180[ig];
	double emass=.511e06;
	double hbar=1973.3e0;
	double alpha, alpha2, emn=0;
	double fmnl, fmnr,gwizz;
//
	double yo = 100e0, zo =4e0;
	double pi = acos(-1.0);
	complex <double> unitz;
	unitz= std::complex<double>(0,1);
	double t=((((hbar*hbar)/emass)/a)/a);
	double efermi=6.0e0;

	double a1x,a1y,a1z,a2x,a2y,a2z;
	complex <double> result;
	double ap;


	for (i=0;i<ig;i++){
		for (j=0;j<ig;j++){
	      gl[i][j]=0;
 	      gr[i][j]=0;
					}
							}

	for(i=0;i<4;i++){
	///
	// Get all the x,y, and z from the atoms 1 and 2
//			FILE* file = fopen("/Users/leond/Documents/Code/Research/positions", "r"); // opens the file
//			ifstream inFile;
//			size_t BUFSIZE = 1000;
//			char buf[BUFSIZE];
//			int linenum=0;
//			for(linenum=0;linenum<isl180[i]*4-4;linenum++) fscanf(file, "%s", buf); //   string<<C
//			if(EOF!=fscanf(file, "%lf %lf %lf %lf ",&a1x,&a1y,&a1z,&ap))
//			{
//	//			printf("%d : %lf %lf %lf %lf\n",++linenum,a1x,a1y,a1z,ap);
//			}
//			yl180[i]=a1y;
//			zl180[i]=a1z;
//			rewind(file);
//			for(linenum=0;linenum<isr180[i]*4-4;linenum++) fscanf(file, "%s", buf);
//			if(EOF!=fscanf(file, "%lf %lf %lf %lf ",&a2x,&a2y,&a2z,&ap))
//			{
//	//			printf("%d : %lf %lf %lf %lf\n",++linenum,a2x,a2y,a2z,ap);
//			}
//			yr180[i]=a2y;
//			zr180[i]=a2z;
//			fclose(file);
	///

	    a1x = coord[isl180[i]-1][0];
	    yl180[i] = coord[isl180[i]-1][1];
	    zl180[i] = coord[isl180[i]-1][2];

	    a2x = coord[isr180[i]-1][0];
	    yr180[i] = coord[isr180[i]-1][1];
	    zr180[i] = coord[isr180[i]-1][2];

	}



	for (i=0;i<ig;i++){
		for (j=0;j<ig;j++){
			for (m=1;m<mtop+1;m++){
				for (n=1;n<ntop+1;n++){

			  emn=(hbar*hbar*pi*pi*(((m*m)/yo)/yo+n*n/zo/zo)/2/emass);
//			  cout<<"emn-"<<emn<<endl;
			  alpha=(e+efermi-emn)/t-1.0e0;
//			  cout<<"alpha-"<<alpha<<endl;
			  alpha2=alpha*alpha;
			  fmnl=sin((m*pi*yl180[i])/yo)*sin((m*pi*yl180[j])/yo)*sin((n*pi*zl180[i])/zo)*sin((n*pi*zl180[j])/zo);
//			  cout<<"fmnl-"<<fmnl<<endl;
			  fmnr=sin((m*pi*yr180[i])/yo)*sin((m*pi*yr180[j])/yo)*sin((n*pi*zr180[i])/zo)*sin((n*pi*zr180[j])/zo);
//			  cout<<"fmnr-"<<fmnr<<endl;

//			  if (fmnr != 0 || fmnl != 0){
////				  if (fmnr < .00000001){fmnr=0;}
////				  if (fmnr < .00000001){fmnr=0;}
//				  cout<<"fmnl-"<<fmnl<<endl;
//				  cout<<"fmnr-"<<fmnr<<endl;
//			  }

			  if(alpha<-1.0e0){
			     zfunc=fabs(alpha)-sqrt(alpha2-1e0);
//				  cout<<"1st zfunc-"<<zfunc<<endl;
			  }

			  if (alpha>-1.0e0 && alpha< 1.0e0){
			     zfunc=-alpha+unitz*sqrt(1e0-alpha2);
//				  cout<<"2nd zfunc-"<<zfunc<<endl;
			  }
			  if (alpha>1.0e0){
			     zfunc=-alpha+sqrt(alpha2-1e0);
				  cout<<"3rd zfunc-"<<zfunc<<endl;
			  }

			  gl[i][j]=fmnl*zfunc+gl[i][j];
			  gr[i][j]=fmnr*zfunc+gr[i][j];
				}
			}
		}
	}

//	for(i=0;i<4;i++){
//	for(j=0;j<4;j++){
//		cout << gl[i][j] << gr[i][j]<<endl;
//	}
//	}

	for (i=0;i<ig;i++){
		for (j=0;j<ig;j++){
			  gl[i][j]=-((((8e0/yo)/zo)/a)/t)*gl[i][j];
 			  gr[i][j]=-(((((8e0)/yo)/zo)/a)/t)*gr[i][j];
		}
	}
//
//	for(i=0;i<4;i++){
//	for(j=0;j<4;j++){
//		cout << gl[i][j] << gr[i][j]<<endl;
//	}
//	cout << endl;
//	}
//
//
//	 gl[1][1] = std::complex <double>(-4.678802688135959E-002,-1.203406422832061E-002)   ;
//     gr[1][1] = std::complex <double>(-4.678802688135929E-002,-1.203406422831929E-002)   ;
//	 gl[1][3] = std::complex <double>(4.678802688135961E-002, 1.203406422832102E-002) ;
//	 gr[1][3] = std::complex <double>(4.678802688135954E-002, 1.203406422832026E-002) ;
//	 gl[3][1] = std::complex <double>(4.678802688135961E-002, 1.203406422832102E-002) ;
//	 gr[3][1] = std::complex <double>( 4.678802688135954E-002, 1.203406422832026E-002)  ;
//	 gl[3][3] = std::complex <double>(-4.678802688135972E-002,-1.203406422832143E-002) ;
//	 gr[3][3] = std::complex <double>(-4.678802688135972E-002,-1.203406422832122E-002)  ;




	}
