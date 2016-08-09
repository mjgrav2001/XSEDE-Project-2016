//This header includes all classes used
#ifndef CLASSES_H_
#define CLASSES_H_
#include <iostream>
#include <complex>
#include <ctime>
#include "petsc.h"
using namespace std;

#define blockDimension 40

class Parameters
{
public:
	Parameters(PetscInt, PetscInt, PetscInt, double, double, double, double, double, double, int, double, double, int, int, int, double);
	//set (used as argument)
	PetscInt rows;
	PetscInt columns;
	PetscInt idim;
	double energyMin;
	double energyMax;
	double energyStep;
	double thop;
	double vsdMax;
	double vhop;
	double BField;
	int leadSize;
	int mtop;
	int ntop;
	time_t invertStart;
	
	//set (non argument)
	double vsdMin;
	double a;
	int vsdSteps;
	
	//derived
	PetscInt ndiag;
	int numEnergySteps;
	double energy;
	int imax1;
	int imax2;
	int imaxd21;
	int imaxd22;
	int imaxd21i;
	int imaxd22i;
}; //end block

Parameters::Parameters(PetscInt inRows, PetscInt inColumns, PetscInt inIdim, double inEnergyMin, double inEnergyMax, double inEnergyStep, double inThop,
					   double inVsdMin, double inVsdMax, int inVsdSteps, double inVhop, double inBField, int inLeadSize, int inMtop, int inNtop, double inA)
{
	rows = inRows;
	columns = inColumns;
	idim = inIdim;
	energyMin = inEnergyMin;
	energyMax = inEnergyMax;
	energyStep = inEnergyStep;
	thop = inThop;
	vsdMin = inVsdMin;
	vsdMax = inVsdMax;
	vsdSteps = inVsdSteps;
	vhop = inVhop;
	BField = inBField;
	leadSize = inLeadSize;
	mtop = inMtop;
	ntop = inNtop;
	a = inA;
	
	energy = energyMin;
	ndiag = rows/idim;
	numEnergySteps = (int) ((energyMax >= 0 && energyMin <= 0) ? (energyMax - energyMin) / energyStep + 0.5 
							: abs(energyMax - energyMin) / energyStep + 1.5);
	
	imax1 = rows - idim + 1;
	imax2 = rows - idim + 2;
	imaxd21 = rows/2 + 1;
	imaxd22 = rows/2 + 2;
	imaxd21i = rows/2 - idim + 1;
	imaxd22i = rows/2 - idim + 2;	
} //end block.block

//This function acts the same way as the bphs function in the old code.
//The main difference is I combined the bphs function with the earg function as they were only used together.
complex<double> bphs(int a1, int a2, double BField, double** disCoordinates)
{
	double a1x, a1y, a1z, a2x, a2y, a2z;
	double phi1, phi2, ri1, ri2, rhoij;
	double cearg, searg, earg = 0;
	complex<double> result;
	
    a1x = disCoordinates[a1-1][0];
    a1y = disCoordinates[a1-1][1];
    a1z = disCoordinates[a1-1][2];
    a2x = disCoordinates[a2-1][0];
    a2y = disCoordinates[a2-1][1];
    a2z = disCoordinates[a2-1][2];
	
    rhoij = (a1x + a2x)*(a1x + a2x) + (a1y + a2y)*(a1y + a2y);
	rhoij = .5 * sqrt(rhoij);
	
	if(rhoij < 0.00001)
		earg = 0.0;
	else
	{
		phi1 = -(a1y + a2y)/rhoij/2.0;
		phi2 = (a1x + a2x)/rhoij/2.0;
		ri1 = a1x - a2x;
		ri2 = a1y - a2y;
		earg = phi1*ri1 + phi2*ri2;
	} //endelse
	
	earg = .5 * BField * rhoij * earg;
	earg = 1.52e-5 * earg;
	searg = sin(earg);
	cearg = cos(earg);
	result = complex<double>(cearg, searg);
	
	return result;
} //end bphs


//Old: Before 09/24/11
double delpos0(int i, int j, double** atomCoordinates, double** disCoordinates)
{
	//const int acc=1.42;
	double del, delnew, delta, deltanew;
	double result;
	double scale;
	
	scale=50.;
	
	del = pow(atomCoordinates[i-1][0]-atomCoordinates[j-1][0],2)
		+ pow(atomCoordinates[i-1][1]-atomCoordinates[j-1][1],2)
		+ pow(atomCoordinates[i-1][2]-atomCoordinates[j-1][2],2);
	delta = sqrt(del);
	
	delnew = pow(disCoordinates[i-1][0]-disCoordinates[j-1][0],2)
		   + pow(disCoordinates[i-1][1]-disCoordinates[j-1][1],2)
		   + pow(disCoordinates[i-1][2]-disCoordinates[j-1][2],2);
	deltanew = sqrt(delnew)/scale;
	
	result = abs(deltanew-delta);

	//cout << " delta " << result << " " << i-1 << " " << j-1 << endl;

	return result;
} //end delpos0


//New: Since 09/24/11
double delpos(int i, int j, double** atomCoordinates, double** disCoordinates)
{
	//const int acc=1.42;
	double del, delnew, delnorm;
	double deltai[3], deltaj[3];
	double result;
	double scale;
	
	scale=50.;
	
	deltai[0] = (disCoordinates[i-1][0]-atomCoordinates[i-1][0])/scale;
	deltai[1] = (disCoordinates[i-1][1]-atomCoordinates[i-1][1])/scale;
	deltai[2] = (disCoordinates[i-1][2]-atomCoordinates[i-1][2])/scale;

	deltaj[0] = (disCoordinates[j-1][0]-atomCoordinates[j-1][0])/scale;
	deltaj[1] = (disCoordinates[j-1][1]-atomCoordinates[j-1][1])/scale;
	deltaj[2] = (disCoordinates[j-1][2]-atomCoordinates[j-1][2])/scale;
	
	delnew = (atomCoordinates[i-1][0]-atomCoordinates[j-1][0])*(deltai[0]-deltaj[0])
	+ (atomCoordinates[i-1][1]-atomCoordinates[j-1][1])*(deltai[1]-deltaj[1])
	+ (atomCoordinates[i-1][2]-atomCoordinates[j-1][2])*(deltai[2]-deltaj[2]);
	
	del = pow(atomCoordinates[i-1][0]-atomCoordinates[j-1][0],2)
	+ pow(atomCoordinates[i-1][1]-atomCoordinates[j-1][1],2)
	+ pow(atomCoordinates[i-1][2]-atomCoordinates[j-1][2],2);
	delnorm = sqrt(del);

	if(i == j)
	{
		result = 0.;
	} 
	else {
		result = delnew / delnorm;
	}
	
	//cout << " delta " << result << " " << i-1 << " " << j-1 << endl;
	
	return result;
} //end delpos


class block
{
	int mi, id;
	complex<double> data[blockDimension][blockDimension];
	char type;
  public:
	void set_values (const PetscInt, const PetscInt, const PetscInt, double, int, int, double, double, char, 
					 double** atomCoordinates, double** disCoordinates);
	void set_data (int, int, complex<double>);
	complex<double> getData(int i, int j){return data[i][j];}
	char getType(){return type;}
	int area () {return (mi*id);}
	void print_block(const PetscInt);
}; //end block


//This is a simple function. Takes the value C and places at a,b in data. Note that a and b much be less than 12
void block::set_data(int a, int b, complex<double> c)
{
	data[a][b] = c;
}  //end block.set_data


//Prints the 12 by 12 data block to the screen.
void block::print_block(const PetscInt idim)
{
	for(int i=0; i<idim; i++)
	{
		for(int j=0; j<idim; j++)
		{
			cout << data[i][j];
			cout << endl;
		} //endfor
	} //endfor
} //end block.print_block


//Sets not only the data, but also the id, and mi. Furthermore, it sets the data to the values determined by its type.
//For instance, the if type = "D" then set_values sets the data along the main diagonal, the two sub-diagonals, and the two corner elements.
void block::set_values(const PetscInt rows, const PetscInt idim, const PetscInt ndiag, double vhop, 
					   int a, int b, double c, double BField, char d, double** atomCoordinates, double** disCoordinates)
{
	mi = a;
	type = d;
	bool check;
	const PetscInt itor = idim/2;
    int ic[itor+1];
	id = mi + b*ndiag;
	
	//Numerical accuracies
	double eta2=5e-9;
	
	complex<double> diaener;
	diaener = complex<double>(c, eta2);

	//Occupation number:
	double occupation=103.0;
	
	//Electron-phonon coupling included:
	double vphon=-5.3;

	//NO Electron-phonon coupling:
	//double vphon=0.0;
	
	double vnew;
	
	vphon = vphon / occupation;
	
	//cout << " id " << id << " b " << b << endl;
	//cout << " idim " << idim << " ndiag " << ndiag << endl;
	
	//Reset all data in block to zero, over riding any previous data.
	for(int j=0; j<idim; j++)
	{
		for(int i=0; i<idim; i++)
		{
			data[i][j] = 0.0;
		} //endfor
	} //endfor

	switch(type)
	{
		case 'D' :
			//Setting the diagonal
			for(int i=0; i<idim; i++)
			{
				data[i][i] = diaener;
			} //endfor

			for(int i=1; i<idim; i++)
			{
				//Above diagonal
				vnew = -vhop + vphon*delpos(idim*mi+i-1+1, idim*mi+i+1, atomCoordinates, disCoordinates);						
				data[i-1][i] = bphs(idim*mi+i-1+1, idim*mi+i+1, BField, disCoordinates) * vnew;
				//Below diagonal
				vnew = -vhop + vphon*delpos(idim*mi+i+1, idim*mi+i-1+1, atomCoordinates, disCoordinates);					
				data[i][i-1] = bphs(idim*mi+i+1, idim*mi+i-1+1, BField, disCoordinates) * vnew;
			} //endfor

			vnew = -vhop + vphon*delpos(1, idim, atomCoordinates, disCoordinates);				
			data[0][idim-1] = vnew;
			vnew = -vhop + vphon*delpos(idim, 1, atomCoordinates, disCoordinates);				
			data[idim-1][0] = vnew;
			break;

		case 'W' :
			//Setting the t-dagger blocks.
			check = true;
			for(int j=1; j<idim; j+=2)
			{
				if(check == true)
				{
					vnew = -vhop + vphon*delpos(mi*idim+j+idim, mi*idim+j+1, atomCoordinates, disCoordinates);					
					data[j-1][j] = bphs(mi*idim+j+idim, mi*idim+j+1, BField, disCoordinates) * vnew; //1,5,9
					check = !check;
				} //endif
				else if(check == false)
				{
					vnew = -vhop + vphon*delpos(mi*idim+j+idim+1, mi*idim+j, atomCoordinates, disCoordinates);		
					data[j][j-1] = bphs(mi*idim+j+idim+1, mi*idim+j, BField, disCoordinates) * vnew; //3,7,11
					check = !check;
				} //endif
			} //endfor
			break;

		case 'T' :
			//Setting the T blocks
			check = true;
			for(int j=1; j<idim; j+=2)
			{
				if(check == true)
				{
					vnew = -vhop + vphon*delpos(mi*idim+j+idim, mi*idim+j+1, atomCoordinates, disCoordinates);	
					data[j][j-1] = conj(bphs(mi*idim+j+idim, mi*idim+j+1, BField, disCoordinates)) * vnew; //3,7,11
					check = !check;
				} //endif
				else if(check == false)
				{
					vnew = -vhop + vphon*delpos(mi*idim+j+idim+1, mi*idim+j, atomCoordinates, disCoordinates);	
					data[j-1][j] = conj(bphs(mi*idim+j+idim+1, mi*idim+j, BField, disCoordinates)) * vnew; //1,5,9
					check = !check;
				} //endif
			} //endfor
			break;

		case 'U' :
			//Setting the U blocks
			check = true;
			for(int j=1; j<=itor; j++)
			{
				if(check == true)
				{
					ic[j] = 2*j-1;
					vnew = -vhop + vphon*delpos(2*j-1, rows-idim+2*(j-1)+2, atomCoordinates, disCoordinates);	
					data[ic[j]-1][ic[j]] = bphs(2*j-1, rows-idim+2*(j-1)+2, BField, disCoordinates) * vnew;
					check = !check;
				} //endif
				else if(check == false)
				{
					ic[j] = 2*(j-1);
					vnew = -vhop + vphon*delpos(2*j, rows-idim+2*(j-1)+1, atomCoordinates, disCoordinates);	
					data[ic[j]+1][ic[j]] = bphs(2*j, rows-idim+2*(j-1)+1, BField, disCoordinates) * vnew;
					check = !check;
				} //endif
			} //endfor
			//data[0][1] = bphs(1, 3590, BField, atomCoordinates) * vnew;
			//data[3][2] = bphs(4, 3591, BField, atomCoordinates) * vnew;
			//data[4][5] = bphs(5, 3594, BField, atomCoordinates) * vnew;
			//data[7][6] = bphs(8, 3595, BField, atomCoordinates) * vnew;
			//data[8][9] = bphs(9, 3598, BField, atomCoordinates) * vnew;
			//data[11][10] = bphs(idim, 3599, BField, atomCoordinates) * vnew;
			break;

		case 'Q' :
			//Setting the u dagger blocks
			check = true;
			for(int j=1; j<=itor; j++)
			{
				if(check == true)
				{
					ic[j] = 2*j-1;
					vnew = -vhop + vphon*delpos(rows-idim+2*(j-1)+2, 2*j-1, atomCoordinates, disCoordinates);	
					data[ic[j]][ic[j]-1] = bphs(rows-idim+2*(j-1)+2, 2*j-1, BField, disCoordinates) * vnew;
					check = !check;
				} //endif
				else if(check == false)
				{
					ic[j] = 2*(j-1);
					vnew = -vhop + vphon*delpos(rows-idim+2*(j-1)+1, 2*j, atomCoordinates, disCoordinates);	
					data[ic[j]][ic[j]+1] = bphs(rows-idim+2*(j-1)+1, 2*j, BField, disCoordinates) * vnew;
					check = !check;
				} //endif
			} //endfor
			//data[1][0] = bphs(3590, 1, BField, atomCoordinates) * vnew;
			//data[2][3] = bphs(3591, 4, BField, atomCoordinates) * vnew;
			//data[5][4] = bphs(3594, 5, BField, atomCoordinates) * vnew;
			//data[6][7] = bphs(3595, 8, BField, atomCoordinates) * vnew;
			//data[9][8] = bphs(3598, 9, BField, atomCoordinates) * vnew;
			//data[10][11] = bphs(3599, idim, BField, atomCoordinates) * vnew;
			break;
	}
} //end block.set_values

#endif /* CLASSES_H_ */
