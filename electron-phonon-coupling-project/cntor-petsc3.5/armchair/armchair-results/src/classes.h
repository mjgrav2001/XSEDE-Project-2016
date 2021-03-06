//This header includes all classes used
//Contains:
//Class Parameters
//Method Parameters
//bphs, delpos, delpos0
//Class block
//Methods set_data, print_block, set_zeros, set_values
//Method AttachLeads
//Method DetachLeads
//wiregf
//
#ifndef CLASSES_H_
#define CLASSES_H_
#include "grid_jobs.h"

#include <petsc.h>
#include <petscksp.h>

#ifndef NOBOOST
#include <omp.h>
#include <boost/program_options.hpp>
namespace opt = boost::program_options;
#endif

#ifndef NOMKL
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#endif

#include <iostream>
#include <complex>
#include <new>
#include <ctime>
#include <cmath>
#include <fstream>

using namespace std;

#define blockDimension 12

class Parameters
{
public:
        Parameters(PetscInt, PetscInt, PetscInt, double, double, double, double, double, double, int, double, double, int, int, int, double, PetscInt, int, int, char[]);
	//set (used as argument)
	PetscInt rows;
	PetscInt columns;
	PetscInt idim;
	double energyMin;
	double energyMax;
	double energyStep;
	double thop;
	double vsdMin;
	double vsdMax;
	int vsdSteps;	
	double vhop;
	double BField;
	int leadSize;
	int mtop;
	int ntop;
	double a;	
	PetscInt groupSize;
	time_t invertStart;
	int useDense;
	int useZh;
	char* fileOpts;
	
	//derived
	PetscInt ndiag;
	int numEnergySteps;
	bool crossedZero;
	double energy;
	int imax1;
	int imax2;
	int imaxd21;
	int imaxd22;
	int imaxd21i;
	int imaxd22i;
	PetscInt worldSize;
	PetscInt myWorldRank;
	int numGroups;
	int myGroupNum;
	int procsInUse;
	int myGroupRank;
	int myGroupSize;
}; //end block


Parameters::Parameters(PetscInt inRows, PetscInt inColumns, PetscInt inIdim, double inEnergyMin, double inEnergyMax, double inEnergyStep, double inThop, double inVsdMin, double inVsdMax, int inVsdSteps, double inVhop, double inBField, int inLeadSize, int inMtop, int inNtop, double inA, PetscInt inGroupSize, int inUseDense, int inUseZh, char inFileOpts[])
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
	groupSize = inGroupSize;
	useDense = inUseDense;
	useZh = inUseZh;
	fileOpts = inFileOpts;
	
	energy = energyMin;
	numEnergySteps = (int) ((energyMax >= 0 && energyMin <= 0) ? (energyMax - energyMin) / energyStep + 0.5 : abs(energyMax - energyMin) / energyStep + 1.5);
	crossedZero = (energyMax >= 0 && energyMin <= 0);
	
	ndiag = rows/idim;
	imax1 = rows - idim + 1;
	imax2 = rows - idim + 2;
	imaxd21 = rows/2 + 1;
	imaxd22 = rows/2 + 2;
	imaxd21i = rows/2 - idim + 1;
	imaxd22i = rows/2 - idim + 2;
	
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
	MPI_Comm_rank(PETSC_COMM_WORLD, &myWorldRank);
	
	numGroups = (worldSize % groupSize == 0) ? worldSize/groupSize : worldSize/groupSize + 1;
	myGroupNum = myWorldRank / groupSize;
	myGroupRank = myWorldRank % groupSize;
	
	if(worldSize % groupSize == 0) {
		numGroups = worldSize / groupSize;
		myGroupSize = groupSize;
	}
	else {
		numGroups = worldSize/groupSize + 1;
		myGroupSize = (myWorldRank < (numGroups-1)*groupSize) ? groupSize : worldSize % groupSize;
	}
} //end block

#ifndef NOBOOST
Parameters applyCmdLine(Parameters d, int argc, char **argv) {
	//Replace defaults with whatever parameters are specified at runtime
	//TODO add flags for skipping or only doing the integration
	try {
		opt::options_description desc("Cntor-specific options");
		desc.add_options()
			("help,h", "Produce help message")
			("begin-energy,b", opt::value<double>(), "Energy at which to begin calculating occupancies")
			("energy-step,s", opt::value<double>(), "Step size through energy values")
			("end-energy,e", opt::value<double>(), "Ending energy")
			("n-atoms,n", opt::value<int>(), "Number of atoms in system")
			("group-size,g", opt::value<int>(), "Size of MPI process group inverting 1 matrix")
			("n-threads,t", opt::value<int>(), "Number of threads for each process")
#ifndef NOMKL
			("dense,d", "Solve with MKL dense solver rather than PETSc")
			("hermitian,z", "Solve with methods specifically for hermitian matrices")
#endif
		;

		opt::variables_map vm;
		opt::store(opt::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
		opt::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            exit(EXIT_SUCCESS);
        }
        if (vm.count("begin-energy")) {
			d.energyMin = vm["begin-energy"].as<double>();
			d.fileOpts = "a"; //Append if not starting at beginning
		}
        if (vm.count("energy-step")) {
			d.energyStep = vm["energy-step"].as<double>();
		}
        if (vm.count("end-energy")) {
			d.energyMax = vm["end-energy"].as<double>();
		}
        if (vm.count("n-atoms")) {
			d.rows = vm["n-atoms"].as<int>();
			d.columns = vm["n-atoms"].as<int>();
		}

        if (vm.count("group-size")) {
			d.groupSize = vm["group-size"].as<int>();
		}
        if (vm.count("n-threads")) {
			omp_set_num_threads(vm["n-threads"].as<int>());
		}

#ifndef NOMKL
        if (vm.count("dense")) {
        	d.useDense = 1;
        }
        if (vm.count("hermitian")) {
        	d.useZh = 1;
        }
#endif

	} catch(exception& e) {
		cerr << "WARNING: " << e.what() << endl;
	}
	//Initialize the return with any modifications from cmd line
	Parameters ret(d.rows,d.columns,d.idim, d.energyMin, d.energyMax, d.energyStep, d.thop,
			   d.vsdMin, d.vsdMax, d.vsdSteps, d.vhop, d.BField, d.leadSize, d.mtop, d.ntop, d.a,
			   d.groupSize, d.useDense, d.useZh, d.fileOpts);
	return ret;
}
#endif

//This function acts the same way as the bphs function in the old code.
//The main difference is I combined the bphs function with the earg function as they were only used together.
complex<double> bphs(int a1, int a2, double BField, double disCoordinates[][3])
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
double delpos0(int i, int j, double atomCoordinates[][3], double disCoordinates[][3])
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
double delpos(int i, int j, double atomCoordinates[][3], double disCoordinates[][3])
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
	void set_zeros (Parameters);
	void set_values (Parameters, int, int, char, double [][3], double [][3]);
	void set_data (int, int, complex<double>);
	complex<double> get_data(int i, int j){return data[i][j];}

	void AttachLeads(Parameters, int, int, int[], int[], int[], int[], complex<double>**, complex<double>**, double[][3], double[][3], int);
	void DetachLeads(Parameters, int, int, int[], int[], int[], int[], complex<double>**, complex<double>**, double[][3], double[][3], int);

	char getType(){return type;}
	int area () {return (mi*id);}
	void print_block(const PetscInt);
}; //end block


//This is a simple function. Takes the value C and places at a,b in data. Note that a and b must be less than idim
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


//Sets all data to zero
void block::set_zeros(Parameters p)
{
	//cout << " idim " << idim << " ndiag " << ndiag << endl;
	
	//Reset all data in block to zero, over riding any previous data.
	for(int j=0; j<p.idim; j++)
	{
		for(int i=0; i<p.idim; i++)
		{
			data[i][j] = complex<double>(0,0);
		} //endfor
	} //endfor
	
} //end block.set_zeros	
	

//Sets not only the data, but also the id, and mi. Furthermore, it sets the data to the values determined by its type.
//For instance, the if type = "D" then set_values sets the data along the main diagonal, the two sub-diagonals, and the two corner elements.
void block::set_values(Parameters p, int a, int b, char d, double atomCoordinates[][3], double disCoordinates[][3])
{
	mi = a;
	type = d;
	bool check;
	const PetscInt itor = p.idim/2;
	int ic[itor+1];
	id = mi + b*p.ndiag;
	
	//Numerical accuracies
	double eta2 = 5e-9;
	
	complex<double> diaener;
	diaener = complex<double>(p.energy, eta2);

	//Occupation number:
	double occupation = 103.0;
	
	//Electron-phonon coupling included:
	double vphon = -5.3;

	//NO Electron-phonon coupling:
	//double vphon = 0.0;
	
	double vnew;
	
	vphon = vphon / occupation;
	
	//cout << " id " << id << " b " << b << endl;
	//cout << " idim " << p.idim << " ndiag " << p.ndiag << endl;
	
	//Reset all data in block to zero, over riding any previous data.
	//for(int j=0; j<p.idim; j++)
	//{
	//	for(int i=0; i<p.idim; i++)
	//	{
	//		data[i][j] = 0.0;
	//	} //endfor
	//} //endfor

	switch(type)
	{
		case 'D' :
			//Setting the diagonal
			for(int i=0; i<p.idim; i++)
			{
				data[i][i] = diaener;
			} //endfor

			for(int i=1; i<p.idim; i++)
			{
				//Above diagonal
				vnew = -p.vhop + vphon*delpos(p.idim*mi+i-1+1, p.idim*mi+i+1, atomCoordinates, disCoordinates);						
				data[i-1][i] = bphs(p.idim*mi+i-1+1, p.idim*mi+i+1, p.BField, disCoordinates) * vnew;
				//Below diagonal
				vnew = -p.vhop + vphon*delpos(p.idim*mi+i+1, p.idim*mi+i-1+1, atomCoordinates, disCoordinates);					
				data[i][i-1] = bphs(p.idim*mi+i+1, p.idim*mi+i-1+1, p.BField, disCoordinates) * vnew;
			} //endfor

			vnew = -p.vhop + vphon*delpos(1, p.idim, atomCoordinates, disCoordinates);				
			data[0][p.idim-1] = vnew;
			vnew = -p.vhop + vphon*delpos(p.idim, 1, atomCoordinates, disCoordinates);				
			data[p.idim-1][0] = vnew;
			break;

		case 'W' :
			//Setting the t-dagger blocks.
			check = true;
			for(int j=1; j<p.idim; j+=2)
			{
				if(check == true)
				{
					vnew = -p.vhop + vphon*delpos(mi*p.idim+j+p.idim, mi*p.idim+j+1, atomCoordinates, disCoordinates);					
					data[j-1][j] = bphs(mi*p.idim+j+p.idim, mi*p.idim+j+1, p.BField, disCoordinates) * vnew; //1,5,9
					check = !check;
				} //endif
				else if(check == false)
				{
					vnew = -p.vhop + vphon*delpos(mi*p.idim+j+p.idim+1, mi*p.idim+j, atomCoordinates, disCoordinates);		
					data[j][j-1] = bphs(mi*p.idim+j+p.idim+1, mi*p.idim+j, p.BField, disCoordinates) * vnew; //3,7,11
					check = !check;
				} //endif
			} //endfor
			break;

		case 'T' :
			//Setting the T blocks
			check = true;
			for(int j=1; j<p.idim; j+=2)
			{
				if(check == true)
				{
					vnew = -p.vhop + vphon*delpos(mi*p.idim+j+p.idim, mi*p.idim+j+1, atomCoordinates, disCoordinates);	
					data[j][j-1] = conj(bphs(mi*p.idim+j+p.idim, mi*p.idim+j+1, p.BField, disCoordinates)) * vnew; //3,7,11
					check = !check;
				} //endif
				else if(check == false)
				{
					vnew = -p.vhop + vphon*delpos(mi*p.idim+j+p.idim+1, mi*p.idim+j, atomCoordinates, disCoordinates);	
					data[j-1][j] = conj(bphs(mi*p.idim+j+p.idim+1, mi*p.idim+j, p.BField, disCoordinates)) * vnew; //1,5,9
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
					vnew = -p.vhop + vphon*delpos(2*j-1, p.rows-p.idim+2*(j-1)+2, atomCoordinates, disCoordinates);	
					data[ic[j]-1][ic[j]] = bphs(2*j-1, p.rows-p.idim+2*(j-1)+2, p.BField, disCoordinates) * vnew;
					check = !check;
				} //endif
				else if(check == false)
				{
					ic[j] = 2*(j-1);
					vnew = -p.vhop + vphon*delpos(2*j, p.rows-p.idim+2*(j-1)+1, atomCoordinates, disCoordinates);	
					data[ic[j]+1][ic[j]] = bphs(2*j, p.rows-p.idim+2*(j-1)+1, p.BField, disCoordinates) * vnew;
					check = !check;
				} //endif
			} //endfor
			break;

		case 'Q' :
			//Setting the u dagger blocks
			check = true;
			for(int j=1; j<=itor; j++)
			{
				if(check == true)
				{
					ic[j] = 2*j-1;
					vnew = -p.vhop + vphon*delpos(p.rows-p.idim+2*(j-1)+2, 2*j-1, atomCoordinates, disCoordinates);	
					data[ic[j]][ic[j]-1] = bphs(p.rows-p.idim+2*(j-1)+2, 2*j-1, p.BField, disCoordinates) * vnew;
					check = !check;
				} //endif
				else if(check == false)
				{
					ic[j] = 2*(j-1);
					vnew = -p.vhop + vphon*delpos(p.rows-p.idim+2*(j-1)+1, 2*j, atomCoordinates, disCoordinates);	
					data[ic[j]][ic[j]+1] = bphs(p.rows-p.idim+2*(j-1)+1, 2*j, p.BField, disCoordinates) * vnew;
					check = !check;
				} //endif
			} //endfor
			break;
	}
} //end block.set_values


void block::AttachLeads(Parameters p, int aibloc, int ajbloc, int *isl180, int *isr180, int *isl, int *isr, complex<double> **gl, complex<double> **gr, double atomCoordinates[][3], double disCoordinates[][3], int check)
{
	int temp1, temp2, temp3, temp4,
	    a1bloc, a2bloc, a3bloc, a4bloc,
	    a1pl, a2pl, a3pl, a4pl;
	
	complex<double> tempcd1, tempcd2;
	complex<double> ham;
	
	if ((check == 1)&&(ajbloc == 0))
	  {
	    aibloc=aibloc+1;
	  }
	
	for(int i=0; i<p.leadSize; i++)
	  {
	    for(int k=0; k<p.leadSize; k++)
	      {
		temp1 = isl180[i]-1;
		temp2 = isl180[k]-1;
		temp3 = isr180[i]-1;
		temp4 = isr180[k]-1;
		
		a1pl  = temp1 % p.idim;
		a2pl  = temp2 % p.idim;
		a3pl  = temp3 % p.idim;
		a4pl  = temp4 % p.idim;
		
		if(temp1 >= temp2)
		  {
		    a1bloc = (int) floor(temp1 / p.idim);
		    a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
		  } //endif
		else
		  {
		    a1bloc = (int) floor(temp2 / p.idim);
		    a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);	
		  } //endelse
		
		if(a2bloc == 3)
		  {
		    //In U
		    //Set a2bloc and a1bloc to U's position in blockHamiltonian
		    a2bloc = 0;
		    a1bloc = 0;
		  } //endif
		
		if(a2bloc == 4)
		  {
		    //In Q
		    a2bloc = 2;
		    a1bloc = p.ndiag-1;
		  } //endif
		
		if ((aibloc == a1bloc)&&(ajbloc == a2bloc))
		  {
		    tempcd1 = get_data(a1pl, a2pl);
		    set_data(a1pl, a2pl, tempcd1 - p.thop*p.thop*gl[i][k]);
		    ham = get_data(a1pl, a2pl);
		    data[a1pl][a2pl] = ham;
		  } //endif
		
		if(temp3 >= temp4)
		  {
		    a3bloc = (int) floor(temp3 / p.idim);
		    a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
		  } //endif
		else
		  {
		    a3bloc = (int) floor(temp4 / p.idim);
		    a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);	
		  } //endelse
		
		if(a4bloc == 3)
		  {
		    //In U
		    a4bloc = 0;
		    a3bloc = 0;
		  } //endif
		
		if(a4bloc == 4)
		  {
		    //In Q
		    a4bloc = 2;
		    a3bloc = p.ndiag-1;
		  } //endif
		
		if ((aibloc == a3bloc)&&(ajbloc == a4bloc))
		  {
		    tempcd2 = get_data(a3pl, a4pl);
		    set_data(a3pl, a4pl, tempcd2 - p.thop*p.thop*gr[i][k]);
		    ham = get_data(a3pl, a4pl);
		    data[a3pl][a4pl] = ham;
		  } //endif	
	      } //endfor
	  } //endfor
} //end AttachLeads


void block::DetachLeads(Parameters p, int aibloc, int ajbloc, int *isl180, int *isr180, int *isl, int *isr, complex<double> **gl, complex<double> **gr, double atomCoordinates[][3], double disCoordinates[][3], int check)
{
	int temp1, temp2, temp3, temp4,
	    a1bloc, a2bloc, a3bloc, a4bloc,
            a1pl, a2pl, a3pl, a4pl;
	
	complex<double> tempcd1, tempcd2;
	complex<double> ham;
	
	if ((check == 1)&&(ajbloc == 0))
	  {
	    aibloc=aibloc+1;
	  }
	
	for(int i=0; i<p.leadSize; i++)
	  {
	    for(int k=0; k<p.leadSize; k++)
	      {
		temp1 = isl180[i]-1;
		temp2 = isl180[k]-1;
		temp3 = isr180[i]-1;
		temp4 = isr180[k]-1;
			
		a1pl  = temp1 % p.idim;
		a2pl  = temp2 % p.idim;
		a3pl  = temp3 % p.idim;
		a4pl  = temp4 % p.idim;
		
		if(temp1 >= temp2)
		  {
		    a1bloc = (int) floor(temp1 / p.idim);
		    a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
		  } //endif
		else
		  {
		    a1bloc = (int) floor(temp2 / p.idim);
		    a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);	
		  } //endelse
		
		if(a2bloc == 3)
		  {
		    //In U
		    //Set a2bloc and a1bloc to U's position in blockHamiltonian
		    a2bloc = 0;
		    a1bloc = 0;
		  } //endif
		
		if(a2bloc == 4)
		  {
		    //In Q
		    a2bloc = 2;
		    a1bloc = p.ndiag-1;
		  } //endif
		
		if ((aibloc == a1bloc)&&(ajbloc == a2bloc))
		  {
		    tempcd1 = get_data(a1pl, a2pl);
		    set_data(a1pl, a2pl, tempcd1 + p.thop*p.thop*gl[i][k]);
		    ham = get_data(a1pl, a2pl);
		    data[a1pl][a2pl] = ham;
		  } //endif
		
		
		if(temp3 >= temp4)
		  {
		    a3bloc = (int) floor(temp3 / p.idim);
		    a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
		  } //endif
		else
		  {
		    a3bloc = (int) floor(temp4 / p.idim);
		    a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);	
		  } //endelse
		
		if(a4bloc == 3)
		  {
		    //In U
		    a4bloc = 0;
		    a3bloc = 0;
		  } //endif
			
		if(a4bloc == 4)
		  {
		    //In Q
		    a4bloc = 2;
		    a3bloc = p.ndiag-1;
		  } //endif
		
		if ((aibloc == a3bloc)&&(ajbloc == a4bloc))
		  {
		    tempcd2 = get_data(a3pl, a4pl);
		    set_data(a3pl, a4pl, tempcd2 + p.thop*p.thop*gr[i][k]);
		    ham = get_data(a3pl, a4pl);
		    data[a3pl][a4pl] = ham;
		  } //endif		
	      } //endfor
	  } //endfor
} //end DetachLeads


#undef __FUNCT__
#define __FUNCT__ "wiregf"
PetscErrorCode wiregf(Parameters p, int *isl180, int *isr180, int *isl, int *isr, complex<double> **gl, complex<double> **gr, double disCoordinates[][3])
{
	int ig = p.leadSize;
	
	complex<double> zfunc, result, unitz;
	
	unitz = complex<double>(0,1);
	
	double yl180[ig], yr180[ig], zl180[ig], zr180[ig];
	double alpha, alpha2, fmnl, fmnr, a1x, a2x;
	
	double emass = .511e06, hbar = 1973.3, emn = 0.0, efermi = 6.0, yo = 100.0, zo = 4.0;
	
	double t = hbar*hbar/emass/p.a/p.a;
	double pi = acos(-1.0);
	//double pi = M_PI;
	
	for(int i=0; i<ig; i++)
	{
		for(int j=0; j<ig; j++)
		{
			gl[i][j] = complex<double>(0,0);
			gr[i][j] = complex<double>(0,0);
		} //endfor
	} //endfor
	
	for(int i=0; i<ig; i++)
	{
	    a1x = disCoordinates[isl180[i]-1][0];
	    yl180[i] = disCoordinates[isl180[i]-1][1];
	    zl180[i] = disCoordinates[isl180[i]-1][2];
		
	    a2x = disCoordinates[isr180[i]-1][0];
	    yr180[i] = disCoordinates[isr180[i]-1][1];
	    zr180[i] = disCoordinates[isr180[i]-1][2];
	} //endfor
	
	for(int i=0; i<ig; i++)
	{
		for(int j=0; j<ig; j++)
		{
			for(int m=1; m<p.mtop+1; m++)
			{
				for(int n=1; n<p.ntop+1; n++)
				{
					emn = hbar*hbar*pi*pi* (m*m/yo/yo + n*n/zo/zo)/2.0/emass;
					alpha = (p.energy + efermi - emn)/t - 1.0;
					alpha2 = alpha * alpha;
					fmnl = sin(m*pi*yl180[i]/yo) * sin(m*pi*yl180[j]/yo) * sin(n*pi*zl180[i]/zo) * sin(n*pi*zl180[j]/zo);
					fmnr = sin(m*pi*yr180[i]/yo) * sin(m*pi*yr180[j]/yo) * sin(n*pi*zr180[i]/zo) * sin(n*pi*zr180[j]/zo);
					
					if(alpha < -1.0)
					{
						zfunc = fabs(alpha) - sqrt(alpha2 - 1.0);
					} //endif
					
					if(alpha >= -1.0 && alpha <= 1.0)
					{
						zfunc = -alpha + unitz*sqrt(1.0 - alpha2);
					} //endif
					
					if(alpha > 1.0)
					{
						zfunc = -alpha + sqrt(alpha2 - 1.0);
						//						cout << "3rd zfunc-" << zfunc << endl;
					} //endif
					
					gl[i][j] += fmnl*zfunc;
					gr[i][j] += fmnr*zfunc;
				} //endfor
			} //endfor
		} //endfor
	} //endfor
	
	for(int i=0; i<ig; i++)
	{
	  for(int j=0; j<ig; j++)
	    {
	      gl[i][j] *= -(8.0/yo/zo/p.a/t);
	      gr[i][j] *= -(8.0/yo/zo/p.a/t);
	    } //endfor
	} //endfor
	
	PetscFunctionReturn(0);
} //end wiregf

#endif /* CLASSES_H_ */
