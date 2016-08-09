//Contains general functions
//ReadInCoordinates
//ReadInDisplaced
//AttachLeads
//DetachLeads
//wiregf
#ifndef GENERAL_H_
#define GENERAL_H_
#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

void ReadInCoordinates(const PetscInt rows, double** atomCoordinates)
{
	ifstream inFile("ARRAY_posr.dat");
	//ifstream inFile("positions");
	int atomNum;

	if(inFile.is_open())
	{
		for(int lineNum=1; (lineNum<=rows && !inFile.eof()); lineNum++)
		{
			inFile >> atomCoordinates[lineNum-1][0] >> atomCoordinates[lineNum-1][1] >> atomCoordinates[lineNum-1][2] >> atomNum;
		} //endfor
		inFile.clear();              // forget we hit the end of file
		inFile.seekg(0, ios::beg); 
		inFile.close();
	} //endif
	else cout << "Error: Positions file not opened!\n";
} //end ReadInCoordinates


void ReadInDisplaced(const PetscInt rows, double** disCoordinates)
{
	ifstream inFile("ARRAY_posr_new.dat");
	//ifstream inFile("positions");
	int disNum;
	
	if(inFile.is_open())
	{
		for(int lineNum=1; (lineNum<=rows && !inFile.eof()); lineNum++)
		{
			inFile >> disCoordinates[lineNum-1][0] >> disCoordinates[lineNum-1][1] >> disCoordinates[lineNum-1][2] >> disNum;
		} //endfor
		inFile.clear();              // forget we hit the end of file
		inFile.seekg(0, ios::beg); 
		inFile.close();
	} //endif
	else cout << "Error: Displaced positions file not opened!\n";
} //end ReadInDisplaced


void AttachLeads(Parameters p, block** stm, int isl180[], int isr180[], int isl[], int isr[], 
				 complex<double> gl[4][4], complex<double> gr[4][4])
{
	int temp1, temp2, temp3, temp4,
	a1bloc, a2bloc, a3bloc, a4bloc,
	a1pl, a2pl, a3pl, a4pl;
	
	complex<double> tempcd1, tempcd2;

	int ig = p.leadSize;
	
	for(int i=0; i<ig; i++)
	{
		for(int j=0; j<ig; j++)
		{
			temp1 = isl180[i]-1;
			temp2 = isl180[j]-1;
			temp3 = isr180[i]-1;
			temp4 = isr180[j]-1;
			
			if(temp1 >= temp2)
			{
				a1bloc = (int) floor(temp1 / p.ndiag);
				a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
			} //endif
			else
			{
				a1bloc = (int) floor(temp2 / p.ndiag);
				a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
			} //endelse
			
			if(temp3 >= temp4)
			{
				a3bloc = (int) floor(temp3 / p.ndiag);
				a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
			} //endif
			else
			{
				a3bloc = (int) floor(temp4 / p.ndiag);
				a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
			} //endelse
			
			if(a2bloc == 3)
			{
				//If in U
				//Set a2bloc and a1bloc to U's position in blockHamiltonian
				a2bloc = 0;
				a1bloc = 0;
			} //endif
			
			if(a2bloc == 4)
			{
				a2bloc = 2;
				a1bloc = p.idim-1;
			} //endif
			
			if(a4bloc == 3)
			{
				a4bloc = 0;
				a3bloc = 0;
			} //endif
			
			if(a4bloc == 4)
			{
				a4bloc = 2;
				a3bloc = p.idim-1;
			} //endif
			
			a1pl  = temp1 % p.ndiag;
			a2pl  = temp2 % p.ndiag;
			a3pl  = temp3 % p.ndiag;
			a4pl  = temp4 % p.ndiag;
			tempcd1 = stm[a1bloc][a2bloc].getData(a1pl, a2pl);
			tempcd2 = stm[a3bloc][a4bloc].getData(a3pl, a4pl);
			
			stm[a1bloc][a2bloc].set_data(a1pl, a2pl, tempcd1 - p.thop*p.thop*gl[i][j]);
			stm[a3bloc][a4bloc].set_data(a3pl, a4pl, tempcd2 - p.thop*p.thop*gr[i][j]);
			
		} //endfor
	} //endfor
} //end AttachLeads


void DetachLeads(Parameters p, block** stm, int isl180[], int isr180[], int isl[], int isr[], 
				 complex<double> gl[4][4], complex<double> gr[4][4])
{
	int temp1, temp2, temp3, temp4,
	a1bloc, a2bloc, a3bloc, a4bloc,
	a1pl, a2pl, a3pl, a4pl;
	
	complex<double> tempcd1, tempcd2;
	
	int ig = p.leadSize;
	
	for(int i=0; i<ig; i++)
	{
		for(int j=0; j<ig; j++)
		{
			temp1 = isl180[i]-1;
			temp2 = isl180[j]-1;
			temp3 = isr180[i]-1;
			temp4 = isr180[j]-1;
			
			if(temp1 >= temp2)
			{
				a1bloc = (int) floor(temp1 / p.ndiag);
				a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
			} //endif
			else
			{
				a1bloc = (int) floor(temp2 / p.ndiag);
				a2bloc = get_atom_connect_block(p.idim, p.ndiag, temp2, temp1);
			} //endelse
			
			if(temp3 >= temp4)
			{
				a3bloc = (int) floor(temp3 / p.ndiag);
				a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
			} //endif
			else
			{
				a3bloc = (int) floor(temp4 / p.ndiag);
				a4bloc = get_atom_connect_block(p.idim, p.ndiag, temp4, temp3);
			} //endelse
			
			if(a2bloc == 3)
			{
				//If in U
				//Set a2bloc and a1bloc to U's position in blockHamiltonian
				a2bloc = 0;
				a1bloc = 0;
			} //endif
			
			if(a2bloc == 4)
			{
				a2bloc = 2;
				a1bloc = p.idim-1;
			} //endif
			
			if(a4bloc == 3)
			{
				a4bloc = 0;
				a3bloc = 0;
			} //endif
			
			if(a4bloc == 4)
			{
				a4bloc = 2;
				a3bloc = p.idim-1;
			} //endif
			
			a1pl  = temp1 % p.ndiag;
			a2pl  = temp2 % p.ndiag;
			a3pl  = temp3 % p.ndiag;
			a4pl  = temp4 % p.ndiag;
			tempcd1 = stm[a1bloc][a2bloc].getData(a1pl, a2pl);
			tempcd2 = stm[a3bloc][a4bloc].getData(a3pl, a4pl);
			
			stm[a1bloc][a2bloc].set_data(a1pl, a2pl, tempcd1 + p.thop*p.thop*gl[i][j]);
			stm[a3bloc][a4bloc].set_data(a3pl, a4pl, tempcd2 + p.thop*p.thop*gr[i][j]);
		} //endfor
	} //endfor
} //end DetachLeads


void wiregf(Parameters p, int isl180[], int isr180[], int isl[], int isr[], complex<double> gl[4][4], complex<double> gr[4][4], double** disCoordinates)
{
	int ig = p.leadSize;
	
	complex<double> zfunc, result, unitz;
	
	unitz = complex<double>(0,1);
	
	double yl180[ig], yr180[ig], zl180[ig], zr180[ig];
	double alpha, alpha2, fmnl, fmnr, a1x, a2x;
	
	double emass = .511e06, hbar = 1973.3, emn = 0.0, efermi = 6.0, yo = 100.0, zo = 4.0;
	
	double t = hbar*hbar/emass/p.a/p.a, pi = acos(-1.0);
	
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
} //end wiregf

#endif /* GENERAL_H_ */
