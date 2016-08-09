#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

const int MAX_LENGTH = 10000;

int main()
{
	ifstream inFile("ARRAY_disr_all.dat");
	//ofstream outFile("ARRAY_posr_new.dat");
	   
	FILE * outfile = fopen("ARRAY_posr_new.dat", "w" );
	
	int i = 0;
	float variable1[MAX_LENGTH], variable2[MAX_LENGTH], variable3[MAX_LENGTH], variable4[MAX_LENGTH], variable5[MAX_LENGTH];
	
	//Read in from file
	if(inFile.is_open())
	{
		while(!inFile.eof())
		{
			inFile >> variable1[i] >> variable2[i] >> variable3[i] >> variable4[i] >> variable5[i];
			i++;
		} //endfor
		inFile.close();
	} //endif
	else cout << "Error: inFile not opened!\n";
	
	//Output to screen
	for(int j=0; j<i; j++)
	{
		cout << variable3[j] << " " << variable4[j] << " " << variable5[j] << endl;
	} //endfor
	cout << endl;
	
	
	//Write to new file
	//if(outFile.is_open())
	//{
	//	outFile.precision(6);
	//	for(int j=0; j<i-1; j++)
	//	{
	//		outFile << fixed << static_cast<float>(variable3[j]) << " " << static_cast<float>(variable4[j])  " " << static_cast<float>(variable5[j]) << endl;		
	//		
	//	}
	//	outFile.close();
	//}
	//else cout << "Error: outFile not opened!\n";

	
	for(int j=0; j<i; j++)
	{
		fprintf(outfile, "%15.6f %14.6f %14.6f %8i \n", variable4[j], variable3[j], variable5[j], j+1);	
	}
	fclose(outfile); 
	
	return 0;
} //end main