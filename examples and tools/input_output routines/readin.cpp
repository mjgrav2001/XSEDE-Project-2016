#include <iostream>
#include <fstream>
using namespace std;

const int MAX_LENGTH = 100;

int main()
{
	ifstream inFile("input.txt");
	ofstream outFile("output.txt");
	int i = 0;
	int variable1[MAX_LENGTH], variable2[MAX_LENGTH], variable3[MAX_LENGTH];
	
	//Read in from file
	if(inFile.is_open())
	{
		while(!inFile.eof())
		{
			inFile >> variable1[i] >> variable2[i] >> variable3[i];
			i++;
		} //endfor
		inFile.close();
	} //endif
	else cout << "Error: inFile not opened!\n";
	
	//Output to screen
	for(int j=0; j<i; j++)
	{
		cout << variable1[j] << " " << variable2[j] << " " << variable3[j] << " ";
	} //endfor
	cout << endl;
	
	//Write to new file
	if(outFile.is_open())
	{
		outFile.precision(5);
		for(int j=0; j<i; j++)
		{
			outFile << fixed << static_cast<float>(variable1[j]) << endl << static_cast<float>(variable2[j])
					<< endl	<< static_cast<float>(variable3[j]) << endl;
		}
		outFile.close();
	}
	else cout << "Error: outFile not opened!\n";
	
	return 0;
} //end main