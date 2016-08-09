#include "block.h"
#include <iostream>
using namespace std;

class block
{
public:
	block(int);
	~block();
	int length;
	//float data[];
	float *data;
}; //end block

block::block(int a)
{
	length = a;
	data = new float[length];
	for(int i=0; i<length; i++) data[i] = i;
} //end block::block

block::~block()
{
	delete [] data;
} //end block::~block

int main(int argc,char **args)
{
	int length = 12;
	block hmod(length);
	for(int i=0; i<length; i++) cout << hmod.data[i] << " ";
	cout << endl;
	
	return 0;
} //end main
