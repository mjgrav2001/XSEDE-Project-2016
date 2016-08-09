/*
 * general.h
 *
 *  Created on: Jul 30, 2010
 *      Author: leond
 */


    void set_coord(double coord[3600][3]){
	            double a1x,a1y,a1z,a2x,a2y,a2z,ap;
				FILE* file = fopen("/Users/leond/Documents/Code/Research/positions", "r"); // opens the file
				ifstream inFile;
				size_t BUFSIZE = 1000;
				char buf[BUFSIZE];
				int linenum=0;
				//for(linenum=0;linenum<3600*4;linenum++){   // fscanf(file, "%s", buf); //   string<<C
				while(EOF!=fscanf(file, "%lf %lf %lf %lf ",&a1x,&a1y,&a1z,&ap))
				{   linenum++;
				    if (linenum < 3601){
                    	coord[linenum-1][0]=a1x;
                    	coord[linenum-1][1]=a1y;
                    	coord[linenum-1][2]=a1z;
  //                  	cout <<  coord[linenum][0] <<"  "<<coord[linenum][1] <<"  "<<coord[linenum][2] <<"  "<< ap <<"  "<< linenum <<endl;
                    }


				}
//				}
			//	cout <<  "Read in done!  " <<endl;
   }
