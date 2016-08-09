/*
 * leads.h
 *
 *  Created on: Jul 28, 2010
 *      Author: leond
 */

void attach_leads(block stm[300][3], int isl180[],int isr180[],int isl[], int isr[], complex <double> gl[4][4],complex <double> gr[4][4]){

	int i,j,temp1,temp2,temp3,temp4,a1bloc,a2bloc,a3bloc,a4bloc,a1pl,a2pl,a3pl,a4pl;
	complex <double> tempcd1, tempcd2;
	double thop;
	thop = -0.25e0;

//	for(i=0;i<4;i++){
//	for(j=0;j<4;j++){
////	cout <<"i j-"<<i << " "<<j<<"  gl[i]-" << gl[i][j] <<"  gr[i]-"<< gr[i][j]<< endl;
//		gl[i][j]= std::complex <double>(0,0);
//	}
//	 gl[1][1] = std::complex <double>(-4.678802688135959E-002,-1.203406422832061E-002);
//   gr[1][1] = std::complex <double>(-4.678802688135929E-002,-1.203406422831929E-002);
//	 gl[1][3] = std::complex <double>(4.678802688135961E-002, 1.203406422832102E-002);
//	 gr[1][3] = std::complex <double>(4.678802688135954E-002, 1.203406422832026E-002);
//	 gl[3][1] = std::complex <double>(4.678802688135961E-002, 1.203406422832102E-002);
//	 gr[3][1] = std::complex <double>( 4.678802688135954E-002, 1.203406422832026E-002);
//	 gl[3][3] = std::complex <double>(-4.678802688135972E-002,-1.203406422832143E-002);
//	 gr[3][3] = std::complex <double>(-4.678802688135972E-002,-1.203406422832122E-002);
//
//	}


	for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			temp1=isl180[i]-1;
			temp2=isl180[j]-1;
			temp3=isr180[i]-1;
			temp4=isr180[j]-1;

			if (temp1>=temp2){
			a1bloc=floor(temp1 / 12);
			a2bloc=get_atom_connect_block(temp2,temp1);
			}
			else{
			a1bloc=floor(temp2 / 12);
			a2bloc=get_atom_connect_block(temp2,temp1);
			}

			if (temp3>=temp4){
			a3bloc=floor(temp3 / 12);
			a4bloc=get_atom_connect_block(temp4,temp3);
			}
			else{
			a3bloc=floor(temp4 / 12);
			a4bloc=get_atom_connect_block(temp4,temp3);
			}

			if (a2bloc == 3) {  /// If in U
				a2bloc = 0;     /// Set a2bloc and a1bloc to U's position in hmod
				a1bloc = 0;
			}

			if (a2bloc == 4) {
					a2bloc = 2;
					a1bloc = 299;
				}

			if (a4bloc == 3) {
				a4bloc = 0;
				a3bloc = 0;
			}

			if (a4bloc == 4) {
					a4bloc = 2;
					a3bloc = 299;
				}

//			cout <<  "Blocks-"<<a1bloc<<" "<<a2bloc<<" "<<a3bloc<<" "<<a4bloc<<" "<<endl;

			a1pl  =temp1%12;
			a2pl  =temp2%12;
			a3pl  =temp3%12;
			a4pl  =temp4%12;
			tempcd1 = stm[a1bloc][a2bloc].getData(a1pl,a2pl);
			tempcd2 = stm[a3bloc][a4bloc].getData(a3pl,a4pl);

//			  cout << "B-Left i-"<<temp1<<"   j-"<<temp2 <<endl;
//			  cout << "B-Right i-"<<temp3<<"   j-"<<temp4<<endl;
//			  cout << " "<<endl;
//	          cout << " html(i,j)"<<tempcd1<< "|||"<< tempcd2 <<endl;
//	          cout << " tp*tp*gl(i,j)  "<<thop*thop*gl[i][j]<< "|||"<< thop*thop*gr[i][j]<<endl;
//            cout << tempcd1-thop*thop*gl[i][j]<< "|||"<< tempcd2-thop*thop*gr[i][j]<<endl;
//			  cout << " "<<endl;

			stm[a1bloc][a2bloc].set_data(a1pl,a2pl,(tempcd1)-(thop*thop*(gl[i][j])));
			stm[a3bloc][a4bloc].set_data(a3pl,a4pl,(tempcd2)-(thop*thop*(gr[i][j])));

//			  cout << stm[a1bloc][a2bloc].getData(a1pl,a2pl)<< "   "<< stm[a3bloc][a4bloc].getData(a3pl,a4pl);
//			  cout <<" "<<endl;
//
		}
	}

}
