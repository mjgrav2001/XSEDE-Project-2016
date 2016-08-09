#include <iostream>
#include <cstdlib>

#include <cmath>
//#define M_E         2.71828182845904523536028747135266250   /* e */
//#define M_LOG2E     1.44269504088896340735992468100189214   /* log 2e */
//#define M_LOG10E    0.434294481903251827651128918916605082  /* log 10e */
//#define M_LN2       0.693147180559945309417232121458176568  /* log e2 */
//#define M_LN10      2.30258509299404568401799145468436421   /* log e10 */
//#define M_PI        3.14159265358979323846264338327950288   /* pi */
//#define M_PI_2      1.57079632679489661923132169163975144   /* pi/2 */
//#define M_PI_4      0.785398163397448309615660845819875721  /* pi/4 */
//#define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi */
//#define M_2_PI      0.636619772367581343075535053490057448  /* 2/pi */
//#define M_2_SQRTPI  1.12837916709551257389615890312154517   /* 2/sqrt(pi) */
//#define M_SQRT2     1.41421356237309504880168872420969808   /* sqrt(2) */
//#define M_SQRT1_2   0.707106781186547524400844362104849039  /* 1/sqrt(2) */
#define M_1_2PI	(M_1_PI/2)	/* 1/2pi */
#define M_2PI	(2*M_PI)	/* 2pi   */

using namespace std;

int wrap;
int stack;
int chn, chm, chj, chi;
float chx, chy;
int tn, _tm, tj, ti;
int gcdtij;
float tx, ty;
int Tn, Tm, Tj, Ti;
float Tx, Ty;

int MINj;
int MAXj;
int MINi;
int MAXi;

//int FIRSTNEIGH[4][2] = {{-1,1},{-1,0},{1,0},{1,-1}};
int FIRSTNEIGHj[4] = { -1, -1, 1,  1};
int FIRSTNEIGHi[4] = {  1,  0, 0, -1};

int N=0;
class Atom
{
public:
	int n,m,j,i,k;
	float x,y,z;
	int iNeigh[3];
	Atom* pNeigh[3];
};
Atom** atoms;

Atom* getAtomAt(int j, int i)
{
	for(int k = 0; k<N; k++)
		if(atoms[k]->i==i && atoms[k]->j==j) return atoms[k];
	return NULL;
}

float ji2x(int j, int i)
{
	return (float)(j+i);
}

float ji2y(int j, int i)
{
	return sqrt(3.0) * (i  - (j&1)/3.0);
}

void ji2Cyl(Atom* a, float rotation)
{
	int i = a->i;
	int j = a->j;
	float x = ji2x(j,i);
	float y = ji2y(j,i);
	float chsq = chx*chx + chy*chy;
	float r = sqrt(chsq) * M_1_2PI;
	float b = (x*chx+y*chy) / chsq * M_2PI + rotation;
	a->x = r*cos(b);
	a->y = r*sin(b);
	a->z = (x*tx+y*ty)/sqrt(tx*tx+ty*ty);
}

void ji2Thor(Atom* a, float rotation)
{
	ji2Cyl(a, rotation);
	float R = sqrt(Tx*Tx + Ty*Ty) * M_1_2PI;
	float b = a->z/R;
	R+=a->x;
	a->z = a->y;
	a->x = R*cos(b);
	a->y = R*sin(b);
}

int gcd(int p, int q)
{
	if(p<0) p=-p;
	if(q<0) q=-q;
	while(p>0)
	{
		int r = q % p;
		q=p;
		p=r;
	}
	return q;
}

bool insider(int j, int i)
{
	int p = (chj+chi)*j + (chj+4*chi)*i - chi*(j&1);
	if(p<0) return false;
	p -= chj*chj + 2*chj*chi + 4*chi*chi;
	if(p>=0) return false;
	p = (Tj+Ti)*j + (Tj+4*Ti)*i - Ti*(j&1);
	if(p<0) return false;
	p -= Tj*Tj + 2*Tj*Ti + 4*Ti*Ti;
	if(p>=0) return false;
	return true;
}

//void readFile(char* path)
//{
//	FILE* file = fopen(path, 'r');
//	
//	fscanf(file, "%d %d %d %d %f %f", &chn, &chm, &chj, &chi, &chx, &chy);
//	fscanf(file, "%d", &stack);
//	fscanf(file, "%d %d %d %d %f %f", &tn, &tm, &tj, &ti, &tx, &ty);
//	fscanf(file, "%d %d %d %d %f %f", &tn, &tm, &tj, &ti, &tx, &ty);
//	fscanf(file, "%d", &stack);
//	fscanf(file, "%d", &stack);
//	
//	
//	fscanf(file, "%d %d %d %d %d %d %d %d %f %f %f", &tn, &tm, &tj, &ti, &tx, &ty);
//}

void printMyFormat(){
	cout << "Chiral Vector (n,m,j,i,x,y)\n";
	cout << chn << " " << chm << " " << chj << " " << chi << " " << chx << " " << chy << endl;
	cout << "GCD\n" << gcdtij << endl;
	cout << "Translational Vector (j,i,x,y)\n";
	cout << tj << " " << tn << " " << _tm << " "  << ti << " " << tx << " " << ty << endl;
	cout << "Unit Cell Count\n";
	cout << stack << endl;
	cout << "Big Translational Vector (j,i,x,y)\n";
	cout << Tj << " " << Tn << " " << Tm << " "<< Ti << " " << Tx << " " << Ty << endl;
	cout << "MIN j,i, MAX j,i, Interval j,i\n";
	cout << MINj << " " << MAXj << " " << MINi << " " << MAXi << " " << (MAXj-MINj+1) << " " << (MAXi-MINi+1) << endl;
	cout << "Nodes Count\n";
	cout << N << endl;
	cout << "Big Radius, Small Radius\n";
	cout << 0 << " " << 0 << endl;
	cout << "NODES [node num, (n,m) , (j,i) , neighbours nums 1,2,3 , (x,y,z) , (phi,theta)]:" << endl;
	
	for(int k=0; k<N; k++)
	{
		Atom* a = atoms[k];
		printf("%d %d %d %d %d", a->k, a->n, a->m, a->j, a->i);
		printf(" %d %d %d", a->iNeigh[0], a->iNeigh[1], a->iNeigh[2]);
		printf(" %f %f %f", a->x, a->y, a->z); 
//		printf(" %f %f", 0.0, 0.0); 
		cout << endl;
	}
}

void printWGraph3D(int includeBasis)
{
	cout << "GraphPlot3D[{" ;
	
	char separ=' ';
	for(int k=0; k<N; k++)
	{
		Atom* a = atoms[k];
		//cout << separ << a->k << "->" << a->k; separ = ',';
		for(int l=0; l<3; l++)
		{
			Atom* b = a->pNeigh[l];
			if(b==NULL) continue;
			if(a->k > b->k) continue;
			cout << separ<< a->k << "->" << b->k;
			separ = ',';
		}
	}
	if(includeBasis){
		int origin = (wrap?N+3:1);
		printf("%cTooltip[%d,\"Chiral Vector\"]->", separ, N+1);
		if(wrap)
			printf("Tooltip[%d,\"Origin\"]", origin);
		else
			printf("%d");
		
		printf("%cTooltip[%d,\"Big Translation Vector\"]->%d", separ, N+2, origin);
	}
	cout << "}";
	
	cout << ", VertexCoordinateRules -> {";
	separ=' ';
	for(int k=0; k<N; k++)
	{
		Atom* a = atoms[k];
		printf("%c%d->{%5.3f,%5.3f,%5.3f}", separ,a->k, a->x, a->y, a->z);
		separ = ',';
	}

	if(includeBasis)
	{
		printf("%c%d->{%5.3f,%5.3f,%5.3f}", separ, N+1, chx, chy, 0.0);
		printf("%c%d->{%5.3f,%5.3f,%5.3f}", separ, N+2, Tx , Ty , 0.0);
		if(wrap)
			printf("%c%d->{%5.3f,%5.3f,%5.3f}", separ, N+3, 0.0, 0.0, 0.0);
	}
	
	cout << "}";
	
	cout << ", SelfLoopStyle->0.0001,\nVertexLabeling -> Tooltip,\n EdgeLabeling -> None" << endl;
//	cout << ",EdgeRenderingFunction -> (Cylinder[#1, .05] &),\n \
//	VertexRenderingFunction -> ({ColorData[\"Atoms\"][12], Sphere[#1, .2]} &),\n \
//	PlotStyle -> Directive[Specularity[White, 20]]\n "
	cout << "]\n\n";
}

int main(int argc, char* argv[])
{
	chn = atoi(argv[1]);
	chm = atoi(argv[2]);
	stack = atoi(argv[3]);
	wrap = atoi(argv[4]);
	
	chi = 2*chm;
	chj = 2*chm-4*chn;
	chx = ji2x(chj, chi);
	chy = ji2y(chj, chi);
	
	tj = -chj - 4*chi;
	ti = chj + chi;
	gcdtij = gcd(tj,ti);
	if((tj/gcdtij)%2) gcdtij /=2;
	tj /=gcdtij;
	ti /=gcdtij;
	
	tx = ji2x(tj,ti);
	ty = ji2y(tj,ti);
	
	Tj = stack * tj;
	Ti = stack * ti;
	Tx = stack * tx;
	Ty = stack * ty;
	
	MINj = Ti+Tj-2;
	MAXj = chj + chi + 2;
	MINi = -1;
	MAXi = chi + Ti + 1;
	
	//	cout << FIRSTNEIGH[0][0] << ' ' << FIRSTNEIGH[0][1] << endl
	//	     << FIRSTNEIGH[1][0] << ' ' << FIRSTNEIGH[1][1] << endl
	//	     << FIRSTNEIGH[2][0] << ' ' << FIRSTNEIGH[2][1] << endl
	//	     << FIRSTNEIGH[3][0] << ' ' << FIRSTNEIGH[3][1] << endl
	//	;
	
	
	N=0;
	for(int i=MINi; i<=MAXi; i++)
	{
		int minj = MINj - i;
		int maxj = MAXj - i;
		for(int j = minj; j<=maxj; j++)
			if(insider(j,i)) N++;
	}
	//cout << N << endl;
	
	atoms = new Atom*[N];
	
	N=0;
	for(int i=MINi; i<=MAXi; i++)
	{
		int minj = MINj - i;
		int maxj = MAXj - i;
		for(int j = minj; j<=maxj; j++)
		{
			if(!insider(j,i)) continue;
			
			Atom* a = atoms[N++] = new Atom();
			a->i=i;
			a->j=j;
			a->k=N;
			switch(wrap)
			{
				case 0:
					a->x=ji2x(j,i);
					a->y=ji2y(j,i);
					a->z=0;
					break;
				case 1:
					ji2Cyl(a, 0);
					break;
				case 2:
					ji2Thor(a, 0);
					break;
			}
			a->n = a->m = 0;
			a->iNeigh[0] = a->iNeigh[1] = a->iNeigh[2] = 0;
			a->pNeigh[0] = a->pNeigh[1] = a->pNeigh[2] = NULL;
		}
	}
	
int N2=0;
	for(int k=0; k<N; k++)
	{
		Atom* a = atoms[k];
		int base = a->j&1;
		for(int l=0; l<3; l++)
		{
			int i2 = a->i + FIRSTNEIGHi[base+l];
			int j2 = a->j + FIRSTNEIGHj[base+l];
			Atom* b = getAtomAt(j2, i2);
			if(wrap==1 || wrap==2)
			{
				if(b == NULL) b=getAtomAt(j2+chj, i2+chi);
				if(b == NULL) b=getAtomAt(j2-chj, i2-chi);
				if(wrap==2)
				{
					if(b == NULL)  b=getAtomAt(j2+Tj, i2+Ti);
					if(b == NULL)  b=getAtomAt(j2-Tj, i2-Ti);
				}
				
			}
			
			if(b!=NULL) a->iNeigh[l] = b->k;
			a->pNeigh[l] = b;
			if(b==NULL)N2++;
		}
	}
	
////////////////////////////////////////////////////////////	

printMyFormat();
//printWGraph3D(0);
	
//    for(int i=0; i<N; i++)
//    {
//		Atom* a = atoms[i];
//	int ka = a->k;
//	int k0 = 0; if(a->pNeigh[0]!=NULL) k0 = a->pNeigh[0]->k;
//	int k1 = 0; if(a->pNeigh[1]!=NULL) k1 = a->pNeigh[1]->k;
//	int k2 = 0; if(a->pNeigh[2]!=NULL) k2 = a->pNeigh[2]->k;
//	printf("%d %d %d %d %d %d %5.3f %5.3f\n",
//	
//	  a->k, a->j, a->i, k0, k1, k2, a->x, a->y, a->z);
//    }
//	

//	cout << endl << "H = {";
//    for(int i=0; i<N; i++)
//    {
//	Atom* a = atoms[i];
//	int ka = a->k;
//	int k0 = -1; if(a->pNeigh[0]!=NULL) k0 = a->pNeigh[0]->k;
//	int k1 = -1; if(a->pNeigh[1]!=NULL) k1 = a->pNeigh[1]->k;
//	int k2 = -1; if(a->pNeigh[2]!=NULL) k2 = a->pNeigh[2]->k;
//
//	if(i>0) cout << ',' << endl;
//	cout << '{';
//	for(int j=0; j<N; j++)
//	{
//	    int k = atoms[j]->k;
//		if(j>0) cout << ',';
//	    cout << ((ka==k)?'0':(k==k0||k==k1||k==k2)?'1':'0');
//	}	
//	cout << '}';
//    }

//	separ=' ';
//	cout << endl << "H = SparseArray[{";
//    for(int k=0; k<N; k++)
//    {
//		Atom* a = atoms[k];
//		int i = a->k;
//		for(int l=0 ; l<3; l++)
//		{
//			if(a->pNeigh[l]==NULL) continue;
//			int j = a->iNeigh[l];
//			cout << separ << "{" << i << "," << j << "}->1";
//			separ=',';
//		}
//		cout << endl;
//	}
//	cout << endl << "}]; MatrixPlot[H]" << endl << "G = Inverse[H]; MatrixPlot[G]" <<endl;
return 0;
}
;	;

//Animate[

//  ViewVector -> {{14*Cos[x], 14*Sin[x], 4*(1 - Cos[x])}, {-22*Sin[x], 
//   22*Cos[x], 0}}]
// , {x, 0, 6.25}, AnimationRunning -> False, AnimationRate -> 0.05]