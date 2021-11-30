#include<string>
#include<vector>
#include<cmath>

// the unit cell vectors
const Coord a1(0.0,0.5,0.5);
const Coord a2(0.5,0.0,0.5);
const Coord a3(0.5,0.5,0.0);

// the reciprocal lattice vectors
const Coord b1(-TWOPI, TWOPI, TWOPI);
const Coord b2( TWOPI,-TWOPI, TWOPI);
const Coord b3( TWOPI, TWOPI,-TWOPI);

const double bsqr=scalarproduct(b1,b1);


bool InFirstBz(Coord q)
{
  double qb1=scalarproduct(q,b1);
  double qb2=scalarproduct(q,b2);
  double qb3=scalarproduct(q,b3);
  
  return !(qb1>bsqr || qb1<-bsqr || qb2<-bsqr || qb2>bsqr || qb3<-bsqr || qb3>bsqr);
}



void SetupQpts(double* par,vector<Coord> qpts)
{

  const int N1=par[N1];
  const int N2=par[N2];
  const int N3=par[N3];

  
  
  for(int i3=-N3; i3<=N3; i3++)
    for(int i2=-N2; i2<=N2; i2++)
	for(int i1=-N1; i1<=N1; i1++)
	  {
	    Coord q=b1*(i1/(1.*N1))*b1+(i2/(1.*N2))*b2+(i3/(1.*N3))*b3;
	    if( InFirstBz(q) ){ qpts.push_back(q);}
	  }
}


class Model
{
 public:
 Model(double J1in,double J2in,double J3in):J1(J1in),J2(J2in),J3(J3in){}
  
  double FindMinEnergyValue(vector<Coord> qpts);
  int FindMinimalSurface(double maxval,vector<Coord> qpts,vector<Coord> minqs, vector<double> minergs);
  double lambdamin(Coord);
  
 private:
  double J1;
  double J2;
  double J3;
};


double Model::lambdamin(Coord q)
{
  double cx=cos(0.5*q.x());
  double cy=cos(0.5*q.y());
  double cz=cos(0.5*q.z());

  double cx2=cx*cx;
  double cy2=cy*cy;
  double cz2=cz*cz;
  
  
  double A=cx +cy +cz;
  double B=cx2+cy2+cz2;
  double F=cx2*cx2+cy2*cy2+cz2*cz2;

  double G=0.5*(A*A-B);
  double K=0.5*(B*B-F);

  double S11=1+G;
  double S13=2*( G*(1+G)-K+2*B-3);
  double S33=G*( 2*G+4*B-7) +6*K-8*B+9;

  return 2*J2*G- sqrt(J1*J1*S11+J1*J3*S13+J3*J3*S33);
}


double Model::FindMinEnergyValue(vector<Coord> qpts)
{
  int N=qpts.size();
  double minerg=9999999;
  for(int i=0; i<N; i++)
    {
      Coord q=qpts[i];
      double erg=lambdamin(q);
      if(erg <= minerg){minerg=erg;}
    }
  return minerg;
}

int Model::FindMinimalSurface(double maxval,vector<Coord> qpts,vector<Coord> minqs, vector<double> minergs)
{
  int N=qpts.size();
  int counter=0;
  
  for(int i=0; i<N; i++)
    {
      Coord q=qpts[i];
      double erg=lambdamin(q);
      if(erg <= maxval)
	{
	  minqs.push_back(q);
	  minergs.push_back(erg);
	  counter++;
	}
    }
  return counter;
}
  
