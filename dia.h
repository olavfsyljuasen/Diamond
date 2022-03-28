#include<string>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

struct Record
{
 public:
  double erg;
  Coord q;
};

bool CompareRecords(Record a, Record b)
{
  return a.erg < b.erg;
}


// the unit cell vectors
const Coord a1(0.0,0.5,0.5);
const Coord a2(0.5,0.0,0.5);
const Coord a3(0.5,0.5,0.0);

// the reciprocal lattice vectors
const Coord b1(-TWOPI, TWOPI, TWOPI);
const Coord b2( TWOPI,-TWOPI, TWOPI);
const Coord b3( TWOPI, TWOPI,-TWOPI);

const Coord X1( TWOPI, 0    , 0    );
const Coord X2( 0    , TWOPI, 0    );
const Coord X3( 0    , 0    , TWOPI);

const double Xsqr=scalarproduct(X1,X1);
const double Xsqri=Xsqr*INNERRADIUS*INNERRADIUS;

const Coord L1( PI, PI, PI);
const Coord L2(-PI, PI, PI);
const Coord L3(-PI,-PI, PI);
const Coord L4( PI,-PI, PI);

const double Lsqr=scalarproduct(L1,L1);
const double Lsqri=Lsqr*INNERRADIUS*INNERRADIUS;



bool Useq(Coord q)
{
  double qX1=scalarproduct(q,X1);
  double qX2=scalarproduct(q,X2);
  double qX3=scalarproduct(q,X3);

  double qL1=scalarproduct(q,L1);
  double qL2=scalarproduct(q,L2);
  double qL3=scalarproduct(q,L3);
  double qL4=scalarproduct(q,L4);

  bool InFirstBZ= !( abs(qX1)>Xsqr || abs(qX2)>Xsqr || abs(qX3)>Xsqr ||
		     abs(qL1)>Lsqr || abs(qL2)>Lsqr || abs(qL3)>Lsqr || abs(qL4)>Lsqr);

  bool OutsideInterior = (abs(qX1)>=Xsqri || abs(qX2)>=Xsqri || abs(qX3)>=Xsqri ||
			  abs(qL1)>=Lsqri || abs(qL2)>=Lsqri || abs(qL3)>=Lsqri || abs(qL4)>=Lsqri);

  return InFirstBZ && OutsideInterior;
}



void SetupQpts(const int n1,const int n2,const int n3,vector<Coord>& qpts)
{  
  for(int i3=-n3; i3<=n3; i3++)
    {
      if(i3%100==0) logfile << "progress: i3=" << i3 << endl;
 
      for(int i2=-n2; i2<=n2; i2++)
	for(int i1=-n1; i1<=n1; i1++)
	  {
	    Coord q= (i1*(1./n1))*b1+(i2*(1./n2))*b2+(i3*(1./n3))*b3;
	    
	    if( Useq(q) ){
	      qpts.push_back(q);}
	  }
    }
}


const double epsilon=1e-9;


inline bool deq(double a,double b){return abs(a-b)<=epsilon;}
//inline bool deq(Coord a,Coord b){return length(a-b)<=epsilon;}
//inline bool deq(Coord a,double b){return abs(sqrt(l2(a))-b)<=epsilon;}


double HexDiff(Coord q)
{
  double dqL1=(Lsqr-abs(scalarproduct(q,L1)))/Lsqr;
  double dqL2=(Lsqr-abs(scalarproduct(q,L2)))/Lsqr;
  double dqL3=(Lsqr-abs(scalarproduct(q,L3)))/Lsqr;
  double dqL4=(Lsqr-abs(scalarproduct(q,L4)))/Lsqr;

  return min(min(dqL1,dqL2),min(dqL3,dqL4));  
}

double SqrDiff(Coord q)
{
  double dqX1=(Xsqr-abs(scalarproduct(q,X1)))/Xsqr;
  double dqX2=(Xsqr-abs(scalarproduct(q,X2)))/Xsqr;
  double dqX3=(Xsqr-abs(scalarproduct(q,X3)))/Xsqr;

  return min(min(dqX1,dqX2),dqX3);  
}




enum phases{                      NOTFOUND  ,GAMMA,  X , K , U , W,   L,   XU,  KU,  XK, GAMMAK,GAMMAX,GAMMAL, KC  };
vector<string> phasedescription={"Not found", "G" , "X","K","U","W", "L", "XU","KU","XK","GK"  , "GX" , "GL" ,"KC" };

int IdentifyPhase(Coord q,double& alpha)
{
  const double qx=q.x();
  const double qy=q.y();
  const double qz=q.z();
  
  const double aqx=abs(qx);
  const double aqy=abs(qy);
  const double aqz=abs(qz);
  
  alpha=0.;
  if( deq(aqx,0.) && deq(aqy,0.) && deq(aqz,0.) ){return GAMMA;}
  
  if( deq(aqx,PI) && deq(aqy,PI) && deq(aqz,PI) ){return L;}
  
  if( (deq(aqx,TWOPI) && deq(aqy,0.) && deq(aqz,0.)) ||
      (deq(aqy,TWOPI) && deq(aqz,0.) && deq(aqx,0.)) ||
      (deq(aqz,TWOPI) && deq(aqx,0.) && deq(aqy,0.))){ return X;}
  
  if( (deq(aqx,0.75*TWOPI) && deq(aqy,0.75*TWOPI) && deq(aqz,0.)) ||
      (deq(aqy,0.75*TWOPI) && deq(aqz,0.75*TWOPI) && deq(aqx,0.)) ||
      (deq(aqz,0.75*TWOPI) && deq(aqx,0.75*TWOPI) && deq(aqy,0.))){ return K;}
  
  if( (deq(aqx,TWOPI) && deq(aqy,0.25*TWOPI) && deq(aqz,0.25*TWOPI)) ||
      (deq(aqy,TWOPI) && deq(aqz,0.25*TWOPI) && deq(aqx,0.25*TWOPI)) ||
      (deq(aqz,TWOPI) && deq(aqx,0.25*TWOPI) && deq(aqy,0.25*TWOPI))){return U;}
  
  if( deq(aqx,TWOPI) && deq(aqy,aqz) ){ alpha=aqy/(0.5*PI); return XU;}
  if( deq(aqy,TWOPI) && deq(aqz,aqx) ){ alpha=aqz/(0.5*PI); return XU;}
  if( deq(aqz,TWOPI) && deq(aqx,aqy) ){ alpha=aqx/(0.5*PI); return XU;}
  
  if( (deq(aqx,0.) && ( (deq(aqy,TWOPI) && deq(aqz,PI)) || (deq(aqz,TWOPI) && deq(aqy,PI)) )) || 
      (deq(aqy,0.) && ( (deq(aqz,TWOPI) && deq(aqx,PI)) || (deq(aqx,TWOPI) && deq(aqz,PI)) )) || 
      (deq(aqz,0.) && ( (deq(aqx,TWOPI) && deq(aqy,PI)) || (deq(aqy,TWOPI) && deq(aqx,PI)) )) ){return W;} 
  
  if( deq(aqx,aqy) && !deq(aqz,aqx) )
    { 
      alpha=aqz/TWOPI; double beta=0.75*TWOPI-PI*alpha; 
      if(deq(aqx,beta)){ return KU;}
    }
    
  if( deq(aqy,aqz) && !deq(aqx,aqy) )
    { 
      alpha=aqx/TWOPI; double beta=0.75*TWOPI-PI*alpha; 
      if(deq(aqy,beta)){ return KU;}
    }
  
  if( deq(aqz,aqx) && !deq(aqy,aqz) )
      { 
	alpha=aqy/TWOPI; double beta=0.75*TWOPI-PI*alpha; 
	if(deq(aqz,beta)){ return KU;}
      }
     
  if( deq(aqx,aqy) && !deq(aqz,aqx) )
    { 
      alpha=aqx/(0.75*TWOPI); double beta=TWOPI-TWOPI*alpha; 
      if(deq(aqz,beta)){ return XK;}
    }

  if( deq(aqy,aqz) && !deq(aqx,aqy) )
    { 
      alpha=aqy/(0.75*TWOPI); double beta=TWOPI-TWOPI*alpha; 
      if(deq(aqx,beta)){ return XK;}
    }

  if( deq(aqz,aqx) && !deq(aqy,aqz) )
    { 
      alpha=aqz/(0.75*TWOPI); double beta=TWOPI-TWOPI*alpha; 
      if(deq(aqy,beta)){ return XK;}
    }

 
  if( deq(aqx,aqy) && deq(aqz,0.) ){ alpha=aqx/(0.75*TWOPI) ; return GAMMAK;} 
  if( deq(aqy,aqz) && deq(aqx,0.) ){ alpha=aqy/(0.75*TWOPI) ; return GAMMAK;} 
  if( deq(aqz,aqx) && deq(aqy,0.) ){ alpha=aqz/(0.75*TWOPI) ; return GAMMAK;} 
  
  if( deq(aqx,aqy) && deq(aqy,aqz) && deq(aqz,aqx) ){alpha=aqx/PI; return GAMMAL;}
  
  if( !deq(aqx,0.) && deq(aqy,0.) && deq(aqz,0.) ){alpha=aqx/TWOPI; return GAMMAX;}
  if( !deq(aqy,0.) && deq(aqz,0.) && deq(aqx,0.) ){alpha=aqy/TWOPI; return GAMMAX;}
  if( !deq(aqz,0.) && deq(aqx,0.) && deq(aqy,0.) ){alpha=aqz/TWOPI; return GAMMAX;}
  
  if( deq(aqx,aqy) ){alpha=aqz/TWOPI; return KC;}
  if( deq(aqy,aqz) ){alpha=aqx/TWOPI; return KC;}
  if( deq(aqz,aqx) ){alpha=aqy/TWOPI; return KC;}
  
  
  return NOTFOUND;
}


class Model
{
 public:
 Model(double J1Ain,double J1Bin,double J2Ain,double J2Bin,double J3in,double J4in):J1A(J1Ain),J1B(J1Bin),J2A(J2Ain),J2B(J2Bin),J3(J3in),J4(J4in){}
  
  double FindMinEnergyValue(vector<Coord>& qpts);
  int FindMinimalSurface(double maxval,vector<Coord>& qpts,vector<Record>& myr);
  double lambdamin(Coord);
  
 private:
  double J1A;
  double J1B;
  double J2A;
  double J2B;
  double J3;
  double J4;
};


double Model::lambdamin(Coord q)
{
  /*
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

  double G4 = 2*B-3.; // the fourth neigbor coupling cos[qx]+cos[qy]+cos[qz];
  
  double S11=1+G;
  double S13=2*( G*(1+G)-K+2*B-3);
  double S33=G*( 2*G+4*B-7) +6*K-8*B+9;

  // must put in K1 etc here.

  return 2*J2*G + J4*G4 - sqrt(J1*J1*S11+J1*J3*S13+J3*J3*S33);
  */

  double cx=cos(0.5*q.x());
  double cy=cos(0.5*q.y());
  double cz=cos(0.5*q.z());

  double Cxy=cx*cy;
  double Cyz=cy*cz;
  double Czx=cz*cx;

  if(J1A != J1B)
    {
      double cypz=cos(0.25*(q.y()+q.z()));
      double cymz=cos(0.25*(q.y()-q.z()));
      
      double J1Acypz=J1A*cypz;
      double J1Bcymz=J1B*cymz;
      
      return 2*J2A*Czx + 2*J2B*(Cxy+Cyz) - sqrt(J1Acypz*J1Acypz+J1Bcymz*J1Bcymz+J1A*J1B*(Cxy+Czx));
    }
  else
    {
      return 2*J2A*Czx + 2*J2B*(Cxy+Cyz) - abs(J1A)*sqrt(1+Cxy+Cyz+Czx);
    }

}


double Model::FindMinEnergyValue(vector<Coord>& qpts)
{
  int N=qpts.size();
  double minerg=1E+10;
  for(int i=0; i<N; i++)
    {
      Coord q=qpts[i];
      double erg=lambdamin(q);
      if(erg <= minerg){minerg=erg;}
    }
  return minerg;
}

int Model::FindMinimalSurface(double maxval,vector<Coord>& qpts,vector<Record>& myr)
{
  int N=qpts.size();
  int counter=0;
  
  for(int i=0; i<N; i++)
    {
      Coord q=qpts[i];
      double erg=lambdamin(q);
      if(erg <= maxval)
	{
	  Record r;
	  r.erg=erg;
	  r.q = q;
	  myr.push_back(r);
	  counter++;
	}
    }

  sort(myr.begin(),myr.end(),CompareRecords);
  return counter;
}
  
