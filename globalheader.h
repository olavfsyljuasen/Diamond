#ifndef GLOBALHEADER_H
#define GLOBALHEADER_H

#include <fstream>
#include <iomanip>
#include<stdio.h>
#include<stdlib.h>
#include <string>
#include<vector>
#include<time.h>
using namespace std;


//enum latticedirections{RIGHTDIR,LEFTDIR,UPDIR,DOWNDIR,INDIR,OUTDIR};


string bitstring(int s,int cycle)
{
  //  string str=(s==0 ? "0":"");
  string str="";
  //  for(int i=0; (s!=0 || i<cycle); s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  for(int i=0; i<cycle; s>>=1, i++){if(s & 1){str='1'+str;}else{str='0'+str;}}
  return str;
}




class Timer{
  friend ostream&
    operator<<(ostream& os, Timer timer)
    { time_t t;
      if( time(&t) == time_t(-1)){exit(1);}
      os << ctime(&t); return os;
    }  
 public:
  Timer(){};
  ~Timer(){if(TRACE) cout << "deleting timer\n";}
  void Start(){t1=time(0);}
  void Stop(){t2=time(0);}
  double GetTimeElapsed(){return difftime(t2,t1);}
 private:
  time_t t1;
  time_t t2;
};

bool iszero(const double v)
{
  const double tiny=1.e-13;
  return (v > -tiny && v < tiny);
}



class Coord
{  friend ostream& operator<<(ostream& os,Coord& c)
    { os << c.x() << " " << c.y() << " " << c.z(); return os; }
  friend ostream& operator<<(ostream& os,const Coord& c)
  { os << c.x() << " " << c.y() << " " << c.z(); return os; }
 public:
  Coord(double xx=0.,double yy=0.,double zz=0.):cx(xx),cy(yy),cz(zz){}
  void x(double a){cx=a;}
  void y(double a){cy=a;}
  void z(double a){cz=a;}
  double x() const {return cx;}
  double y() const {return cy;}
  double z() const {return cz;}
 private:
  double cx,cy,cz;
};

inline bool operator==(const Coord& a,const Coord& b)
{ return bool((a.x() == b.x()) && (a.y()==b.y()) && (a.z()==b.z()));}

Coord operator-(const Coord& l,const Coord& r){
  return Coord(l.x()-r.x(),l.y()-r.y(),l.z()-r.z());
}

Coord operator+(const Coord& l,const Coord& r){
  return Coord(l.x()+r.x(),l.y()+r.y(),l.z()+r.z());
}

double operator*(const Coord& l,const Coord& r){
  return double(l.x()*r.x()+l.y()*r.y()+l.z()*r.z());
}

Coord operator*(const int& l,const Coord& r){
  return Coord(l*r.x(),l*r.y(),l*r.z());
}

Coord operator*(const double& l,const Coord& r){
  return Coord(l*r.x(),l*r.y(),l*r.z());
}

Coord operator*(const Coord& l,const int& r){
  return Coord(l.x()*r,l.y()*r,l.z()*r);
}

Coord operator*(const Coord& l,const double& r){
  return Coord(l.x()*r,l.y()*r,l.z()*r);
}

Coord& operator+=(Coord& l,const Coord& r){
  l.x(l.x()+r.x());  l.y(l.y()+r.y()); l.z(l.z()+r.z());
  return l;
}

Coord& operator-=(Coord& l,const Coord& r){
  l.x(l.x()-r.x());  l.y(l.y()-r.y()); l.z(l.z()-r.z());
  return l;
}


double l2(const Coord& a)
{
  return a.x()*a.x()+a.y()*a.y()+a.z()*a.z();
}

double scalarproduct(const Coord& l,const Coord& r){
  return double(l.x()*r.x()+l.y()*r.y()+l.z()*r.z());
}


Coord crossproduct(const Coord& a,const Coord& b)
{
  Coord result(a.y()*b.z()-a.z()*b.y(),a.z()*b.x()-a.x()*b.z(),a.x()*b.y()-a.y()*b.x()); 
  return result;
}

// a routine that converts integers to strings
const int MAXDIGITS = 20;
const int ASCII0 = 48;
string int2string(int a)
{
  int digit[MAXDIGITS];
  int d=0;
  string s("");
  if (a == 0){ s+=ASCII0; return s;}
  while( a > 0){ digit[d++] = a-10*(a/10); a /=10 ;}
  for(int i=d-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}


// this routine also handles leading zeros
string int2string(int a,int nsiffer)
{
  int digit[MAXDIGITS];
  string s("");
  for(int d=0; d<MAXDIGITS; d++){ digit[d] = a-10*(a/10); a /=10 ;}
  for(int i=nsiffer-1; i>=0; i--) s+= digit[i]+ASCII0;
  return s;
}


class Observables;
Observables operator*(Observables& left,const double& factor);

class Observables{
  friend ostream& operator<<(ostream& os, Observables& t)
  { 
    for(int i=0; i<t.GetN(); i++) os << t.Get(i) << "  ";
    return os;
  }
  friend istream& operator>>(istream& is, Observables& t)
  { 
    for(int i=0; i<t.GetN(); i++) 
      { double temp; is >> temp; t.Set(i,temp);}
    return is;
  }
 public:
  Observables(int n=NSTATICOBS):N(n),nm(0),reweight(0.),obs(n)
    {
      if(TRACE) cout << "n: " << n << endl;
      Init(); 
      if(TRACE) cout << "creating Observables: " << N << " elements\n";
    }
  Observables(const Observables& old):N(old.GetN()),nm(old.GetNm()),obs(N)
    { 
      //      if(TRACE) cout << "calling copy-constructor.........................\n";
      for(int i=0; i<N; i++) obs[i]=old.Get(i);
    }             
  ~Observables()
    { 
      if(TRACE) cout << "deleting Observables: " << N 
		     << " elements " << *this << endl;
    }
  Observables& operator=(const Observables& old)
  {
    if(TRACE) cout << "calling assignment" << endl; 
    if( this == &old) return *this;
    N  = old.GetN();
    nm = old.GetNm();
    reweight= old.GetReweight();
    for(int i=0; i<N; i++) obs[i]=old.Get(i);
    return *this;
  }
  void Init(){for(int i=0; i<N; i++) obs[i]=0.; nm=0; reweight=0.;}
  vector<double>& Vec(){return obs;}
  double* Start(){return &obs[0];}
  void Set(const int i,const double d){obs[i] = d;}
  double Get(const int i) const {return obs[i];}
  int GetN() const {return N;}
  void SetNm(const int i){nm=i;}
  int GetNm() const {return nm;}
  void SetReweight(const double a){reweight=a;}
  double GetReweight() const {return reweight;}
  void AddReweight(const double ar){reweight+=ar;}
  void Increment(){nm++;}
  void AddNm(const int i){nm+=i;}
  void AddTo(const int i,const double d){obs[i]+=d;}
  /*
  ostream& WriteAll(ostream& outfile)
    {
      outfile << nm << endl;
      outfile << reweight << endl;
      for(int i=0; i<N; i++) outfile << setprecision(16) << obs[i] << " ";
      return outfile;
    }
  void ReadAll(istream& is)
    {
      is >> nm;
      is >> reweight;
      for(int i=0; i<N; i++) is >> obs[i];
    }
  */
  /*
  void Save(string name)
    {
      if(TRACE) cout << "Saving observables to " << name << endl;
      ofstream outfile("tempobservables",ios::out);
      outfile << nm << endl;
      outfile << reweight << endl;
      for(int i=0; i<N; i++) outfile << setprecision(16) << obs[i] << " ";
      outfile.close();
      std::rename("tempobservables",name.c_str());
      if(TRACE) cout << "Saving done " << endl;
    }
  
  void Load(string name)
    {
      if(TRACE) cout << "Loading observables from " << name << endl;

      ifstream infile(name.c_str(),ios::in);

      if(!infile)
	{
	  if(TRACE) cout << "Cannot open " << name << " ,initializing cumulative obs" << endl;
	  Init(); 
	  return;
	}
      else
	{ 
	  if(TRACE) cout << "File " << name << " exists" << endl;
	  infile >> nm;
	  infile >> reweight;
	  if( infile )
	    {
	      for(int i=0; i<N; i++){ infile >> obs[i];}
	      if(TRACE) cout << "loading successful " << endl;
	      return;
	    }
	  else
	    if(TRACE) cout << "Error on reading nm" << endl;
	}    
    }
  */

  void Save(string name)
    {
      if(TRACE) cout << "Saving observables to " << name << endl;
      ofstream outfile("tempobservables",ios::out);
      outfile.write((char*) &nm,sizeof(nm));
      outfile.write((char*) &reweight,sizeof(reweight));
      outfile.write((char*) &obs[0],N*sizeof(obs[0]));
      outfile.close();
      std::rename("tempobservables",name.c_str());
      if(TRACE) cout << "Saving done " << endl;
    }
  
  void Load(string name)
    {
      if(TRACE) cout << "Loading observables from " << name << endl;

      ifstream infile(name.c_str(),ios::in);

      if(!infile)
	{
	  if(TRACE) cout << "Cannot open " << name << " ,initializing cumulative obs" << endl;
	  Init(); 
	  return;
	}
      else
	{ 
	  if(TRACE) cout << "File " << name << " exists" << endl;
	  infile.read((char*) &nm,sizeof(nm));
	  infile.read((char*) &reweight,sizeof(reweight));
	  infile.read((char*) &obs[0],N*sizeof(obs[0]));
	}    
      infile.close();
    }

#ifdef BINARYOUTPUTFILES  
  void PrintAve(string name,int i=0)
    {
      if(TRACE)    cout << "Printing averages to binary file " << name << endl;
      ofstream file(name.c_str(),ios::app);

      double factor=(reweight == 0. ? (this->GetNm() != 0 ? 1./this->GetNm():0.) : 1./reweight);
   
      Observables outobs = (*this)*factor;

      file.write( (char*) &i,sizeof(i));
      file.write( (char*) &N,sizeof(N));
      file.write( (char*) outobs.Start(),N*sizeof(double));
    }
  void PrintAve(string name,int i,vector<int>& indx,int length)
    {
      if(TRACE)    cout << "Printing averages to binary file " << name << endl;
      ofstream file(name.c_str(),ios::app);

      double factor=(reweight == 0. ? (this->GetNm() != 0 ? 1./this->GetNm():0.) : 1./reweight);
   
      Observables outobs = (*this)*factor;

      double* vobs=outobs.Start();

      file.write( (char*) &i,sizeof(i));
      file.write( (char*) &length,sizeof(int));
      for(int j=0; j<length; j++) file.write( (char*) &vobs[indx[j]],sizeof(double));
    }
  void PrintObs(string name)
    {
      if(TRACE) cout << "Printing observables to binary file " << name << endl;
      ofstream file(name.c_str(),ios::app);

      int i=0; // Dummy for no line id
      file.write( (char*) &i,sizeof(i));
      file.write( (char*) &N,sizeof(N));
      file.write( (char*) &obs[0],N*sizeof(obs[0]));
    } 
#else
  void PrintAve(string name,int i=0)
    {
      if(TRACE)    cout << "Printing averages to " << name << endl;
      ofstream file(name.c_str(),ios::app);
      
      //      cout << "N: " << this->GetNm() << endl;

      double factor=(reweight == 0. ? (this->GetNm() != 0 ? 1./this->GetNm():0.) : 1./reweight);
      Observables outobs = (*this)*factor;
      if(i!=0) file << " " << i << "    "; 
      file << " " << setprecision(16) << outobs << endl;
    }
  void PrintAve(string name,int i,vector<int>& indx,int length)
    {
      if(TRACE)    cout << "Printing averages to " << name << endl;
      ofstream file(name.c_str(),ios::app);
      double factor=(reweight == 0. ? (this->GetNm() != 0 ? 1./this->GetNm():0.) : 1./reweight);

      Observables outobs = (*this)*factor;
      double* vobs=outobs.Start();
      if(i!=0) file << " " << i << "    "; 
      for(int j=0; j<length; j++) file << " " << setprecision(16) << vobs[indx[j]];
      file << endl;
    }

  void PrintObs(string name)
    {
      if(TRACE) cout << "Printing observables to " << name << endl;
      ofstream file(name.c_str(),ios::app);
      Observables outobs = (*this);
      file << "    " << setprecision(16) << outobs << endl;
    } 
#endif
private:
  int N;
  int nm;
  double reweight; 
  vector<double> obs;
};

Observables& operator+=(Observables& left,const Observables& right){
  for(int i=0; i<left.GetN(); i++){left.Set(i,left.Get(i)+right.Get(i));}
  return left;
}

Observables& operator*=(Observables& left,const Observables& right){
  for(int i=0; i<left.GetN(); i++){left.Set(i,left.Get(i)*right.Get(i));}
  return left;
}

Observables& operator*=(Observables& left,const double& right){
  for(int i=0; i<left.GetN(); i++){left.Set(i,left.Get(i)*right);}
  return left;
}

Observables& operator/=(Observables& left,const double& right){
  double factor= 1./right;
  for(int i=0; i<left.GetN(); i++){ left.Set(i,left.Get(i)*factor);}
  return left;
}

Observables operator/(const Observables& left,int right){
  double factor= 1./right;
  Observables newobs(left.GetN());
  for(int i=0; i<left.GetN(); i++){ newobs.Set(i,left.Get(i)*factor);}
  newobs.SetNm(left.GetNm());
  return newobs;
}

Observables operator*(const Observables& left,const double& factor){
  Observables newobs(left.GetN());
  for(int i=0; i<left.GetN(); i++){ newobs.Set(i,left.Get(i)*factor);}
  return newobs;
}

Observables operator*(Observables& left,const double& factor){
  Observables newobs(left.GetN());
  for(int i=0; i<left.GetN(); i++){ newobs.Set(i,left.Get(i)*factor);}
  return newobs;
}

Observables operator*(const double& factor,const Observables& right){
  Observables newobs(right.GetN());
  for(int i=0; i<right.GetN(); i++){ newobs.Set(i,right.Get(i)*factor);}
  return newobs;
}

Observables operator*(const Observables& left,const Observables& right){
  Observables newobs(left.GetN());
  for(int i=0; i<left.GetN(); i++){ newobs.Set(i,left.Get(i)*right.Get(i));}
  return newobs;
}



#endif
