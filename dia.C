#include<string>
#include<iostream>
#include<cassert>
#include<fstream>
#include<iomanip>
#include<vector>


using namespace std;

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#ifndef M_SQRTTWO 
#define M_SQRTTWO 1.414213562373095048801689
#endif
#ifndef M_SQRTTHREE
#define M_SQRTTHREE 1.732050807568877293527446
#endif
#ifndef M_ONEHALFSQRTTHREE
#define M_ONEHALFSQRTTHREE 0.8660254037844386467637232
#endif

// what model to consider:


// universal definitions
const double PI=M_PI;
const double TWOPI=2.*PI;
const double FOURPI=4.*PI;
//const int NIL=0;
enum directions{XDIR,YDIR,ZDIR};

bool TRACE = false;

const int NPARAMS=4;
enum inparams{N1,N2,N3,TOLERANCE};

const int MAXNQ=1000000;

// bogus parameters
const int NSTATICOBS=0;
const int MAXRUNS=1;
const int MAXPAR=100;
const int NBINS=1;
const int EQFLAG=1;


// filenames
const string RESFILE= "res.dat";
const string INPUTFILE = "J1J2J3.in";

const string READIN   = "read.in";




#include <fstream>
ofstream logfile("log.txt",ios::app);

#include "globalheader.h"

#include "RunParameter.h"

#include "dia.h"

// for use with the PBS system when signal is received, 
// process raises the signal_caught_exit_now flag which is handled by the program.

volatile __sig_atomic_t signal_caught_exit_now = 0;

#include <stdlib.h>
#include <signal.h>
void sighandler(int a)
{
  //  logfile << "Sighandler: Exiting gracefully." << endl;
  logfile << "Sighandler: Raising exit flag" << endl;
  signal_caught_exit_now = 1;
}

int main()
{
  if(TRACE) cout << "Starting Main()" << endl;
  Timer mytimer;  
  logfile << "\nStarting at: " << mytimer << endl;
  mytimer.Start();

  // for use with the PBS system, when time runs out
  // the exit signal SIGTERM is sent to the program
  // which is handled so that it exits gracefully

  signal(SIGTERM, sighandler);

  
  RunParameters rp;
  logfile << "parameters: " << rp;

  double* par = rp.GetPars(0);

  vector<Coord> qpts;

  SetupQpts(par,qpts);

  logfile << "The number of q pts in BZ= " << qpts.size() << endl;
  
  ifstream infile(INPUTFILE.c_str());
  ofstream outfile(RESFILE.c_str(),ios::app);
  
  while(infile)
    {
      double J1, J2, J3;
      infile >> J1;
      infile >> J2;
      infile >> J3;
      if(!infile){continue;}

      Model mymodel(J1,J2,J3);
      
      double minerg=mymodel.FindMinEnergyValue(qpts);
      double maxerg=minerg+par[TOLERANCE];
      
      vector<Coord> minqs;
      vector<double> ergs;
      
      int N=mymodel.FindMinimalSurface(maxerg,qpts,minqs,ergs);
      
      logfile << "J1=" << J1 << "J2=" << J2 << "J3=" << J3 << " minerg=" << minerg << " maxerg=" << maxerg << " #= " << N << " " << endl;
      
      outfile << J1 << " " << J2 << " " << J3 << " ";
      for(int i=0; i<N; i++)
	outfile << minqs[i] << " " << ergs[i] << " ";
      outfile << endl;
    }
  outfile.close();
  infile.close();
    
  
  mytimer.Stop();
  double time_spent = mytimer.GetTimeElapsed();
  logfile << "Time elapsed:  " << time_spent << " sec. :" 
	  << time_spent/3600. << " hours." << endl;
  logfile << "Ending at: " << mytimer;
  logfile << "---Done--- " << endl;
  logfile.close();

  ofstream terminationsign("end.exec");
  terminationsign << "execution ended at " << mytimer;
  terminationsign.close();

  cout << "program ended \n";
  exit(0);
}









