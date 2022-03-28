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
#ifdef BZSHELL
const double INNERRADIUS=0.999; /* only consider q pts closer to the BZ boundary */
#else
const double INNERRADIUS=0.; /* pick all points in the BZ */
#endif 

// universal definitions
const double PI=M_PI;
const double TWOPI=2.*PI;
const double FOURPI=4.*PI;
//const int NIL=0;
enum directions{XDIR,YDIR,ZDIR};

bool TRACE = false;

const int NPARAMS=27;
enum inparams{LINEID,J1ASTART,J1ASTOPP,NJ1A,J1BSTART,J1BSTOPP,NJ1B,RATIOFLAG1,INVRATIOFLAG1,J2ASTART,J2ASTOPP,NJ2A,J2BSTART,J2BSTOPP,NJ2B,RATIOFLAG2,INVRATIOFLAG2,J3START,J3STOPP,NJ3,J4START,J4STOPP,NJ4,N1,N2,N3,TOLERANCE};

//The RATIOFLAG1 variable is set to 1 when J1BSTART etc is interpreted as (J1B/J1A)START etc.
//The INVRATIOFLAG1 variable is set to 1 when J1ASTART etc is interpreted as (J1A/J1B)START etc.


const int MAXNQ=1000000;

// bogus parameters
const int NSTATICOBS=0;
const int MAXRUNS=1;
const int MAXPAR=100;
const int NBINS=1;
const int EQFLAG=1;


// filenames
const string RESFILE= "res.dat";
const string INFOFILE= "info.dat";
const string INPUTFILE = "J1J2J3J4.in";

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

  double* par = rp.GetPars(1);

  vector<Coord> qpts;

  const int n1=par[N1];
  const int n2=par[N2];
  const int n3=par[N3];

  logfile << "N1=" << n1 << " N2=" << n2 << " N3=" << n3 << " Tolerance=" << par[TOLERANCE] << endl;

  SetupQpts(n1,n2,n3,qpts);

  logfile << "The number of q pts in BZ= " << qpts.size() << "( INNERRADIUS= " << INNERRADIUS << ")" << endl;

  ofstream outfile(RESFILE.c_str(),ios::app);
  ofstream infofile(INFOFILE.c_str(),ios::app);



  ifstream infile(INPUTFILE.c_str());

  int nj1a=par[NJ1a];
  int nj1b=par[NJ1b];
  int nj2a=par[NJ2a];
  int nj2b=par[NJ2b];
  int nj3=par[NJ3];
  int nj4=par[NJ4];

  bool RF1=(par[RATIOFLAG1]==1.?);
  bool RF2=(par[RATIOFLAG2]==1.?);
  bool IRF1=(par[INVRATIOFLAG1]==1.?);
  bool IRF2=(par[INVRATIOFLAG2]==1.?);
  
  if(RF1==1 && IRF1==1){cout << "RATIOFLAG1 and INVRATIOFLAG1 cannot both be 1" << endl; exit(1);}
  if(RF2==1 && IRF2==1){cout << "RATIOFLAG2 and INVRATIOFLAG2 cannot both be 1" << endl; exit(1);}

  double DELTAJ1A=(nj1a==1 ? 0 : (par[J1ASTOPP]-par[J1ASTART])/(nj1a-1));
  double DELTAJ1B=(nj1b==1 ? 0 : (par[J1BSTOPP]-par[J1BSTART])/(nj1b-1));
  double DELTAJ2A=(nj2a==1 ? 0 : (par[J2ASTOPP]-par[J2ASTART])/(nj2a-1));
  double DELTAJ2B=(nj2b==1 ? 0 : (par[J2BSTOPP]-par[J2BSTART])/(nj2b-1));

  double DELTAJ3=(nj3==1 ? 0 : (par[J3STOPP]-par[J3START])/(nj3-1));
  double DELTAJ4=(nj4==1 ? 0 : (par[J4STOPP]-par[J4START])/(nj4-1));

  if(IRF1)
    logfile << "J1A/JB1 in [" << par[J1ASTART] << ":" << par[J1ASTOPP] << "] in steps of " << DELTAJ1A << " (" << nj1a << ") values" << endl;
  else
    logfile << "J1A in [" << par[J1ASTART] << ":" << par[J1ASTOPP] << "] in steps of " << DELTAJ1A << " (" << nj1a << ") values" << endl;

  if(RF1)
    logfile << "J1B/J1A in [" << par[J1BSTART] << ":" << par[J1BSTOPP] << "] in steps of " << DELTAJ1B << " (" << nj1b << ") values" << endl;
  else
    logfile << "J1B in [" << par[J1BSTART] << ":" << par[J1BSTOPP] << "] in steps of " << DELTAJ1B << " (" << nj1b << ") values" << endl;

  if(IRF2)
    logfile << "J2A/J2B in [" << par[J2ASTART] << ":" << par[J2ASTOPP] << "] in steps of " << DELTAJ2A << " (" << nj2a << ") values" << endl;
  else
    logfile << "J2A in [" << par[J2ASTART] << ":" << par[J2ASTOPP] << "] in steps of " << DELTAJ2A << " (" << nj2a << ") values" << endl;

  if(RF2)
    logfile << "J2B/J2A in [" << par[J2BSTART] << ":" << par[J2BSTOPP] << "] in steps of " << DELTAJ2B << " (" << nj2b << ") values" << endl;
  else
    logfile << "J2B in [" << par[J2BSTART] << ":" << par[J2BSTOPP] << "] in steps of " << DELTAJ2B << " (" << nj2b << ") values" << endl;

  logfile << "J3 in [" << par[J3START] << ":" << par[J3STOPP] << "] in steps of " << DELTAJ3 << " (" << nj3 << ") values" << endl;
  logfile << "J4 in [" << par[J4START] << ":" << par[J4STOPP] << "] in steps of " << DELTAJ4 << " (" << nj4 << ") values" << endl;


  
  
  for(int j4=0; j4<nj4; j4++)
      for(int j3=0; j3<nj3; j3++)
	for(int j2b=0; j2b<nj2b; j2b++)
	  for(int j2a=0; j2a<nj2a; j2a++)
	    for(int j1b=0; j1b<nj1b; j1b++)
	      for(int j1a=0; j1a<nj1a; j1a++)
		  {
		    double J1A=par[J1ASTART]+DELTAJ1A*j1a;
		    double J1B=par[J1BSTART]+DELTAJ1B*j1b;
		    if(RF1) J1B*=J1A;
		    if(IRF1) J1A*=J1B;
		    double J2A=par[J2ASTART]+DELTAJ2A*j2a;
		    double J2B=par[J2BSTART]+DELTAJ2B*j2b;
		    if(RF2) J2B*=J2A;
		    if(IRF2) J2A*=J2B;
		    double J3=par[J3START]+DELTAJ3*j3;
		    double J4=par[J4START]+DELTAJ4*j4;

		    
		    if(infile) // override for loop if infile exists
		      {
			logfile << "Infile exists, reading it: " << INPUTFILE << endl;
	      
			infile >> J1A;
			infile >> J1B;
			infile >> J2A;
			infile >> J2B;
			infile >> J3;
			infile >> J4;

			if(!infile){continue;}
		      }
		    
		    Model mymodel(J1A,J1B,J2A,J2B,J3,J4);
	
		    double minerg=mymodel.FindMinEnergyValue(qpts);
		    double maxerg=minerg+par[TOLERANCE];
	  
		    vector<Record> myr;
		    
		    int N=mymodel.FindMinimalSurface(maxerg,qpts,myr);
	  
		    logfile << "J1A=" << J1A  << "J1B=" << J1B << " J2A=" << J2A << " J2B=" << J2B 
			    << " J3=" << J3 << " J4= " << J4 << " minerg=" << minerg << " maxerg=" << maxerg << " #= " << N << " " << endl;
	  
		    outfile << J1A << " " << J1B << " " << J2A << " " << J2B << " " << J3 << " " << J4 << " ";
		      for(int i=0; i<N; i++)
			outfile << setprecision(10) << myr[i].q << " " << myr[i].erg << " ";
		    outfile << endl;

		    infofile << endl;
		    infofile << J1A << " " << J1B << " " << J2A << " " << J2B << " " << J3 << " " << J4 << " ";			  
		    infofile << "min erg=" << setprecision(10) << myr[0].erg << " N=" << N << endl;
		    for(int i=0; i<N; i++)
		      {
			Coord q=myr[i].q;
			
			double alpha=0.;
			int qstatus=IdentifyPhase(q,alpha);
			logfile << " qstatus=" << qstatus << " " << phasedescription[qstatus] << endl;

			double dHex=sqrt(HexDiff(q));
			double dSqr=sqrt(SqrDiff(q));

			infofile << "               q=[" << setprecision(10) << myr[i].q << "] (E=" << myr[i].erg << ") ";
			infofile << "(dHex= " << dHex << ",dSqr=" << dSqr << ") ";
			infofile << phasedescription[qstatus] << " ";
			if(alpha !=0.) infofile << "(a=" << alpha << ") ";
			infofile << endl;
		      }
		    infofile << endl;
		  }

  outfile.close();
  infofile.close();
  
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









