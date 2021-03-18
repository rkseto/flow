// figure out C{6} with cumulants
// now do differential flow

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "TComplex.h"
//#include "TBenchMark.h"
#include <fstream>
#include <sstream>
#include <string>

#include "TMath.h"

#define PI 3.14159
#define TWOPI 6.28318

void read() {
  /*
Now this adds pt dependence and differential flow
This adds the direct calculations of the Q's to save time

Parameters of fit should be vn^2
3 particle correlation = vn*vn*v2n
4 particle correlation should be -vn^4  
   */


  // we will define our histograms here

  TH1D *pttrackh = new TH1D("pttrackh","pttrackh",100,0,4);
  TH1D *phitrackh = new TH1D("phitrackh","phitrackh",100,0,2*PI);
  TH1D *ytrackh = new TH1D("ytrackh","ytrackh",100,-1.2,1.2);
  

  std::ifstream infile("amptrunning.dat");

      int iev, testnum, nparticles, npart_proj, npart_targ, nprojelas, nprojinelas, ntargelas, ntarginelas;
    double b, zero;

  // loop over events
  std::string line;
  int nevents=100;
  int ncount=10;
  while(std::getline(infile,line)){
    cout<<line<<endl;
    std::istringstream iss(line);
    if(!(iss >> iev>> testnum>> nparticles >> b>>npart_targ>>npart_proj>>nprojelas>>nprojinelas>>ntargelas>>ntarginelas>>zero)){cout<<" problem"<<endl; break;} //error
    if(iev%ncount==0)cout<<" iev="<<iev<<endl;
    cout<<" nparticles="<<nparticles<<endl;
    for(int ipart=0; ipart<nparticles; ipart++){
      std::getline(infile,line);
      std::istringstream iss(line);
      //      cout<<line<<endl;
   }
  }
}
