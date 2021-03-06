// figure out C{6} with cumulants
// now do differential flow

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "TComplex.h"
#include "TMath.h"

#define PI 3.14159
#define TWOPI 6.28318

int nparticles=0;
double phitrack[10000];
double pttrack[10000];
double ytrack[10000];
int nparticles_previous=0;
double phitrack_previous[10000];
double pttrack_previous[10000];
double ytrack_previous[10000];

const int nparticlesmin=6;

// define a function with parameters
Double_t fitf(Double_t *xx,Double_t *par);
//TF1 *flowf = new TF1("flowf",fitf,0,TWOPI,5);

int throwevent(int ievent);

void flow8work() {
  /*
Now this adds pt dependence and differential flow
This adds the direct calculations of the Q's to save time

Parameters of fit should be vn^2
3 particle correlation = vn*vn*v2n
4 particle correlation should be -vn^4  
   */


  // we will define our histograms here
  gROOT->Time();
  TStopwatch stopwatch;
  stopwatch.Start();
  
  TH1D *pttrackh = new TH1D("pttrackh","pttrackh",100,0,4);
  TH1D *phitrackh = new TH1D("phitrackh","phitrackh",100,0,2*PI);
  TH1D *ytrackh = new TH1D("ytrackh","ytrackh",100,-1.2,1.2);
  
  TH1D *v2h = new TH1D("v2h","v2 contribution track by track",200,-1.2,1.2);
  TH1D *reseventplaneh = new TH1D("reseventplaneh","event plane reslution event",2000,0.0,1.00);
   TH1D *phitrackseventh = new TH1D("phitrackseventh","phitrackseventh",100,0,TWOPI);

   TH1D *dphirealh = new TH1D("dphirealh","dphirealh",100,-PI,PI);
   TH1D *dphimeash = new TH1D("dphimeash","dphimeash",100,-PI,PI);
   
   TH1D *dphipairs_same = new TH1D("dphipairs_same","dphipairs_same",100,-PI,PI);
   TH1D *dphipairs_mixed = new TH1D("dphipairs_mixed","dphipairs_mixed",100,-PI,PI);

   TH1D *d4phih = new TH1D("d4phih","d4phi",200,-4*PI,4*PI);

   TH1D *PSI2h = new TH1D("PSI2h","event plane angle",100,-PI,PI);
   TH1D *PSI2PSI2foundh = new TH1D("PSI2PSI2foundh","event plane angle real-found",100,-PI,PI);
   TH1D *PHI2resh = new TH1D("PHI2resh","event plane resoulution 2 subevents",100,-PI,PI);   

   int nharmonic=2;  //2
   int nevents=4000; //1000, for nharmonic=4 need 4000

   int inputtype=0;
   cout<<" input type 0=throw 1=ampt 2=rqmd 4=therminator: ";
   cin>>inputtype;
   cout<<" inputtype set to "<<inputtype<<endl;
   if(inputtype==0)throwevent(0);

   bool doloop=false;
   bool doloop2=false;
   bool doloop3=false;
   bool doloop4=false;
   bool ptdepv=false;

   cout<<" nharmonic(2): ";
   cin>>nharmonic;
   cout<<" nharmonic set to "<<nharmonic<<endl;

   cout<<" nevents(1000): ";
   cin>>nevents;
   cout<<" nevents set to "<<nevents<<endl;
   
   cout<<" do loops 2 (needed to do a fit to the 2particle CF) 0-no 1-yes(1): ";
   cin>>doloop2;
   cout<<" doloop2 set to "<<doloop2<<endl;
   cout<<" do loops 3 0-no 1-yes: ";
   cin>>doloop3;
   cout<<" doloop3 set to "<<doloop3<<endl;
   cout<<" do loops 4 0-no 1-yes: ";
   cin>>doloop4;
   cout<<" doloop4 set to "<<doloop4<<endl;

   int ncount=10;
   if(nevents>=10000)ncount=10000;
   //******** NOW START WORK ***************
   // OK now lets throw events.
   
   double v2obs=0;
   int nv2count=0;

   double res_eventplane=0;
   int nrescount=0;
   
   double cdphi3sumtot=0;
   double cdphisumtot=0;
   double cdphi4sumtot=0;

   double CF2=0; // cf from <<2>>
   double CF4=0; // cf from <<4>> <<2>>
   double CF6=0; // cf from <<6>> <<4>> <<2>>
   double bb2bb=0; // <<2>>
   double bb4bb=0; // <<4>>
   double bb6bb=0; // <<6>>
   double W2sum=0;   //<
   double W4sum=0;   //<
   double W6sum=0;  //<
   double W3sum=0;  //<
   
   TH1D *hb2pb = new TH1D("b2pb","b2pb",20,0,2);
   hb2pb->Sumw2(true);
   
   // loop over events
   for(int ievent=0; ievent<nevents; ievent++){
     if(ievent%ncount==0)cout<<" ievent="<<ievent<<endl;
     
     // make an event **************** set V2 and PHI here
     if(inputtype==0)throwevent(1);

     // analyze event
     int ncdphisum=0;
     double cdphisum=0;
     int ncdphi3sum=0;
     double cdphi3sum=0;
     int ncdphi4sum=0;
     double cdphi4sum=0;

     double Q2x=0;      // used to find the reaction plane
     double Q2y=0;

     double Q2xsub1=0;
     double Q2ysub1=0;
     double Q2xsub2=0;
     double Q2ysub2=0;
     int nhalf=nparticles/2;

     for(int i=0; i<nparticles;i++){     
       pttrackh->Fill(pttrack[i]);
       phitrackh->Fill(phitrack[i]);
       ytrackh->Fill(ytrack[i]);
     }
     //Finding reaction plane
     Double Weight=1;
     for(int i=0; i<nparticles;i++){
       phitrackseventh->Fill(phitrack[i]);
       Q2x=Q2x+Weight*cos(nharmonic*phitrack[i]);
       Q2y=Q2y+Weight*sin(nharmonic*phitrack[i]);
       
       if(i<nhalf){
	 Q2xsub1=Q2xsub1+Weight*cos(nharmonic*phitrack[i]);
	 Q2ysub1=Q2ysub1+Weight*sin(nharmonic*phitrack[i]);
       }else if(i<2*nhalf){
	 Q2xsub2=Q2xsub2+Weight*cos(nharmonic*phitrack[i]);
	 Q2ysub2=Q2ysub2+Weight*sin(nharmonic*phitrack[i]);
       }else{
       }
     }
     double PSI2found=atan2(Q2y,Q2x)/nharmonic;
     double dpsieventfound=0-PSI2found; // CHECK THIS and plot
     if(dpsieventfound>PI){
       dpsieventfound=TWOPI-dpsieventfound;
     }else if(dpsieventfound<-PI){
       dpsieventfound=-TWOPI-dpsieventfound;
     }else{
     }
     
     if(dpsieventfound>PI/2){
       dpsieventfound=PI-dpsieventfound;
     }else if(dpsieventfound<-PI/2){
       dpsieventfound=-PI-dpsieventfound;
     }else{
     }
     
     
     PSI2PSI2foundh->Fill(dpsieventfound);     
     double PSI2foundsub1=atan2(Q2ysub1,Q2xsub1)/nharmonic;
     double PSI2foundsub2=atan2(Q2ysub2,Q2xsub2)/nharmonic;
     double PSI2restmp=0;     
     double dpsi=PSI2foundsub1-PSI2foundsub2; // CHECK THIS and plot
     if(dpsi>PI){
       dpsi=TWOPI-dpsi;
     }else if(dpsi<-PI){
       dpsi=-TWOPI-dpsi;
     }else{
     }
     if(dpsi>PI/2){
       dpsi=PI-dpsi;
     }else if(dpsi<-PI/2){
       dpsi=-PI-dpsi;
     }else{
     }
     PHI2resh->Fill(dpsi);
     
     PSI2restmp=cos(nharmonic*(dpsi)); 
     //     cout<<" PSI2found="<<PSI2found<<" subevent1="<<PSI2foundsub1<<" subevent2="<<PSI2foundsub2<<" restmp="<<PSI2restmp <<endl;
     
     reseventplaneh->Fill(PSI2restmp);
     //finish of reaction plane
     
     
     res_eventplane+=PSI2restmp;
     nrescount++;
     
     TComplex Qnharmonic[12]={0,0,0,0,0,0,0,0,0,0,0,0}; //<<   
	  
     // now that we have the reaction plane calculate v2 and Q's
     for(int i=0; i<nparticles;i++){
       double dphi=phitrack[i]-PSI2found-PI;  // CHECK THIS
       if(dphi>PI){
	 dphi=TWOPI-dphi;
       }else if(dphi<-PI){
	 dphi=-TWOPI-dphi;
       }else{
       }
       
       double dphireal=phitrack[i]-0-PI;  // CHECK THIS
       if(dphireal>PI){
	 dphireal=TWOPI-dphireal;
       }else if(dphireal<-PI){
	 dphireal=-TWOPI-dphireal;
       }else{
       }
       
       dphirealh->Fill(dphireal);
       dphimeash->Fill(dphi);
       
       
       if(dphi>PI){
	 dphi=TWOPI-dphi;
       }else if(dphi<-PI){
	 dphi=-TWOPI-dphi;
       }else{
       }
       //-------------------
       
       v2obs+=cos(nharmonic*(dphi)); 
       v2h -> Fill(cos(nharmonic*(dphi)));
       nv2count++;
       
       // calculate Qvectors
       for(int j=0; j<12; j++){
	 Qnharmonic[j]+=TComplex(cos((j+1)*phitrack[i]),sin((j+1)*phitrack[i]),false);     
       }
       
     }// end loop over particles
     
     // see PRC 83 044913
     
     double qn2=Qnharmonic[nharmonic-1].Rho2();// Q_nharmonic^2
     double b2b=qn2;
     b2b=(b2b-nparticles)/nparticles/(nparticles-1);
     
     TComplex tmpc;
     
     double qn4=qn2*qn2;// Q_nharmonic^4 
     double q2n2=Qnharmonic[2*nharmonic-1].Rho2();// Q_2*nharmonic^2
     tmpc=Qnharmonic[2*nharmonic-1]*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1]);
     double q2nnstarnstar=tmpc.Re();
     double b4b = qn4+q2n2-2*q2nnstarnstar-2*2*(nparticles-2)*qn2+2*nparticles*(nparticles-3);
     b4b=b4b/(nparticles*(nparticles-1)*(nparticles-2)*(nparticles-3));
     
     double qn6=qn4*qn2; //Q_n^6
     double q3n2=Qnharmonic[3*nharmonic-1].Rho2(); //Q3n^2
     double q2n2n2=q2n2*qn2;
     tmpc=Qnharmonic[2*nharmonic-1]*Qnharmonic[nharmonic-1]*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1]);
     double q2nnnstarnstarnstar = tmpc.Re();
     tmpc=Qnharmonic[3*nharmonic-1]*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1]);
     double q3nnstarnstarnstar = tmpc.Re();
     tmpc=Qnharmonic[3*nharmonic-1]*TComplex::Conjugate(Qnharmonic[2*nharmonic-1])*TComplex::Conjugate(Qnharmonic[nharmonic-1]);
     double q3n2nstarnstar = tmpc.Re();
     double t1=qn6+9*q2n2*qn2-6*q2nnnstarnstarnstar;
     t1=t1/nparticles/(nparticles-1)/(nparticles-2)/(nparticles-3)/(nparticles-4)/(nparticles-5);
     double t2=q3nnstarnstarnstar-3*q3n2nstarnstar;
     t2=t2/nparticles/(nparticles-1)/(nparticles-2)/(nparticles-3)/(nparticles-4)/(nparticles-5);
     double t3=9*(nparticles-4)*q2nnstarnstar+2*q3n2;
     t3=t3/nparticles/(nparticles-1)/(nparticles-2)/(nparticles-3)/(nparticles-4)/(nparticles-5);
     double t4=qn4+q2n2;
     t4=t4/nparticles/(nparticles-1)/(nparticles-2)/(nparticles-3)/(nparticles-5);
     double t5=qn2;
     t5=t5/nparticles/(nparticles-1)/(nparticles-3)/(nparticles-4);
     double t6=6/(nparticles-1)/(nparticles-2)/(nparticles-3);
     double b6b=t1+4*t2+2*t3-9*t4+18*t5-t6;
     b6b=b6b/(nparticles*(nparticles-1)*(nparticles-2)*(nparticles-3)*(nparticles-4)*(nparticles-5));  //<<
       
     //     cout<<" ts "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<" "<<t5<<" "<<t6<<endl;
     
     double W2=nparticles*(nparticles-1);   //<
     double W4=nparticles*(nparticles-1)*(nparticles-2)*(nparticles-3);   //<
     double W6=nparticles*(nparticles-1)*(nparticles-2)*(nparticles-3)*(nparticles-4)*(nparticles-5);  //<
     double W3=nparticles*(nparticles-1)*(nparticles-2);   //<
     //    W2=1; W3=1; W4=1; W6=1;

     bb2bb+=W2*b2b;
     bb4bb+=W4*b4b;
     bb6bb+=W6*b6b;
     W2sum+=W2;   //<
     W4sum+=W4;   //<
     W6sum+=W6;  //<
     W3sum+=W3;   //<
   
     // pairs (cumulants calculations - by looping)
     for(int i=0; i<nparticles;i++){
       if(!doloop2)break;
       for(int j=0; j<nparticles;j++){   //<<
	 if(i==j)continue;  //<<
	 double dphi=phitrack[i]-phitrack[j];
	 double cdphi=cos(nharmonic*(dphi));
	 ncdphisum++;
	 cdphisum+=cdphi;
	 
	 if(dphi>PI){
	   dphi=TWOPI-dphi;
	 }else if(dphi<-PI){
	   dphi=-TWOPI-dphi;
	 }else{
	 }
	 dphipairs_same->Fill(dphi);
	 
	 //-------------3 and 4 particle correlations------------------	 
	 
	 for(int k=0; k<nparticles; k++){ //<<
	   if(!doloop3) break;
	   if(i==k)continue;  //<<
	   if(j==k)continue;  //<<	   

	   
	   for(int l=0; l<nparticles; l++){    //<<
	   if(!doloop4) break;
	     if(i==l)continue;  //<<
	     if(j==l)continue;  //<<
	     if(k==l)continue;  //<<	   
	     double d4phi=phitrack[i]+phitrack[j]-phitrack[k]-phitrack[l];
	     ncdphi4sum++;
	     cdphi4sum+=cos(nharmonic*d4phi);
	     d4phih->Fill(d4phi);
	   } // end l
	   
	   double d3phi=phitrack[i]+phitrack[j]-2*phitrack[k];
	   ncdphi3sum++;
	   cdphi3sum+=cos(nharmonic*d3phi);
	 } // end k
	 
	 //----------------------------------------------------
	 
       }// end j same
       
       
       if(ievent>0){
	 for(int j=i+1; j<nparticles_previous;j++){
	   double dphi=phitrack[i]-phitrack_previous[j];
	   
	   if(dphi>PI){
	     dphi=TWOPI-dphi;
	   }else if(dphi<-PI){
	     dphi=-TWOPI-dphi;
	   }else{
	   }
	   
	   double cdphi_mixed=cos(nharmonic*(dphi));
	   dphipairs_mixed->Fill(dphi);
	 } // end j mixed
       } 
     } // end i - all particles
     
     // note- in the following the W's cancel if W=nparticles*(nparticles-1)... - see eqn 5 and 6 from 1010.0233
     cdphisumtot+=W2*cdphisum/(nparticles*(nparticles-1)); //<
     cdphi3sumtot+=W3*cdphi3sum/(nparticles*(nparticles-1)*(nparticles-2));  //<
     cdphi4sumtot+=W4*cdphi4sum/(nparticles*(nparticles-1)*(nparticles-2)*(nparticles-3));  //<
     
     // finally save the event
     nparticles_previous=nparticles;
     for(int i=0; i<nparticles;i++){
       phitrack_previous[i]=phitrack[i];
       pttrack_previous[i]=pttrack[i];
       ytrack_previous[i]=ytrack[i]; 
     }
     
     //------------------------------------------------------------
     // now figure out differential flow
     TH1D *ptcosh = new TH1D("ptcosh","ptcosh",20,0,2);
     TH1D *ptsinh = new TH1D("ptsinh","ptsinh",20,0,2);
     TH1D *ptcounth = new TH1D("ptcounth","ptcpunth",20,0,2);
     ptcosh->Sumw2(true);
     ptsinh->Sumw2(true);
     ptcounth->Sumw2(true);
     
     int npartref=0;
     double cosref=0;
     double sinref=0;
     double yrefpoi=0.;  // y< yrefpoi = refparticle, y> yrefpoi = particle of interest
     for (int i=1; i<nparticles; i++){
       if(ytrack[i]>yrefpoi){
	 ptcounth->Fill(pttrack[i],1);
	 ptcosh->Fill(pttrack[i],cos(nharmonic*phitrack[i]));
	 ptsinh->Fill(pttrack[i],sin(nharmonic*phitrack[i]));
       }else{
	 npartref++;
	 cosref+=cos(nharmonic*phitrack[i]);
	 sinref+=sin(nharmonic*phitrack[i]);
       }
     } //  i
     
     /*    
	   for(int ii=0; ii<ptcounth->GetNbinsX();ii++){
	   cout<<" ii="<<ii<<" n="<<ptcounth->GetBinContent(ii)
	   <<" cos="<<ptcosh->GetBinContent(ii)<<" sin="<<ptsinh->GetBinContent(ii)<<endl;
	   }
     */
     
     for(int ii=0; ii<ptcounth->GetNbinsX();ii++){
       if(ptcounth->GetBinContent(ii)==0)ptcounth->SetBinContent(ii,1);
     }
     
     ptcosh->Divide(ptcounth);
     ptsinh->Divide(ptcounth);
     cosref/=double(npartref);
     sinref/=double(npartref);
     hb2pb->Add(ptcosh,cosref);
     hb2pb->Add(ptsinh,sinref);
     
     
     delete ptcosh;
     delete ptsinh;
     delete ptcounth;
     
     
   } // end loop over events

   
     //------------------------------------------------------------   
   // CF for differential flow
   hb2pb->Scale(1./double(nevents));
   
   // calculate CFs
   //   bb2bb=bb2bb/float(nevents);
   bb2bb=bb2bb/W2sum; //<<
   
   hb2pb->Scale(1./sqrt(bb2bb));   //this is v2(pt)
   //   cout<<" bb4bb="<<bb4bb<<endl;
   //   bb4bb=bb4bb/float(nevents);
   //   bb6bb=bb6bb/float(nevents);
   bb4bb=bb4bb/W4sum; //<<
   bb6bb=bb6bb/W6sum; //<<
   CF2=bb2bb;
   CF4=bb4bb-2*bb2bb*bb2bb;
   CF6=bb6bb-9*bb2bb*bb4bb+12*bb2bb*bb2bb*bb2bb;
   
   // calcualte eventplane resolution
   res_eventplane/=nrescount;
   
   //calculate a final v2
   v2obs/=nv2count;
   
   // DIVIDE
   TH1D *dphipairs_divided = (TH1D*)dphipairs_same->Clone("dphipairs_divided");
   dphipairs_divided->Sumw2();
   dphipairs_divided->Divide(dphipairs_mixed);
   
   // now do a fit
   // the parameters should be vn^2
   TF1 *myfit = new TF1("myfit","[0]*(1.0+2*[1]*cos(1*x) + 2*[2]*cos(2*x)+ 2*[3]*cos(3*x)+ 2*[4]*cos(4*x)+ 2*[5]*cos(5*x) )", -PI, PI);
   dphipairs_divided->Fit("myfit");
   double fitnorm=myfit->GetParameter(0);
   double v1sq=myfit->GetParameter(1);
   double v2sq=myfit->GetParameter(2);
   double v3sq=myfit->GetParameter(3);
   double v4sq=myfit->GetParameter(4);
   double v5sq=myfit->GetParameter(5);
   
   //    myfit->SetParameters(fitnorm,v1sq,v2sq,v3sq,v4sq,0);
   
   myfit->SetParameters(fitnorm,v1sq,0,0,0,0);
   myfit->SetLineColor(3);
   myfit->DrawCopy("SAME");
   myfit->SetParameters(fitnorm,0,v2sq,0,0,0);
   myfit->SetLineColor(4);
   myfit->DrawCopy("SAME");
   myfit->SetParameters(fitnorm,0,0,v3sq,0,0);
   myfit->SetLineColor(5);
   myfit->DrawCopy("SAME");
   myfit->SetParameters(fitnorm,0,0,0,v4sq,0);
   myfit->SetLineColor(6);
   myfit->DrawCopy("SAME");
   
   
   //--- now print out results ---------------------    
   cout<<endl;
   cout<<" nevents="<<nevents<<" nharmonic="<<nharmonic<<" doloop3="<<doloop3<<" doloop4="<<doloop4<<endl;
   cout<<endl;
   
   cout<<" CF2="<<CF2<<" sqrt(CF2):  v"<<nharmonic<<"{2}="<<sqrt(CF2)<<endl;
   cout<<" CF4="<<CF4<<" sqrt(sqrt(-CF4)):  v"<<nharmonic<<"{4}="<<sqrt(sqrt(-CF4))<<endl;
   cout<<" CF6="<<CF6<<" 6throot((CF6)/4.):  v"<<nharmonic<<"{6}="<<std::cbrt(  sqrt( (CF6/4.) )  )<<endl<<endl;
   
   cout<<"FINAL values from event plane" << endl << endl;
   cout<<" FINAL event plane resolution="<<res_eventplane<<endl;
   cout<<" FINAL v"<<nharmonic<<" from Event plane="<<v2obs<<" nharmonic="<<nharmonic<<endl;
   
   cout<<" 2 particle correlation fit: v1="<<sqrt(v1sq)
       <<"  v2="<<sqrt(v2sq)
       <<"  v3="<<sqrt(v3sq)
       <<"  v4="<<sqrt(v4sq)
       <<"  v5sq="<<v5sq
       <<endl<<endl;
   
   cout<<"if loops done FINAL values from particle cumulants" << endl<<endl;
   
   //   double cumulant2 = 1.*cdphisumtot/nevents; 
   //   double cumulant3 = 1.*cdphi3sumtot/nevents; 
   //   double cumulant4 = 1.*cdphi4sumtot/nevents;
   double cumulant2 = 1.*cdphisumtot/W2sum; //<<
   double cumulant3 = 1.*cdphi3sumtot/W3sum; //<<
   double cumulant4 = 1.*cdphi4sumtot/W4sum-2*cumulant2*cumulant2; //<<
   if(!doloop3)cumulant3=0;
   if(!doloop4)cumulant4=0;
   
   cout<<" 2 particle cumulant = "<<cumulant2<<endl;
   cout<<" 3 particle cumulant = "<<cumulant3<<endl;
   cout<<" 4 particle cumulant = "<<cumulant4<<endl;
   cout<<endl;
   
   cout<<"if loops done FINAL flow calculations from particle cumulants" << endl<<endl;
   cout<<" from 2 particle cumulant: v"<<nharmonic<<"="<<sqrt(cumulant2)<<endl;
   
   //   double cumulant4minus2 = 1.*(cumulant4-(2*(cumulant2*cumulant2)));   //<<
   cout<<" from 4 particle cumulant: v"<<nharmonic<<"="<<sqrt(sqrt(-cumulant4))<<endl<<endl; //<<
   double cum3=0;

   //   if(nharmonic==1)cum3=v1*v1*v2;
   //   if(nharmonic==2)cum3=v2*v2*v4;
   //   cout<<" nharmonic="<<nharmonic<<", 3 particle cumulant="<<cdphi3sumtot/W3sum<<" we expect "<<cum3<<" i.e. v"<<nharmonic<<"*v"<<nharmonic<<"*v"<<2*nharmonic<<endl<<endl;
   cout<<" note: W3sum="<<W3sum<<" nevents="<<nevents<<endl;
   //---------------------
   
   // Open a ROOT file and save the formula, function and histogram
   //
   TFile *myfile = new TFile("flow8.root","RECREATE");
   
   pttrackh->Write();
   phitrackh->Write();
   ytrackh->Write();
   
   //   flowf->Write();
   //   h1f->Write();
   reseventplaneh->Write();
   v2h->Write();
   phitrackseventh->Write();
   dphipairs_same->Write();
   dphipairs_mixed->Write();
   dphipairs_divided->Write();
   d4phih->Write();
   PSI2h->Write();
   PSI2PSI2foundh->Write();
   PHI2resh->Write();
   dphirealh->Write();
   dphimeash->Write();
   hb2pb->Write();
   
   stopwatch.Stop();
   stopwatch.Print();
   
}

Double_t fitf(Double_t *xx,Double_t *par) {
  double v1=par[1];
  double v2=par[2];
  double v3=par[3];
  double v4=par[4];
  double PSI=par[0];
  
  double x=xx[0];  

  double  fitval = 1./TWOPI*(
			     1.+2*v1*cos(1*(x-PSI))+2*v2*cos(2*(x-PSI))+2*v3*cos(3*(x-PSI))+2*v4*cos(4*(x-PSI))
			     );  
  return fitval;
}

int throwevent(int ido=0){

  // zero for init
  // 1  for throw event

  static bool iprintthrowevent=true;  
  static TF1 *flowf = new TF1("flowf",fitf,0,TWOPI,5);
  static  TRandom3 *ran0 = new TRandom3();
  static TRandom3 *ran1 = new TRandom3();   
  
  static int nparticlesav=200; //100
  static int ptdepv=0;
  static double v1=0.10; //0.02   // 0.05; // 0.1
  static double v2=0.20; //0.05   //0.15;  // 0.2
  static double v3=0.15; //0.01  //0.07;  //0.15
  static double v4=0.05; //0.02  //0.03;  //0.05


  
  if(ido==0){
    cout<<" inputtype: throwevents INIT"<<endl;
    cout<<"Initial values: v1="<<v1<<"  v2="<<v2<<"  v3="<<v3<<"  v4="<<v4
	<<endl<<endl;
    
    cout<<" nparticles/event average(100) max is 10000, for pt dependence need at least 500 for T=300 MeV: ";
    cin>>nparticlesav;
    if(nparticlesav>10000)nparticlesav=10000;
    cout<<" nparticles set to "<<nparticlesav<<endl;
    
    cout<<" pt dependent v's 0-no 1-yes: ";
    cin>>ptdepv;
    cout<<" ptdepv set to "<<ptdepv<<endl;
    double PSI2=0.;
    
    flowf->SetParNames("PSI2","v1","v2","v3","v4");   
    flowf->SetParameters(PSI2,v1,v2,v3,v4);

    TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,10,700,900);
    c1->SetFillColor(18);
    
    TPad *pad1 = new TPad("pad1","The pad with the function",0.05,0.50,0.95,0.95,21);
    TPad *pad2 = new TPad("pad2","The pad with the histogram",0.05,0.05,0.95,0.45,21);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->GetFrame()->SetFillColor(42);
    pad1->GetFrame()->SetBorderMode(-1);
    pad1->GetFrame()->SetBorderSize(5);
    flowf->SetLineColor(4);
    flowf->SetLineWidth(6);
    flowf->Draw();
    TPaveLabel *lfunction = new TPaveLabel(5,39,9.8,46,"v2 function");
    lfunction->SetFillColor(41);
    lfunction->Draw();
    c1->Update();
    //
    // Create a one dimensional histogram (one float per bin)
    // and fill it following the distribution in function flowf.
    //
    TH1D *h1f = new TH1D("h1f","Test of an Observed Distribution",200,0,TWOPI);
    
    pad2->cd();
    pad2->GetFrame()->SetFillColor(42);
    pad2->GetFrame()->SetBorderMode(-1);
    pad2->GetFrame()->SetBorderSize(5);
    h1f->SetFillColor(45);
    h1f->FillRandom("flowf",1000000);
    h1f->Draw();
    c1->Update();
    
    // the previous stuff was just for illustration
  }  //init


  if(ido!=0){

    if(iprintthrowevent)cout<<" inputtype: throwevents THOWEVENT"<<endl;       
    iprintthrowevent=false;
    //   nparticles=nparticlesav;
   nparticles=ran0->Poisson(nparticlesav);
   if(nparticles<nparticlesmin)cout<<" nparticles small "<<nparticles<<endl;

     double v2Event=v2; //potentially can wiggle around (fluctuate)
     double PSI2Event=0.0;
     PSI2Event=ran0->Rndm()*TWOPI-PI;

     //     PSI2h->Fill(PSI2Event);
     //     PSI2PSI2foundh->Fill(PSI2Event);
       

     flowf->SetParameters(PSI2Event,v1,v2,v3,v4); // for testing &&&
       
     for(int i=0; i<nparticles;i++){
       pttrack[i]=gRandom->Exp(0.300);  // .250 GeV tau
       //       pttrack[i]=gRandom->Rndm()*2.;  // flat in pt
       pttrack[i]+=0.180;
       double vv1=v1;
       double vv2=v2;
       double vv3=v3;
       double vv4=v4;
       
       double ppt=pttrack[i];
       if(ppt>1.)ppt=1;
       vv1=ppt*v1;
       vv2=ppt*v2;
       vv3=ppt*v3;
       vv4=ppt*v4;
       
       if(ptdepv)flowf->SetParameters(PSI2Event,vv1,vv2,vv3,vv4);      
       //       if(!ptdepv)flowf->SetParameters(PSI2Event,v1,v2,v3,v4); // for testing &&&
       
       phitrack[i]=flowf->GetRandom();
       ytrack[i]=(gRandom->Rndm()-.5)*2; // 
       //      ytrack[i]=1.; // set all y to 1 now
        //       if(phitrack[i]>1. && phitrack[i]<1.9)i--; // mess up acceptance
      }

  }
  
  return nparticles;  
}
