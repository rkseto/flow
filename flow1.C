#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TFormula.h"
#include "TPaveLabel.h"
#include "TFile.h"
//#include "TBenchMark.h"

#include "TMath.h"

#define PI 3.14159
#define TWOPI 6.28318

// define a function with parameters
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


void flow1() {
  /*
Parameters of fit should be vn^2
3 particle correlation = vn*vn*v2n
4 particle correlation should be -vn^4    (there seems to be a sign wrong?)
   */


  // we will define our histograms here
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
   
   TCanvas *c1 = new TCanvas("c1","The FillRandom example",200,10,700,900);
   c1->SetFillColor(18);

   TPad *pad1 = new TPad("pad1","The pad with the function",0.05,0.50,0.95,0.95,21);
   TPad *pad2 = new TPad("pad2","The pad with the histogram",0.05,0.05,0.95,0.45,21);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();

   //   gBenchmark->Start("flow1");
   //
   double v1=0.1; // 0.05;   // 0.1
   double v2=0.2; //0.15;    // 0.2
   double v3=0.15; //0.07;   //0.15
   double v4=0.05; //0.03;    //0.05

   int nharmonic=2;  //2
   int nparticlesav=100;
   int nparticlesmin=6;
   int nevents=1000; //4000
   int ncount=100;
   if(nevents>=10000)ncount=1000;
   double PSI2=0.;
   
   TF1 *flowf = new TF1("flowf",fitf,0,TWOPI,5);
   flowf->SetParNames("PSI2","v1","v2","v3","v4");   

   flowf->SetParameters(PSI2,v1,v2,v3,v4);
      
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
   pad2->cd();
   pad2->GetFrame()->SetFillColor(42);
   pad2->GetFrame()->SetBorderMode(-1);
   pad2->GetFrame()->SetBorderSize(5);
   TH1D *h1f = new TH1D("h1f","Test of an Observed Distribution",200,0,TWOPI);
   h1f->SetFillColor(45);
   h1f->FillRandom("flowf",100000);
   h1f->Draw();
   c1->Update();
   // the previous stuff was just for illustration

   //******** NOW START WORK ***************
   // OK now lets throw events.
   int nparticles=0;
   double phitrack[10000];
   double pttrack[10000];
   double ytrack[10000];
   int nparticles_previous=0;
   double phitrack_previous[10000];
   double pttrack_previous[10000];
   double ytrack_previous[10000];

   double v2obs=0;
   int nv2count=0;

   double res_eventplane=0;
   int nrescount=0;

   TRandom *ran0 = new TRandom();
   TRandom *ran1 = new TRandom();   

   double cdphi3sumtot=0;
   double cdphisumtot=0;
   double cdphi4sumtot=0;

   nparticles=nparticlesav;

   // loop over events
   for(int ievent=0; ievent<nevents; ievent++){
     if(ievent%ncount==0)cout<<" ievent="<<ievent<<endl;
   // make an event **************** set V2 and PHI here
     nparticles=ran0->Poisson(nparticlesav);
     if(nparticles<nparticlesmin)continue;
     //     nparticlesav=nparticlesav;
     double v2Event=v2; //potentially can wiggle around (fluctuate)
     double PSI2Event=0.0;
     //     PSI2Event=ran0->Rndm()*TWOPI;
     PSI2Event=ran0->Rndm()*TWOPI-PI;
     PSI2h->Fill(PSI2Event);
     PSI2PSI2foundh->Fill(PSI2Event);
       
     //     flowf->SetParameters(v2Event,PSI2Event,float(nharmonic));
     flowf->SetParameters(PSI2Event,v1,v2,v3,v4);

     for(int i=0; i<nparticles;i++){
       phitrack[i]=flowf->GetRandom();
       pttrack[i]=ran1->Poisson(1.); // set all pt to 1 GeV now
       ytrack[i]=1.; // set all pt to 1 GeV now

       //       if(phitrack[i]>1. && phitrack[i]<1.9)i--; // mess up acceptance
       
     }
     
     // analyze event
     int ncdphisum=0;
     double cdphisum=0;
     int ncdphi3sum=0;
     double cdphi3sum=0;
     int ncdphi4sum=0;
     double cdphi4sum=0;
     
     double Q2x=0;
     double Q2y=0;

     double Q2xsub1=0;
     double Q2ysub1=0;
     double Q2xsub2=0;
     double Q2ysub2=0;
     int nhalf=nparticles/2;
     //Finding reaction plane
     for(int i=0; i<nparticles;i++){
       phitrackseventh->Fill(phitrack[i]);
       Q2x=Q2x+pttrack[i]*cos(nharmonic*phitrack[i]);
       Q2y=Q2y+pttrack[i]*sin(nharmonic*phitrack[i]);
       
       if(i<nhalf){
	 Q2xsub1=Q2xsub1+pttrack[i]*cos(nharmonic*phitrack[i]);
	 Q2ysub1=Q2ysub1+pttrack[i]*sin(nharmonic*phitrack[i]);
       }else if(i<2*nhalf){
	 Q2xsub2=Q2xsub2+pttrack[i]*cos(nharmonic*phitrack[i]);
	 Q2ysub2=Q2ysub2+pttrack[i]*sin(nharmonic*phitrack[i]);
       }else{
       }
     }
     double PSI2found=atan2(Q2y,Q2x)/nharmonic;
     double dpsieventfound=PSI2Event-PSI2found; // CHECK THIS and plot
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
          
     // now that we have the reaction plane calculate v2
     for(int i=0; i<nparticles;i++){
       double dphi=phitrack[i]-PSI2found-PI;  // CHECK THIS
       if(dphi>PI){
	   dphi=TWOPI-dphi;
	     }else if(dphi<-PI){
	   dphi=-TWOPI-dphi;
	     }else{
	 }

       double dphireal=phitrack[i]-PSI2Event-PI;  // CHECK THIS
       if(dphireal>PI){
	   dphireal=TWOPI-dphireal;
	     }else if(dphireal<-PI){
	   dphireal=-TWOPI-dphireal;
	     }else{
	 }

       dphirealh->Fill(dphi);
       dphimeash->Fill(dphireal);

       
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

     }// end loop over particles

     // pairs (cumulants calculations)
     for(int i=0; i<nparticles;i++){
       for(int j=i+1; j<nparticles;j++){
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

	 for(int k=j+1; k<nparticles; k++){

	   for(int l=k+1; l<nparticles; l++){
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

     cdphisumtot+=cdphisum/ncdphisum;
     cdphi3sumtot+=cdphi3sum/ncdphi3sum;
     cdphi4sumtot+=cdphi4sum/ncdphi4sum;

     // finally save the event
     nparticles_previous=nparticles;
     for(int i=0; i<nparticles;i++){
       phitrack_previous[i]=phitrack[i];
       pttrack_previous[i]=pttrack[i];
       ytrack_previous[i]=ytrack[i]; 
     }

   } // end loop over events

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
    cout<<" nevents="<<nevents<<"  nparticlesav="<<nparticlesav<<" nharmonic="<<nharmonic<<endl;
    cout<<endl<<endl<<endl;
    cout<<"Initial values: v1="<<v1
      	<<"  v2="<<v2
      	<<"  v3="<<v3
      	<<"  v4="<<v4
	<<endl<<endl;
    
    cout<<"FINAL values from event plane" << endl << endl;
    cout<<" FINAL event plane resolution="<<res_eventplane<<endl;
    cout<<" FINAL v"<<nharmonic<<" from Event plane="<<v2obs<<" nharmonic="<<nharmonic<<endl;
    
    cout<<" 2 particle correlation fit: v1="<<sqrt(v1sq)
	<<"  v2="<<sqrt(v2sq)
      	<<"  v3="<<sqrt(v3sq)
      	<<"  v4="<<sqrt(v4sq)
	<<"  v5sq="<<v5sq
	<<endl<<endl;

    cout<<"FINAL values from particle cumulants" << endl<<endl;

    double cumulant2 = 1.*cdphisumtot/nevents;
    double cumulant3 = 1.*cdphi3sumtot/nevents;
    double cumulant4 = 1.*cdphi4sumtot/nevents;

    cout<<" 2 particle cumulants = "<<cumulant2<<endl;
    cout<<" 3 particle cumulants = "<<cumulant3<<endl;
    cout<<" 4 particle cumulants = "<<cumulant4<<endl;
    cout<<endl;

    cout<<"FINAL calculations from particle cumulants" << endl<<endl;
    cout<<" from 2 particle cumulant: v"<<nharmonic<<"="<<sqrt(cumulant2)<<endl;
    
    double cumulant4minus2 = 1.*(cumulant4-(2*(cumulant2*cumulant2)));
    cout<<" from 4 particle cumulant: v"<<nharmonic<<"="<<sqrt(sqrt(-cumulant4minus2))<<endl<<endl;
    double cum3=0;
    if(nharmonic==1)cum3=v1*v1*v2;
    if(nharmonic==2)cum3=v2*v2*v4;
    cout<<" nharmonic="<<nharmonic<<", 3 particle cumulant="<<cdphi3sumtot/nevents<<" we expect "<<cum3<<" i.e. v"<<nharmonic<<"*v"<<nharmonic<<"*v"<<2*nharmonic<<endl<<endl;
    //---------------------
    
   // Open a ROOT file and save the formula, function and histogram
   //
   TFile *myfile = new TFile("flow1.root","RECREATE");

   flowf->Write();
   h1f->Write();
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

   
   //   gBenchmark->Show("flow1");
}

