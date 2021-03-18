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
#include "TBenchMark.h"

#include "TMath.h"

#define PI 3.14159
#define TWOPI 6.28318

void tester()
{
  
  int counter=0;
  int sum = 0;
  double temp_av = 0;
  double holder_av = 0;

  int counter_cont=0;
  int sum_cont = 0;
  double temp_av_cont = 0;

  int jk = 0;
  TRandom *ran0 = new TRandom();

  double answer_no_reset = 0;
  double answer_reset = 0;



  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 4; j++){
      for (int k = j +1; k < 4; k++){
	counter++;
	counter_cont++;
	jk = ran0->Poisson(10);
	sum+=jk;
	sum_cont+=jk;

	cout << "for i="<< i << " and j=" << j << " and k=" << k << " then counter=" << counter << " jk=" << jk << " and sum=" << sum << endl;
      }
    }
 
    
    //without reset
    temp_av_cont = sum_cont/counter_cont;
    cout << "temp_av_cont = " << temp_av_cont << endl;
    answer_no_reset+=temp_av_cont;

    //with reset
    temp_av = sum/counter;
    cout << "temp_av = " << temp_av << endl;
    answer_reset+=temp_av;

    //reseting
    counter=0;
    sum=0;   
  }


  cout << "answer_no_reset = " << answer_no_reset/5 << endl;
  cout << "answer_reset = " << answer_reset/5 << endl;

}
