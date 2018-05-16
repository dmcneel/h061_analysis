#define Monitors_mg26_cxx
//test
#include "Monitors_mg26.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>

#define NUMPRINT 20 //>0
ULong64_t NUMSORT=100000000;
ULong64_t NumEntries = 0;
ULong64_t ProcessedEntries = 0;
Float_t Frac = 0.1; //Progress bar
TStopwatch StpWatch;

TH2F* hxfxn[24];
TH2F* heVx[24]; 
TH2F* hecalVxcal[24];
TH2F* htacVxcal[24];
TH2F* htacVxcalg[24];
TH2F* hRecoiltacVxcal[24];
TH2F* hecalVz;
TH2F* hecalVzR;
TH2F* hecalVzR2;
TH2F* hecalVzRcuta;
TH2F* hecalVzRcutb;
TH2F* hecalVzRnotcutb;
TH2F* hrdt[4];
TH2F* helum[2];
TH2F* he0dee;//ezero
TH2F* he0det;
TH2F* he0et;
TH1F* h0detet;
TH1F* h0dettact;
TH1F* h0ettact;
TH1F* h0de;
TH1F* h0e;
TH1F* h0tac;
TH1F* hBIC;
TH2F* htacVelum;
TH1F* helum1d;
TH1I* hrtac[4]; //timing between recoil detectors
TH2F* hrg[4]; //gated recoils
TH1D* projExcuta;
TH1D* projExcutb;
TH1I* hrmult;

TH1F* hExcuta;
TH1F* hExnotcutb;
TH1F* hExcuta_bkg;
TH1F* hExnotcutb_bkg;
TH1F* hExcuta_sub;

TH1F* hExcuta_evz;
TH1F* hExcuta_evz_bkg;
TH1F* hExcuta_evz_sub;

TFile* outfile;
    
TH2I* hRecoilChvArrayCh;
TH2I* hTACArrayRecoil;
TH2I* hTACArrayRecoilg;
TH2F* hexVzRcuta ;
TH2F* hexVzRnotcutb ;
TH2F* hERecoilVex;
TH2F* hDERecoilVex;
TH2F* hKERecoilVex;
TCutG* cut0a;
TCutG* cut0b;
TCutG* cut1a;
TCutG* cut1b;
TCutG* cut2a;
TCutG* cut2b;
TCutG* cut3a;
TCutG* cut3b;

TCutG* open0;
TCutG* open1;
TCutG* open2;
TCutG* open3;


TCutG* tp0;
TCutG* tp1;
TCutG* tp2;
TCutG* tp3;

TCutG* cutb2_0;
TCutG* cutb2_1;
TCutG* cutb2_2;
TCutG* cutb2_3;

TCutG* evzcut;

TH1F* htacE;
//TH1F* hex;
TH1F* hexR;
TH1I* htacArray[24];
TH1F* htac[4];//0 array-rdt, 1 elum-rdt

Float_t x[24],z[24];
Float_t xcal[24],ecal[24],xfcal[24],xncal[24],ecrr[24],ezero[4];
Int_t tacA[24];
Float_t z_array_pos[6] = {35.868,29.987,24.111,18.248,12.412,6.676};//in cm

//Float_t z_off=26.5;//to the physical end of array from the target (-20 to si start)
Float_t z_off=46.5;//to the physical end of array from the target (-20 to si start)
Float_t xcoef[24] = {0.949029,1.01802,0.971783,0.927059,0.991241,0.881107,
		     0.906056,0.926736,0.989037, 1.00974,0.950697, 0,
		     1.03673,1.00772,1.00793,0.970457,0.955211,1.12365,
		     0.965497,0.977822, 0.889425,0.983675,1.08844,1.09201};
Float_t ecoef[24] = {0.895836,0.879174, 1.05728,1.10927,0.894856,0.964865,
		     1.04468,1.00183,0.97847,0.949798,0.927987,0,
		     0.831529,0.927974,0.922803,0.96689,1.04279,0.892727,
		     0.896227,0.958843,1.05903,0.899419,1.10259,0.88497};
/*Float_t kcoef[24] = {0.003444976, 0.003457902, 0.003132563, 0.003341067, 0.003563254, 0.003654822, 
		     0.003366944, 0.003534284, 0.003526171, 0.003309869, 0.003735409, 0.0,
		     0.003791968, 0.003526171, 0.00364557, 0.003384254, 0.003353712, 0.003434705,
		     0.003505249, 0.003378299, 0.003289549, 0.003448586, 0.002994541, 0.003576529};
Float_t bcoef[24] = {0.453923445, 0.417753264, 0.501465942, 0.538881671, 0.449696876, 0.546873096,
		     0.436930294, 0.433740451, 0.448512397, 0.595999138, 0.389035019, 0.0,
		     0.369876234, 0.470022039, 0.533858228, 0.589323149, 0.387687336, 0.52972093,
		     0.494175262, 0.41855132, 0.473385494, 0.584180512, 0.49499818, 0.584009314};*/
Float_t kcoef[24] = {0.00377855,
		     0.004016195,
		     0.003470516,
		     0.003675252,
		     0.004003605,
		     0.00408688,
		     0.003707257,
		     0.003818087,
		     0.003829535,
		     0.003745308,
		     0.004016195,
		     0.0,
		     0.004217104,
		     0.003911639,
		     0.003991094,
		     0.003806706,
		     0.003739824,
		     0.003876024,
		     0.003889303,
		     0.003737089,
		     0.003669971,
		     0.003917638,
		     0.003330248,
		     0.004080351};
Float_t bcoef[24] = {0.096813683,
		     0.065432704,
		     0.075887908,
		     0.061711223,
		     0.067195298,
		     0.07188432,
		     0.049367925,
		     0.051168909,
		     0.04944018,
		     0.054667889,
		     0.069448899,
		     0.0,
		     0.11279158,
		     0.04486585,
		     0.072937969,
		     0.075727571,
		     -0.004330307,
		     -0.06998437,
		     0.069584964,
		     0.059662546,
		     0.029494684,
		     0.047889571,
		     0.038246154,
		     0.048290096};
  
Float_t tempTime=-1000;
Long64_t tempTimeLong=10001;

void Monitors_mg26::Begin(TTree *tree)
{
  TString option = GetOption();
  NumEntries = tree->GetEntries();


  //get cuts
  TFile * cuts = new TFile("cuts_tp_323.root");
  cut0a = (TCutG*)cuts->Get("cut0a");
  cut0b = (TCutG*)cuts->Get("cut0b");
  cut1a = (TCutG*)cuts->Get("cut1a");
  cut1b = (TCutG*)cuts->Get("cut1b");
  cut2a = (TCutG*)cuts->Get("cut2a");
  cut2b = (TCutG*)cuts->Get("cut2b");
  cut3a = (TCutG*)cuts->Get("cut3a");
  cut3b = (TCutG*)cuts->Get("cut3b");

  TFile * cuts2 = new TFile("cuts325.root");
  tp0=(TCutG*)cuts2->Get("tp0");
  tp1=(TCutG*)cuts2->Get("tp1");
  tp2=(TCutG*)cuts2->Get("tp2");
  tp3=(TCutG*)cuts2->Get("tp3");

  TFile * cuts3 = new TFile("evzcut_gspos.root");
  evzcut=(TCutG*)cuts3->Get("evzcut");

  TFile * cuts4 = new TFile("tp_cuts_326.root");
  open0 = (TCutG*)cuts4->Get("open0");
  open1 = (TCutG*)cuts4->Get("open1");
  open2 = (TCutG*)cuts4->Get("open2");
  open3 = (TCutG*)cuts4->Get("open3");

  TFile * cuts5= new TFile("newcutb_326.root");
  cutb2_0=(TCutG*)cuts5->Get("newcutb0");
  cutb2_1=(TCutG*)cuts5->Get("newcutb1");
  cutb2_2=(TCutG*)cuts5->Get("newcutb2");
  cutb2_3=(TCutG*)cuts5->Get("newcutb3");
  
  outfile = new TFile("position2.root","RECREATE");

  
  //Generate all of the histograms needed for drawing later on

  for (Int_t i=0;i<24;i++) {//array loop
    hxfxn[i] = new TH2F(Form("hxfxn%d",i),
			Form("Raw PSD XF vs. XN (ch=%d);XF (channel);XN (channel)",i),
			500,0,4000,500,0,4000);
    heVx[i] = new TH2F(Form("heVx%d",i),
		       Form("Raw PSD E vs. X (ch=%d);X (channel);E (channel)",i),
		       500,-0.1,1.1,500,0,4000);
    hecalVxcal[i] = new TH2F(Form("hecalVxcal%d",i),
			     Form("Cal PSD E vs. X (ch=%d);X (channel);E (channel)",i),
			     500,-0.1,1.1,500,0,12);
    htacVxcal[i] = new TH2F(Form("htacVxcal%d",i),
			    Form("Cal tac vs. X (ch=%d);X (channel);tac (channel)",i),
			    200,-0.1,1.1,200,2000,3500);
    htacVxcalg[i] = new TH2F(Form("htacVxcalg%d",i),
			    Form("Cal tac vs. X gated on array-rf(ch=%d);X (channel);tac (channel)",i),
			    200,-0.1,1.1,200,2000,3500);
    hRecoiltacVxcal[i] = new TH2F(Form("hRecoiltacVxcal%d",i),
			    Form("Cal tac vs. X (ch=%d);X (channel);tac (channel)",i),
			    200,-0.1,1.1,400,-200,200);
    
  }//array loop
  hecalVz = new TH2F("hecalVz","E vs. Z;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hecalVzR = new TH2F("hecalVzR","E vs. Z gated;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hecalVzR2 = new TH2F("hecalVzR2","E vs. Z gated;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hecalVzRcuta = new TH2F("hecalVzcuta","E vs. Z gated;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hecalVzRcutb = new TH2F("hecalVzcutb","E vs. Z gated;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hecalVzRnotcutb = new TH2F("hecalVznotcutb","E vs. Z gated;Z (cm);E (MeV)",450,-95,-10,600,0,12);
  hexVzRcuta = new TH2F("hexVzRcuta","Ex vs. Z gated;Z (cm);Ex (MeV)",450,-95,-10,700,-2,12);
  hexVzRnotcutb = new TH2F("hexVzRnotcutb","Ex vs. Z gated;Z (cm);Ex (MeV)",450,-75,-10,700,-2,12);
  hERecoilVex= new TH2F("hERecoilVex","Ex vs. RE gated;Ex (MeV);E (MeV)",700,-2,12,2000,0,6000);
  hDERecoilVex= new TH2F("hDERecoilVex","Ex vs. RDE gated;Ex (MeV);E (MeV)",700,-2,12,2000,0,6000);
  hKERecoilVex= new TH2F("hKERecoilVex","Ex vs. DE+E;Ex (MeV);E (MeV)",500,-2,12,500,0,300);

  //Recoils
  for (Int_t i=0;i<4;i++) {
    hrdt[i] = new TH2F(Form("hrdt%d",i),
		       Form("Raw Recoil DE vs Eres (ch=%d); Eres (channel); DE (channel)",i),
		       500,0,8000,500,0,8000);
    hrg[i] = new TH2F(Form("hrg%d",i),
		       Form("Gated Recoil DE vs Eres (ch=%d); Eres (channel); DE (channel)",i),
		       500,0,8000,500,0,8000);
      hrtac[i] = new TH1I(Form("hrtac%d",i),
		       Form("TAC between recoil pairs (ch=%d); t-diff (channel);",i),
		       500,-1000,1000);
    
  }
    hrmult=new TH1I("hrmult","recoil multiplicity",20,0,20);

  //ELUM
  helum[0] = new TH2F("helum0","Elum Ring Energies; E (channels); Ring Number",
		      500,100,4000,16,0,16);
  helum[1] = new TH2F("helum1","Elum Wedge Energies; E (channels); Ring Number",
		      500,100,4000,16,0,16);
  htacVelum = new TH2F("htacVelum","ELUM vs. TAC; TAC (channels); E (channels)",
		       500,2000,3500,500,50,7550);
  hBIC = new TH1F("hBIC","BIC counter; BIC; Counts",
		  10,0,10);
  helum1d = new TH1F("helum1d","Elum;",
		  500,10,8000);

  hRecoilChvArrayCh = new TH2I("recoilchvarraych","recoilchvarraych",24,0,24,8,0,8);
  hTACArrayRecoil = new TH2I("tacarrayrecoil","tacarrayrecoil",24,0,24,600,-300,300);
  hTACArrayRecoilg = new TH2I("tacarrayrecoilg","tacarrayrecoil gated evz",24,0,24,600,-300,300);
  
  //TAC
  htac[0] = new TH1F("htac0","Array-RDT0 TAC; DT [clock ticks]; Counts",6,0,6);
  htac[1] = new TH1F("htac1","Array-RDT1 TAC; DT [clock ticks]; Counts",6,0,6);
  htac[2] = new TH1F("htac2","Array-RDT2 TAC; DT [clock ticks]; Counts",6,0,6);
  htac[3] = new TH1F("htac3","Array-RDT3 TAC; DT [clock ticks]; Counts",6,0,6);

  htacE = new TH1F("htacE","Elum-RDT TAC; DT [clock ticks]; Counts",4,0,4);

  // hex = new TH1F("hex","excitation spectrum",200,-2,10);
  hexR = new TH1F("hexR","excitation spectrum with Recoil",500,-2,10);
  hExcuta = new TH1F("hexcuta","excitation spectrum with Recoil",150,-2,10);
  hExcuta_bkg = new TH1F("hexcuta_bkg","excitation spectrum with Recoil",150,-2,10);
  hExcuta_sub = new TH1F("hexcuta_sub","excitation spectrum with Recoil",150,-2,10);

  hExcuta_evz = new TH1F("hexcuta_evz","excitation spectrum with Recoil",150,-2,10);
  hExcuta_evz_bkg = new TH1F("hexcuta_evz_bkg","excitation spectrum with Recoil",150,-2,10);
  hExcuta_evz_sub = new TH1F("hexcuta_evz_sub","excitation spectrum with Recoil",150,-2,10);

  
  hExnotcutb = new TH1F("hexnotcutb","excitation spectrum with Recoil",150,-2,10);
  hExnotcutb_bkg = new TH1F("hexnotcutb_bkg","excitation spectrum with Recoil",150,-2,10);

  for (Int_t i=0;i<24;i++) {
    htacArray[i] = new TH1I(Form("htacArray%d",i), Form("Array-RDT TAC for ch%d",i), 600, -300,300);
  }

  //EZERO

  he0dee = new TH2F("he0dee","EZERO DE-E; E [ch]; DE [ch]",500,0,8000,500,0,8000);//ezero
  he0det = new TH2F("he0det","EZERO DE-RF; RF [ch]; DE [ch]",500,2000,3500,500,0,8000);//
  he0et = new TH2F("he0et","EZERO E-RF; RF [ch]; DE [ch]",500,2000,3500,500,0,8000);//
  h0detet = new TH1F("h0detet","EZERO DE Time - E Time; DET-ET [ch]",500,-250,250);//
  h0dettact = new TH1F("h0dettact","EZERO DE Time - TAC Time; DET-ET [ch]",2000,-1000,1000);//
  h0ettact = new TH1F("h0ettact","EZERO E Time - TAC Time; DET-ET [ch]",2000,-1000,1000);//
  h0de = new TH1F("h0de","EZERO DE ; DE [ch]",500,50,4050);//
  h0e = new TH1F("h0e","EZERO - E; E [ch]",500,50,4050);//
  h0tac = new TH1F("h0tac","EZERO RF; RF [ch]",500,50,4050);//
  StpWatch.Start();
}

void Monitors_mg26::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

}

Bool_t Monitors_mg26::Process(Long64_t entry)
{
  ProcessedEntries++;
  if (ProcessedEntries<NUMSORT) {
    
    if (ProcessedEntries>NumEntries*Frac-1) {
      printf(" %3.0f%% (%llu/%llu k) processed in %6.1f seconds\n",
	     Frac*100,ProcessedEntries/1000,NumEntries/1000,StpWatch.RealTime());
      StpWatch.Start(kFALSE);
      Frac+=0.1;
    }

    b_Energy->GetEntry(entry);
    b_XF->GetEntry(entry);
    b_XN->GetEntry(entry);
    b_RDT->GetEntry(entry);
    b_TAC->GetEntry(entry);
    b_ELUM->GetEntry(entry);
    b_EZERO->GetEntry(entry);
    b_EnergyTimestamp->GetEntry(entry);
    b_RDTTimestamp->GetEntry(entry);
    b_TACTimestamp->GetEntry(entry);
    b_ELUMTimestamp->GetEntry(entry);
    b_EZEROTimestamp->GetEntry(entry);
    b_RDTTimestampCorrection->GetEntry(entry);
    b_EnergyTimestampCorrection->GetEntry(entry);

    //Do calculations and fill histograms
    //Array calcs first
    for (Int_t i=0;i<24;i++) {
      //Calibrations go here
      xfcal[i] = xf[i];
      xncal[i] = xn[i]*xcoef[i];
      ecal[i] = e[i]/ecoef[i];
      ecrr[i] = e[i]*kcoef[i]+bcoef[i];

      if (xf[i]>0 || xn[i]>0) {
	x[i] = 0.5*((xf[i]-xn[i]) / (xf[i]+xn[i]))+0.5;//for downstream?
	//x[i] = 0.5*((xf[i]-xn[i]) / (xf[i]+xn[i]))+0.5;

	if (xfcal[i]>xncal[i]) {
	  xcal[i] = xfcal[i]/ecal[i];
	} else if (xncal[i]>=xfcal[i]) {
	  xcal[i] = 1.0 - xncal[i]/ecal[i];
	}
	//	z[i] = 5.0*(xcal[i]-0.5) + z_off + z_array_pos[i%6];//for downstream?
	z[i] = 5.0*(xcal[i]-0.5) - z_off - z_array_pos[i%6];
      }

      //Array fill next
      hxfxn[i]->Fill(xf[i],xn[i]);
      
      if (e[i]>50&&(xn[i]>0||xf[i]>0)) {
	heVx[i]->Fill(x[i],e[i]);
	hecalVxcal[i]->Fill(xcal[i],ecrr[i]);
	hecalVz->Fill(z[i],ecrr[i]);
	htacVxcal[i]->Fill(xcal[i],tac[0]);
      }

    }//array loop

    /* ELUM */
    if (elum[0]!=0) {
      // htacVelum->Fill(tac[0],TMath::Abs(elum[0]));
      //helum1->Fill(TMath::Abs(elum[0]));
    }
    for (Int_t i=0;i<16;i++) {
      helum[0]->Fill(elum[i],i);
      helum[1]->Fill(elum[i+16],i);
      for(Int_t j=0;j<4;j++){
	if(rdt[j]-elum_t[i]>-10&&rdt[j]-elum[i]<10)htacE->Fill(0.5+j);
      }
    }

    //EZERO
   
    he0dee->Fill(ezero[1],ezero[0]);
    he0det->Fill(TMath::Abs(tac[0]),ezero[0]);
    he0et->Fill(TMath::Abs(tac[0]),ezero[1]);
    h0detet->Fill((int)(ezero_t[0]-ezero_t[1]));
    h0dettact->Fill((int)(ezero_t[0]-tac_t[0]));
    h0ettact->Fill((int)(ezero_t[1]-tac_t[0]));
    h0de->Fill(ezero[0]);
    h0e->Fill(ezero[1]);
    h0tac->Fill(TMath::Abs(tac[0]));

    //TACs
    //if (tac[1]>0) {
    //  for (Int_t mm=0;mm<tac[1];mm++)
    //	hBIC->Fill(5);
    // }

    for(int i=0;i<24;i++){
      if(e[i]>50&&(xn[i]>0||xf[i]>0)){
	for(int j=0;j<8;j++){
	  if(rdt[j]>250){
	    hRecoilChvArrayCh->Fill(i,j);
	    int relt = (int)(rdt_t[j]-e_t[i]);
	    hTACArrayRecoil->Fill(i,(int)(rdt_t[j]-e_t[i]));
	  }
	}
      }
      
    }


    double ex = 0;


    
    for (Int_t i=0;i<4;i++){
      for(Int_t j=0;j<6;j++){
	if(i==0&&e[i*6+j]>100){
	  tacA[i*6+j]= (int)(rdt_t[0]-e_t[i*6+j]);
	  htacArray[i*6+j]->Fill(tacA[i*6+j]);
	  if(tacA[i*6+j]>-10&&tacA[i*6+j]<10&&rdt[0]>1000){hecalVzR->Fill(z[i*6+j],ecrr[i*6+j]);
	    ex = (z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137;
	    hexR->Fill((z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137);
	    htac[0]->Fill(j+0.5);}
	}
	  if(i==1&&e[i*6+j]>100){
	    tacA[i*6+j]= (int)(rdt_t[6]-e_t[i*6+j]);
	  htacArray[i*6+j]->Fill(tacA[i*6+j]);
	  if(tacA[i*6+j]>-10&&tacA[i*6+j]<10&&rdt[6]>1000){hecalVzR->Fill(z[i*6+j],ecrr[i*6+j]);
	    ex = (z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137;
	    hexR->Fill((z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137);
	    htac[3]->Fill(j+0.5);}
	  }
	  if(i==2&&e[i*6+j]>100){
	    	  tacA[i*6+j]= (int)(rdt_t[4]-e_t[i*6+j]);
	  htacArray[i*6+j]->Fill(tacA[i*6+j]);
	  if(tacA[i*6+j]>-10&&tacA[i*6+j]<10&&rdt[4]>1000){hecalVzR->Fill(z[i*6+j],ecrr[i*6+j]);
	    ex = (z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137;
	    hexR->Fill((z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137);
	  htac[2]->Fill(j+0.5);}
	  }
	  if(i==3&&e[i*6+j]>100){
	    tacA[i*6+j]= (int)(rdt_t[2]-e_t[i*6+j]);
	  htacArray[i*6+j]->Fill(tacA[i*6+j]);
	  if(tacA[i*6+j]>-10&&tacA[i*6+j]<10&&rdt[2]>1000){hecalVzR->Fill(z[i*6+j],ecrr[i*6+j]);
	    ex = (z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137;
	    hexR->Fill((z[i*6+j]*0.14028+13.638-ecrr[i*6+j])*1.1068-0.3137);
	    htac[1]->Fill(j+0.5);}
	  }
      }
    }
    
    int cutaflag = 0;
    int cutbflag = 0;
    int evzcutflag = 0;
    /*
    if(cut0a && cut0a->IsInside(rdt[1],rdt[0])){cutaflag = 1;}
    
    if(cut1a && cut1a->IsInside(rdt[3],rdt[2])){cutaflag = 1;}

    if(cut2a && cut2a->IsInside(rdt[5],rdt[4])){cutaflag = 1;}

    if(cut3a && cut3a->IsInside(rdt[7],rdt[6])){cutaflag = 1;}
    
    
    if(cut0b && cut0b->IsInside(rdt[1],rdt[0])){cutbflag = 1;}
    if(cut1b && cut1b->IsInside(rdt[3],rdt[2])){cutbflag = 1;}
    if(cut2b && cut2b->IsInside(rdt[5],rdt[4])){cutbflag = 1;}
    if(cut3b && cut3b->IsInside(rdt[7],rdt[6])){cutbflag = 1;}
    
     
    if(tp0 && tp0->IsInside(rdt[1],rdt[0])){cutaflag = 1;}
    if(tp1 && tp1->IsInside(rdt[3],rdt[2])){cutaflag = 1;}
    if(tp2 && tp2->IsInside(rdt[5],rdt[4])){cutaflag = 1;}
    if(tp3 && tp3->IsInside(rdt[7],rdt[6])){cutaflag = 1;}
    */
    if(cutb2_0 && cutb2_0->IsInside(rdt[1],rdt[0])){cutbflag = 1;}
    if(cutb2_1 && cutb2_1->IsInside(rdt[3],rdt[2])){cutbflag = 1;}
    if(cutb2_2 && cutb2_2->IsInside(rdt[5],rdt[4])){cutbflag = 1;}
    if(cutb2_3 && cutb2_3->IsInside(rdt[7],rdt[6])){cutbflag = 1;}

    if(open0 && open0->IsInside(rdt[1],rdt[0])){cutaflag = 1;}
    if(open1 && open1->IsInside(rdt[3],rdt[2])){cutaflag = 1;}
    if(open2 && open2->IsInside(rdt[5],rdt[4])){cutaflag = 1;}
    if(open3 && open3->IsInside(rdt[7],rdt[6])){cutaflag = 1;}
    
    double ex_st = -888;
    Int_t rmult=0;
    for (Int_t i=0;i<4;i++){

      
      Float_t de_r=rdt[2*i];
      Float_t e_r=rdt[2*i+1];

      
      if(de_r>0&&e_r>0) rmult=rmult+1;
    }
    /* RECOILS */
    for (Int_t i=0;i<4;i++){
      Float_t de_r=rdt[2*i];
    
      Float_t e_r=rdt[2*i+1];
      Float_t degain[4]={1,2.166,1.095,2.148};
      Float_t egain[4]={2.062,1.185,2.101,1};
      Float_t dec=de_r*degain[i];
      Float_t ec=e_r*egain[i];
      Float_t decal=dec*0.0194;
      Float_t ecal=ec*0.0182;

 
   


	  if(de_r>1000&&e_r>200&&rmult==1){

	hrdt[i]->Fill(e_r,de_r);
	hrtac[i]->Fill(rdt_t[2*i+1]-rdt_t[2*i]);
	for(Int_t j=0; j<24;j++){
	  int time_rel=(int)(e_t[j]-rdt_t[i*2]);
	  int weight=0;
	  hRecoiltacVxcal[j]->Fill(x[j],time_rel); 
	  if(e[j]>100){
	    //ex_st = z[j]*0.13778 +9.99862 - ecrr[j]; //fit straighten on d,p ground state
	    //ex_st = z[j]*0.1363 +10.04 - ecrr[j];//heliomatic
	    //ex_st = z[j]*0.12744 +17.33 - ecrr[j];//heliomatic for t,p
	    ex_st = z[j]*0.128823 +17.29 - ecrr[j];//scaling t,p heliomatic by ratio of d,p fit:heliomatic
	    //ex_st = z[j]*0.1297757 + 17.47 - ecrr[j];//scaling t,p heliomatic by ratio of d,p fit:heliomatic
	    if(time_rel>-10&&time_rel<10) weight=1;
	    if(time_rel<100&&time_rel>80) weight=-1;
	    htacVxcalg[j]->Fill(xcal[j],tac[0],weight);
	    hrg[i]->Fill(ec,dec,weight);
	    hecalVzR2->Fill(z[j],ecrr[j],weight);
	    hTACArrayRecoilg->Fill(j,time_rel);
	    if(cutaflag==1&&weight==1){
	     
	      hecalVzRcuta->Fill(z[j],ecrr[j],weight);
	      hexVzRcuta->Fill(z[j],ex_st,weight);
	      hERecoilVex->Fill(ex_st,ec,weight); // these arent currently gainmatched so not useful
	      hDERecoilVex->Fill(ex_st,dec,weight);
	      hKERecoilVex->Fill(ex_st,decal+ecal,weight);
	      if(time_rel>-10&&time_rel<10){
		hExcuta->Fill(ex_st);
		if(evzcut && evzcut->IsInside(z[j],ecrr[j])){
		  hExcuta_evz->Fill(ex_st);
	
		}
	      }
	      if(time_rel>80&&time_rel<100){
		hExcuta_bkg->Fill(ex_st);
		if(evzcut && evzcut->IsInside(z[j],ecrr[j])){
		  hExcuta_evz_bkg->Fill(ex_st);
		}
	      }
	     
	      
	    }
	    if(cutbflag==1){hecalVzRcutb->Fill(z[j],ecrr[j],weight);}
	    else if(cutbflag==0){
	      if(time_rel>-10&&time_rel<10){
		hExnotcutb->Fill(ex_st);

	     }
	      if(time_rel>80&&time_rel<100){
		hExnotcutb_bkg->Fill(ex_st);
	      }
	      
	      hexVzRnotcutb->Fill(z[j],ex_st,weight);
	      hecalVzRnotcutb->Fill(z[j],ecrr[j],weight);}
	    
	  }
	}
      }
      
    }
    hrmult->Fill(rmult);
  }
  return kTRUE;
}

void Monitors_mg26::SlaveTerminate()
{

}

void Monitors_mg26::Terminate()
{
   
  TCanvas *cxfxn = new TCanvas("cxfxn","XFXN",1200,800);
  cxfxn->Clear(); cxfxn->Divide(6,4);
  TCanvas *ceVx = new TCanvas("ceVx","EVX",1200,800);
  ceVx->Clear(); ceVx->Divide(6,4);
  TCanvas *cecalVxcal = new TCanvas("cecalVxcal","ECALVXCAL",1200,800);
  cecalVxcal->Clear(); cecalVxcal->Divide(6,4);
  TCanvas *ctac = new TCanvas("ctac","TAC",1200,800);
  ctac->Clear(); ctac->Divide(6,4);
  TCanvas *crecoiltac = new TCanvas("crecoiltac","recoilTAC",1200,800);
  crecoiltac->Clear(); crecoiltac->Divide(6,4);

  for (Int_t i=0;i<24;i++) {
    cxfxn->cd(i+1); hxfxn[i]->Draw("col");
    ceVx->cd(i+1); heVx[i]->Draw("col");
    cecalVxcal->cd(i+1); hecalVxcal[i]->Draw("col");
    htacVxcalg[i]->SetMinimum(0);
    ctac->cd(i+1); htacVxcalg[i]->Draw("col");
    crecoiltac->cd(i+1); hRecoiltacVxcal[i]->Draw("col");
   
  }

    TCanvas *celum = new TCanvas("celum","ELUM",1200,800);
  celum->Clear(); //celum->Divide(1,3);
  helum[0]->Draw("col");
  
  


  TCanvas *crdt = new TCanvas("crdt","RDT",1000,1000);
  crdt->Clear();crdt->Divide(2,2);
  for (Int_t i=0;i<4;i++) {
    crdt->cd(i+1); hrdt[i]->Draw("col");
  }
  TCanvas *cecalVz = new TCanvas("cevalVz","ECALVZ",1000,650);
  cecalVz->Clear();cecalVz->Divide(2,2);
  cecalVz->cd(1);hecalVz->Draw("col");
  cecalVz->cd(2);hecalVzR->Draw("col");
  cecalVz->cd(3);hRecoilChvArrayCh->Draw("col");
  cecalVz->cd(4);hTACArrayRecoil->Draw("col");
  
  /*
    projEx = 
  */


  TCanvas *cex_sub = new TCanvas("cex_sub","cexsub",1800,1200);
  cex_sub->Clear();cex_sub->Divide(2,2);
  cex_sub->cd(1);
  hExcuta_sub->Add(hExcuta);
  hExcuta_sub->Add(hExcuta_bkg,-1);
  hExcuta_sub->Draw();
  cex_sub->cd(2);
  hExcuta_evz_sub->Add(hExcuta_evz);
  hExcuta_evz_sub->Add(hExcuta_evz_bkg,-1);
  hExcuta_evz_sub->Draw();
  cex_sub->cd(3);
  hecalVzRcuta->SetMinimum(0);
  hecalVzRcuta->Draw("col");
  cex_sub->cd(4);
  hrg[3]->SetMinimum(0);
  hrg[3]->Draw("col");

  
  /*
  TCanvas *cEB = new TCanvas("cEB","Elum and BiC",1000,800);
  cEB->Clear(); cEB->Divide(1,3);
  cEB->cd(1);hBIC->Draw();
  cEB->cd(2);helum1->Draw();
  cEB->cd(3);htacVelum->Draw("col");
 
  TCanvas *cecalVzR = new TCanvas("cevalVzR","ECALVZ Gated",1000,650);
  cecalVzR->Clear();hecalVzR->Draw("col");
  
  TCanvas *ctac = new TCanvas("ctac","TAC",1200,800);
  ctac->Clear(); ctac->Divide(6,4);
  for (Int_t i=0;i<24;i++) {
    ctac->cd(i+1); htacArray[i]->Draw("");
  }
 
  TCanvas *ctacA = new TCanvas("ctacA","TACA",1200,800);
  ctacA->Clear(); ctacA->Divide(2,2);
  for (Int_t i=0;i<4;i++) {
    ctacA->cd(i+1); htac[i]->Draw("");
  }
 
 TCanvas *cex = new TCanvas("cex","Excitation Spectrum",1200,800);
  cex->Clear(); cex->Divide(1,2);
  //  cex->cd(1); hex->Draw("");
  cex->cd(2); hexR->Draw("");
 
  TCanvas *celum = new TCanvas("celum","ELUM",1200,800);
  celum->Clear(); celum->Divide(1,3);
  celum->cd(1); helum[0]->Draw("col");
  celum->cd(2); helum[1]->Draw("col");
  celum->cd(3); htacE->Draw("");

 
 TCanvas *ce0 = new TCanvas("ce0","EZERO - Energy - RF",1200,800);
 ce0->Clear(); ce0->Divide(2,2);
 ce0->cd(1); he0dee->Draw("col");
 ce0->cd(2); he0det->Draw("col");
 ce0->cd(3); he0et->Draw("col");
 ce0->Modified(); ce0->Update();

 TCanvas *ce00 = new TCanvas("ce00","EZERO - DeltaTs",1200,800);
 ce00->Clear(); ce00->Divide(2,2);
 ce00->cd(1); h0de->Draw("");
 ce00->cd(2); h0e->Draw("");
 ce00->cd(3); h0tac->Draw("");
 ce00->Modified(); ce00->Update();
 // ce00->WaitPrimitive();
 */
  if (ProcessedEntries>=NUMSORT)
    printf("Sorted only %llu\n",NUMSORT);
  printf("Total time for sort: %3.1f\n",StpWatch.RealTime());
  printf("Which is a rate of: %3.1f k/s\n",(Float_t)ProcessedEntries/StpWatch.RealTime()/1000.);
  StpWatch.Start(kFALSE);
  outfile->Write();
}
