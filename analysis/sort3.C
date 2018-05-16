#define sort3_cxx
// The class definition in sort3.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("sort3.C")
// Root > T->Process("sort3.C","some options")
// Root > T->Process("sort3.C+")
//

#include "sort3.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TStopwatch.h>



TStopwatch StpWatch;
Float_t Frac = 0.1; //Progress bar
ULong64_t ProcessedEntries = 0;
ULong64_t SkippedEntries = 0;
ULong64_t FilledEntries=0;
ULong64_t NumEntries = 0;
#define maxhits 10
#define M 100
Int_t idDetMap[160] = {300,-1,-1,-1,-1,-1,-1,-1,-1,-1,//RF TIMING STUFF
		       200,201,202,203,204,205,206,207,-1,-1,//ELUM
		       208,209,210,211,212,213,214,215,-1,-1,//ELUM
		       100,100,101,101,102,102,103,103,-1,-1,//Recoil Pairs have same detid
		       1,0,5,4,3,2,1,0,-1,-1,/*1*/
		       3,2,1,0,5,4,3,2,-1,-1,/*2*/
		       11,10,9,8,7,6,5,4,-1,-1,/*3*/
		       -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,//Empty
		       7,6,11,10,9,8,7,6,-1,-1,/*4*/
		       15,14,13,12,11,10,9,8,-1,-1,/*5*/
		       17,16,15,14,13,12,17,16,-1,-1,/*6*/
		       -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,//Empty
		       19,18,17,16,15,14,13,12,-1,-1,/*7*/
		       21,20,19,18,23,22,21,20,-1,-1,/*8*/
		       23,22,21,20,19,18,23,22,-1,-1,/*9*/      
		       -1,-2,-3,-4,-5,-6,-7,-8,-9,-10};///

Int_t idKindMap[160] = {6,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			5,5,5,5,5,5,5,5,-1,-1,
			5,5,5,5,5,5,5,5,-1,-1,
			3,4,3,4,3,4,3,4,-1,-1,//3=DE,4=E
			1,1,0,0,0,0,0,0,-1,-1,//1 0==energy 1==xf 2==xn
			2,2,2,2,1,1,1,1,-1,-1,//2
			0,0,0,0,0,0,2,2,-1,-1,//3
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			2,2,1,1,1,1,1,1,-1,-1,//4
			0,0,0,0,2,2,2,2,-1,-1,//5
			2,2,2,2,2,2,0,0,-1,-1,//6
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
			0,0,1,1,1,1,1,1,-1,-1,//7
			1,1,1,1,0,0,0,0,-1,-1,//8
			2,2,2,2,2,2,1,1,-1,-1,//9
			-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

Int_t fitrangeA[24]={120,140,120,140,140,140,
		     80,120,100,120,120,200,
		     80,80,80,120,120,120,
		     100,120,120,120,80,80};

Int_t fitrangeR[4]={80,80,100,80};
Bool_t traceprocessing=1;
TFile *oFile;
TTree *gen_tree;
TCanvas *cc1;
TCanvas *cc2;
TCanvas *cc3;
TCanvas *cta;
TCanvas *ct0;
TCanvas *ct1;
TCanvas *ct2;
TCanvas *ct3;
TCanvas *ct4;
TCanvas *ct5;
TCanvas *ct6;
TF1 *f1;
TH2F *hrise0;
TH2F *hrise3;
TH2I *hmult1;
TH2I *hmult2;
TH2I *hmult3;
TH2I *hmult4;
TH1I *evttype;
TH1I *hfill;
TH1F *htrace[6];
TH2F *htrace0[24];
// TH2I *htrace2[24];
TH2F *htrace3[4];
TH2F *htrace4[4];
// TH2I *htrace5;
// TH2I *htrace6;

//Outfile struct
typedef struct {
  Int_t EventType;
  Float_t Energy[24];
  Float_t XF[24];
  Float_t XN[24];
  Float_t RDT[8];
  Float_t TAC[10];
  Float_t ELUM[16];
  Float_t EZERO[4];//0,1 - DE0,E0  

  ULong64_t EnergyTimestamp[24];
  Float_t ETC[24];
  ULong64_t XFTimestamp[24];
  ULong64_t XNTimestamp[24];
  ULong64_t RDTTimestamp[8];
  Float_t RTC[4];
  ULong64_t TACTimestamp[10];
  ULong64_t ELUMTimestamp[16];
  ULong64_t EZEROTimestamp[4];
   
   
} OUTPUT;

typedef struct {

  Double_t trace[203];
  ULong64_t counter;

} TRACE;

TRACE average[7];
TRACE current[7];

OUTPUT psd;
void sort3::Begin(TTree * tree)
{
  StpWatch.Start();
  NumEntries = tree->GetEntries();
  f1=new TF1("f1","pol4",50,100);
  evttype=new TH1I("evttype","type of event",11,-1,10);
  hfill=new TH1I("hfill","detector type",11,-1,10);
  hmult1=new TH2I("hmult1","array mult vs recoil mult",100,0,100,100,0,100);
  hmult2=new TH2I("hmult2","array mult vs tacmult",100,0,100,100,0,100);
  hmult3=new TH2I("hmult3","elum vs recoil mult",100,0,100,100,0,100);
  hrise0=new TH2F("hrise0","rise time vs channel id",24,0,24,300,500,800);
  hrise3=new TH2F("hrise3","rise time vs channel id",4,0,4,300,500,800);
  for(int i=0;i<7;i++){
    
    htrace[i]=new TH1F(Form("htrace%d",i),Form("trace for det kind %d",i),203,0,203);
  }


  for(int i=0;i<24;i++){
    
    htrace0[i]=new TH2F(Form("htrace0%d",i),Form("average trace for det kind 0 det id %d",i),203,0,203,1600,0,32000);
  }
  for(int i=0;i<4;i++){
    
    htrace3[i]=new TH2F(Form("htrace3%d",i),Form("average trace for det kind 3 det id %d",i),203,0,203,1600,0,32000);
  }
  for(int i=0;i<4;i++){
    
    htrace4[i]=new TH2F(Form("htrace4%d",i),Form("average trace for det kind 4 det id %d",i),203,0,203,1600,0,32000);
  }


  for(int i=0;i<6;i++){
    for(int j=0; j<203;j++){
      average[i].trace[j]=0;
      if(i==0) average[i].counter=0;


    }
  }
 

  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
 
  oFile = new TFile("gen.root","RECREATE");
  gen_tree = new TTree("gen_tree","PSD Tree");

  gen_tree->Branch("e",psd.Energy,"Energy[24]/F");
  gen_tree->Branch("e_t",psd.EnergyTimestamp,"EnergyTimestamp[24]/l");
  gen_tree->Branch("e_tc",psd.ETC,"EnergyTimestampCorrection[24]/F");

  gen_tree->Branch("xf",psd.XF,"XF[24]/F");
  gen_tree->Branch("xf_t",psd.XFTimestamp,"XFTimestamp[24]/l");
 
 
  gen_tree->Branch("xn",psd.XN,"XN[24]/F");
  gen_tree->Branch("xn_t",psd.XNTimestamp,"XNTimestamp[24]/l"); 
 

  gen_tree->Branch("rdt",psd.RDT,"RDT[8]/F");
  gen_tree->Branch("rdt_t",psd.RDTTimestamp,"RDTTimestamp[8]/l"); 
  gen_tree->Branch("rdt_tc",psd.RTC,"RDTTimestampCorrection[24]/F");
		   

  gen_tree->Branch("tac",psd.TAC,"TAC[10]/F");
  gen_tree->Branch("tac_t",psd.TACTimestamp,"TACTimestamp/l"); 
  
  gen_tree->Branch("elum",psd.ELUM,"ELUM[16]/F");
  gen_tree->Branch("elum_t",psd.ELUMTimestamp,"ELUMTimestamp[16]/l");
 
  gen_tree->Branch("ezero",psd.EZERO,"EZERO[4]/F");
  gen_tree->Branch("ezero_t",psd.EZEROTimestamp,"EZEROTimestamp[4]/l"); 
 
  TString option = GetOption();
}

void sort3::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

Bool_t sort3::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either sort3::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  //zero struct
  for (Int_t i=0;i<24;i++) {//num dets
     
    psd.Energy[i]=0;
    psd.XF[i]=0;
    psd.XN[i]=0;
       
    psd.EnergyTimestamp[i]=0;
    psd.ETC[i]=0;
    psd.XFTimestamp[i]=0;
    psd.XNTimestamp[i]=0;

    

    if (i<16){
      psd.ELUM[i]=0;
      psd.ELUMTimestamp[i]=0;

    }
    if (i<8){
      psd.RDT[i]=0;
      psd.RDTTimestamp[i]=0;
    }
    if(i<4){
      psd.RTC[i]=0;
    }
    if (i==0){
      psd.TAC[0]=0;
      psd.EventType=0;
    }
  }
  for(int i=0;i<6;i++){
    for(int j=0;j<203;j++){
      current[i].trace[j]=0;
    }
  }
  //Begin doing stuff here.
  //Start with the multiplicity check

 
  b_NumHits->GetEntry(entry);
  
  Bool_t goodN=0;
  Int_t event_type=-1;

  if(NumHits<maxhits) goodN=1;
  // std::cout<<NumHits<<endl;
 
  if(goodN){
    ProcessedEntries++;
    //since there's no way to break the Process loop, I guess I will build everything that's good in this loop.
    //Next, find out WHAT the hits were. If they are good, proceed.
    b_id->GetEntry(entry);
    b_pre_rise_energy->GetEntry(entry);
    b_post_rise_energy->GetEntry(entry);
    b_event_timestamp->GetEntry(entry);
    b_sync_error_flag->GetEntry(entry);
    b_cfd_valid_flag->GetEntry(entry);
    b_trace_length->GetEntry(entry);
    b_m2_end_sample->GetEntry(entry);



    /////////////////////////////IDENTIFICATION BLOCK/////////////////////////

    Int_t idConst=1010;
    Int_t channum[maxhits]; //channum goes from 0 to 159 for the 16 digitizers
    Int_t detid[maxhits]; //detid goes from 0 to 5
    Int_t detkind[maxhits];
    for(Int_t i=0;i<maxhits;i++){
      channum[i]=-1;
      detid[i]=-1;
      detkind[i]=-1;
    }

    for(Int_t i=0;i<NumHits;i++){
      channum[i]=id[i]-idConst;
      detid[i]=idDetMap[id[i]-idConst];
      detkind[i]=idKindMap[id[i]-idConst];

      // std::cout<<"channel number "<<channum[i]<<" detid "<<detid[i]<<endl;
    }//end identification and values loop

    /////////////////////////End Identification block////////////////

    /////////////////////////Begin multiplicity block////////////////
    Int_t arraymult=0;
    Int_t recoilmult=0;
    Int_t elummult=0;
    Int_t tacmult=0;
    for(Int_t i=0;i<NumHits;i++){
      if(detkind[i]==6) tacmult++;
      for(Int_t j=0;j<i;j++){
	if(detid[i]==detid[j]&&detkind[i]!=detkind[j]){
	  //	  std::cout<<" match found "<<channum[i]<<" "<<channum[j]<<" det kinds "<<detkind[i]<<" and "<<detkind[j]<<" det number "<<detid[i]<<" "<<detid[j]<<endl;
	  
	  if(detkind[i]==0||detkind[i]==1||detkind[i]==2) arraymult++;
	  if(detkind[i]==3||detkind[i]==4) recoilmult++;
	}
	if(detid[i]!=detid[j]&&detkind[i]==detkind[j]){
	  if(detkind[i]==5) elummult++;
	}
	if(detid[i]==detid[j]&&detkind[i]==detkind[j]){
	  //	  std::cout<<" double hit "<<detid[i]<<" "<<detid[j]<<" "<<detkind[i]<<" "<<detkind[j]<<endl;
	}
      }
    }//end multiplicity loop
    //   std::cout<<" array mult "<<arraymult<<" recoil mult "<<recoilmult<<" elum mult " << elummult<< " tac mult "<<tacmult<<endl;
 
    hmult1->Fill(recoilmult,arraymult);
    hmult2->Fill(tacmult,arraymult);
    hmult3->Fill(recoilmult,elummult);

    //////////////END MULTIPLICITY BLOCK/////////////////

    ///////////Fill output structure block////////////////


    Float_t ee[maxhits];
    ULong64_t ts[maxhits];
    

    for(Int_t i=0;i<NumHits;i++){
      ee[i]=((Float_t)post_rise_energy[i]-(Float_t)pre_rise_energy[i])/M;
      ts[i]=event_timestamp[i];
      hfill->Fill(detkind[i]);
      if(detkind[i]==0){ 
	psd.Energy[detid[i]]=ee[i];

	psd.EnergyTimestamp[detid[i]]=ts[i];
	  }
      if(detkind[i]==1){ 
	psd.XF[detid[i]]=ee[i];
	psd.XFTimestamp[detid[i]]=ts[i];
	  }
      if(detkind[i]==2){ 
	psd.XN[detid[i]]=ee[i];
	psd.XNTimestamp[detid[i]]=ts[i];
	  }
      if(detkind[i]==3){ //DE
	if(ee[i]<0){
	  ee[i]=TMath::Abs(ee[i]);
	  psd.RDT[(detid[i]-100)*2]=ee[i];
	  psd.RDTTimestamp[(detid[i]-100)*2]=ts[i];
	}
      }
      if(detkind[i]==4){ //E
	if(ee[i]<0){
	  ee[i]=TMath::Abs(ee[i]);	  
	  psd.RDT[(detid[i]-100)*2+1]=ee[i];
	  psd.RDTTimestamp[(detid[i]-100)*2+1]=ts[i];
	}
      }
      if(detkind[i]==5){ 
	psd.ELUM[detid[i]-200]=ee[i];
	psd.EnergyTimestamp[detid[i]-200]=ts[i];
	  }   
      if(detkind[i]==6){ 
	psd.TAC[0]=ee[i];
	psd.TACTimestamp[0]=ts[i];
	  }

    }//end fill
    gen_tree->Fill();
    FilledEntries++;

    ///////////End Fill output structure block////////////


    ////////////////event type determination/////////////
    if(arraymult==3){
      event_type=0;
      evttype->Fill(event_type);
    }
    if(arraymult==3&&recoilmult==1&&tacmult==1){
      event_type=1;
      evttype->Fill(event_type);
    } //3 array pairs (3 signals), 1 recoil pair, 1 TAC. Good event.
    if(arraymult==1&&recoilmult==1&&tacmult==1){
      event_type=2;
      evttype->Fill(event_type);
    } //1 array pair. Reconstruct 3rd signal. Recoil and TAC present.
    if(arraymult==3&&recoilmult==1&&tacmult==0){
      event_type=3; //otherwise good event, missing TAC.
      evttype->Fill(event_type);
    }
    if(arraymult==0&&recoilmult==1){
      event_type=4;   
      evttype->Fill(event_type);
    }          //recoil only event
    if(elummult==1){
      evttype->Fill(event_type);
      event_type=5;       
    }      //elum event
    else evttype->Fill(event_type);
      //other event types can be added if necessary; possible decomposition of events with multple recoil hits but otherwise good array events; multiple array hits;
    psd.EventType=event_type;

    //////////////////////////END EVENT TYPE DETERMINATION////////////////
   

    //////////////////////////////TRACE PROCESSING///////////////////////
   Int_t rangelookup;
	
   if(traceprocessing){

	b_trace->GetEntry(entry);
	for(Int_t i=0;i<NumHits;i++){
	  if(detkind[i]>-1){
	  
	    average[detkind[i]].counter++;
	    for(Int_t j=0;j<203;j++){
	      Float_t norm=(trace[i][j]-m2_end_sample[i])*8000/ee[i];
	      // average[detkind[i]].trace[j]+=(Double_t)(trace[i][j]&0x3fff);
	      // htrace[detkind[i]]->SetBinContent(j,average[detkind[i]].trace[j]);
	      htrace[detkind[i]]->SetBinContent(j,norm);
	      // htrace[detkind[i]]->Fill(j,(trace[i][j]-m2_end_sample[i])*8000/ee[i]);
		  if(detkind[i]==0){
		    htrace0[detid[i]]->Fill(j,norm);
		    rangelookup=fitrangeA[detid[i]];
		  }
		  if(detkind[i]==3){
		    htrace3[detid[i]-100]->Fill(j,norm);
		    rangelookup=fitrangeR[detid[i]-100];
		  
		  }
		  if(detkind[i]==4){
		    htrace4[detid[i]-100]->Fill(j,norm);
		  }
	    }
	 
	    switch(detkind[i]){
	    case 0:
	    
	      f1->SetRange(50,rangelookup);
	      htrace0[detid[i]]->Fit(f1,"NQR");
	      psd.ETC[detid[i]]=(Float_t)f1->GetX(4000,50,rangelookup);
	      hrise0->Fill(detid[i],psd.ETC[detid[i]]*10);
	      break;
	    case 3:
	     
	      f1->SetRange(50,rangelookup);
	      htrace[3]->Fit(f1,"NQR");
	      if(f1->GetX(4000,50,rangelookup)>0) psd.RTC[detid[i]-100]=(Float_t)f1->GetX(4000,50,rangelookup);
	      else psd.RTC[detid[i]]=0;
	      hrise3->Fill(detid[i]-100,psd.RTC[detid[i]-100]*10);
	      break;
	    default: break;
	    }

	    // switch(detkind[i]){
	    // case 0:
	    //   ct0->cd();
	    //   htrace[0]->Draw("colz");
	    // case 1:
	    //   ct1->cd();
	    //   htrace[1]->Draw("colz");
	    // case 2:
	    //   ct2->cd();
	    //   htrace[2]->Draw("colz");
	    // case 3:
	    //   ct3->cd();
	    //   htrace[3]->Draw("colz");
	    // case 4:
	    //   ct4->cd();
	    //   htrace[4]->Draw("colz");
	    // case 5:
	    //   ct5->cd();
	    //   htrace[5]->Draw("colz");
	    // case 6:
	    //   ct6->cd();
	    //   htrace[6]->Draw("colz");
	    // default: break;
	    // }
	    

	  }
	}
   }
	////////////////////////////END TRACE PROCESSING//////////////////////
    
  
















         
  }//goodN loop
  
 



























  if(!goodN) SkippedEntries++;


  if (ProcessedEntries+SkippedEntries>NumEntries*Frac-1) {
    printf(" %3.0f%% (%llu/%llu k) processed in %6.1f seconds\n",Frac*100,(ProcessedEntries+SkippedEntries)/1000,NumEntries/1000,StpWatch.RealTime());
    StpWatch.Start(kFALSE);
    Frac+=0.1;

  }
 

  return kTRUE;
}

void sort3::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void sort3::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  // cta=new TCanvas("cta");
    // ct1->Divide(1,7);
    // ct1->cd(1);
 
  // for(int i=0;i<7;i++){
  //   // cout<<"counter "<<average[i].counter<<endl;
  //   for(int j=0;j<203;j++){
  //       average[i].trace[j]/=average[i].counter;
  //     htrace[i]->SetBinContent(j,average[i].trace[j]);
  //   }

  //   cta->cd(i+1);
  //   htrace[i]->Draw();
  // }
  // ct0=new TCanvas("ct0");
  // ct1=new TCanvas("ct1");
  // ct2=new TCanvas("ct2");
  // ct3=new TCanvas("ct3");
  // ct4=new TCanvas("ct4");
  // ct5=new TCanvas("ct5");
  // ct6=new TCanvas("ct6");

  // ct0->cd();
  // htrace[0]->Draw("colz");
  // ct1->cd();
  // htrace[1]->Draw("colz");
  // ct2->cd();
  // htrace[2]->Draw("colz");
  // ct3->cd();
  // htrace[3]->Draw("colz");
  // ct4->cd();
  // htrace[4]->Draw("colz");
  // ct5->cd();
  // htrace[5]->Draw("colz");
  // ct6->cd();
  // htrace[6]->Draw("colz");
  // cc1=new TCanvas("cc1","cc1",800,600);
  // cc1->Clear();
  // hmult1->Draw("colz");
  // cc2=new TCanvas("cc2","cc2",800,600);
  // cc2->Clear();
  // hmult2->Draw("colz");
  // cc3=new TCanvas("cc3","cc3",800,600);
  // cc3->Clear();
  // hmult3->Draw("colz");

  gen_tree->Write();
  oFile->Close();
  
  //  cc0 = new TCanvas("cc0","cc0",800,600);
  // cc0->Clear(); hEvents->Draw();  
  
  printf("Total processed entries : %3.1f k\n",(ProcessedEntries+SkippedEntries)/1000.0);
  printf("Total skipped entries : %3.1f k\n",(SkippedEntries)/1000.0);
  printf("Total filled entries : %3.1f k\n",(FilledEntries)/1000.0);
  printf("Total time for sort: %3.1f\n",StpWatch.RealTime());
  printf("Rate for sort: %3.1f k/s\n",(Float_t)ProcessedEntries/StpWatch.RealTime()/1000.0);
  // StpWatch.Start(kFALSE);

}
