//*******************************************************************************
//*                                                                             *
//*  simple code to make histograms for the Q^2 bin efficiencies                *
//*                                                                             *
//*  HxReco : Hist of the reconstructed events                                  *
//*  HxGene : Hist of the generated     events                                  *
//*  HEffi3D: HxReco/HxGene                                                     *
//*                                                                             *
//*                                                                             *
//*  P.Dini fecit, Anno Domini MMXVIII                                          *
//*                                                                             *
//*******************************************************************************
//g++ -O3 -o testEffi3DB0-2016-makeHisto testEffi3DB0-2016-makeHisto.cc `root-config --cflags --libs` -lProof  -lRooFit


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h> 
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>   
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TEfficiency.h>
#include  <TStopwatch.h>
using namespace std;

void CreateInputHistoFile();  

//********************************************************************************************************
// these files are the hadd root command outputs of the MC RECO and GENE events
//
  char InputFileNameMCReco[300]   = "2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root";
  char InputFileNameMCGene[300]   = "2016MC_GEN_LMNR_double_add2.root";
//
//********************************************************************************************************
char InputRecoB0TreeName[10]	= "ntuple";
char InputGeneB0TreeName[10]	= "ntuple";
char eosRecoDir[100] = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/";
char eosGeneDir[100] = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN/";
char RecoDir[100]    =  ".";
char GeneDir[100]    =  ".";
double Q2Min = 0.; 
double Q2Max = 0.; 
int    Q2Bin = -1;
int    xCosLHBin =  0;
int    xCosKHBin =  0;
int    xPhiHBin  =  0;
double XMinCosThetaL = -1.;
double XMaxCosThetaL = 1.;
double XMinCosThetaK = -1.;
double XMaxCosThetaK = 1.;
double XMinPhi       =-TMath::Pi();
double XMaxPhi       = TMath::Pi();
double BinWCosThetaL=(XMaxCosThetaL-XMinCosThetaL)/double(xCosLHBin);
double BinWCosThetaK=(XMaxCosThetaK-XMinCosThetaK)/double(xCosKHBin);
double BinWPhi      =(XMaxPhi-XMinPhi)/double(xPhiHBin);
double xMinQMuMu = 1.;
double xMaxQMuMu = 19.;
double XMinSign = 5.0;
double XMaxSign = 5.6;
double B0Mass	= 5.27962;
double B0Sigma  = 0.030;
double NSigma  = 3.;
double XMinSignW = B0Mass - NSigma*B0Sigma;
double XMaxSignW = B0Mass + NSigma*B0Sigma;
double XStepSign = 0.001;
float xMassHBin  = (XMaxSign -XMinSign)/XStepSign;
float xQ2HBin	 = (xMaxQMuMu -xMinQMuMu)/0.1;

char OutFileNameInputHisto[300] =  "";
TFile*OutFileInputHisto=0;

// Pointer to Histograms
  TH3D *HxReco=0;
  TH3D *HwReco=0;
  TH3D *Hw1Reco=0;
  TH3D *Hw2Reco=0;
  TH3D *HxGene=0;
  TEfficiency* pEff       = 0	;
  TEfficiency* pEffCosL   = 0	;
  TEfficiency* pEffCosK   = 0	;
  TEfficiency* pEffPhi    = 0	;
  TEfficiency* pEffCosLK  = 0	;    
  float xMassHBin2  =  xMassHBin /5; // plot only!
  TH1D* HxMass    = new TH1D( "HxMass"     , "B^{0} Mass",	xMassHBin2, XMinSign,  XMaxSign);
  TH1D* HxRecoMass= new TH1D( "HxRecoMass" , "B^{0} Mass",	xMassHBin2, XMinSign,  XMaxSign);
  TH1D* HxGeneMass= new TH1D( "HxGeneMass" , "B^{0} Mass",	xMassHBin2, XMinSign,  XMaxSign);

  TH1D* HxRecoQ2     = new TH1D( "HxRecoQ2"    , "Q^{2} Bin correct tagged",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* HwRecoQ2     = new TH1D( "HwRecoQ2"    , "Q^{2} Bin wrong   tagged",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* Hw1RecoQ2    = new TH1D( "Hw1RecoQ2"   , "Q^{2} Bin wrong   tagged 1",   xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* Hw2RecoQ2    = new TH1D( "Hw2RecoQ2"   , "Q^{2} Bin wrong   tagged 2",   xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* HxRecoTotQ2  = new TH1D( "HxRecoTotQ2" , "Q^{2} Tot correct tagged ",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* HwRecoTotQ2  = new TH1D( "HwRecoTotQ2" , "Q^{2} Tot wrong   tagged ",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* Hw1RecoTotQ2 = new TH1D( "Hw1RecoTotQ2", "Q^{2} Tot wrong   tagged 1",   xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* Hw2RecoTotQ2 = new TH1D( "Hw2RecoTotQ2", "Q^{2} Tot wrong   tagged 2",   xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* HxGeneQ2     = new TH1D( "HxGeneQ2"    , "Q^{2} Bin ",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);
  TH1D* HxGeneTotQ2  = new TH1D( "HxGeneTotQ2" , "Q^{2} Tot ",	 xQ2HBin   , xMinQMuMu, xMaxQMuMu);

  char PDFNameRecoHisto[50]	  =  "";
  char PDFNameGeneHisto[50]	  =  "";
  char PDFNameEffiHisto[50]       =  "";
  float MarkerSize    = 0.35;

int main (int argc, char** argv) {
if (argc<=1 ){
    cout<<"Q2Bin not set"<<endl;
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8]\n"<<endl;
    exit(1);
}   


 
switch ( *argv[1] ) {

  case '0' : 
   Q2Min = 1.; 
   Q2Max = 2.;
   Q2Bin = 0;
      xCosLHBin =  25;
      xCosKHBin =  25;
      xPhiHBin  =  25;
    break;
  case '1' : 
   Q2Min = 2.; 
   Q2Max = 4.3; 
   Q2Bin = 1;
   xCosLHBin =  25;
   xCosKHBin =  25;
   xPhiHBin  =  25;
    break;
  case '2' : 
   Q2Min = 4.3; 
   Q2Max = 6.; 
   Q2Bin = 2;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '3' : 
   Q2Min = 6.; 
   Q2Max = 8.68; 
   Q2Bin = 3;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '4' : 
   Q2Min = 8.68; 
   Q2Max = 10.09; 
   Q2Bin = 4;
   xCosLHBin =  25;
   xCosKHBin =  25;
   xPhiHBin  =  25;
    break;
  case '5' : 
   Q2Min = 10.09; 
   Q2Max = 12.86; 
   Q2Bin = 5;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '6' : 
   Q2Min = 12.86; 
   Q2Max = 14.18; 
   Q2Bin = 6;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '7' : 
   Q2Min = 14.18; 
   Q2Max = 16.; 
   Q2Bin = 7;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '8' : 
   Q2Min = 16; 
   Q2Max = 19.; 
   Q2Bin = 8;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;

  default : 
    // Process for all other cases.
    cout<<"Q2Bin not set correctly!!!"<<endl;
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8] dir [where dir=e if you want to read files from Sara's directories]\n"<<endl;
    exit(1);

}

  std::cout<<"--------------------------------------------\n"<<endl;
  std::cout<<" Setting selection for Q^2 bin: "<<*argv[1]<<" ==> "<<Q2Min<<"<Q^2<"<<Q2Max<<std::endl;
  std::cout<<"--------------------------------------------\n"<<endl;
  std::cout<<"--------------------------------------------\n"<<endl;
if (argc>2 && (strcmp(argv[2],"e") == 0) ){

  std::cout<<" Setting EOS ad directory files  \n"<<endl;
  memcpy ( GeneDir, eosGeneDir, strlen(eosGeneDir)+1 );
  memcpy ( RecoDir, eosRecoDir, strlen(eosRecoDir)+1 );

}
  std::cout<<" Reading MC Reconstructed  files in directory: "<<RecoDir<<std::endl;
  std::cout<<" Reading MC Generated      file: in directory: "<<GeneDir<<std::endl;
  std::cout<<"--------------------------------------------\n"<<endl;

  sprintf(PDFNameRecoHisto,"B0-RecoHist-2016-OnlyPlot-Q2Bin-%d-Bins-%d-%d-%d.pdf", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin); 
  sprintf(PDFNameGeneHisto,"B0-GeneHist-2016-OnlyPlot-Q2Bin-%d-Bins-%d-%d-%d.pdf", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin); 
  sprintf(PDFNameEffiHisto,"B0-EffiHist-2016-OnlyPlot-Q2Bin-%d-Bins-%d-%d-%d.pdf", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin); 


  sprintf(OutFileNameInputHisto,"testGoofitEffi3DB0-2016-InputHisto-Q2Bin-%d-Bins-%d-%d-%d.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin); 
  TStopwatch TimeWatch;
  TimeWatch.Start();

  if (!TFile::Open(OutFileNameInputHisto,"READ"))
  {
    std::cout<<"File:"<<OutFileNameInputHisto<<" not found!!! create..."<<std::endl;
  }else{
    std::cout<<"File:"<<OutFileNameInputHisto<<" found!!! I will save the previous file in "<<OutFileNameInputHisto<<".tmp"<<std::endl;
    gSystem->Exec(Form("mv %s %s.tmp",OutFileNameInputHisto,OutFileNameInputHisto));
  }
  CreateInputHistoFile();
  TimeWatch.Stop();
  TimeWatch.Print();
  return 0 ;
}  
void CreateInputHistoFile(){   

  TCanvas* c3 = new TCanvas("c3","Reco Histograms",200,10,900,780);
  TCanvas* c4 = new TCanvas("c4","Gene Histograms",200,10,900,780);
  TCanvas* c5 = new TCanvas("c5","Effi Histograms",200,10,900,780);
  c3->Divide(2,2);  
  c4->Divide(2,2);  
  c5->Divide(2,2);  
  
  TChain* RecoB0Tree = new TChain();

  int nfile = 0;
  
  nfile = RecoB0Tree->Add(Form("%s/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/%s",RecoDir,InputRecoB0TreeName));  
  if(nfile==0 ||  !RecoB0Tree->GetFile() ){
    cout<<"Error:  no Reco files found!!!\n"<<endl;
    exit(1);
  }else{
    printf("Try to open %s/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/%s \n",RecoDir,InputRecoB0TreeName);
    cout<<"Opening "<<nfile <<" Reco files found!!!\n"<<endl;
  }  
  if(!RecoB0Tree ){
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" not found!!!\n"<<endl;
    exit(1);
  }else{
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" OK FOUND!!!\n"<<endl;
  }  
  nfile = 0;
  TChain* GeneB0Tree = new TChain();  
  nfile = GeneB0Tree->Add(Form("%s/2016MC_GEN_LMNR_double_sub*_p*.root/%s",GeneDir,InputGeneB0TreeName));
  if( nfile==0 ||  !GeneB0Tree->GetFile() ){
    cout<<"Error:  no Gene files found!!!\n"<<endl;
    exit(1);
  }else{
    printf("Try to open %s/2016MC_GEN_LMNR_double_sub*_p*.root/%s \n",GeneDir,InputGeneB0TreeName);
    cout<<"Opening "<<nfile <<" Gene files found!!!\n"<<endl;
  }  
  if(!GeneB0Tree ){
    cout<<"TTree Gene Data: "<< InputGeneB0TreeName <<" not found!!!\n"<<endl;
    exit(1);
  }else{
    cout<<"TTree Gene Data: "<< InputGeneB0TreeName <<" OK FOUND!!!\n"<<endl;
  }  

  printf("(Mass Window     : xB0Mass>%8f && xB0Mass<%8f \n",XMinSign,XMaxSign);







  TH3D* HxGene = new   TH3D( "HxGene"    , "B^{0} Gene",		 xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
  TH3D* HxReco = new   TH3D( "HxReco"    , "B^{0} Reco correct tagged",  xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
	
  TH3D* HwReco = new   TH3D( "HwReco"    , "B^{0} Reco wrong   tagged",  xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
						
  TH3D* Hw1Reco = new  TH3D( "Hw1Reco"   , "B^{0} 1 Reco wrong   tagged",xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
						
  TH3D* Hw2Reco = new  TH3D( "Hw2Reco"   , "B^{0} 2 Reco wrong   tagged",xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );
						
  std::cout<<"Hist Mass Bin = "<<xMassHBin<<std::endl;
 

//======================================================================
//======================================================================
//======================================================================
//
//			      RECONSTRUCTED EVENTS
//
//======================================================================
//======================================================================
//======================================================================
  float  tagged_massF	;
  float  mumuMassF	;
  double tagged_mass	;
  double cos_theta_l	;
  double cos_theta_k	;
  double phi_kst_mumu	;
  double mumuMass	;
  double recQ2          ;
  float   truthMatchMum             ;
  float   truthMatchMup             ;
  float   truthMatchTrkm            ;
  float   truthMatchTrkp            ;
  float   genSignal                 ;
  float   tagB0                     ;
  RecoB0Tree->SetBranchAddress("tagged_mass"   ,&tagged_massF);
  RecoB0Tree->SetBranchAddress("cos_theta_l"   ,&cos_theta_l);
  RecoB0Tree->SetBranchAddress("cos_theta_k"   ,&cos_theta_k);
  RecoB0Tree->SetBranchAddress("phi_kst_mumu"  ,&phi_kst_mumu);
  RecoB0Tree->SetBranchAddress("mumuMass"      ,&mumuMassF);
  RecoB0Tree->SetBranchAddress("truthMatchMum" ,&truthMatchMum);
  RecoB0Tree->SetBranchAddress("truthMatchMup" ,&truthMatchMup);
  RecoB0Tree->SetBranchAddress("truthMatchTrkm",&truthMatchTrkm);
  RecoB0Tree->SetBranchAddress("truthMatchTrkp",&truthMatchTrkp);
  RecoB0Tree->SetBranchAddress("tagB0"         ,&tagB0);
  RecoB0Tree->SetBranchAddress("genSignal"     ,&genSignal);
  
  int nentries = (int)RecoB0Tree->GetEntries();
  
  for (Int_t i=0;i<nentries;i++) { 
   RecoB0Tree->GetEntry(i);
   tagged_mass   =  tagged_massF       ;
   mumuMass	 =  mumuMassF          ;  
   recQ2         =  mumuMass*mumuMass  ;

   if(truthMatchMum>0.&&truthMatchMup>0.&&truthMatchTrkm>0.&&truthMatchTrkp>0.){
    HxRecoTotQ2->Fill(recQ2);
    HxMass  ->Fill(tagged_mass);
// folding theta_l e phi     
//     cos_theta_l=fabs(cos_theta_l);
//     cos_theta_k=fabs(cos_theta_k);
//       phi_kst_mumu=fabs(phi_kst_mumu);

    if(tagged_mass>=XMinSignW&&tagged_mass<=XMaxSignW){
// first Q^2 Bin
     if(recQ2>Q2Min&&recQ2<Q2Max){
      if(cos_theta_l>=XMinCosThetaL&&cos_theta_l<=XMaxCosThetaL&&cos_theta_k>=XMinCosThetaK&&cos_theta_k<=XMaxCosThetaK){
 	  HxRecoMass->Fill(tagged_mass);
        if( (tagB0 == 1 && genSignal == 1) || ( tagB0 == 0 && genSignal == 2)){  //(correctly tagged)
 	  HxRecoQ2->Fill(recQ2);
 	  HxReco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
	}  
	if( (tagB0 == 1 && genSignal == 2) || ( tagB0 == 0 && genSignal == 1)){  //(wrong     tagged)
 	  HwRecoQ2->Fill(recQ2);
 	  HwReco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
        }
	if( (tagB0 == 1 && genSignal == 2) ){  //(wrong B0Bar?    tagged)
 	  Hw2RecoQ2->Fill(recQ2);
 	  Hw2Reco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
        }
	if( (tagB0 == 0 && genSignal == 1) ){  //(wrong B0?       tagged)
 	  Hw1RecoQ2->Fill(recQ2);
 	  Hw1Reco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
        }
      }
     }
    }
   }
  }
  char TXT[200];
  sprintf(TXT,"Mass Reco   Entries = %7f",HxRecoMass->GetEntries());
  cout<<"***********************************"<<endl;
  cout<<"***** RECONSTRUCTED EVENTS ********\n"<<endl;
  cout<<"***********************************\n"<<endl;
  cout<<"RecoB0Tree   Entries      = "<<nentries<<endl;
  cout<<"RecoB0Tree Wrong Tagged   = "<<HwReco->GetEntries()<<endl;
  cout<<"RecoB0Tree Wrong Tagged 1 = "<<Hw1Reco->GetEntries()<<endl;
  cout<<"RecoB0Tree Wrong Tagged 2 = "<<Hw2Reco->GetEntries()<<endl;
  cout<<TXT<<endl;
  cout<<"\n***********************************"<<endl;
  cout<<"***********************************"<<endl;
//======================================================================
//======================================================================
//======================================================================
//
//			      GENERATED EVENTS
//
//======================================================================
//======================================================================
//======================================================================
  double  genB0Mass	   ;
  double  gen_cos_theta_l   ;
  double  gen_cos_theta_k   ;
  double  gen_phi_kst_mumu  ;
  double  genQ2	           ;


  GeneB0Tree->SetBranchAddress("genB0Mass"	   ,&genB0Mass);
  GeneB0Tree->SetBranchAddress("gen_cos_theta_l"   ,&gen_cos_theta_l);
  GeneB0Tree->SetBranchAddress("gen_cos_theta_k"   ,&gen_cos_theta_k);
  GeneB0Tree->SetBranchAddress("gen_phi_kst_mumu"  ,&gen_phi_kst_mumu);
  GeneB0Tree->SetBranchAddress("genq2"  	   ,&genQ2);

  int ngentries = (int)GeneB0Tree->GetEntries();
  for (Int_t i=0;i<ngentries;i++) { 
     GeneB0Tree->GetEntry(i);
     genQ2=genQ2*genQ2;
     HxGeneTotQ2->Fill(genQ2);
// folding theta_l e phi     
//      gen_cos_theta_l = fabs(gen_cos_theta_l);
//      gen_cos_theta_k = fabs(gen_cos_theta_k);
//        gen_phi_kst_mumu= fabs(gen_phi_kst_mumu);
//
//  Q^2 Bin
     if(genQ2>Q2Min&&genQ2<Q2Max){
      if(gen_cos_theta_l>=XMinCosThetaL&&gen_cos_theta_l<=XMaxCosThetaL&&gen_cos_theta_k>=XMinCosThetaK&&gen_cos_theta_k<=XMaxCosThetaK){
 	 HxGeneQ2->Fill(genQ2);
 	 HxGeneMass->Fill(genB0Mass);
 	 HxGene->Fill(gen_cos_theta_l,gen_cos_theta_k,gen_phi_kst_mumu);
      } 
     } 
  }

  sprintf(TXT,"Mass Gen   Entries = %7f",HxGeneMass->GetEntries());
  cout<<"***********************************\n"<<endl;
  cout<<"***** GENERATED EVENTS   **********\n"<<endl;
  cout<<"***********************************\n"<<endl;
  cout<<"GenB0Tree   Entries = "<<ngentries<<endl;
  cout<<TXT<<endl;
  cout<<"\n***********************************"<<endl;
  cout<<"***********************************"<<endl;


   pEff = new TEfficiency(*HxReco,*HxGene);
   TH3D *HEffi3D   = (TH3D*)HxReco ->Clone(); HEffi3D  ->SetName("HEffi3D"  );HEffi3D  ->Sumw2(); HEffi3D  ->Divide(HxGene);
   TH3D *HEffi3Dw  = (TH3D*)HwReco ->Clone(); HEffi3Dw ->SetName("HEffi3Dw" );HEffi3Dw ->Sumw2(); HEffi3Dw ->Divide(HxGene);
   TH3D *HEffi3Dw1 = (TH3D*)Hw1Reco->Clone(); HEffi3Dw1->SetName("HEffi3Dw1");HEffi3Dw1->Sumw2(); HEffi3Dw1->Divide(HxGene);
   TH3D *HEffi3Dw2 = (TH3D*)Hw2Reco->Clone(); HEffi3Dw2->SetName("HEffi3Dw2");HEffi3Dw2->Sumw2(); HEffi3Dw2->Divide(HxGene);
 

  OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"RECREATE");
  TH1D* HEffi3DX = (TH1D*) HEffi3D->ProjectionX("HEffi3DX",1,HEffi3D->GetNbinsY(),1,HEffi3D->GetNbinsZ());HEffi3DX->SetTitle("HEffi3D Projection Cos#theta_{L}");
  TH1D* HEffi3DY = (TH1D*) HEffi3D->ProjectionY("HEffi3DY",1,HEffi3D->GetNbinsX(),1,HEffi3D->GetNbinsZ());HEffi3DY->SetTitle("HEffi3D Projection Cos#theta_{K}");
  TH1D* HEffi3DZ = (TH1D*) HEffi3D->ProjectionZ("HEffi3DZ",1,HEffi3D->GetNbinsX(),1,HEffi3D->GetNbinsY());HEffi3DZ->SetTitle("HEffi3D Projection #phi");

  TH1D* HEffi3DwX = (TH1D*) HEffi3Dw->ProjectionX("HEffi3DwX",1,HEffi3Dw->GetNbinsY(),1,HEffi3Dw->GetNbinsZ());HEffi3DwX->SetTitle("HEffi3Dw Projection Cos#theta_{L}");
  TH1D* HEffi3DwY = (TH1D*) HEffi3Dw->ProjectionY("HEffi3DwY",1,HEffi3Dw->GetNbinsX(),1,HEffi3Dw->GetNbinsZ());HEffi3DwY->SetTitle("HEffi3Dw Projection Cos#theta_{K}");
  TH1D* HEffi3DwZ = (TH1D*) HEffi3Dw->ProjectionZ("HEffi3DwZ",1,HEffi3Dw->GetNbinsX(),1,HEffi3Dw->GetNbinsY());HEffi3DwZ->SetTitle("HEffi3Dw Projection #phi");
  
  TH1D* HEffi3Dw1X = (TH1D*) HEffi3Dw1->ProjectionX("HEffi3Dw1X",1,HEffi3Dw1->GetNbinsY(),1,HEffi3Dw1->GetNbinsZ());HEffi3Dw1X->SetTitle("HEffi3Dw1 Projection Cos#theta_{L}");
  TH1D* HEffi3Dw1Y = (TH1D*) HEffi3Dw1->ProjectionY("HEffi3Dw1Y",1,HEffi3Dw1->GetNbinsX(),1,HEffi3Dw1->GetNbinsZ());HEffi3Dw1Y->SetTitle("HEffi3Dw1 Projection Cos#theta_{K}");
  TH1D* HEffi3Dw1Z = (TH1D*) HEffi3Dw1->ProjectionZ("HEffi3Dw1Z",1,HEffi3Dw1->GetNbinsX(),1,HEffi3Dw1->GetNbinsY());HEffi3Dw1Z->SetTitle("HEffi3Dw1 Projection #phi");

  TH1D* HEffi3Dw2X = (TH1D*) HEffi3Dw2->ProjectionX("HEffi3Dw2X",1,HEffi3Dw2->GetNbinsY(),1,HEffi3Dw2->GetNbinsZ());HEffi3Dw2X->SetTitle("HEffi3Dw2 Projection Cos#theta_{L}");
  TH1D* HEffi3Dw2Y = (TH1D*) HEffi3Dw2->ProjectionY("HEffi3Dw2Y",1,HEffi3Dw2->GetNbinsX(),1,HEffi3Dw2->GetNbinsZ());HEffi3Dw2Y->SetTitle("HEffi3Dw2 Projection Cos#theta_{K}");
  TH1D* HEffi3Dw2Z = (TH1D*) HEffi3Dw2->ProjectionZ("HEffi3Dw2Z",1,HEffi3Dw2->GetNbinsX(),1,HEffi3Dw2->GetNbinsY());HEffi3Dw2Z->SetTitle("HEffi3Dw2 Projection #phi");
  
  
  TH1D* HxRecoX  = (TH1D*) HxReco->ProjectionX("HxRecoX",1,HxReco->GetNbinsY(),1,HxReco->GetNbinsZ());HxRecoX->SetTitle("HxReco Projection Cos#theta_{L}");
  TH1D* HxRecoY  = (TH1D*) HxReco->ProjectionY("HxRecoY",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsZ());HxRecoY->SetTitle("HxReco Projection Cos#theta_{K}");
  TH1D* HxRecoZ  = (TH1D*) HxReco->ProjectionZ("HxRecoZ",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsY());HxRecoZ->SetTitle("HxReco Projection #phi");
  TH1D* HxGeneX  = (TH1D*) HxGene->ProjectionX("HxGeneX",1,HxGene->GetNbinsY(),1,HxGene->GetNbinsZ());HxGeneX->SetTitle("HxGene Projection Cos#theta_{L}");
  TH1D* HxGeneY  = (TH1D*) HxGene->ProjectionY("HxGeneY",1,HxGene->GetNbinsX(),1,HxGene->GetNbinsZ());HxGeneY->SetTitle("HxGene Projection Cos#theta_{K}");
  TH1D* HxGeneZ  = (TH1D*) HxGene->ProjectionZ("HxGeneZ",1,HxGene->GetNbinsX(),1,HxGene->GetNbinsY());HxGeneZ->SetTitle("HxGene Projection #phi");
  
  int count = 0;
  BinWCosThetaL=(XMaxCosThetaL-XMinCosThetaL)/double(xCosLHBin);
  BinWCosThetaK=(XMaxCosThetaK-XMinCosThetaK)/double(xCosKHBin);
  BinWPhi      =(XMaxPhi-XMinPhi)/double(xPhiHBin);
  if(TEfficiency::CheckConsistency(*HxReco,*HxGene))
  {
   pEff = new TEfficiency(*HxReco,*HxGene);
  }else{
   std::cout<<"NO Check consistency HxReco,HxGene!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      for(Int_t i = 0; i < xCosLHBin+2 ; ++i) {
       for(Int_t j = 0; j < xCosKHBin+2 ; ++j) {
  	for(Int_t k = 0; k < xPhiHBin+2 ; ++k) {
         if(HxReco->GetBinContent(i,j,k) > HxGene->GetBinContent(i,j,k)) {
          count++;
          std::cout<<"============================================="<<std::endl;
          std::cout<<"Check consistency false bin i,j,k= "<<i<<","<<j<<","<<k<<" Pass = "<<HxReco->GetBinContent(i,j,k)<<" Gen = "<<HxGene->GetBinContent(i,j,k)<<std::endl;
         }
        }
       }
      }
  }
  
  HxMass->Write();
  HxGeneMass->Write();
  HxRecoMass->Write();
  HxGene->Write();
  HxReco->Write(); 
  HwReco->Write(); 
  Hw1Reco->Write(); 
  Hw2Reco->Write(); 
  HxGeneQ2->Write();
  HxRecoQ2->Write(); 
  HwRecoQ2->Write(); 
  HxGeneTotQ2->Write();
  HxRecoTotQ2->Write(); 
  HwRecoTotQ2->Write(); 
  Hw1RecoTotQ2->Write(); 
  Hw2RecoTotQ2->Write(); 
  HEffi3D->Write();
  HEffi3DX->Write();
  HEffi3DY->Write();
  HEffi3DZ->Write();
  HEffi3DwX->Write();
  HEffi3DwY->Write();
  HEffi3DwZ->Write();
  HEffi3Dw1X->Write();
  HEffi3Dw1Y->Write();
  HEffi3Dw1Z->Write();
  HEffi3Dw2X->Write();
  HEffi3Dw2Y->Write();
  HEffi3Dw2Z->Write();
  HxRecoX->Write(); 
  HxRecoY->Write(); 
  HxRecoZ->Write(); 
  HxGeneX->Write();
  HxGeneY->Write();
  HxGeneZ->Write();
  
  gSystem->Exec(Form("mv %s %s.tmp",PDFNameRecoHisto,PDFNameRecoHisto));
  c3->cd(1);HxRecoX->Draw();
  c3->cd(2);HxRecoY->Draw();
  c3->cd(3);HxRecoZ->Draw();
  c3->Print(PDFNameRecoHisto);
  gSystem->Exec(Form("mv %s %s.tmp",PDFNameGeneHisto,PDFNameGeneHisto));
  c4->cd(1);HxGeneX->Draw();
  c4->cd(2);HxGeneY->Draw();
  c4->cd(3);HxGeneZ->Draw();
  c4->Print(PDFNameGeneHisto);
  gSystem->Exec(Form("mv %s %s.tmp",PDFNameEffiHisto,PDFNameEffiHisto));
  c5->cd(1);HEffi3DX->Draw();
  c5->cd(2);HEffi3DY->Draw();
  c5->cd(3);HEffi3DZ->Draw();
  c5->Print(PDFNameEffiHisto);
  OutFileInputHisto->Close();
  cout<<"**********************************************************************\n"<<endl;
  cout<<"save HxGene HxRec  in "<<OutFileNameInputHisto<<"\n"<<endl;
  cout<<"**********************************************************************\n"<<endl;
  if(count!=0){
   std::cout<<"TEfficiency::CheckConsistency ==> NO!!!"<<std::endl;
   std::cout<<"Tot false bins = "<<count<<std::endl;    
   }else{
   std::cout<<"TEfficiency::CheckConsistency OK!!!"<<std::endl;
  } 
  return ;  
}

