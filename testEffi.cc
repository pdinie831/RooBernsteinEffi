//*******************************************************************************
//*                                                                             *
//*  test code to compare efficiency projection and fit parametrization         *
//*  using RooBernsteinClass.                                                   *
//*                                                                             *
//*  P.Dini fecit, Anno Domini MMXVIII                                          *
//*                                                                             *
//*******************************************************************************


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
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooRealVar.h>
#include "RooBernsteinEffi.h"
using namespace std;


int main (int argc, char** argv) {

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
bool   wrongTagged = false;
char OutFileNameInputFit[300] =  "testGoofitEffi3DB0-2016-Q2Bin-3-Bins-25-25-25-BernDeg-5-5-4-integraBin.root";
char OutFileName[300]         =  "testEffi-Bin3.root";
char ListParName[300] 	      =  "ListParValues-2016-Q2Bin-3-Bins-25-25-25-BernDeg-5-5-4-integraBin.plo";

//============================
// maxDegree START
//============================
  int maxDegree1 =5;
  int maxDegree2 =5;
  int maxDegree3 =4;
//============================
// maxDegree END
//============================

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
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8]\n"<<endl;
    exit(1);

}
  std::cout<<"--------------------------------------------\n"<<endl;
  std::cout<<" Setting selection for Q^2 bin: "<<*argv[1]<<" ==> "<<Q2Min<<"<Q^2<"<<Q2Max<<std::endl;
  std::cout<<" Num Param CosL:  "<<maxDegree1<<std::endl;
  std::cout<<" Num Param CosK:  "<<maxDegree2<<std::endl;
  std::cout<<" Num Param Phi :  "<<maxDegree3<<std::endl;
  std::cout<<"--------------------------------------------\n"<<endl;
if (argc>2 && (strcmp(argv[2],"w") == 0) ){
  std::cout<<"====================================="<<endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<"====================================="<<endl;
  sprintf(OutFileNameInputFit,"testGoofitEffi3DB0-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(OutFileName,"testEffiRooFit-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(ListParName,"ListParValues-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.plo", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  wrongTagged = true;
}else{
  std::cout<<"======================================="<<endl;
  std::cout<<"======== DEFAULT: TAGGED =============="<<std::endl;
  std::cout<<"======================================="<<endl;
  sprintf(OutFileNameInputFit,"testGoofitEffi3DB0-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(OutFileName,"testEffiRooFit-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(ListParName,"ListParValues-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.plo", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  wrongTagged = false;
}
if (argc>3 && (strcmp(argv[3],"c") == 0) ){
  std::cout<<"=========================================================================="<<endl;
  std::cout<<"Setting the option: Efficiency function evaluated in the CENTER of the bin"<<std::endl;
  std::cout<<"=========================================================================="<<endl;
//  
  std::string str(OutFileNameInputFit) ;
  str.replace(str.end()-5,str.end(),"-centerBin.root");
  sprintf(OutFileNameInputFit,str.c_str());
//  

  std::string plo(ListParName) ;
  plo.replace(plo.end()-4,plo.end(),"-centerBin.plo");
  sprintf(ListParName,plo.c_str());
}else{
  std::cout<<"======================================================================="<<endl;
  std::cout<<"======== DEFAULT: INTEGRAL of the Efficiency in the bin  =============="<<std::endl;
  std::cout<<"======================================================================="<<endl;
//  
  std::string str(OutFileNameInputFit) ;
  str.replace(str.end()-5,str.end(),"-integraBin.root");
  sprintf(OutFileNameInputFit,str.c_str());
//  
  std::string plo(ListParName) ;
  plo.replace(plo.end()-4,plo.end(),"-integraBin.plo");
  sprintf(ListParName,plo.c_str());
}
  

  TCanvas* ProjEffiPlots = new TCanvas("ProjEffiPlots","X  Efficiency",200,10,900,780);
  ProjEffiPlots->Divide(2,2);  

  TCanvas* ClosurePlots = new TCanvas("ClosurePlots","X  Efficiency",200,10,900,780);
  ClosurePlots->Divide(2,2);  

  int numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);

  TFile*OutFileInputFit=0;
  TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
  TH3D *HxReco=0;
  TH3D *HxGene=0;
  TH3D *HEffi3D=0;
  TH1D *HEffi3DX=0;
  TH1D *HEffi3DY=0;
  TH1D *HEffi3DZ=0;
  RooRealVar x("x", "x",  XMinCosThetaL,XMaxCosThetaL);
  RooRealVar y("y", "y",  XMinCosThetaK,XMaxCosThetaK);
  RooRealVar z("z", "z",  XMinPhi,XMaxPhi);
  double xBinw = (XMaxCosThetaL-XMinCosThetaL)/xCosLHBin;
  double yBinw = (XMaxCosThetaK-XMinCosThetaK)/xCosKHBin;
  double zBinw = (XMaxPhi-XMinPhi)/xPhiHBin;

  if (TFile::Open(OutFileNameInputFit,"READ"))
  {
   OutFileInputFit = TFile::Open(OutFileNameInputFit,"READ");
   cout<<"File:"<<OutFileNameInputFit<<" FOUND !!!"<<endl;
  }else{
    cout<<"File:"<<OutFileNameInputFit<<" not found!!! Exit..."<<endl;
    exit(1);
  }
  if(wrongTagged){ 
   HxReco    = (TH3D*)OutFileInputFit->Get("HwReco");
   cout<<"Setting wrong tagged => events = "<<HxReco->GetEntries()<<endl;
  }else{ 
   HxReco    = (TH3D*)OutFileInputFit->Get("HxReco");
   cout<<"Setting  tagged => events = "<<HxReco->GetEntries()<<endl;
  } 
  if(!HxReco ){
    cout<<"HxReco Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HxReco Histo: OK FOUND!!! Entries: "<<HxReco->GetEntries()<<endl;
  }
  HxGene    = (TH3D*)OutFileInputFit->Get("HxGene");
  if(!HxGene ){
    cout<<"HxGene Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HxGene Histo: OK FOUND!!! Entries: "<<HxGene->GetEntries()<<endl;
  }  
  HEffi3D    = (TH3D*)OutFileInputFit->Get("HEffi3D");
  if(!HEffi3D ){
    cout<<"HEffi3D Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HEffi3D Histo: OK FOUND!!! Entries: "<<HEffi3D->GetEntries()<<endl;
  }  
  HEffi3DX    = (TH1D*)OutFileInputFit->Get("HEffi3DX");
  if(!HEffi3DX ){
    cout<<"HEffi3DX Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HEffi3DX Histo: OK FOUND!!! Entries: "<<HEffi3DX->GetEntries()<<endl;
  }  
  HEffi3DY    = (TH1D*)OutFileInputFit->Get("HEffi3DY");
  if(!HEffi3DY ){
    cout<<"HEffi3DY Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HEffi3DY Histo: OK FOUND!!! Entries: "<<HEffi3DY->GetEntries()<<endl;
  }  
  HEffi3DZ    = (TH1D*)OutFileInputFit->Get("HEffi3DZ");
  if(!HEffi3DZ ){
    cout<<"HEffi3DZ Histo: not found!!! Exit..."<<endl;
    exit(1);
  }else{
    cout<<"HEffi3DZ Histo: OK FOUND!!! Entries: "<<HEffi3DZ->GetEntries()<<endl;
  }  

  RooArgList *coefLis = new RooArgList();
  std::vector<RooRealVar> parLis; 

 
  char varName[10];
  double parIni=0.0;
  std::cout<<"Try to open list of parameters :"<< ListParName <<std::endl;
  std::fstream *parListInput = new std::fstream(ListParName,std::ifstream::in);
  if(parListInput->is_open()){
     std::cout<<"List of initial parameters :"<< ListParName <<" FOUND!!!"<<std::endl;
     for (int i=0;i<numParameters;++i){
      	    sprintf(varName, "p%d", i);
 	    *parListInput >> parIni;
     			 parLis.emplace_back(varName,varName, parIni,-100.0,100.0);
			 std::cout<<"Coeff ("<<varName<<") = "<<parLis[i].getVal()<<std::endl;
     }
    parListInput->close();
  }else{
     std::cout<<"List of initial parameters "<< ListParName <<" not found!"<<std::endl;
     exit(1);
  }   
  for (int i=0;i<numParameters;++i){
   coefLis->add(parLis[i]);
  }	      


  RooBernsteinEffi * Effi = new RooBernsteinEffi("Effi","Effi",x,y,z,*coefLis,maxDegree1,maxDegree2,maxDegree3);
  
  TH3D* HTestEffi = new   TH3D( "HTestEffi"    , "HTestEffi",	xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 								xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
								xPhiHBin , XMinPhi, XMaxPhi );
  TH3D* HtestReco = new   TH3D( "HtestReco"    , "HtestReco",	xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 								xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
								xPhiHBin , XMinPhi, XMaxPhi );
  HTestEffi->Sumw2();  
  std::cout<<"Number of bins for CosL = "<<xCosLHBin<<std::endl;
  std::cout<<"Number of bins for CosK = "<<xCosKHBin<<std::endl;
  std::cout<<"Number of bins for Phi  = "<<xPhiHBin <<std::endl;
  std::cout<<"Bin Width for CosL      = "<<xBinw<<std::endl;
  std::cout<<"Bin Width for CosK      = "<<yBinw<<std::endl;
  std::cout<<"Bin Width for Phi       = "<<zBinw<<std::endl;
  for(int i = 1; i <= xCosLHBin ; ++i) {
   double xi = HxGene->GetXaxis()->GetBinCenter(i);
   for(int j = 1 ; j <= xCosKHBin ; ++j) {
    double yj = HxGene->GetYaxis()->GetBinCenter(j);
    for(int k = 1 ; k <= xPhiHBin  ; ++k) {
    double zk = HxGene->GetZaxis()->GetBinCenter(k);
     x.setVal(xi);
     y.setVal(yj);
     z.setVal(zk);
	
//     
     double evalGene = HxGene->GetBinContent(i,j,k);
     double evalEffi = Effi->getVal();
//     double evalEffi = Effi->evaluateInt(xBinw,yBinw,zBinw );
     HTestEffi->Fill(xi,yj,zk,evalEffi);      
     HtestReco->Fill(xi,yj,zk,evalEffi*evalGene);      
     } 
    }
   }
   TH1D* HTestEffiX = (TH1D*) HTestEffi->ProjectionX("HTestEffiX",1,HTestEffi->GetNbinsY(),1,HTestEffi->GetNbinsZ());HTestEffiX->SetTitle("Effi Cos#theta_{L}");
   TH1D* HTestEffiY = (TH1D*) HTestEffi->ProjectionY("HTestEffiY",1,HTestEffi->GetNbinsX(),1,HTestEffi->GetNbinsZ());HTestEffiY->SetTitle("Effi Cos#theta_{K}");
   TH1D* HTestEffiZ = (TH1D*) HTestEffi->ProjectionZ("HTestEffiZ",1,HTestEffi->GetNbinsX(),1,HTestEffi->GetNbinsY());HTestEffiZ->SetTitle("Effi #phi");
   TH1D* HtestRecoX = (TH1D*) HtestReco->ProjectionX("HtestRecoX",1,HtestReco->GetNbinsY(),1,HtestReco->GetNbinsZ());HtestRecoX->SetTitle("Reco Cos#theta_{L}");
   TH1D* HtestRecoY = (TH1D*) HtestReco->ProjectionY("HtestRecoY",1,HtestReco->GetNbinsX(),1,HtestReco->GetNbinsZ());HtestRecoY->SetTitle("Reco Cos#theta_{K}");
   TH1D* HtestRecoZ = (TH1D*) HtestReco->ProjectionZ("HtestRecoZ",1,HtestReco->GetNbinsX(),1,HtestReco->GetNbinsY());HtestRecoZ->SetTitle("Reco #phi");
   TH1D* HxRecoX    = (TH1D*) HxReco->ProjectionX("HxRecoX",1,HxReco->GetNbinsY(),1,HxReco->GetNbinsZ());HxRecoX->SetTitle("Reco Cos#theta_{L}");
   TH1D* HxRecoY    = (TH1D*) HxReco->ProjectionY("HxRecoY",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsZ());HxRecoY->SetTitle("Reco Cos#theta_{K}");
   TH1D* HxRecoZ    = (TH1D*) HxReco->ProjectionZ("HxRecoZ",1,HxReco->GetNbinsX(),1,HxReco->GetNbinsY());HxRecoZ->SetTitle("Reco #phi");
     std::cout<<"HTestEffi Integral = "<<HTestEffi->Integral() <<std::endl;
     std::cout<<"HEffi3D   Integral = "<<HEffi3D->Integral()   <<std::endl;
   ProjEffiPlots->cd(1);
   HEffi3DX->Draw("E1");
   HTestEffiX->SetLineColor(kRed);
   HTestEffiX->Draw("same,HIST C");
   ProjEffiPlots->cd(2);
   HEffi3DY->Draw("E1");
   HTestEffiY->SetLineColor(kRed);
   HTestEffiY->Draw("same,HIST C");
   ProjEffiPlots->cd(3);
   HEffi3DZ->Draw("E1");
   HTestEffiZ->SetLineColor(kRed);
   HTestEffiZ->Draw("same,HIST C");
   ClosurePlots->cd(1);
   HxRecoX->Draw("E1");
   HtestRecoX->SetLineColor(kRed);
   HtestRecoX->Draw("same,HIST C");
   ClosurePlots->cd(2);
   HxRecoY->Draw("E1");
   HtestRecoY->SetLineColor(kRed);
   HtestRecoY->Draw("same,HIST C");
   ClosurePlots->cd(3);
   HxRecoZ->Draw("E1");
   HtestRecoZ->SetLineColor(kRed);
   HtestRecoZ->Draw("same,HIST C");
   OutFile->cd();
   HTestEffi->Write();
   HTestEffiX->Write();
   HTestEffiY->Write();
   HTestEffiZ->Write();
   HEffi3D->Write();
   HEffi3DX->Write();
   HEffi3DY->Write();
   HEffi3DZ->Write();
   ClosurePlots->Write();
   ProjEffiPlots->Write();
   OutFile->Close();
}
