#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TFrame.h"

#include "TLatex.h"
//#include "TASImage.h"

//#include "mt2_bisect.h"

TString lumi_13TeV = "2.11 fb^{-1}";

using namespace std;

int _multiWP[3] = {2,1,0};//Medium - ele, Loose - mu

const int numSUSY1 = 1; //2 - T1tttt1200, 1 - T1tttt1500 
const int numSUSY2 = 1;
     
const int signalSUSY = 2;//2 - T1tttt1200, 3 - T1tttt1500 
   
const int SRNumberAll = 72;

bool comp (const pair<double, int> i, const pair<double, int> j) { return (i.first>j.first); }

void CMS_lumi( TPad* pad, int iPeriod, int iPosX ); 

double uncCalc(double a, double b, double c);
    
int SRID (double met, double HT, int njets, int nbjets, int ptBin, double mt2ll) {
    int index = 0;
    
    //Jan's SR
    //int njetsbin = int(njets > 4);
    int nbjetsbin = int(nbjets > 0) +  int(nbjets > 1) + int(nbjets > 2);
    //int nbjetsbin = int(nbjets > 0);
    int njetsbin = int(njets > 4);
    //int htbin = int (HT > 500);
    int htbin = int (HT > 400) + int (HT > 600);
    //int htbin = int (HT > 300) + int (HT > 1600);
    //int metbin = int(met > 150);
    //int metbin = int(met > 150) + int(met > 400);
    int metbin = int(met > 150) + int(met > 300);
    //int metbin = int(met > 400);

    //int mt2llbin = int(mt2ll > 120);

    
    // MET > 500 
    
    if(met > 300){
        htbin = 0;    
        njetsbin = 0;    
        nbjetsbin = 0;
    }

    if(HT > 600 && met < 300){
        njetsbin = 0;    
        nbjetsbin = 0;
        metbin = 0;
    }

    if(nbjets > 2 && HT < 600 && met < 300){
        htbin = 0;
        metbin = 0;
        njetsbin = 0;
    }

    if(HT < 400){
        njetsbin = 0;
    }

    if(HT > 400 && HT < 600){
        njetsbin = 0;
    }
    
    
    // 200 < MET < 500
    /*
    if(nbjetsbin == 2 && metbin == 1){
        htbin = 0;    
        njetsbin = 0;    
        if(mt2llbin == 1 && met < 200)
            nbjetsbin = 1;
    
    }
    
    if(nbjetsbin == 0 && metbin == 1){
        if(mt2llbin == 1){
            htbin = 0;    
            njetsbin = 0;    
        }
        if(htbin == 0 && mt2llbin == 0){
            njetsbin = 0;    
        }
        if(htbin > 0 && njetsbin < 2){
            njetsbin = 0;
        }
    
    }

    if(nbjetsbin == 1 && metbin == 1){
        if(htbin == 0){
            mt2llbin = 0;    
            njetsbin = 0;    
        }
        if(htbin > 0 && njetsbin < 2){
            mt2llbin = 0;    
            njetsbin = 0;    
        }
        if(htbin == 2 && njetsbin == 0){
            nbjetsbin = 0; 
        }
    
    }

    // MET < 200
    if(nbjetsbin == 2 && metbin == 0){
        if(mt2llbin == 1){
            htbin = 0;    
            njetsbin = 0;    
        }
        if(mt2llbin == 0)
            htbin = 0;
        if(mt2llbin == 1 && met < 200)
            nbjetsbin = 1;
    
    }

    if(nbjetsbin == 0 && met < 200 && mt2llbin == 1){
        htbin = 0;    
        njetsbin = 0;    
    
    }

    if(nbjetsbin == 1 && met < 200){
        if (mt2llbin == 1){
            htbin = 0;    
            njetsbin = 0;    
        }

        if (mt2llbin == 0 && htbin == 0 && njetsbin == 2)
            nbjetsbin = 0;
    
    }

    if(nbjetsbin == 0 && met < 200 && htbin == 2){
        njetsbin = 0;    
    }

    if(nbjetsbin == 1 && met < 200 && htbin == 2){
        njetsbin = 0;    
    }

    // Merging all bins where signal >> bkg
    if( metbin == 1 && nbjetsbin == 1 && htbin > 0 && njetsbin == 2){
        metbin = 2;
        mt2llbin = 1;
        njetsbin = 0;
        nbjetsbin = 0;
        htbin = 0;
    }
    
    if( metbin == 1 && nbjetsbin == 2 && htbin == 0 && njetsbin == 0){
        metbin = 2;
        mt2llbin = 1;
        njetsbin = 0;
        nbjetsbin = 0;
        htbin = 0;
    }
    */


    // Index calculating

    //index = 54 * metbin + 18 * htbin + 6 * njetsbin + 3 * mt2llbin + nbjetsbin; // for offZ region
    //index = 27 * metbin + 9 * htbin + 3 * njetsbin + nbjetsbin; // for offZ region
    index = 24 * metbin + 8 * htbin + 4 * njetsbin + nbjetsbin; 

    return index;

}

void showHist(TVirtualPad *, TH1D *, TH1D *, THStack *, string, string, string, double, TLegend*);

void readTree7414()
{
    
    const int nLeptonsMax = 10;
    TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.18, "X");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.15, "X");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

    // List of branches
    TBranch        *b__eventNb;   //!
    TBranch        *b__runNb;   //!
    TBranch        *b__lumiBlock;   //!
    TBranch        *b__weight;   //!
    TBranch        *b__genqpt;   //!
    TBranch        *b__nLeptons;   //!
    TBranch        *b__lPt;   //!
    TBranch        *b__lEta;   //!
    TBranch        *b__lPhi;   //!
    TBranch        *b__lE;   //!
    TBranch        *b__lPtmc;   //!
    TBranch        *b__lEtamc;   //!
    TBranch        *b__lPhimc;   //!
    TBranch        *b__lEmc;   //!
    TBranch        *b__nuPtmc;   //!
    TBranch        *b__nuEtamc;   //!
    TBranch        *b__nuPhimc;   //!
    TBranch        *b__nuEmc;   //!
    TBranch        *b__mtmc;   //!
    TBranch        *b__nEle;   //!
    TBranch        *b__nMu;   //!
    TBranch        *b__nTau;   //!
    TBranch        *b__flavors;   //!
    TBranch        *b__charges;   //!
    TBranch        *b__indeces;   //!
    TBranch        *b__isolation;   //!
    TBranch        *b__miniisolation;   //!
    TBranch        *b__multiisolation;   //!
    TBranch        *b__ptrel;   //!
    TBranch        *b__ptratio;   //!
    TBranch        *b__origin;   //!
    TBranch        *b__originReduced;   //!
    TBranch        *b__PVchi2;   //!
    TBranch        *b__PVerr;   //!
    TBranch        *b__ipPV;   //!
    TBranch        *b__ipPVerr;   //!
    TBranch        *b__ipZPV;   //!
    TBranch        *b__ipZPVerr;   //!
    TBranch        *b__ipPVmc;   //!
    TBranch        *b__3dIP;   //!
    TBranch        *b__3dIPerr;   //!
    TBranch        *b__3dIPsig;   //!
    TBranch        *b__mt;   //!
    TBranch        *b__mt2;   //!
    TBranch        *b__mt2ll;   //!
    TBranch        *b__mt2blbl;   //!
    TBranch        *b__isloose;   //!
    TBranch        *b__istight;   //!
    TBranch        *b__istightID;   //!
    TBranch        *b__closeJetPtAll;   //!
    TBranch        *b__closeJetAngAll;   //!
    TBranch        *b__ptRelAll;   //!
    TBranch        *b__closeJetPtAllMC;   //!
    TBranch        *b__closeJetPtAllstatus;   //!
    TBranch        *b__partonIdMatched;   //!
    TBranch        *b__sameParton;   //!
    TBranch        *b__n_PV;   //!
    TBranch        *b__met;   //!
    TBranch        *b__met_phi;   //!
    TBranch        *b_HT;   //!
    TBranch        *b_HT40;   //!
    TBranch        *b__genmet;   //!
    TBranch        *b__genmet_phi;   //!
    TBranch        *b__mompt;   //!
    TBranch        *b__momphi;   //!
    TBranch        *b__mometa;   //!
    TBranch        *b__mompdg;   //!
    TBranch        *b__n_bJets;   //!
    TBranch        *b__n_Jets;   //!
    TBranch        *b__bTagged;   //!
    TBranch        *b__jetE;
    TBranch        *b__jetEta;   //!
    TBranch        *b__jetPhi;   //!
    TBranch        *b__jetPt;   //!
    TBranch        *b__csv;   //!
    TBranch        *b__jetDeltaR;   //!
    TBranch        *b__jetDeltaRloose;   //!

    TBranch        *b__trigDiMuIso;
    TBranch        *b__trigMu8Ele23Iso;
    TBranch        *b__trigMu23Ele12Iso;
    TBranch        *b__trigEle23Ele12Iso;

    
 
    
    // Declaration of leaf types
    ULong64_t       _eventNb;
    ULong64_t       _runNb;
    ULong64_t       _lumiBlock;
    Double_t        _weight;
    Double_t        _genqpt;
    Int_t           _nLeptons;
    Double_t        _lPt[10];
    Double_t        _lEta[10];
    Double_t        _lPhi[10];
    Double_t        _lE[10];
    Double_t        _lPtmc[10];
    Double_t        _lEtamc[10];
    Double_t        _lPhimc[10];
    Double_t        _lEmc[10];
    Double_t        _nuPtmc[10];
    Double_t        _nuEtamc[10];
    Double_t        _nuPhimc[10];
    Double_t        _nuEmc[10];
    Double_t        _mtmc[10];
    Int_t           _nEle;
    Int_t           _nMu;
    Int_t           _nTau;
    Int_t           _flavors[10];
    Double_t        _charges[10];
    Int_t           _indeces[10];
    Double_t        _isolation[10];
    Double_t        _miniisolation[10][2];
    Bool_t          _multiisolation[10][5];
    Double_t        _ptrel[10];
    Double_t        _ptratio[10];
    Int_t           _origin[10];
    Int_t           _originReduced[10];
    Double_t        _PVchi2;
    Double_t        _PVerr[3];
    Double_t        _ipPV[10];
    Double_t        _ipPVerr[10];
    Double_t        _ipZPV[10];
    Double_t        _ipZPVerr[10];
    Double_t        _ipPVmc[10];
    Double_t        _3dIP[10];
    Double_t        _3dIPerr[10];
    Double_t        _3dIPsig[10];
    Double_t        _mt[10];
    Double_t        _mt2[10][10];
    Double_t        _mt2ll;
    Double_t        _mt2blbl;
    Bool_t          _isloose[10];
    Bool_t          _istight[10];
    Bool_t          _istightID[10];
    Double_t        _closeJetPtAll[10];
    Double_t        _closeJetAngAll[10];
    Double_t        _ptRelAll[10];
    Double_t        _closeJetPtAllMC[10];
    Double_t        _closeJetPtAllstatus[10];
    Int_t           _partonIdMatched[10];
    Bool_t          _sameParton[10];
    Int_t           _n_PV;
    Double_t        _met;
    Double_t        _met_phi;
    Double_t        HT;
    Double_t        HT40;
    Double_t        _genmet;
    Double_t        _genmet_phi;
    Double_t        _mompt[10];
    Double_t        _momphi[10];
    Double_t        _mometa[10];
    Int_t           _mompdg[10];
    Int_t           _n_bJets;
    Int_t           _n_Jets;
    Bool_t          _bTagged[20];
    Double_t        _jetE[20];
    Double_t        _jetEta[20];
    Double_t        _jetPhi[20];
    Double_t        _jetPt[20];
    Double_t        _csv[20];
    Double_t        _jetDeltaR[20][10];
    Double_t        _jetDeltaRloose[20];

    bool _trigEmulator[nLeptonsMax];
    bool _isotrigEmulator[nLeptonsMax];

    Double_t _mvaValue[nLeptonsMax];

    Bool_t          _trigDiMuIso;
    Bool_t          _trigMu8Ele23Iso;
    Bool_t          _trigMu23Ele12Iso;
    Bool_t          _trigEle23Ele12Iso;

    
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    
    const int nSamples = 28;
    
    TString fileList[nSamples] = {"ttbar_2l_13.root", 
        "dy550_13.root", "dy550_1_13.root", "dy550_2_13.root", "dy550_3_13.root", "dy550_4_13.root", 
        "dy50_13.root", "dy50_1_13.root", "dy50_2_13.root", "dy50_3_13.root", "dy50_4_13.root", 
        "wjets_13.root", "wjets_1_13.root", "wjets_2_13.root", "wjets_3_13.root", "wjets_4_13.root", "wjets_5_13.root", "wjets_6_13.root", "wjets_7_13.root",
        "ttW13.root", "ttZ13.root", "ttH13.root", "ttG13.root",
        "WZ13.root",
        //"wzjets13.root",
        "ZZ13.root", "HZZ13.root",
        "tZq13.root",
        "data_combine.root"
    };


    double xSections[nSamples] = {
        87.35,
        71310, 224.2, 37.2, 3.581, 1.124,
        6025.2, 139.4, 42.75, 5.497, 2.21,
        61526, 1347, 360, 48.9, 12.8, 5.26, 1.33, 0.03089,
        0.2043, 0.2529, 0.5085 * 0.4, 3.697,
        4.42965, 
        //5.26,
        1.256, 1.1595*0.1*0.1,
        0.0758,
        1.
    };

    Color_t colsStack[nSamples] = {18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, kGreen+3, kGreen-6, kGreen-6, kViolet+2, kOrange,  kMagenta-7,  kMagenta-7,  kMagenta-7,
        kBlack
    };
    double weights[nSamples];
    
    
    TFile *hfile[nSamples];
    TTree* inputTreeFake[nSamples];
    
    TString names[nSamples] = {"fakes","DY","dy_M-10_50_1","dy_M-10_50_2","dy_M-10_50_3","dy_M-10_50_4","dy_M-50","dy_M-50_1","dy_M-50_2","dy_M-50_3","dy_M-50_4","WJets", "WJets_1_13", "WJets_2_13", "WJets_3_13", "WJets_4_13", "WJets_5_13", "WJets_6_13", "WJets_7_13","t#bar{t}W","t#bar{t}Z/H","t#bar{t}H","X+#gamma","WZ","rare","HZZ","tZq","data"};
    

    TH1D* h_leadpt[nSamples];
    TH1D* h_2ndleadpt[nSamples];
    TH1D* h_trailpt[nSamples];
    TH1D* h_sumpt[nSamples];
    TH1D* h_mll[nSamples];
    TH1D* h_njets[nSamples];
    TH1D* h_SR[nSamples];

    TH1D* h_trilep[nSamples];
    
    
    for (int sample=0; sample!=nSamples; ++sample)  {
        weights[sample] = 0;

        TString name;

        name = Form("h_each_trilepton_%d",sample);
        h_trilep[sample] = new TH1D(name, "Tri-lepton events"+names[sample]+";trilep; events / 1", 5, -0.5, 4.5);
        h_trilep[sample]->SetLineColor(colsStack[sample]);
        h_trilep[sample]->SetFillColor(colsStack[sample]);
        h_trilep[sample]->SetFillStyle(1001);
        h_trilep[sample]->SetMarkerColor(colsStack[sample]);
        h_trilep[sample]->Sumw2();

        h_trilep[sample]->GetXaxis()->SetBinLabel(1, "Total");
        h_trilep[sample]->GetXaxis()->SetBinLabel(2, "#mu#mu#mu");
        h_trilep[sample]->GetXaxis()->SetBinLabel(3, "#mu#mue");
        h_trilep[sample]->GetXaxis()->SetBinLabel(4, "#muee");
        h_trilep[sample]->GetXaxis()->SetBinLabel(5, "eee");

        h_trilep[sample]->GetXaxis()->SetTitleSize(0.15);
        h_trilep[sample]->GetXaxis()->SetLabelSize(0.25);
        h_trilep[sample]->GetXaxis()->SetLabelOffset(0.02);


        name = Form("lead_leptons_%d",sample);
        h_leadpt[sample] = new TH1D(name, "Leading lepton p_{T} "+names[sample]+";Leading lepton p_{T} [GeV]; events / 10 GeV", 20, 0, 200);
        h_leadpt[sample]->SetLineColor(colsStack[sample]);
        h_leadpt[sample]->SetFillColor(colsStack[sample]);
        h_leadpt[sample]->SetFillStyle(1001);
        h_leadpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_leadpt[sample]->SetMarkerStyle(20+sample*5);
        h_leadpt[sample]->Sumw2();
        
        name = Form("2ndlead_leptons_%d",sample);
        h_2ndleadpt[sample] = new TH1D(name, "2nd Leading lepton p_{T} "+names[sample]+";2nd Leading lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_2ndleadpt[sample]->SetLineColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillStyle(1001);
        h_2ndleadpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_leadpt[sample]->SetMarkerStyle(20+sample*5);
        h_2ndleadpt[sample]->Sumw2();
        
        name = Form("trail_leptons_%d",sample);
        h_trailpt[sample] = new TH1D(name, "Trailing lepton p_{T} "+names[sample]+";Trailing lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_trailpt[sample]->SetLineColor(colsStack[sample]);
        h_trailpt[sample]->SetFillColor(colsStack[sample]);
        h_trailpt[sample]->SetFillStyle(1001);
        h_trailpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_trailpt[sample]->Sumw2();
        
        name = Form("sumpt_leptons_%d",sample);
        h_sumpt[sample] = new TH1D(name, "Sum lepton p_{T} "+names[sample]+";Sum lepton p_{T} [GeV]; events / 10 GeV", 40, 0, 400);
        h_sumpt[sample]->SetLineColor(colsStack[sample]);
        h_sumpt[sample]->SetFillColor(colsStack[sample]);
        h_sumpt[sample]->SetFillStyle(1001);
        h_sumpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_sumpt[sample]->Sumw2();
        
        name = Form("mll_%d",sample);
        h_mll[sample] = new TH1D(name, "Invariant mass of 2 lepton "+names[sample]+";Invariant mass of 2 lepton [GeV]; events / 10 GeV", 20, 70, 110);
        h_mll[sample]->SetLineColor(colsStack[sample]);
        h_mll[sample]->SetFillColor(colsStack[sample]);
        h_mll[sample]->SetFillStyle(1001);
        h_mll[sample]->SetMarkerColor(colsStack[sample]);
        h_mll[sample]->Sumw2();
        
        
        name = Form("h_njets_%d",sample);
        h_njets[sample] = new TH1D(name, "N_{jets} "+names[sample]+";N_{jets}; events / 1", 10, 0, 10);
        h_njets[sample]->SetLineColor(colsStack[sample]);
        h_njets[sample]->SetFillColor(colsStack[sample]);
        h_njets[sample]->SetFillStyle(1001);
        h_njets[sample]->SetMarkerColor(colsStack[sample]);
        //h_njets[sample]->SetMarkerStyle(20+sample*5);
        h_njets[sample]->Sumw2();
        
        name = Form("h_SR_%d",sample);
        h_SR[sample] = new TH1D(name, "N_{jets} "+names[sample]+";SR; events / 1", 400, 0, 400);
        h_SR[sample]->SetLineColor(colsStack[sample]);
        h_SR[sample]->SetFillColor(colsStack[sample]);
        h_SR[sample]->SetFillStyle(1001);
        h_SR[sample]->SetMarkerColor(colsStack[sample]);
        h_SR[sample]->Sumw2();
    }
    
    THStack* st_leadpt = new THStack("st_leadpt","Leading lepton p_{T}");
    THStack* st_2ndleadpt = new THStack("st_2ndleadpt","2nd Leading lepton p_{T}");
    THStack* st_trailpt = new THStack("st_trailpt","Trailing lepton p_{T}");
    THStack* st_sumpt = new THStack("st_sumpt","Sum lepton p_{T}");
    THStack* st_mll = new THStack("st_mll","Invariant mass of 2 leptons");
    THStack* st_njets = new THStack("st_njets","N_{jets}");
    THStack* st_SR = new THStack("st_SR","SR");
    THStack* st_trilep = new THStack("st_trilep","Di-lepton events");
    
    
    for (int i=nSamples-2; i!=-1; --i) {
        st_leadpt->Add(h_leadpt[i]);
        st_2ndleadpt->Add(h_2ndleadpt[i]);
        st_trailpt->Add(h_trailpt[i]);
        st_sumpt->Add(h_sumpt[i]);
        st_mll->Add(h_mll[i]);
        st_njets->Add(h_njets[i]);
        st_SR->Add(h_SR[i]);
        st_trilep->Add(h_trilep[i]);
    }
    
    const int nVars  = 17;
    TH1D* distribs[nVars][nSamples];
    THStack* distribsST[nVars];
    TString varN[nVars] = {"p_{T}^{leading} [GeV]","p_{T}^{trailing} [GeV]",
        "|#eta|^{max}","M_{ll}(OSSF)","M_{T} [GeV]","R_{0.3}","R_{mini}",
        "E_{T}^{miss}","H_{T}","N_{jets}","N_{jets40}","N_{b jets}","M_{T2}","M_{T2}^{true}","#DeltaM(OSSF,Z)", "M_{T2ll}", "M_{T2blbl}"};
    double varMin[nVars] = {0,0,
        0,70,0,0,0,
        30,0,0,0,0,0,0,0, 0, 0};
    
    double varMax[nVars] = {200,100,
        2.4,120,200,10,10,
        100,200,4,10,2,200,200,100, 200, 400
    };
    
    int nBins[nVars] = {40,20,
        12,10,20,50,50,
        7,5,4,10,2,40,40,20, 20, 40
    };
    
    for (int i=0; i!=nVars;++i) {
        TString name = Form("varST_%d",i);
        distribsST[i] = new THStack(name,varN[i]);
        for (int j=nSamples-1; j!=-1; --j) {
            name = Form("var_%d_%d",i,j);
            distribs[i][j] = new TH1D(name,name+";"+varN[i],nBins[i],varMin[i],varMax[i]);
            distribs[i][j]->SetLineColor(colsStack[j]);
            if (j < nSamples-1)
                distribs[i][j]->SetFillColor(colsStack[j]);
            distribs[i][j]->SetMarkerColor(colsStack[j]);
            distribs[i][j]->SetMarkerStyle(20);
            distribs[i][j]->SetMarkerSize(0.5);
            distribs[i][j]->SetLineWidth(1);
            distribs[i][j]->Sumw2();
            if (j < nSamples-1)
                distribsST[i]->Add(distribs[i][j]);
        }
        
    }

    
    double _multiConst[3][3];
    
    //loose
    _multiConst[0][0] = 0.22;
    _multiConst[0][1] = 0.63;
    _multiConst[0][2] = 6;
    
    //medium
    _multiConst[1][0] = 0.14;
    _multiConst[1][1] = 0.68;
    _multiConst[1][2] = 6.7;
    
    //tight
    _multiConst[2][0] = 0.10;
    _multiConst[2][1] = 0.70;
    _multiConst[2][2] = 7;
    
    int nLoc = 0;
    int nLocEle = 0;
    int nLocMu = 0;
    int nJLoc = 0;
    int nBLoc = 0;
    double HTLoc = 0;
    int leptInd[3];
    
    TLorentzVector l0p4, l1p4;
    
    int dilep_event_bkg[2] = {0};
    int trilep_event_bkg[2] = {0};
    int both_event_bkg[2] = {0};
    int total_event_bkg[2] = {0};

    int dilep_event_sig[2][6] = {0};
    int trilep_event_sig[2][6]= {0};
    int both_event_sig[2][6] = {0};
    int total_event_sig[2][6] = {0};

    ofstream outFileEventNumber;
    outFileEventNumber.open("tableEventNb.txt");
    outFileEventNumber.precision(3);

    //0 - vloose, 1 - vloose FOIDEmu, 2 - vlooseFOIDISOEMU, 3 - tight
    double valuesMVA[3][4];
    valuesMVA[0][0] = -0.16;
    valuesMVA[1][0] = -0.65;
    valuesMVA[2][0] = -0.74;

    valuesMVA[0][1] = -0.70;
    valuesMVA[1][1] = -0.83; 
    valuesMVA[2][1] = -0.92; 

    valuesMVA[0][2] = -0.155;
    valuesMVA[1][2] = -0.56;
    valuesMVA[2][2] = -0.76; 

    valuesMVA[0][3] = 0.87;
    valuesMVA[1][3] = 0.60;
    valuesMVA[2][3] = 0.17;

    double nEventsNumber[nSamples] = {0.};
 
    for (int sam = 0; sam != nSamples; ++sam) {
        //if (sam != 12) continue;
        //if (sam == 2) continue;
        //if (sam != 27) continue;
        //if ((sam > 0) && (sam < 6)) continue;
        //if ((sam > 10) && (sam < 19)) continue;
        //hfile[sam] = new TFile("data/rootfiles/"+fileList[sam],"read");
        hfile[sam] = new TFile("/Users/ikhvastu/Desktop/CERN/MCsamples/"+fileList[sam],"read");
        
        hfile[sam]->cd("FakeElectrons");
        inputTreeFake[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        
        inputTreeFake[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
        inputTreeFake[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
        inputTreeFake[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
        inputTreeFake[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
        
        inputTreeFake[sam]->SetBranchAddress("_nLeptons", &_nLeptons, &b__nLeptons);
        inputTreeFake[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
        inputTreeFake[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
        inputTreeFake[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
        inputTreeFake[sam]->SetBranchAddress("_lE", _lE, &b__lE);
        inputTreeFake[sam]->SetBranchAddress("_lPtmc", _lPtmc, &b__lPtmc);
        inputTreeFake[sam]->SetBranchAddress("_lEtamc", _lEtamc, &b__lEtamc);
        inputTreeFake[sam]->SetBranchAddress("_lPhimc", _lPhimc, &b__lPhimc);
        inputTreeFake[sam]->SetBranchAddress("_lEmc", _lEmc, &b__lEmc);
        inputTreeFake[sam]->SetBranchAddress("_nuPtmc", _nuPtmc, &b__nuPtmc);
        inputTreeFake[sam]->SetBranchAddress("_nuEtamc", _nuEtamc, &b__nuEtamc);
        inputTreeFake[sam]->SetBranchAddress("_nuPhimc", _nuPhimc, &b__nuPhimc);
        inputTreeFake[sam]->SetBranchAddress("_nuEmc", _nuEmc, &b__nuEmc);
        //inputTreeFake[sam]->SetBranchAddress("_mtmc", _mtmc, &b__mtmc);
        //inputTreeFake[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
        //inputTreeFake[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
        //inputTreeFake[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
        inputTreeFake[sam]->SetBranchAddress("_flavors", _flavors, &b__flavors);
        inputTreeFake[sam]->SetBranchAddress("_charges", _charges, &b__charges);
        //inputTreeFake[sam]->SetBranchAddress("_indeces", _indeces, &b__indeces);
        inputTreeFake[sam]->SetBranchAddress("_isolation", _isolation, &b__isolation);
        inputTreeFake[sam]->SetBranchAddress("_miniisolation", _miniisolation, &b__miniisolation);
        inputTreeFake[sam]->SetBranchAddress("_multiisolation", _multiisolation, &b__multiisolation);
        inputTreeFake[sam]->SetBranchAddress("_ptrel", _ptrel, &b__ptrel);
        inputTreeFake[sam]->SetBranchAddress("_ptratio", _ptratio, &b__ptratio);
        //inputTreeFake[sam]->SetBranchAddress("_origin", _origin, &b__origin);
        //inputTreeFake[sam]->SetBranchAddress("_originReduced", _originReduced, &b__originReduced);
        inputTreeFake[sam]->SetBranchAddress("_PVchi2", &_PVchi2, &b__PVchi2);
        inputTreeFake[sam]->SetBranchAddress("_PVerr", _PVerr, &b__PVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPV", _ipPV, &b__ipPV);
        inputTreeFake[sam]->SetBranchAddress("_ipPVerr", _ipPVerr, &b__ipPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipZPV", _ipZPV, &b__ipZPV);
        inputTreeFake[sam]->SetBranchAddress("_ipZPVerr", _ipZPVerr, &b__ipZPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPVmc", _ipPVmc, &b__ipPVmc);
        inputTreeFake[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
        inputTreeFake[sam]->SetBranchAddress("_3dIPerr", _3dIPerr, &b__3dIPerr);
        inputTreeFake[sam]->SetBranchAddress("_3dIPsig", _3dIPsig, &b__3dIPsig);
        inputTreeFake[sam]->SetBranchAddress("_mt", _mt, &b__mt);
        //inputTreeFake[sam]->SetBranchAddress("_mt2", _mt2, &b__mt2);
        //inputTreeFake[sam]->SetBranchAddress("_mt2ll", &_mt2ll, &b__mt2ll);
        //inputTreeFake[sam]->SetBranchAddress("_mt2blbl", &_mt2blbl, &b__mt2blbl);
        //inputTreeFake[sam]->SetBranchAddress("_isloose", _isloose, &b__isloose);
        inputTreeFake[sam]->SetBranchAddress("_istight", _istight, &b__istight);
        //inputTreeFake[sam]->SetBranchAddress("_istightID", _istightID, &b__istightID);
        //inputTreeFake[sam]->SetBranchAddress("_closeJetPtAll", _closeJetPtAll, &b__closeJetPtAll);
        //inputTreeFake[sam]->SetBranchAddress("_closeJetAngAll", _closeJetAngAll, &b__closeJetAngAll);
        //inputTreeFake[sam]->SetBranchAddress("_ptRelAll", _ptRelAll, &b__ptRelAll);
        //inputTreeFake[sam]->SetBranchAddress("_closeJetPtAllMC", _closeJetPtAllMC, &b__closeJetPtAllMC);
        //inputTreeFake[sam]->SetBranchAddress("_closeJetPtAllstatus", _closeJetPtAllstatus, &b__closeJetPtAllstatus);
        //inputTreeFake[sam]->SetBranchAddress("_partonIdMatched", _partonIdMatched, &b__partonIdMatched);
        //inputTreeFake[sam]->SetBranchAddress("_sameParton", _sameParton, &b__sameParton);
        inputTreeFake[sam]->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
        inputTreeFake[sam]->SetBranchAddress("_met", &_met, &b__met);
        inputTreeFake[sam]->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);
        //inputTreeFake[sam]->SetBranchAddress("HT", &HT, &b_HT);
        //inputTreeFake[sam]->SetBranchAddress("HT40", &HT40, &b_HT40);
        inputTreeFake[sam]->SetBranchAddress("_genmet", &_genmet, &b__genmet);
        inputTreeFake[sam]->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);
        //inputTreeFake[sam]->SetBranchAddress("_mompt", _mompt, &b__mompt);
        //inputTreeFake[sam]->SetBranchAddress("_momphi", _momphi, &b__momphi);
        //inputTreeFake[sam]->SetBranchAddress("_mometa", _mometa, &b__mometa);
        //inputTreeFake[sam]->SetBranchAddress("_mompdg", _mompdg, &b__mompdg);
        inputTreeFake[sam]->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
        inputTreeFake[sam]->SetBranchAddress("_n_Jets", &_n_Jets, &b__n_Jets);
        //inputTreeFake[sam]->SetBranchAddress("_bTagged", _bTagged, &b__bTagged);
        inputTreeFake[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
        inputTreeFake[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
        inputTreeFake[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
        inputTreeFake[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
        inputTreeFake[sam]->SetBranchAddress("_csv", _csv, &b__csv);
        inputTreeFake[sam]->SetBranchAddress("_jetDeltaR", _jetDeltaR, &b__jetDeltaR);

        inputTreeFake[sam]->SetBranchAddress("_mvaValue", &_mvaValue);
        //inputTreeFake[sam]->SetBranchAddress("_jetDeltaRloose", _jetDeltaRloose, &b__jetDeltaRloose);
        
        if(sam != 27){
            inputTreeFake[sam]->SetBranchAddress("_genqpt", &_genqpt, &b__genqpt);
        }
        
        inputTreeFake[sam]->SetBranchAddress("_trigDiMuIso", &_trigDiMuIso , &b__trigDiMuIso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigMu8Ele23Iso", &_trigMu8Ele23Iso , &b__trigMu8Ele23Iso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigMu23Ele12Iso", &_trigMu23Ele12Iso , &b__trigMu23Ele12Iso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigEle23Ele12Iso", &_trigEle23Ele12Iso , &b__trigEle23Ele12Iso ); // Di-muon isolation trigger
        

        _hCounter->Read("hCounter");
        Double_t scale = xSections[sam]*2110/(_hCounter->GetBinContent(1));
        
        Long64_t nEntries = inputTreeFake[sam]->GetEntries();
        std::cout<<"Entries in "<<fileList[sam]<<" "<<nEntries<<std::endl;
        std::cout<<xSections[sam]<<" "<<_hCounter->GetBinContent(1)<<" "<<scale<<std::endl;

        int allEvents = 0;
        int negEvents = 0;

        
        for (Long64_t it=0; it!=nEntries; ++it) {

            inputTreeFake[sam]->GetEntry(it);
            if(sam == 27){

                scale = 1.;
                _weight = 1.;
            }

            //if (!((_eventNb ==109887547) && (_lumiBlock ==1098876))) continue;
            //std::cout << "Ev info: " << _eventNb << " " << _lumiBlock << std::endl;
            if (it%100000 == 0)
                cout<<'.'<<flush;

            //if( it == 10000) break;
 
            if (sam == 1 && _genqpt > 100) continue;
            if (sam == 6 && _genqpt > 100) continue;
            if (sam == 11 && _genqpt > 100) continue;
            if (_nLeptons < 3) continue;
            //if (_met < 50) continue;
            //std::cout << "What is going on?" << std::endl;
            nLoc = 0;
            nLocEle = 0;
            nLocMu = 0;

            //std::cout << "Lepton Number: " << _nLeptons << std::endl;

            for (int i=0; i!=_nLeptons; ++i) {
                if (_flavors[i] > 1) continue;
                if(_lPt[i] < 10) continue;
                
                //std::cout << " Lepton info " << _lPt[i] << " " << _mvaValue[i] << " " << _flavors[i] << " " << _istight[i] << " " << _miniisolation[i][0] << " " << _ptratio[i] << " " << _ptrel[i] << " ";
                //if(_flavors[i] == 0)
                //    std::cout << _isotrigEmulator[i];
                //std::cout << std::endl;

                if(!_istight[i]) continue;

                if(TMath::Abs(_3dIPsig[i]) > 4) continue;
                
                    if(_flavors[i] == 0 && ((_miniisolation[i][0] < 0.16) && ((_ptratio[i] > 0.76) || (_ptrel[i] > 7.2)))) { 
                    
                    /*  
                    bool passedMVA = false;
                    if (TMath::Abs(_lEta[i]) < 0.8 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[0][3];
                    } else if (TMath::Abs(_lEta[i]) < 1.479 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[1][3];
                    } else {
                        passedMVA = _mvaValue[i]> valuesMVA[2][3];
                    }

                    if(!passedMVA) continue;
                    */

                    leptInd[nLoc] = i;
                    nLoc++;
                    nLocEle++;
                }
                if(_flavors[i] == 1 && ((_miniisolation[i][0] < 0.20) && ((_ptratio[i] > 0.69) || (_ptrel[i] > 6.0))) ){
                    leptInd[nLoc] = i;
                    nLoc++;
                    nLocMu++;
                }
                //}
                
            }

            //std::cout << "Lepton Numbers: " << nLoc << " " << nLocMu << " " << nLocEle << std::endl;
            
            if (nLoc != 3) continue;
            //if (nLocEle != 3) continue;
            //if (nLocMu > 1) continue;
            //if ((nLocMu != 1 && nLocEle != 2) && nLocEle != 3) continue;
            //if (nLocMu != 3) continue;

            //if (sam < 27)
            if((!_trigDiMuIso) && (!_trigMu8Ele23Iso) && (!_trigMu23Ele12Iso) && (!_trigEle23Ele12Iso)) continue; 


            //std::cout << "Pass 3 lepton cut" << std::endl;
            
            double maxPt = 0.;
            double max2ndPt = 0.;
            double minPt = 0.;

            vector<std::pair<double,int> > leptonPt;

            /*
            for(int i = 0; i < nLoc; i++)
                std::cout << "Leptons pt and index: " << _lPt[leptInd[i]] << " " << leptInd[i] << std::endl;
            */
            
            for(int i = 0; i < nLoc; i++)
                leptonPt.push_back(std::make_pair(_lPt[leptInd[i]], leptInd[i]));

            sort (leptonPt.begin(), leptonPt.end(), comp); 

            maxPt = leptonPt[0].first;
            max2ndPt = leptonPt[1].first;
            minPt = leptonPt[2].first;

            int maxPtInd = leptonPt[0].second;
            int max2ndPtInd = leptonPt[1].second;
            int minPtInd = leptonPt[2].second;


            //3e triggers
            /*
            if(nLocEle > 2){
                if(sam < 12){
                    total_event_bkg[0]++;
                    if((maxPt > 18) && (max2ndPt > 15) && (minPt > 10))            
                        trilep_event_bkg[0]++;
                    if((maxPt > 25) && (max2ndPt > 14) && (minPt > 7))            
                        dilep_event_bkg[0]++;
                    if(((maxPt > 18) && (max2ndPt > 15) && (minPt > 10)) || ((maxPt > 25) && (max2ndPt > 14) && (minPt > 7))) 
                        both_event_bkg[0]++;
                }
                else{
                    total_event_sig[0][sam-12]++;
                    if((maxPt > 18) && (max2ndPt > 15) && (minPt > 10))            
                        trilep_event_sig[0][sam-12]++;
                    if((maxPt > 25) && (max2ndPt > 14) && (minPt > 7))            
                        dilep_event_sig[0][sam-12]++;
                    if(((maxPt > 18) && (max2ndPt > 15) && (minPt > 10)) || ((maxPt > 25) && (max2ndPt > 14) && (minPt > 7))) 
                        both_event_sig[0][sam-12]++;
                }
            }

            //3mu triggers
            if(nLocEle > 2){
                if(sam < 12){
                    total_event_bkg[1]++;
                    if((maxPt > 13) && (max2ndPt > 11) && (minPt > 6))            
                        trilep_event_bkg[1]++;
                    if((maxPt > 18) && (max2ndPt > 9) && (minPt > 5))            
                        dilep_event_bkg[1]++;
                    if(((maxPt > 13) && (max2ndPt > 11) && (minPt > 6)) || ((maxPt > 18) && (max2ndPt > 9) && (minPt > 5))) 
                        both_event_bkg[1]++;
                }
                else{
                    total_event_sig[1][sam-12]++;
                    if((maxPt > 13) && (max2ndPt > 11) && (minPt > 6))            
                        trilep_event_sig[1][sam-12]++;
                    if((maxPt > 18) && (max2ndPt > 9) && (minPt > 5))            
                        dilep_event_sig[1][sam-12]++;
                    if(((maxPt > 13) && (max2ndPt > 11) && (minPt > 6)) || ((maxPt > 18) && (max2ndPt > 9) && (minPt > 5))) 
                        both_event_sig[1][sam-12]++;
                }
            }
            */


            /*
            if (maxPt < 10) continue;
            if (max2ndPt < 10) continue;
            if ((_flavors[minPtID] == 0) && (minPt < 5)) continue;
            if ((_flavors[minPtID] == 1) && (minPt < 7)) continue;
            */

            /*
            if (maxPt > 35) continue;
            if (max2ndPt > 35) continue;
            if (minPt > 35) continue;
            */

            if (maxPt < 20) continue;
            if (max2ndPt < 10) continue;
            if (minPt < 10) continue;

            int intPtBin;

            //if ((maxPt < 35) && (max2ndPt < 35) && (minPt < 35))
            //if ((maxPt < 35) && (max2ndPt < 35))
            //if (maxPt < 35)
            //if (maxPt + max2ndPt + minPt < 100)
            //    intPtBin = 0;
            //else if ((maxPt > 25) && (max2ndPt > 25) && (minPt > 25))
            //    intPtBin = 1;
            //else 
            //    intPtBin = 1;
            //if (maxPt < 20) continue;
            
            nJLoc = 0;
            nBLoc = 0;
            HTLoc = 0;
            
            for (int i=0; i!=_n_Jets; ++i) {
                //std::cout << _csv[i] << std::endl; 
                bool clean = true;

                TLorentzVector jet;
                jet.SetPtEtaPhiE(_jetPt[i], _jetEta[i], _jetPhi[i], _jetE[i]);

                //std::cout << _jetDeltaRloose[i] << std::endl;
                
                for (int j=0; j!=_nLeptons; ++j) {
                    if(TMath::Abs(_3dIPsig[j] > 4)) continue;
                    
                    TLorentzVector lep;
                    //lep.SetPtEtaPhiE(_lPt[leptInd[j]], _lEta[leptInd[j]], _lPhi[leptInd[j]], _lE[leptInd[j]]);
                    lep.SetPtEtaPhiE(_lPt[j], _lEta[j], _lPhi[j], _lE[j]);
                   
                    //std::cout << "Pt of jets: " << _jetPt[i] << " ; Delta R: " << _jetDeltaR[i][leptInd[j]] << std::endl;
                    //if (_jetDeltaR[i][leptInd[j]] < 0.4) {
                    if(jet.DeltaR(lep) < 0.4){
                            //removeLep = true;
                            //_jetPt[i]-=_lPt[leptInd[j]];
                        clean = false;
                        break;
                    } else clean = true;
                }
                
                if (clean && _jetPt[i] > 30 && TMath::Abs(_jetEta[i]) < 2.4) {
                    nJLoc++;
                    HTLoc+=_jetPt[i];
                    //std::cout << _bTagged[i] << std::endl;
                    //if (_bTagged[i])
                    /*
                    if(_csv[i] > 0.605)
                        nBLoc++;
                        */
                    if(_csv[i] > 0.89){
                        nBLoc++;
                    }
                }
            }

            //std::cout << nJLoc << " " << nBLoc << std::endl;

            
            if(_met < 30) continue;
            if(_met > 100) continue;
            if(nBLoc > 0) continue;
            if(nJLoc > 1) continue;
            //if (HTLoc < 60) continue;
            //if(nJLoc < 2) continue;
            //if (HTLoc < 80) continue;

            double deltaMZ = 999999.;
            //double deltaM = 999999.;
            double minll = 10000.;
            int index1 = -1;
            int index2 = -1;
            int index3 = leptInd[0];
            double mll_true = 9999.;
            
            for (int l0 = 0; l0<nLoc; ++l0) {
                l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);
                for (int l1 = l0+1; l1<nLoc; ++l1) {
                    if (_charges[leptInd[l0]] != _charges[leptInd[l1]]) {
                        l1p4.SetPtEtaPhiE(_lPt[leptInd[l1]],_lEta[leptInd[l1]],_lPhi[leptInd[l1]],_lE[leptInd[l1]]);
                        l1p4+=l0p4;
                        double mdiL = l1p4.M();
                        if (_flavors[leptInd[l0]] == _flavors[leptInd[l1]] ) {
                            if (mdiL < minll) minll = mdiL;
                            if (fabs(mdiL - 91) < deltaMZ) {
                                deltaMZ = fabs(mdiL - 91);
                                mll_true = mdiL;
                                index3 = leptInd[3 - l0 - l1];
                                index1 = leptInd[l0];
                                index2 = leptInd[l1];
                            }
                        }
                    }
                }
            }

            
            //std::cout << "3rd lepton flavor: " << _flavors[index3] << std::endl;

            //if(_mt[index3] < 50 ) continue;

            /*
            if((_charges[leptInd[0]] == _charges[leptInd[1]]) && (_charges[leptInd[1]] == _charges[leptInd[2]])){
                if(nLocEle > 1) {
                    for (int l0 = 0; l0<3; ++l0) {
                        if(_flavors[leptInd[l0]] == 1) continue;
                        l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);
                        for (int l1 = l0+1; l1<3; ++l1) {
                            if(_flavors[leptInd[l1]] == 1) continue;
                            l1p4.SetPtEtaPhiE(_lPt[leptInd[l1]],_lEta[leptInd[l1]],_lPhi[leptInd[l1]],_lE[leptInd[l1]]);
                            l1p4+=l0p4;
                            double mdiL = l1p4.M();
                            if (mdiL < minll) minll = mdiL;
                            if (fabs(mdiL - 91) < deltaM) {
                                deltaMZ = fabs(mdiL - 91);
                                index3 = leptInd[3 - l0 - l1];
                            }
                        }
                    }                
                }
                else if(nLocEle == 1){
                    for (int l0 = 0; l0<3; ++l0) {
                        l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);
                        for (int l1 = l0+1; l1<3; ++l1) {
                            l1p4.SetPtEtaPhiE(_lPt[leptInd[l1]],_lEta[leptInd[l1]],_lPhi[leptInd[l1]],_lE[leptInd[l1]]);
                            l1p4+=l0p4;
                            double mdiL = l1p4.M();
                            if (mdiL < minll) minll = mdiL;
                            if (fabs(mdiL - 50) < deltaM) {
                                deltaM = fabs(mdiL - 50);
                                index3 = leptInd[3 - l0 - l1];
                            }
                        }
                    }                
                }
                else{
                    double minll3mu = -999.;
                    for (int l0 = 0; l0<3; ++l0) {
                        l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);
                        for (int l1 = l0+1; l1<3; ++l1) {
                            l1p4.SetPtEtaPhiE(_lPt[leptInd[l1]],_lEta[leptInd[l1]],_lPhi[leptInd[l1]],_lE[leptInd[l1]]);
                            l1p4+=l0p4;
                            double mdiL = l1p4.M();
                            if (mdiL > minll3mu) 
                                minll3mu = mdiL;
                        }
                    }                
                    minll = minll3mu;
                }
            }
            */

            if (deltaMZ > 15) continue;
            if (minll < 12) continue;

            /*
            if(_flavors[index1] > 0 || _flavors[index2] > 0)
                std::cout << "SOMETHING WRONG" << std::endl;
                */

            //std::cout << "1st and 2nd lepton flavor: " <<  _flavors[index1] << " " << _flavors[index2] << std::endl;
            
            distribs[9][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);

            if(nJLoc < 2){
            h_trailpt[sam]->Fill(TMath::Min(minPt, 99.),scale*_weight);
            h_leadpt[sam]->Fill(TMath::Min(maxPt, 199.),scale*_weight);
            h_2ndleadpt[sam]->Fill(TMath::Min(max2ndPt, 199.),scale*_weight);
            h_sumpt[sam]->Fill(TMath::Min(minPt+maxPt+max2ndPt, 399.),scale*_weight);
            h_mll[sam]->Fill(TMath::Min(mll_true, 119.),scale*_weight);
           
            h_njets[sam]->Fill(TMath::Min(nJLoc, 9),scale*_weight);
            
            h_SR[sam]->Fill(SRID(_met,HTLoc,nJLoc,nBLoc,intPtBin,_mt2ll),scale*_weight);

            h_trilep[sam]->Fill(0.,scale*_weight);
            
            if(nLocMu == 3 && nLocEle == 0)
                h_trilep[sam]->Fill(1.,scale*_weight);
            if(nLocMu == 2 && nLocEle == 1)
                h_trilep[sam]->Fill(2.,scale*_weight);
            if(nLocMu == 1 && nLocEle == 2)
                h_trilep[sam]->Fill(3.,scale*_weight);
            if(nLocMu == 0 && nLocEle == 3)
                h_trilep[sam]->Fill(4.,scale*_weight);

            //h_SR[sam]->Fill(SRID(_met,HTLoc,nJLoc,nBLoc,intPtBin,_mt2ll),1);
            
           
            distribs[0][sam]->Fill(TMath::Min(maxPt,varMax[0]-0.1),scale*_weight);
            distribs[1][sam]->Fill(TMath::Min(minPt,varMax[1]-0.1),scale*_weight);
            distribs[2][sam]->Fill(TMath::Max(TMath::Max(fabs(_lEta[leptInd[0]]),fabs(_lEta[leptInd[1]])), fabs(_lEta[leptInd[2]])),scale*_weight);
            distribs[3][sam]->Fill(TMath::Min(mll_true, 119.),scale*_weight);
            //if (deltaM < 15)
            distribs[4][sam]->Fill(TMath::Min(_mt[index3], varMax[4]-0.1),scale*_weight);
            distribs[5][sam]->Fill(TMath::Min(TMath::Min(_isolation[0],_isolation[1]),varMax[5]-0.1),scale*_weight);
            distribs[6][sam]->Fill(TMath::Min(TMath::Min(_miniisolation[0][0],_miniisolation[1][0]),varMax[6]-0.1),scale*_weight);
            //if (deltaM < 15)
            distribs[7][sam]->Fill(TMath::Min(_met,varMax[7]-0.1),scale*_weight);
            distribs[8][sam]->Fill(TMath::Min(HTLoc,varMax[8]-0.1),scale*_weight);

            distribs[10][sam]->Fill(TMath::Min(double(_n_Jets),varMax[10]-1),scale*_weight);
            
            distribs[11][sam]->Fill(TMath::Min(double(nBLoc),varMax[11]-1),scale*_weight);
            distribs[12][sam]->Fill(TMath::Min(_mt2[0][1],varMax[12]-1),scale*_weight);
            if (_originReduced[leptInd[0]] == 0 && _originReduced[leptInd[1]] == 0)
                distribs[13][sam]->Fill(TMath::Min(_mt2[leptInd[0]][leptInd[1]],varMax[12]-1),scale*_weight);
            else if (_originReduced[leptInd[0]] == 0 && _originReduced[leptInd[2]] == 0)
            distribs[13][sam]->Fill(TMath::Min(_mt2[leptInd[0]][leptInd[2]],varMax[12]-1),scale*_weight);
            else
                distribs[13][sam]->Fill(TMath::Min(_mt2[leptInd[1]][leptInd[2]],varMax[12]-1),scale*_weight);

            distribs[14][sam]->Fill(deltaMZ,scale*_weight);
            distribs[15][sam]->Fill(TMath::Min(_mt2ll,varMax[15]-1),scale*_weight);
            distribs[16][sam]->Fill(TMath::Min(_mt2blbl,varMax[16]-1),scale*_weight);
            }

            allEvents++;
            if(_weight < 0)
                negEvents++;
            
        }
        std::cout << "Total nEv: " << allEvents << " ; negEvents: " << negEvents << std::endl;
        cout<<endl;
        outFileEventNumber << fileList[sam]<<" & " << scale*_weight << " &" << allEvents*scale*TMath::Abs(_weight) << " & "<< negEvents*scale*TMath::Abs(_weight) << "\n";
        weights[sam] = scale*TMath::Abs(_weight);
        std::cout << "Weights: " << weights[sam] << std::endl;

        nEventsNumber[sam] = (allEvents - 2 * negEvents) * scale * TMath::Abs(_weight);

        /*
        ofstream outFile;
        outFile.precision(3);

        outFile.open("tableTrig.txt");
        outFile<<"\\begin{tabular}{|c|cccc|}\\hline\\hline\n";
        outFile<<" Trigger & 3 Mu & 2Mu 1Ele & 1Mu 2Ele & 3Ele \\\\\\hline\\hline\n";
    

        for(int dec = 0; dec < 7; dec++){

            //outFile << triggerDec_str[dec] << " & " << goodevents3MuTrig[dec]  * 1./goodevents3Mu[dec] << " & " << goodevents2Mu1EleTrig[dec] * 1./goodevents2Mu1Ele[dec] << " & " << goodevents1Mu2EleTrig[dec] * 1./goodevents1Mu2Ele[dec] << " & " << goodevents3EleTrig[dec]  * 1./goodevents3Ele[dec];
            outFile << triggerDec_str[dec] << " & " << goodevents3MuTrig[dec] << " & " << goodevents2Mu1EleTrig[dec] << " & " << goodevents1Mu2EleTrig[dec] << " & " << goodevents3EleTrig[dec] ;

            outFile<<"\\\\\n";
        }

        outFile<<"\\hline\\hline\n";
        outFile<<"\\end{tabular}";

        outFile.close();
        */
    }
    cout<<endl;
    outFileEventNumber.close();
    /*
    for(int j = 0; j < 2; j++){
        std::cout << "Trilep/Dilep/Both BKG: "  << trilep_event_bkg[j] << " " << dilep_event_bkg[j] << " " << both_event_bkg[j] << ", TOTAL: " << total_event_bkg[j] << std::endl; 
        for(int i = 0; i < 6; i++)
            std::cout << "Trilep/Dilep/Both SIG" << i << ": "  << trilep_event_sig[j][i] << " " << dilep_event_sig[j][i] << " " << both_event_sig[j][i] << ", TOTAL: " << total_event_sig[j][i] << std::endl; 
        std::cout << "-----------------------------------------------------------------------------------------------------" << std::endl;
    }
    */

    // SF calculation
    double nEventTotalBKG = 0.;
    for(int i = 0; i < 28; i++){
        if(i == 23)
            continue;
        if(i == 27)
            continue;
        nEventTotalBKG += nEventsNumber[i];

    }
    double SF = (nEventsNumber[27] - nEventTotalBKG) / nEventsNumber[23];
    double SFerror = uncCalc(nEventsNumber[27], nEventTotalBKG, nEventsNumber[23]);

    std::cout << "SF = " << SF << "+- " << SFerror <<std::endl; 
    std::cout << "Purity: " << nEventsNumber[23] / (nEventsNumber[23] + nEventTotalBKG) << std::endl;
    
    std::cout<<"Done"<<std::endl;
    
    TLegend* mtleg = new TLegend(0.75,0.88,0.95,0.40);
    mtleg->SetFillColor(0);
    mtleg->SetFillStyle(0);
    mtleg->SetBorderSize(0);
    for (int i=0; i!=nSamples-1; ++i) {
        if(i > 0 && i < 19)
            continue;
        if(i == 21)
            continue;
        if(i == 25)
            continue;
        if(i == 26)
            continue;
        mtleg->AddEntry(h_leadpt[i],names[i],"f");
    }

    mtleg->AddEntry(h_leadpt[nSamples-1],names[nSamples-1],"lep");
    //mtleg->AddEntry(h_leadpt[nSamples-numSUSY2],names[nSamples-numSUSY2],"l");

    TCanvas* allPlots = new TCanvas("allPlots","allPlots",1600,1200);
    allPlots->Divide(5,4,0.001,0.001);
    for (int i=0; i!=nVars;++i) {
        allPlots->cd(i+1);
        //distribsST[i]->SetMaximum(distribsST[i]->GetMaximum() * 10);
        distribsST[i]->Draw("hist");
        if ((i == 15) || (i == 16))
            allPlots->cd(i+1)->SetLogy();
        distribs[i][nSamples-2]->Draw("hist same");
        distribs[i][nSamples-1]->Draw("axis same");
        /*
        if (i == 0){
            mtleg->AddEntry(distribs[i][nSamples-2],names[nSamples-2],"l");
            mtleg->AddEntry(distribs[i][nSamples-1],names[nSamples-1],"l");
        } 
        */
    
        mtleg->Draw("same");

    }

    double scale_num = 1.6;

    TCanvas* plot = new TCanvas("ptLep","ptLep",1300,350);
    plot->Divide(4,1);
    
    plot->cd(1);
    plot->cd(1)->SetLogy();
    showHist(plot->cd(1),distribs[8][nSamples-numSUSY1],distribs[8][nSamples-numSUSY2],distribsST[8],"","H_{T}(GeV)","Events / " + std::to_string(int((varMax[8] - varMin[8])/nBins[8])) + " GeV",scale_num, mtleg);

    plot->cd(2);
    plot->cd(2)->SetLogy();
    showHist(plot->cd(2),distribs[9][nSamples-numSUSY1],distribs[9][nSamples-numSUSY2],distribsST[9],"","N_{jets}","Events / " + std::to_string(int((varMax[9] - varMin[9])/nBins[9])),scale_num, mtleg);

    plot->cd(3);
    plot->cd(3)->SetLogy();
    showHist(plot->cd(3),distribs[3][nSamples-numSUSY1],distribs[3][nSamples-numSUSY2],distribsST[3],"","M_{ll}(GeV)","Events / " + std::to_string(int((varMax[3] - varMin[3])/nBins[3])) + " GeV",scale_num, mtleg);

    plot->cd(4);
    showHist(plot->cd(4),distribs[11][nSamples-numSUSY1],distribs[11][nSamples-numSUSY2],distribsST[11],"","N_{bjets}","Events / " + std::to_string(int((varMax[11] - varMin[11])/nBins[11])),scale_num, mtleg);

    TCanvas* plot1 = new TCanvas("ptLep1","ptLep1",975,700);
    plot1->Divide(3,2);

    plot1->cd(1);
    plot1->cd(1)->SetLogy();
    showHist(plot1->cd(1),h_leadpt[nSamples-numSUSY1],h_leadpt[nSamples-numSUSY2],st_leadpt,"","Leading lepton p_{T}(GeV)","Events / 10 GeV" ,scale_num, mtleg);
    //mtleg->Draw();
    
    plot1->cd(2);
    plot1->cd(2)->SetLogy();
    showHist(plot1->cd(2),h_2ndleadpt[nSamples-numSUSY1],h_2ndleadpt[nSamples-numSUSY2],st_2ndleadpt,"","Subleading lepton p_{T}(GeV)","Events / 10 GeV",scale_num, mtleg);

    plot1->cd(3);
    plot1->cd(3)->SetLogy();
    showHist(plot1->cd(3),h_trailpt[nSamples-numSUSY1],h_trailpt[nSamples-numSUSY2],st_trailpt,"","Trailing lepton p_{T}(GeV)","Events / 10 GeV",scale_num, mtleg);
    

    plot1->cd(4);
    plot1->cd(4)->SetLogy();
    showHist(plot1->cd(4),distribs[7][nSamples-numSUSY1],distribs[7][nSamples-numSUSY2],distribsST[7],"","#slash{E}_{T}(GeV)","Events / " + std::to_string(int((varMax[7] - varMin[7])/nBins[7])) + " GeV",scale_num, mtleg);

    
    plot1->cd(5);
    plot1->cd(5)->SetLogy();
    showHist(plot1->cd(5),h_trilep[nSamples-numSUSY1],h_trilep[nSamples-numSUSY2],st_trilep,"","","Events",scale_num, mtleg);

    
    plot1->cd(6);
    showHist(plot->cd(6),distribs[4][nSamples-numSUSY1],distribs[4][nSamples-numSUSY2],distribsST[4],"","M_{T}(GeV)","Events / " + std::to_string(int((varMax[4] - varMin[4])/nBins[4])) + " GeV",scale_num, mtleg);

    /*
    plot->cd(9)->SetLogy();
    showHist(distribs[15][nSamples-numSUSY1],distribs[15][nSamples-numSUSY2],distribsST[15],"M_{T2ll}","GeV","Events / 10 GeV",10);
    */

    //plot->cd(11);
    
    /*
    plot->cd(10)->SetLogy();
    showHist(distribs[16][nSamples-numSUSY1],distribs[16][nSamples-numSUSY2],distribsST[16],"M_{T2blbl}","GeV","Events / 10 GeV",10);
    */

    //plot->cd(12);
   
 
    std::cout<<"Done 2"<<std::endl;
    
//____________________________________________________________________________________________________

    int indSR[SRNumberAll];
    int SRNumber = 0;

    for(int i = 0; i < SRNumberAll ; i++){
        double allBkg = 1.;
        //for(int j = 0; j < nSamples - 2; j++)
        for(int j = 0; j < nSamples; j++)
            allBkg += h_SR[j]->GetBinContent(i+1);
        //allBkg += h_SR[10+signalSUSY]->GetBinContent(i+1);
        if(allBkg != 0.){
        //if(i < upperLimit){
            indSR[SRNumber] = i;
            SRNumber++; 
        }
    }

    double yields[4][SRNumberAll ] = {0.};//fakes, rare, s1, s2 
    double yieldsErrs[4][SRNumberAll ] = {0.};//fakes, rare, s1, s2
    int yInd[5] = {0,19,27,27,27};

    for (int i=0; i!=4; ++i) {
        for (int k=0; k!=SRNumber; ++k) {
            yields[i][k] = 0;
            yieldsErrs[i][k] = 0;
        }
    }

   
    
    /*
    TString names3[nSamples] = {"t$\\bar{\\text{t}}$","DY","DY1","DY2","DY3","DY4","t$\\bar{\\text{t}}$W", "t$\\bar{\\text{t}}$Z","t$\\bar{\\text{t}}$H","WZ","ZZ","H\\toZZ","T54q_315","T54q_325"};

    ofstream outFileAll;
    outFileAll.precision(3);
    
    outFileAll.open("tableAll.txt");
    outFileAll<<"\\resizebox*{1\\textwidth}{!}{\n";

    outFileAll<<"\\begin{tabular}{|c|c|ccccc|cccccc|}\\hline\\hline\n";
    outFileAll<<"SR";
    for (int i=0; i!=12; ++i) {
        outFileAll<<" & "<<names3[i];
    }
    outFileAll<<"\\\\\n";

    for (int k=0; k!=SRNumber; ++k) {
        outFileAll<<k;
        for (int i=0; i!=12; ++i) {
            outFileAll<<" & "<<h_SR[i]->GetBinContent(indSR[k]+1);
        }
        outFileAll<<"\\\\\n";

    }
    outFileAll<<"\\hline\\hline\n";
    outFileAll<<"\\end{tabular}}";
    outFileAll.close();
    */

    for (int i=0; i!=4; ++i) {
        for (int j=yInd[i]; j!=yInd[i+1]; ++j) {
            for (int k=0; k!=SRNumber; ++k) {
                //if (h_SR[j]->GetBinContent(indSR[k]+1) >=0) {
                    yields[i][k]+=h_SR[j]->GetBinContent(indSR[k]+1);
                    if (h_SR[j]->GetBinContent(indSR[k]+1) == 0)
                        yieldsErrs[i][k]+=weights[j]*weights[j];
                    else
                        yieldsErrs[i][k]+=(h_SR[j]->GetBinError(indSR[k]+1))*(h_SR[j]->GetBinError(indSR[k]+1));
                //}
            }
        }
    }

    TH1F* h_SR_yield[4];
    for (int i =0; i!=4; ++i) {
        TString name = Form("h_SR_yield_%d",i);
        h_SR_yield[i] = new TH1F(name,name,400,0,400);
        h_SR_yield[i]->Sumw2();
        for (int k=0; k!=SRNumber; ++k) {
            h_SR_yield[i]->SetBinContent(k+1,yields[i][k]);
            h_SR_yield[i]->SetBinError(k+1,sqrt(yieldsErrs[i][k]));
        }
    }
    
    std::cout<<"Done 3"<<std::endl;
//_____________________________________________________________________________________________________

    const int lptNumb = 3;
    const int bjNumb = 4;
    const int jNumb = 2;
    const int metNumb = 3;
    const int htNumb = 3;
    const int mt2llNumb = 2;

    double METtest[3] = {100., 250., 550.};
    double HTtest[3] = {100., 400., 800.};
    int NJtest[3] = {2, 4, 5};
    int NBJtest[3] = {1, 2, 3};
    double MT2lltest[2] = {50., 150.};
    

    //string lptl[lptNumb] = {"1 lepton $pt < 35 $", "other"};
    string lptl[lptNumb] = {"3 leptons $pt < 35 $","3 leptons $pt > 35 $", "other"};
    //string nbjl[bjNumb] = {"1 b-tags","2 b-tags","$\\geq 3 b-tags$"};
    string nbjl[bjNumb] = {"0 b-tags","1 b-tags","2 b-tags","$\\geq 3 b-tags$"};
    //string njl[jNumb] = {"2-3 jets", "4 jets", "$\\geq$ 5 jets"};
    string njl[jNumb] = {"2-4 jets", "$\\geq$ 5 jets"};
    string metl[metNumb] = {"$ 50 < E_{T}^{miss} < 150 $ ","$ 150 < E_{T}^{miss} < 300 $ ", " $ E_{T}^{miss} > 300 $ "};
    string htl[htNumb] = {"$ 60 < H_{T} < 400 $ ","$ 400 < H_{T} < 600 $ ", " $ H_{T} > 600 $ "};
    string mt2lll[mt2llNumb] = {"$ m_{T2ll} < 120 GeV$ ", "$ m_{T2ll} > 120 GeV$ "};

    

    ofstream tableBkg;
    tableBkg.open("tableBkg.txt");
    //tableBkg << " \\documentclass{article}\n " ;
    //tableBkg << " \\usepackage{multirow}\n " ;
    //tableBkg << " \\usepackage{adjustbox}\n " ;
    //tableBkg << " \\begin{document}\n " ;
    tableBkg<<"\\begin{table}\n";
    tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    tableBkg<<"\\begin{tabular}{ccc|c|c||c|c||c|c||c|c||}\\cline{4-11}\n"; // for 3 bjet regions, which uses in offZ
    tableBkg << std::fixed << setprecision(2) << "\n";
    tableBkg << " & & " << "&\\multicolumn{2}{|c||}{" << nbjl[0] << "}&\\multicolumn{2}{|c||}{" << nbjl[1] << "}&\\multicolumn{2}{|c||}{" << nbjl[2] << "}&\\multicolumn{2}{|c||}{" << nbjl[3] << "}\\\\ \\cline{4-11}\n"; // for offZ
    tableBkg << " & & " << "&\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c||}{signal(1200)}" << " & " << "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c||}{signal(1200)}" << "  & " << "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c||}{signal(1200)}" << " & " << "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c||}{signal(1200)}" << "  \\\\ \\hline\n"; // for offZ

    int num = 0;

    for(int l = 0; l < metNumb; l++){
        tableBkg << "\\multicolumn{1}{|c}{\\multirow{ " << jNumb * htNumb << "}{*}{ " << metl[l] << "}} ";

        for(int j = 0; j < htNumb; j++){
            if (j != 0)
                tableBkg << "\\multicolumn{1}{|c}{}";
            tableBkg << " & " << " \\multicolumn{1}{|c}{\\multirow{"  << jNumb << "}{*}{ " << htl[j] << " }} ";

            for(int i = 0; i < jNumb; i++){
                if (i != 0)
                    tableBkg << "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                tableBkg << " & \\multicolumn{1}{|c|}{" << njl[i] << "} ";

                    for(int ii = 0; ii < bjNumb; ii++){
                        double err = sqrt(yieldsErrs[0][num] + yieldsErrs[1][num]) ;
                        tableBkg << " & $" << yields[0][num] + yields[1][num] << "\\pm" << err << "$ & $" << yields[2][num]  << "\\pm" << sqrt(yieldsErrs[2][num]) << "$  " ; 
                        num++;
                    }

                tableBkg << "\\\\ \\cline{3-11}\n"; // offZ 
            }

            tableBkg << " \\cline{2-11}\n "; //offZ
        }

        tableBkg << "\\hline \\hline \n";
    }


    tableBkg <<"\\end{tabular}\n";
    tableBkg <<"\\end{adjustbox}\n";
    tableBkg <<"\\end{table}\n";
    //tableBkg <<"\\end{document}";
    tableBkg .close();
    

    // table for split in mt2ll
    /*
    ofstream tableBkg;
    tableBkg.open("tableBkg.txt");
    //tableBkg << " \\documentclass{article}\n " ;
    //tableBkg << " \\usepackage{multirow}\n " ;
    //tableBkg << " \\usepackage{adjustbox}\n " ;
    //tableBkg << " \\begin{document}\n " ;
    tableBkg<<"\\begin{table}\n";
    tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    tableBkg<<"\\begin{tabular}{cccc|c|c|c||c|c|c||c|c|c||}\\cline{5-13}\n"; // for 3 bjet regions, which uses in offZ
    tableBkg << setprecision(3) << "\n";
    tableBkg << " & & &" << "&\\multicolumn{3}{|c||}{" << nbjl[0] << "}&\\multicolumn{3}{|c||}{" << nbjl[1] << "}&\\multicolumn{3}{|c||}{" << nbjl[2] << "}\\\\ \\cline{5-13}\n"; // for offZ
    tableBkg << " & & &" <<  "&\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}" <<  "& " << "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}"  <<  "& " << "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}" << "\\\\ \\hline\n"; // for offZ

    int num = 0;

    for(int l = 0; l < metNumb; l++){
        tableBkg << "\\multicolumn{1}{|c}{\\multirow{ " << jNumb * htNumb * mt2llNumb << "}{*}{ " << metl[l] << "}} ";

        for(int j = 0; j < htNumb; j++){
            if (j != 0)
                tableBkg << "\\multicolumn{1}{|c}{}";
            tableBkg << " & " << " \\multicolumn{1}{|c}{\\multirow{"  << jNumb * mt2llNumb << "}{*}{ " << htl[j] << " }} ";

            for(int i = 0; i < jNumb; i++){
                if (i != 0)
                    tableBkg << "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                tableBkg << " & \\multicolumn{1}{|c|}{\\multirow{" << mt2llNumb << "}{*}{ " << njl[i] << "}} ";

                    for(int jj = 0; jj < mt2llNumb; jj++){
                        if (jj != 0)
                            tableBkg << "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                        tableBkg << " & \\multicolumn{1}{|c|}{" << mt2lll[jj] << "} ";

                        for(int ii = 0; ii < bjNumb; ii++){
                            double err = sqrt(yieldsErrs[0][num] + yieldsErrs[1][num]) ;
                            tableBkg << " & $" << yields[0][num] + yields[1][num] << "\\pm" << err << "$ & $" << yields[2][num]  << "\\pm" << sqrt(yieldsErrs[2][num]) << "$ & $" << yields[3][num]  << "\\pm" << sqrt(yieldsErrs[3][num]) << "$ " ; 
                            num++;
                        }

                        tableBkg << "\\\\ \\cline{4-13}\n"; // offZ 
                    }

                tableBkg << "\\cline{3-13}\n"; // offZ 
            }

            tableBkg << " \\cline{2-13}\n "; //offZ
        }

        tableBkg << "\\hline \\hline \n";
    }


    tableBkg <<"\\end{tabular}\n";
    tableBkg <<"\\end{adjustbox}\n";
    tableBkg <<"\\end{table}\n";
    //tableBkg <<"\\end{document}";
    tableBkg .close();

    // table for SR
    std::map<int,int> SR;
    std::map<int,int>::iterator it_SR;
    int nSR = 1;

    SR.clear();
    ofstream table;
    table.open("tableSR.txt");
    table<<"\\begin{table}\n";
    table<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    table<<"\\begin{tabular}{cccc|c|c|c||}\\cline{5-7}\n"; // for 3 bjet regions, which uses in offZ
    table<< setprecision(3) << "\n";
    table<< " & & &" << "& "<< nbjl[0] << "&" << nbjl[1] << "&" << nbjl[2] << "\\\\ \\cline{5-7}\n"; // for offZ

    for(int l = 0; l < metNumb; l++){
        table<< "\\multicolumn{1}{|c}{\\multirow{ " << jNumb * htNumb * mt2llNumb << "}{*}{ " << metl[l] << "}} ";

        for(int j = 0; j < htNumb; j++){
            if (j != 0)
                table<< "\\multicolumn{1}{|c}{}";
            table<< " & " << " \\multicolumn{1}{|c}{\\multirow{"  << jNumb * mt2llNumb << "}{*}{ " << htl[j] << " }} ";

            for(int i = 0; i < jNumb; i++){
                if (i != 0)
                    table<< "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                table<< " & \\multicolumn{1}{|c|}{\\multirow{" << mt2llNumb << "}{*}{ " << njl[i] << "}} ";

                    for(int jj = 0; jj < mt2llNumb; jj++){
                        if (jj != 0)
                            table<< "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                        table<< " & \\multicolumn{1}{|c|}{" << mt2lll[jj] << "} ";

                        for(int ii = 0; ii < bjNumb; ii++){
                            int temp = SRID(METtest[l], HTtest[j], NJtest[i], NBJtest[ii], 0, MT2lltest[jj]);
                            it_SR = SR.find(temp);
                            if (it_SR == SR.end()){
                                SR.insert(std::pair<int,int>(temp,nSR));
                                    nSR++;
                            }
                            table<< " & SR" << SR[temp];
                        }

                        table<< "\\\\ \\cline{4-7}\n"; // offZ 
                    }

                table<< "\\cline{3-7}\n"; // offZ 
            }

            table<< " \\cline{2-7}\n "; //offZ
        }

        table<< "\\hline \\hline \n";
    }


    table<<"\\end{tabular}\n";
    table<<"\\end{adjustbox}\n";
    table<<"\\end{table}\n";
    table.close();
    */

    //table for spliting in lepton PT
/*
    ofstream tableBkg;
    tableBkg.open("tableBkg.tex");
    tableBkg << " \\documentclass{article}\n " ;
    tableBkg << " \\usepackage{multirow}\n " ;
    tableBkg << " \\usepackage{adjustbox}\n " ;
    tableBkg << " \\begin{document}\n " ;
    tableBkg<<"\\begin{table}\n";
    tableBkg<<"\\begin{adjustbox}{width=1\\textwidth}\n";
    tableBkg<<"\\begin{tabular}{cccc|c|c|c||c|c|c||c|c|c||}\\cline{5-13}\n"; // for 3 bjet regions, which uses in offZ
    tableBkg << setprecision(3) << "\n";
    tableBkg << " & & & &" << "\\multicolumn{3}{|c||}{" << nbjl[0] << "} & \\multicolumn{3}{|c||}{" << nbjl[1] << "} & \\multicolumn{3}{|c||}{" << nbjl[2]  << "}\\\\ \\cline{5-13}\n"; // for offZ
    tableBkg << " & & & &"  <<  "\\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}" <<  " & \\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}"  <<  " & \\multicolumn{1}{|c}{bkg}" << "  & " <<  "\\multicolumn{1}{|c}{signal(1200)}" << "  & "<<  "\\multicolumn{1}{|c||}{signal(1500)}" << "\\\\ \\hline\n"; // for offZ

    int num = 0;

    for(int l = 0; l < metNumb; l++){
        tableBkg << "\\multicolumn{1}{|c}{\\multirow{ " << jNumb * htNumb * lptNumb << "}{*}{ " << metl[l] << "}} ";

        for(int j = 0; j < htNumb; j++){
            if (j != 0)
                tableBkg << "\\multicolumn{1}{|c}{}";
            tableBkg << " & " << " \\multicolumn{1}{|c}{\\multirow{"  << jNumb * lptNumb  << "}{*}{ " << htl[j] << " }} ";

            for(int i = 0; i < jNumb; i++){
                if (i != 0)
                    tableBkg << "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                tableBkg << " & \\multicolumn{1}{|c|}{\\multirow{"  << lptNumb << "}{*}{ " << njl[i] << "}} ";

                for(int k = 0; k < lptNumb; k++){
                    if (k != 0)
                        tableBkg << "\\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{} & \\multicolumn{1}{|c}{}";
                    tableBkg << " & \\multicolumn{1}{|c|}{" << lptl[k] << "} ";
                        

                    for(int ii = 0; ii < bjNumb; ii++){
                        double err = sqrt(yieldsErrs[0][num] + yieldsErrs[1][num]) ;
                        tableBkg << " & $" << yields[0][num] + yields[1][num] << "\\pm" << err << "$ & $" << yields[2][num]  << "\\pm" << sqrt(yieldsErrs[2][num]) << "$ & $" << yields[3][num]  << "\\pm" << sqrt(yieldsErrs[3][num]) << "$ " ; 
                        num++;
                    }

                    tableBkg << "\\\\ \\cline{4-13}\n"; // offZ 

                }

                tableBkg << "\\cline{3-13}\n"; // offZ 
            }

            tableBkg << " \\cline{2-13}\n "; //offZ
        }

        tableBkg << "\\hline \\hline \n";
    }


    tableBkg <<"\\end{tabular}\n";
    tableBkg <<"\\end{adjustbox}\n";
    tableBkg <<"\\end{table}\n";
    tableBkg <<"\\end{document}";
    tableBkg .close();
    */
    
//___________________________________________________________________________________________________________________

    /*
    ofstream outFile;
    outFile.precision(3);

    outFile.open("table.txt");
    outFile<<"\\begin{tabular}{|c|cc|c|ccc|}\\hline\\hline\n";
    outFile<<"SR & Non-prompt & MC & Total & T14t (1.2) & T14t (1.5)& T5WZ (1.5) \\\\\\hline\\hline\n";
    
    for (int k=0; k!=SRNumber; ++k) {
        outFile<<k;
        for (int i=0; i!=2; ++i) {
            outFile<<" & "<<yields[i][k];
        }
        outFile<<" & "<<yields[0][k]+yields[1][k];
        for (int i=2; i!=5; ++i) {
            outFile<<" & "<<yields[i][k];
        }
        outFile<<"\\\\\n";
    }
    outFile<<"\\hline\\hline\n";
    outFile<<"\\end{tabular}";
    outFile.close();
    */
    
    std::cout<<"Done 4"<<std::endl;

//_____________________________________________________________________________________________________    

    // DATACARD creation
    //gSystem->Exec("rm datacard.txt"); // delete previous tex file
    ofstream datacard;
    datacard.open("datacard.txt"); // create a new datacard file
    datacard << fixed << showpoint << setprecision(4);
    datacard << "Date: May, 26, 2015" << endl;
    datacard << "Description: SUSY search, first probe" << endl;
    datacard << "lumi 10 pb^-1" << endl;
    const int Nback = 2;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "imax " <<  SRNumber << " number of channels" << endl;
    datacard << "jmax " <<  Nback << " number of backgrounds" << endl;
    datacard << "kmax " <<  (Nback + 1) * SRNumber + 3 << " number of nuisance parameters" << endl;
    datacard << "shapes * * FAKE" << endl;
    //datacard << "kmax " << 3 << " number of nuisance parameters" << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    int iBin = 1;
    datacard << "Observation ";
    for (int i = 0; i < SRNumber; i++) {
        //datacard << 0 << " ";
        if(h_SR_yield[0]->GetBinContent(iBin) +  h_SR_yield[1]->GetBinContent(iBin) <= 0)
            datacard << 0.0001 << " ";
        else
            datacard << h_SR_yield[0]->GetBinContent(iBin) +  h_SR_yield[1]->GetBinContent(iBin) << " ";
        iBin++;
    }

    datacard << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "bin ";
    for (int i = 0; i < SRNumber ; i++)
        datacard << i + 1 << " " << i + 1 << " " << i + 1 << " " ;
    datacard << endl;
    
    datacard << "process " ;
    for (int i = 0; i < SRNumber ; i++)
        datacard << "S B1 B2 ";
    datacard << endl;
    
    datacard << "process ";
    for (int i = 0; i < SRNumber ; i++)
        datacard << "0 1 2 ";
    datacard << endl;
    
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "rate ";
    iBin = 1;
    for (int i = 0; i < SRNumber; i++) {
        
        if ( h_SR_yield[signalSUSY]->GetBinContent(iBin) <= 0.)
            datacard << 0.0001 << " ";
        else
            datacard << h_SR_yield[signalSUSY]->GetBinContent(iBin) << " ";
        
        for (int numB = 0; numB < 2; numB++){
            if (h_SR_yield[numB]->GetBinContent(iBin) <= 0.)
                datacard << 0.0001 << " ";
            else
                datacard << h_SR_yield[numB]->GetBinContent(iBin) << " ";
        }
        iBin++;
    }
    datacard << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    
    iBin = 1;
    
    for (int i = 0; i < SRNumber; i++) {
                    
        datacard << (Nback + 1) * i+1 << " lnN ";
                    
        for (int num = 0; num < (Nback + 1) * i; num++)
        datacard << "1.00 ";
                    
        if (h_SR_yield[signalSUSY]->GetBinContent(iBin) <= 0)
            datacard << "1.0";
        else
            datacard << 1 +  TMath::Abs(h_SR_yield[signalSUSY]->GetBinError(iBin)/h_SR_yield[signalSUSY]->GetBinContent(iBin));
                    
        for (int num = (Nback + 1) * i + 1; num < (Nback + 1) * SRNumber; num++)
            datacard << " 1.00";
        datacard << endl;
                    
        //_________________________________________________________________________
                    
        datacard << (Nback + 1) * i + 2 << " lnN ";
                    
        for (int num = 0; num < (Nback + 1) * i + 1; num++)
             datacard << "1.00 ";
                    
        if ( h_SR_yield[0]->GetBinContent(iBin) <= 0)
             datacard << "1.0";
        else
             datacard << 1 +  TMath::Abs(1/sqrt(10 * h_SR_yield[0]->GetBinContent(iBin))) ;
                    
        for (int num = (Nback + 1) * i + 2; num < (Nback + 1) * SRNumber; num++)
            datacard << " 1.00";
        datacard << endl;
                    
        //_________________________________________________________________________
                    
        datacard << (Nback + 1) * i + 3 << " lnN ";
                    
        for (int num = 0; num < (Nback + 1) * i + 2; num++)
            datacard << "1.00 ";
                    
        if (h_SR_yield[1]->GetBinContent(iBin) <= 0)
             datacard << "1.0";
        else
             datacard << 1  + TMath::Abs(h_SR_yield[1]->GetBinError(iBin)/h_SR_yield[1]->GetBinContent(iBin)) ;
                    
        for (int num =  (Nback + 1) * i + 3; num < (Nback + 1) * SRNumber; num++)
             datacard << " 1.00";
                    
        datacard << endl;
                    
        iBin++;
    }
    
    //_________________________________________________________________________
    
    datacard << (Nback + 1) * SRNumber + 1 << " lnN "; // without split in MT2
    
    for (int num = 0; num < SRNumber ; num++)
        datacard << "1.00 1.50 1.00 ";
    
    datacard << endl;
    
    datacard << (Nback + 1) * SRNumber + 2 << " lnN "; // without split in MT2
    
    for (int num = 0; num < SRNumber; num++)
        datacard << "1.00 1.00 1.20 ";
    
    datacard << endl;
    
    datacard << (Nback + 1) * SRNumber + 3 << " lnN "; // without split in MT2
    
    for (int num = 0; num < SRNumber; num++)
        datacard << "1.10 1.00 1.00 ";
    
    datacard << endl;
    
    datacard.close();
    
    std::cout << "datacard DONE" << std::endl;

    
    //h_SR_yield[0]  - TH1F, which has 30 bins, - non prompt backgrund
    //h_SR_yield[1] - prompt background
    //h_SR_yield[signal] - SUSY signal

    
    //return 0;
    
}


void showHist(TVirtualPad* c1, TH1D *hist, TH1D *hist2, THStack *stack, string title, string titleX, string titleY, double num, TLegend *leg){   

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();
    //pad1->DrawFrame(xmin, -0.1, xmax, 1.1);
    
    hist->SetTitle(title.c_str());
    //hist->GetXaxis()->SetTitle(titleX.c_str());
    hist->GetYaxis()->SetTitle(titleY.c_str());
    hist->SetMaximum(hist->GetMaximum() * num);
    hist->SetLineWidth(1);
    hist->SetFillColor(0);
    hist->SetLineColor(1);
    hist->SetLineStyle(1);
    hist->SetMarkerSize(0.5);
    //hist->GetXaxis()->SetLabelSize(0.15);
    //hist->GetXaxis()->SetTitleSize(0.08);
    hist->GetYaxis()->SetLabelSize(0.08);
    hist->GetYaxis()->SetTitleSize(0.08);
    hist->GetYaxis()->SetTitleOffset(1.);
    hist->SetMarkerStyle(1);
    hist->Draw("e");
    stack->Draw("histsame");
    hist->Draw("esame");
    pad1->cd();
    pad1->Update();
    pad1->RedrawAxis();
    pad1->GetFrame()->Draw();
    leg->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);  

    latex.SetTextFont(42);
    latex.SetTextAlign(31); 
    latex.SetTextSize(0.07);    
    latex.DrawLatex(0.98, 0.93,lumi_13TeV + "(13 TeV)");

    TLatex cmsText;
    cmsText.SetNDC();
    cmsText.SetTextAngle(0);
    cmsText.SetTextColor(kBlack);  

    cmsText.SetTextFont(61);
    cmsText.SetTextAlign(31); 
    cmsText.SetTextSize(0.08);  
    cmsText.DrawLatex(0.25, 0.93,"CMS");


    TLatex extraText;
    extraText.SetNDC();
    extraText.SetTextAngle(0);
    extraText.SetTextColor(kBlack);  

    extraText.SetTextFont(52);
    extraText.SetTextAlign(31); 
    extraText.SetTextSize(0.07);  
    extraText.DrawLatex(0.5, 0.93,"Preliminary");

    

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();
    TH1D * histcopy = (TH1D*)hist->Clone("histcopy");
    //TH1D * histo_stack = (TH1D*)stack->GetHistogram();
    //histcopy->Sumw2();
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    histcopy->SetStats(0);

    histcopy->SetTitle("");
    histcopy->GetXaxis()->SetTitle(titleX.c_str());
    histcopy->GetYaxis()->SetTitle("data/pred");

    //histcopy->GetXaxis()->SetLabelSize(0.15);
    histcopy->GetYaxis()->SetLabelSize(0.12);
    //histcopy->GetXaxis()->SetTitleSize(0.15);
    histcopy->GetYaxis()->SetTitleSize(0.15);

    histcopy->GetYaxis()->SetTitleOffset(0.4);

    histcopy->SetMaximum(1.5);
    histcopy->SetMinimum(0.5);
    //
    histcopy->SetMarkerStyle(1);

    TList *histos = stack->GetHists();
    TIter next(histos);
    TH1D *histo_stack = (TH1D*)stack->GetHistogram();
    while ((hist = (TH1D*)next())) {
      //cout << "Adding " << hist->GetName() << endl;
      histo_stack->Add(hist);
    }

    histcopy->Divide(histo_stack);

    histcopy->Draw("e");
    line->Draw("same");
    /*
    hist2->SetLineWidth(2.);
    hist2->SetFillColor(0);
    hist2->SetLineColor(1);
    hist2->SetLineStyle(3);
    hist2->Draw("hist same");
    hist2->Draw("axis same");
    */

}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);

    readTree7414();

    rootapp->Run();

    return 0;
}

double uncCalc(double a, double b, double c){

    return TMath::Sqrt(TMath::Power((TMath::Sqrt(a) + TMath::Sqrt(b))/c,2) + TMath::Power((a - b)/c/c*TMath::Sqrt(c),2));

}
