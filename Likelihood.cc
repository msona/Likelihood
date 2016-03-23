/*
Information zu Likelihood.cc:

Das Programm liest Events im Format *.csv ein, schreibt die Daten als tree in eine *.root Datei
und schreibt in eine weitere *.root Datei (tree1.root) Histogramme der einzelnen Variablen wie Masse M,
Transversalimpuls pt usw.

Dann werden Cuts in pt durchgef�hrt und damit besondere Events wie Z0, Charmonium, Bottonium selektiert. Die selektierten
Events werden in einem Stack-Plot mit dem Background dargestellt (der Stack-Plot wird ebenfalls in tree1.root gespeichert).

Im vorliegenden Code ist die *.csv Datei "dielectron100k.csv" zum Einlesen hardgecoded, diese befindet sich im gleichen Verzeichnis
wie die ausf�hrbare Datei "Likelihood".

In einem weiteren Schritt wurde die Maximum-Likelihood-Methode (ab Zeile 503) implementiert, sowie eine Gittersuchmethode (also ein Minimierer)
zum Maximieren der Log-Likelihoodfunktion. Damit wird dann ein Fit der Gau�verteilung an den Z0-Peak 
durchgef�hrt. Im Detail wird die Variable "M" als Stichproben-Zufallsvariable verwendet, dabei werden
nur die Events mit (Pt1 und Pt2 >20GeV und M >75GeV verwendet.

Das Ergebnis des Fits wird in Form der Fitparameter im Terminal angezeigt und als Plot in der Datei "CANVASpicture.jpg" abgespeichert.

Der selbst geschriebene Minimierer wird dann mit dem Minimierer Minuit2 verglichen, dabei wird wieder die LogLikelihoodfunktion minimiert.
 
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "vector"
#include "iostream"
#include "iomanip"
#include "cmath"
#include "cassert"
#include "Cintex/Cintex.h"


//neue includes f�r Einlesen von *.csv Dateien:
#include "Riostream.h"
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TSystem.h"   // Compilerfehler ist weg sobald diese Datei included wird!!!
#include "TTree.h"
#include "THStack.h"   // Histogramme aufeinanderlegen (z.b. Background und dann Signal darauf)
#include "TF2.h"


//#include "Math/Minimizer.h"
#include "Math/Functor.h"   // aus C-Funktion einen Functor machen f�r Root-Minimierer

#include "TRandom2.h"
#include "Math/Factory.h"
#include "TError.h"
// Root-Minimierer:
#include "TMinuit.h" 
#include "TMinuitMinimizer.h" 
#include "Minuit2/Minuit2Minimizer.h"
#include "TROOT.h"
#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"

// f�r Breit-Wigner-Verteilung:
#include "Math/SpecFuncMathCore.h"
#include "TStyle.h"

#include "TStopwatch.h"

#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "TError.h"


// f�r Zufallszahlen
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TStopwatch.h>

using namespace std;

Double_t gauss(Double_t *x,Double_t *par);
Double_t fitf(Double_t *x,Double_t *par);
Double_t fitf2(Double_t *x,Double_t *par);
Double_t funktion1(Double_t *x,Double_t *par);

double RosenBrock(const double *xx );
double Rastrigin(const Double_t *x,const Double_t *par);
double Rastrigin12(const Double_t *x,const Double_t *par);
double Rastrigin01(const Double_t *x);

Double_t likelihood(Double_t *x,Double_t *par);
Double_t likelihood1d(const Double_t *x,const Double_t *par);
Double_t likelihoodxy(const Double_t *x);

int NumericalMinimization();
int NumericalMinimization2();
int NumericalMinimization3();
int NumericalMinimization4();


Double_t breitwigner(Double_t* x, Double_t* par);
Double_t breitwigner2(const Double_t* x,const Double_t* par);
Double_t breitwignerRoot(const Double_t *x, const Double_t *par);
Double_t breitwignerRoot2(const Double_t *x, const Double_t *par);

Double_t breitwignerRootGaus(const Double_t *x, const Double_t *par);
Double_t breitgausfun(Double_t *x, Double_t *par);


Double_t likelihoodBreitWigner(const Double_t *x);
Double_t likelihoodBreitWignerMinusDraw(const Double_t *x,const Double_t *par);

// Breit-Wigner-Funktion im Format f�r Minuit2  (par[0] und par[1] )
Double_t likelihoodBreitWignerMinus(const Double_t *x);
Double_t breitwignerMinuit2(const Double_t *x, const Double_t *par);

//globale Variablen (damit bei Funktionsaufrufen nicht weitere Parameter�bergaben n�tig sind)
int AnzahlEvents=0;
Double_t xwerte[7000];
Double_t xwert [1] = {0.0};

// Beginn der Main-Funktion
int main(int argc, char **argv)
{

  if(argc!=2)
  {
    cout << "program needs one argument" << endl;
    return 1;
  }

  TStopwatch sw;
  ROOT::Cintex::Cintex::Enable();

  TString pathName;

  pathName="/home/nebus/blatt1/";    // arbeite lokal auf x230 Laptop
  pathName+=argv[1];

  
  //# Einlesen der *.csv Datei Zee.csv mit TFile->Write und TTree->Readfile #########################################################
  
  
   TString dir = gSystem->UnixPathName(__FILE__);
   cout << "   Dateiname zum Einlesen : " << gSystem->UnixPathName(__FILE__) << "   ";
   //pathName.ReplaceAll("UScsvToRoot.C","");
   //pathName.ReplaceAll("/./","/");

   //Einlesen und konvertieren von *.csv Datei nach *.root
   TFile *f = new TFile("Ergebnisse/dielectron100k.root","RECREATE");
   TTree *tree = new TTree("ntuple","Daten aus einer *.csv Datei");
   //tree->ReadFile("dielectron100k.csv","Type/C:Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');
   tree->ReadFile("Daten/dielectron100k.csv","Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');
   
   f->Write();
  
   double_t px1, py1, pz1,pt1,M, E1,E2,pt2;
   
  TH1F *hM   = new TH1F("hM","M Masse distribution",100,-3,3);
  TH2F *hpxpy = new TH2F("hpxpy","py vs px",30,-3,3,30,-3,3);
  
  
  hpxpy->GetXaxis()->SetTitle(" x-Achse hat den Titel: px1 in GeV");
  hpxpy->GetYaxis()->SetTitle(" y-Achse hat den Titel: py1 in GeV");
  
  
  Int_t nentries = (Int_t)tree->GetEntries(); 
   
  tree->SetBranchAddress("pt1",&pt1);   //Setze die Einzulesende Variable pt1 f�r tree->GetEntry(i) 
 for (Int_t i = 0; i<nentries; i++) {
 tree->GetEntry(i); // Code l�uft durch Compiler mit tree->GetEntry
 //cout << " pt1: " << pt1 << " px1: "<<px1;
 }
   
 tree->SetBranchAddress("px1",&px1);
  tree->SetBranchAddress("py1",&py1);
   tree->SetBranchAddress("pz1",&pz1);
    tree->SetBranchAddress("M",&M);
  for (Int_t i = 0; i<nentries; i++) {
 tree->GetEntry(i); // Code l�uft durch Compiler mit tree->GetEntry  ; Die Werte von pt1 stimmen mit den echten aus den Daten �berein!
 
    hM->Fill(M);
    hpxpy->Fill(px1,py1);
    
 //cout << " M: " << M;
 }  
 
 
 
 // Kopieren der Variable M aus dem Tree tree nach Tree "t1" mit Variable "M1"  funktioniert, die neue Datei heisst tree1.root
 // Bauen eines neuen Trees t1 mit neuen Variablen, die aus den Variablen des alten Trees aus Zmumu.root zusammengesetzt sind
 double_t phi1,phi2,Ht;
 
 tree->SetBranchAddress("M",&M);
 tree->SetBranchAddress("E1",&E1);
 tree->SetBranchAddress("E2",&E2);
 tree->SetBranchAddress("phi1",&phi1);
 tree->SetBranchAddress("phi2",&phi2);
 
 tree->SetBranchAddress("pt1",&pt1);
 tree->SetBranchAddress("pt2",&pt2);
 
 TFile f1("Ergebnisse/tree1.root","recreate");
 TTree t1("t1","a simple Tree with simple variables");
 
 TH1F *hHt   = new TH1F("hHt"," Transversale Impulssumme Ht=pt1+pt2",500,0,310);
   hHt->GetXaxis()->SetTitle(" Transversale Impulssumme Ht = pt1+pt2 in GeV");
   hHt->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hHt->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hHt->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hHt->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hHt->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hHt->SetLineColor(kRed);
   hHt->SetFillColor(kRed);
 
  TH1F *hsumE1E2   = new TH1F("hsumE1E2"," Summe E1+E2",500,0,310);
   hsumE1E2->GetXaxis()->SetTitle(" Summe der Teilchenenergien E1 und E2 in GeV");
   hsumE1E2->GetYaxis()->SetTitle(" Anzahl Ereignisse"); 
   
  TH1F *hPt1berechnet   = new TH1F("hPt1berechnet"," Pt1=Wurzel(px1^2+py1^2)",400,0,120);
   hPt1berechnet->GetXaxis()->SetTitle(" Pt1 in GeV");
   hPt1berechnet->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt1berechnet->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1berechnet->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1berechnet->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1berechnet->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1berechnet->SetLineColor(4);
   hPt1berechnet->SetFillColor(4);
   
  TH1F *hPt1   = new TH1F("hPt1"," Pt1 aus den Daten",400,0,120);
   hPt1->GetXaxis()->SetTitle(" Pt1 in GeV");
   hPt1->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt1->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1->SetLineColor(4);
   hPt1->SetFillColor(4);
   
     TH1F *hPt2   = new TH1F("hPt2"," Pt2 aus den Daten",400,0,120);
   hPt2->GetXaxis()->SetTitle(" Pt2 in GeV");
   hPt2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt2->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt2->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt2->SetLineColor(4);
   hPt2->SetFillColor(4);
   
  TH1F *hMasse   = new TH1F("hMasse"," Masse M",400,0,115);
   hMasse->GetXaxis()->SetTitle(" M in GeV");
   hMasse->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasse->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasse->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasse->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasse->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasse->SetLineColor(kBlue);
   hMasse->SetFillColor(kBlue);
   
   
 Double_t M1, sumE1E2, Pt1berechnet;
 t1.Branch("M1",&M1,"M1/D");  
 t1.Branch("E1",&E1,"E1/D");
 t1.Branch("E2",&E2,"E2/D");  
 t1.Branch("phi1",&phi1,"phi1/D"); 
 t1.Branch("phi2",&phi2,"phi2/D"); 
 t1.Branch("sumE1E2",&sumE1E2,"sumE1E2/D"); 
 
 // transversale Impulssumme Ht = pt1+pt2  berechnen
 t1.Branch("Ht",&Ht,"Ht/D"); 
 //fill the tree
   for (Int_t i=0; i<nentries; i++) {
       tree->GetEntry(i);
       M1 = M;
       sumE1E2 =E1+E2;
       Ht=pt1+pt2;
       
       Pt1berechnet=sqrt(px1*px1 +py1*py1);

       hsumE1E2->Fill(sumE1E2);
       hHt->Fill(Ht);
       hPt1berechnet->Fill(Pt1berechnet);
       hPt1->Fill(pt1);
       hPt2->Fill(pt2);
       
       
       hMasse->Fill(M);
       
 t1.Fill();
   }
   
hpxpy->Write();   // Histogramm hpxpy in die Datei tree1.root schreiben!
hsumE1E2->Write();
hHt->Write();
hPt1berechnet->Write();
hPt1->Write();
hPt2->Write();
hMasse->Write();
t1.Write();
   // Ende Kopieren der Variable
   

//Beginn: Versuche aus pp->ee durch Cuts die Z0 Events herauszufiltern:
  TH1F *hMasseCut   = new TH1F("hMasseCut","  Masse M fuer Cut: pt1>20GeV pt2>20GeV und 110>M>70GeV ",1000,0,130);
   hMasseCut->GetXaxis()->SetTitle(" M in GeV");
   hMasseCut->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCut->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCut->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCut->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCut->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCut->SetLineColor(kBlue);
   hMasseCut->SetFillColor(kBlue);
   
   
   TH1F *hMasseRandom   = new TH1F("hMasseRandom","  Masse M fuer Z0-Peak: Zufallszahlen ",1000,0,130);
   hMasseRandom->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom->SetLineColor(kBlue);
   hMasseRandom->SetFillColor(kBlue);
   
      TH1F *hMasseRandom2   = new TH1F("hMasseRandom2","  simulierter Z0-Peak: Zufallszahlen, BreitWignerVerteilung ",1000,0,130);
   hMasseRandom2->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom2->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom2->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom2->SetLineColor(kBlue);
   hMasseRandom2->SetFillColor(kBlue);
   
      TH1F *hMasseRandom3   = new TH1F("hMasseRandom3","  Masse M fuer Z0-Peak: Zufallszahlen natuerliche Breite als BreitWigner und Detektoreffekte als Gau� ",1000,0,130);
   hMasseRandom3->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom3->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom3->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom3->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandom3->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom3->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom3->SetLineColor(kBlue);
   hMasseRandom3->SetFillColor(kBlue);
   
    TH1F *hMasseRandomVoigt   = new TH1F("hMasseRandomVoigt","  Masse M fuer Z0-Peak: Zufallszahlen per Voigtfunktion ",1000,0,130);
   hMasseRandomVoigt->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandomVoigt->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandomVoigt->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandomVoigt->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseRandomVoigt->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandomVoigt->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandomVoigt->SetLineColor(kBlue);
   hMasseRandomVoigt->SetFillColor(kBlue);
   
   TH1F *hMasseCut2   = new TH1F("hMasseCut2"," pp->ee Masse M fuer Cut: pt1>20GeV pt2>20GeV und 110>M>70GeV ",1000,0,130);
   hMasseCut2->GetXaxis()->SetTitle(" M in GeV");
   hMasseCut2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCut2->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCut2->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCut2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCut2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCut2->SetLineColor(kBlue);
   hMasseCut2->SetFillColor(kBlue);
   
  TH1F *hPt1Cut   = new TH1F("hPt1Cut"," Pt1 nach Cut: pt1>20GeV und pt2>20GeV ",1000,0,130);
   hPt1Cut->GetXaxis()->SetTitle(" Pt1 in GeV");
   hPt1Cut->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt1Cut->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1Cut->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hPt1Cut->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1Cut->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1Cut->SetLineColor(4);
   hPt1Cut->SetFillColor(4);
   
 TH1F *hMasseCutBackground   = new TH1F("hMasseCutBackground"," pp->ee Masse M fuer Cut: Background (kein Z0)  ",1000,0,130);
   hMasseCutBackground->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutBackground->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutBackground->GetXaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCutBackground->GetYaxis()->SetLabelSize(0.02);   // Gr��e der Zahlen des Histogramms
   hMasseCutBackground->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCutBackground->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCutBackground->SetLineColor(kRed);
   hMasseCutBackground->SetFillColor(kRed);
   hMasseCutBackground->SetLineColor(kYellow);
   hMasseCutBackground->SetFillColor(kYellow);
         
    double_t MasseCut;
    
   THStack *hStack = new THStack("hStack","Stacked 1D histograms: Z0 Background (gelb) und Z0 Events (rot); Cut pt1&pt2>20GeV"); // Histogramme aufeinanderlegen
  THStack *hStackJpsi = new THStack("hStackJpsi","Stacked 1D histograms: Background (gelb) und Z0 Events (Z0 rot) und JPsi blau"); // Histogramme aufeinanderlegen
   
   
   for (Int_t i=0; i<nentries; i++) {
       tree->GetEntry(i);
       
      //if(pt2>20&&pt1>20&&E1>20&&E2>20){    // 6675 Events erf�llen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
        if(pt2>20&&pt1>20&&M>70){      // 6675 Events erf�llen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
       hMasseCut->Fill(M);
       hMasseCut2->Fill(M);
       hPt1Cut->Fill(pt1);
      }else{//hMasseCutBackground->Fill(M);}
      }
      
      if( !(pt2>20&&pt1>20&&M>70)){
      hMasseCutBackground->Fill(M);    
      }
   }
   

   hStack->Add(hMasseCutBackground);
   hStack->Add(hMasseCut);
   
   
   //Beginn: Cuts f�r das J/Psi
 
   
   //Ende JPsi Cuts
   
   
      //Beginn: Cuts f�r das Bottonium
   
  
t1.Write();



//Ende: Versuche aus pp->ee durch Cuts die Z0 Events herauszufiltern. --------------------------------------------------------
   
   // obiger Code f�hrt zu einem Einlesen der ersten 10 Spalten der csv-Datei-Tabelle!!
    // der Inhalt ist im ROOT-TBrowser darstellbar!
   // Type	Run	Event	E1	px1	py1	pz1	pt1	eta1	phi1	Q1	E2	px2	py2	pz2	pt2	eta2	phi2	Q2	M

  // Die Anzahl der Events (663 St�ck) stimmt mit der *.csv Datei �berein!
   // Plausibilit�tstest: Die Variable M (invariante Masse) hat ihr Maximum bei 90 GeV, der Masse des Z0 aus dem die beiden Leptonen erzeugt wurden   
// Root  Version   5.34/19       9 July 2014   *
// g++ Version: 4:4.9.2-2ubuntu2      gcc Version: 4.9.2-2ubuntu2
  
 

// Beginn: Fitten mit der Maximum-Likelihood-Methode (diese soll selbst implementiert werden)
// Aufstellen der Likelihood-Funktion und Maximieren dieser Funktion (oder -Log-Likelihood und minimieren)
// d.h. der Parameter a der Wahrscheinlichkeitsdichtefunktion wird so gew�hlt, dass
// die Likelihoodfunktion maximal wird.

  //hMasseCut->Fit("funktion1","R");
 
  //hMasseCut->Fit("gaus");
  

  
  hMasseCut->Write();

 
  // par[3] ist ein Array mit den 3 Parametern der Gaussfunktion
  double par [3] = {500.0,hMasseCut->GetMean(),hMasseCut->GetRMS()};
  

// Beginn Code f�r MaximumLikelihood-Methode

// zum Beschleunigen der Funktionsauswertungen werden die Stichprobenwerte X_i (Variable M) in einem Array xwerte[7000] zwischengespeichert

   AnzahlEvents=0;
  for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
       
      if(pt2>20&&pt1>20&&M>75){
      xwerte[AnzahlEvents] =M;AnzahlEvents++;}     // lese die Daten ein (M ist die Masse; Datensatz ist dielectron100k
  }

//Setze Anfangswerte der Parameter par[1] und par[2]:
double par1=15.0;
double par2=85.0;
par[1]=par1;
par[2]=par2;

//setze die x_i in die Likelihoodfunktion ein und gebe die Funktionswerte mit cout aus:
   //  fitf(xwert,par) ist bereits ein Teil der Likelihoodfunktion f�r den x_i Wert xwert[0]
   Double_t Likelihood=-1000000.0;
   long counter =0;
   
   Double_t MaxLikelihood1=-1000000.0;
   double parameter1=0.0,parameter2=0.0,parameter0=0.0;
   
// erster Teil der Maximierung (Suche Intervalle nach dem Maximum ab und gebe das gr��te Intervall an den n�chsten Maximierungsschritt weiter) 
// 2 for Schleifen f�r die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   



// erster Teil der Maximierung (Intervallhalbierungsverfahren/Gittersuchmethode): feines Gitter (damit wird das Finden von Nebenmaxima vermieden) ---------------------------------------------------------
double intervallgroesse=14.0;
    
for(int n=0;n<1;n++){  //Schleife: f�r jedes n wird ein noch kleineres Intervall gew�hlt

// 2 for Schleifen f�r die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   
for (double j2=-intervallgroesse; j2<intervallgroesse; j2=j2+intervallgroesse/10) {     
   if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}   
   par[2] = j2+par2;
for (double j=-intervallgroesse; j<intervallgroesse; j=j+intervallgroesse/10) {  
    if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}
    //Likelihood f�r verschiedene Parameterwerte auswerten:
    par[1] = j+par1;
    Likelihood=0.0;
    
  // Schleife �ber alle Events f�r den Likelihood-Fit
  //Likelihood=likelihood(par,par);
  //  Likelihood= -Rastrigin12(par,par);
    Likelihood = likelihoodBreitWigner(par);
  
   counter++;
cout << "  Parameter "<<par[1]<<" "<<par[2]<<" Likelihoodfunktion= "<<  Likelihood<<" Anzahl Funktionsauswertungen "<< counter<<endl;
cout << "   intervallgroesse = " << intervallgroesse<<endl;
}
}
intervallgroesse=intervallgroesse/2; // erreicht e^-5 (Minimum 0 der Rastriginfunktion) nach 90 Funktionsauswertungen Parameter -3.05176e-05 -3.05176e-05 Likelihoodfunktion= -3.69534e-07
//  intervall/8 ergibt ein Nebenmaximum!

par2=parameter2;
par1=parameter1;
}

//zweiter Teil Maximierung: verkleinere die Gitterseitenl�nge um Faktor 1/2 und w�hle gr�beren Abtastschritt als in erstem Teil
for(int n=0;n<20;n++){  //Schleife: f�r jedes n wird ein noch kleineres Intervall gew�hlt

// 2 for Schleifen f�r die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   
for (double j2=-intervallgroesse; j2<intervallgroesse; j2=j2+intervallgroesse/2) {     
   if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}   
   par[2] = j2+par2;
for (double j=-intervallgroesse; j<intervallgroesse; j=j+intervallgroesse/2) {  
    if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}
    //Likelihood f�r verschiedene Parameterwerte auswerten:
    par[1] = j+par1;
    Likelihood=0.0;
    
  // Schleife �ber alle Events f�r den Likelihood-Fit
  //Likelihood=likelihood(par,par);
  //  Likelihood= -Rastrigin12(par,par);
   Likelihood = likelihoodBreitWigner(par);
  
   counter++;
cout << "  Parameter "<<par[1]<<" "<<par[2]<<" Likelihoodfunktion= "<<  Likelihood<<" Anzahl Funktionsauswertungen "<< counter<<endl;
cout << "   intervallgroesse = " << intervallgroesse<<endl;
}
}
intervallgroesse=intervallgroesse/2; // erreicht e^-5 (Minimum 0 der Rastriginfunktion) nach 90 Funktionsauswertungen Parameter -3.05176e-05 -3.05176e-05 Likelihoodfunktion= -3.69534e-07
//  intervall/8 ergibt ein Nebenmaximum!

par2=parameter2;
par1=parameter1;
}
// Ende zweiter Teil der Maximierung (Intervallhalbierungsverfahren)---------------------------------------------------------



cout << endl;
cout << endl;
cout << " eigene Gittersuchmethode: Minimum der LogLikelihoodfunktion gefunden fuer Parameter "<<" p1=Gamma= "<<parameter1<<" p2=M= "<<parameter2<<"  Likelihoodfunktion= "<<  MaxLikelihood1<<endl;;
cout<< "  Anzahl Funktionsauswertungen bis das Minimum erreicht ist :"<<counter<<endl;

//cout << " GetMean =  "<< hMasseCut->GetMean()<<" GetRMS = " <<hMasseCut->GetRMS()<< "  Anzahl Events: "<< AnzahlEvents <<endl;
cout <<endl;


//int bmin = hMasseCut->GetXaxis()->FindBin(0)
//int bmax = hMasseCut->GetXaxis()->FindBin(110)
//double integral = hMasseCut->Integral(bmin,bmax);
//cout<< "  Integral ueber Histogramm hMasseCut :   bmin="<<bmin<<" bmax="<<bmax<<"  Integral="<<integral<<endl;

par[1]=parameter1;
par[2]=parameter2;
double xwert[1];
double integralPDF=0.0;

for(double x=0;x<110;x=x+0.01){
    xwert[0]=x;
//cout<<breitwigner2(xwert,par)<<"  x= "<<x<<endl;
integralPDF = integralPDF + 0.01*breitwigner2(xwert,par);
}
cout << " Integral ueber Breit-Wigner-Funktion = " << integralPDF<<endl;
cout << " Integral ueber Histogramm hMasseCut = " << hMasseCut->Integral(0,1001)<<endl;
// Verh�ltnis 6,401188647 


//Minimierung der gleichen Funktion mittels TMinuit2-------------------------------------------------------------
// Zur Info zum Compilerfehler wegen Minuit2:  zus�tzlich zum Minuit2-Includefile muss auch die Bibliothek -lMinuit2 und der Pfad zu dieser Bibliothek im Makefile angegeben werden, sonst kann man nicht
// auf Minuit2 zugreifen und es gibt den Compilerfehler dem Sinn nach "Minuit2 nicht definiert".

ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

//ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kSimplex );

min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
   
ROOT::Math::Functor f5(&likelihoodBreitWignerMinus,2); //Minuit2 findet hier das richtige Minimum
//   ROOT::Math::Functor f5(&Rastrigin01,2);  // Minuit2 findet nicht das richtige Minimum
   double step[2] = {0.01,0.01};
   double variable[2] = {15.0,85.0};
 
   min.SetFunction(f5);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
   //min.SetVariable(2,"z",variable[2], step[2]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout << "---Minuit2:  Minimum: f( Gamma=" << xs[0] << ", M=" << xs[1] << ")  =   " << likelihoodBreitWignerMinus(xs) << endl << endl<<endl;
   
   cout<< "Minuit2 Hesse():  "<<min.Hesse()<<endl;
   
  
   min.PrintResults();   //print result of minimization. Quelle: Minuit2-Class-Reference, Members der Minuit2-Klasse
   cout<<endl;
   cout<< "Minuit2 Correlation between Variable 0 and 1:   "<<min.Correlation(0,1)<<endl;
   
   //cout<< "Minuit2 min.Errors: "<<min.Errors()<<endl;
   cout<< "Minuit2 ProvidesError:  "<<min.ProvidesError()<<endl;
   cout<< "Minuit2 min.NCalls:  number of function calls to reach the minimum :"<<min.NCalls()<<endl;
   
   cout<<endl;
   
//     cout << "   Minuit2:  Minimum: f(" << xs[0] << "," << xs[1] << ")  =   " << Rastrigin01(xs) << endl << endl<<endl;   

//--------------------------------------------------------------------------
   

   
// Beginn-Minimierung per Root (GSLMinimizer)-------------------------------------------------------------------

 


// Ende GSLMinimizer-Minimierung-------------------------------------------------------------------
//   
//NumericalMinimization();   

//NumericalMinimization2();

//NumericalMinimization3();

//NumericalMinimization4();
//-----------------------------------------------------------------------------------------   

par[1]=parameter1;
par[2]=parameter2;

double para2 [3];;
para2[0]=parameter1;
para2[1]=parameter2;



// Wiso ist p2=0 ???

//Baue eine TF1 Funktion der PDF gauss() mit oben per MaxLikelihood bestimmten Parametern
par[1]=parameter1;
par[2]=parameter2;


// BreitWigner Fitten an Z0 Peak

gStyle->SetOptFit(1111);
gStyle->SetOptStat(1111);

   TF1 *gauss1 = new TF1("gauss1",gauss,0,150,3);
   gauss1->SetParameters(500,parameter1,parameter2);
   gauss1->SetParNames ("Constant","Mean_value","Sigma");
   gauss1->Write();
   
   TF1 *bw = new TF1("breitwigner3",breitwigner2,0,130,3);
   bw->SetParameters(25.0,parameter1,parameter2);
   bw->SetParNames ("Normierungskonstante","Gamma","M");
   
   

   
TCanvas *cMasseCutBreitWigner = new TCanvas("cMasseCutBreitWigner","c1",1600,1000);
hMasseCut->Draw();
bw->Draw("same");
cMasseCutBreitWigner->SaveAs("cMasseCutBreitWigner.jpg");
cMasseCutBreitWigner->Write();   

// Fitten Breit-Wigner an Z0-Peak mit Root ##############################################################################################

TF1 *func = new TF1("breitwignerRoot",breitwignerRoot,0, 130,3);
   func->SetParameters(788.1,5.543,90.56);
   func->SetParNames ("Normierungskonstante","Gamma","M");

   
hMasseCut->Fit("breitwignerRoot","IL");
//hMasseCut->Chisquare("breitwignerRoot");

//hMasseCut->Chi2Test(hMasseCut2,"UW P ",res)

TCanvas *cMasseCutBreitWigner2 = new TCanvas("cMasseCutbreitwignerRoot","c1",1600,1000);
hMasseCut->Draw();
cMasseCutBreitWigner2->SaveAs("cMasseCutbreitwignerRootOffsetPlus3.jpg");
cMasseCutBreitWigner2->Write(); 

func->Write();



TF1 *func2 = new TF1("breitwignerRoot2",breitwignerRoot2,0, 130,3);
   func2->SetParameters(1.0,5.0,95.0,3.0);
   func2->SetParNames ("Normierungskonstante","Gamma","M","Offset");
   
TF1 *func3 = new TF1("breitwignerRootGaus1",breitwignerRootGaus,0, 130,4);
   func3->SetParameters(788.0,2.49,90.56,3.0585);
   func3->SetParNames ("Normierungskonstante","Gamma","M","sigma"); 
func3->Write();   



//TF1 *vfit = new TF1("vfit1", "gaus", 0., 130.);    Fit damit funktioniert   fit vfit1

TF1 *vfit = new TF1("vfit1", "TMath::Voigt", 0, 130);


//vfit->SetParameters(1,1);
vfit->Write();

//   TF1 *f2 = new TF1("f2", "TMath::Voigt(x,[0],[1],4)",-4,4);
//   f2->SetParameters(1,1);




TCanvas *c9 = new TCanvas("c1","c1",1600,1000);
hMasseCut2->Draw();
func3->Draw("same");
c9->SaveAs("Voigfunktion_und_Z0Massenpeak.jpg");
c9->Write();


//hMasseCut2->Fit("breitwignerRoot2","L");
//TCanvas *cMasseCutBreitWigner3 = new TCanvas("cMasseCutbreitwignerRoot2_Normierungskonstante","c1",1600,1000);
//hMasseCut2->Draw();
//cMasseCutBreitWigner3->SaveAs("cMasseCutbreitwignerRoot2Offset.jpg");
//cMasseCutBreitWigner3->Write(); 
//func2->Write();

// Ende Fitten Breit-Wigner an Z0-Peak mit Root ##############################################################################################
   
// f�llen von hMasseRandom mit Zufallszahlen +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TRandom *r0 = new TRandom();

//hMasseRandom->FillRandom("gaus",10000);

func->SetParameters(788.1,5.543,90.56);
double random0,random1;

// Histogramm mit Breit-Wigner verteilten Zufallszahlen f�llen:
  for (Int_t i=0; i<58680; i++) {
  //random0 = r0->Rndm(i);
  random0 = func->GetRandom();
  
  hMasseRandom2->Fill(random0);
  }


hMasseRandom2->Fit("breitwignerRoot","IL");
hMasseRandom2->Write();

TCanvas *crandom2 = new TCanvas("crandom2","c1",1600,1000);
hMasseRandom2->Draw();
crandom2->SaveAs("hMasseRandom2_nurBreitWigner.jpg");

//func->SetParameters(788.1,2.49,91.18);
func->SetParameters(788.1,2.49,90.56);

  gRandom = r0;
  sw.Start();
  
  
// Histogramm mit gaussverteilten Zufallszahlen f�llen:  
    for (Int_t i=0; i<58680; i++) {
     // Mean=90.6378 Sigma=4.54337
     random1 = r0->Gaus(90.6378,4.54337);
     hMasseRandom->Fill(random1);
    
  }
  
// Histogramm mit Voigtfunktion- verteilten Zufallszahlen f�llen:  
   for (Int_t i=0; i<58680; i++) {
     random0 = func->GetRandom();
     random1 = r0->Gaus(0,3.0585);
     
     hMasseRandom3->Fill(random0+random1);
  }
  


  //Voigtfunktion:
TF1 *voigt1 = new TF1("voigt1",breitgausfun, 0, 130 ,4);
//voigt1->SetParameters(4.067,90.29,2513,9.266);
voigt1->SetParameters(3.503,90.72,826.1,3.578);

 voigt1->SetParNames ("Gamma","Mz","Normierungskonstante","sigma"); 
voigt1->Write();  
  
  func3->SetParameters(788.0,2.49,90.56,3.0585);
       for (Int_t i=0; i<58680; i++) {
     
     random1 = voigt1->GetRandom();
     hMasseRandomVoigt->Fill(random1);
    
  }
  hMasseRandomVoigt->Write();
  
  


//gStyle->SetOptFit(1111);  
hMasseRandom3->Fit("voigt1","IL");
hMasseRandomVoigt->Fit("voigt1","IL");

hMasseRandom->Fit("breitwignerRoot","IL");

cout << " Fit der Voigtfunktion (also Faltung BreitWigner mit Gauss) an den echten Z0-Peak: "<<endl;
//hMasseCut2->Fit("breitwignerRootGaus1","IL");

hMasseCut2->Fit("voigt1","IL");
hMasseCut2->Write();





TCanvas *crandom = new TCanvas("crandom","c1",1600,1000);
hMasseRandom3->Draw();
crandom->SaveAs("hMasseRandom3_BreitWigner_mit_Gauss.jpg");
  
hMasseRandom->Write();  
hMasseRandom3->Write();  


Double_t res[1000];


   cout << " Chi-Quadrat-Test fuer Vergleich Gauss-Verteilung mit Z0-Peak: (Mean=90.6378 Sigma=4.54337) "<<endl;
hMasseRandom->Chi2Test(hMasseCut,"UW P ",res);
   cout << " Chi-Quadrat-Test fuer Vergleich Breit-Wigner-Verteilung mit Z0-Peak: (Gamma=5.543 M=90.56) "<<endl;
hMasseRandom2->Chi2Test(hMasseCut,"UW P ",res);
 cout << " Chi-Quadrat-Test fuer Vergleich Breit-Wigner-Verteilung (Gamma=2.49 M=90.56) und Gauss-Addition (Sigma=3.0585,Mean=0.0) mit Z0-Peak: "<<endl;
hMasseRandom3->Chi2Test(hMasseCut,"UW P ",res);


cout << " Chi-Quadrat-Test fuer Vergleich Voigtfunktion mit realem Z0-Peak: "<<endl;
hMasseRandomVoigt->Chi2Test(hMasseCut,"UW P ",res);
voigt1->SetParameters(4.067,90.29,2513,9.266);
cout << " Chi-Quadrat-Test fuer Vergleich Voigtfunktion mit Voigt-Zufallszahlen-Peak: "<<endl;
hMasseRandomVoigt->Chi2Test(hMasseRandom3,"UW P ",res);


// Chi-Quadrat-Test zwischen Histogramm und Funktion:
cout << endl<<" Chisqare-Test  Z-Peak mit Voigt : "<<hMasseCut->Chisquare(voigt1)<<endl;
cout << endl<<" Chisqare-Test  Z-Peak mit Breit-Wigner : "<<hMasseCut->Chisquare(func)<<endl;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //hMasseCut->Fit("breitwigner3","L");
   bw->Write();  
   hMasseCut->Write();
   
   
   //hMasseCut2->Fit("breitwigner3");
   hMasseCut2->Write();
   
   //hMasseCut->Fit("gauss1","L");
   
   
   TF2 *f2 = new TF2("f2",Rastrigin,-5,5,-5,5,2);
   f2->Write();
   
   TF2 *f3 = new TF2("f3",Rastrigin12,-5,5,-5,5,2); // xmin,xmax,ymin,ymax,n_par,n_dimensions
   f3->Write();
   
   TF2 *Likelihood1 = new TF2("likelihood1d",likelihood1d,85,95,4,20,0);
   Likelihood1->Write();
   
   TF2 *breitwigner4 = new TF2("likelihoodBreitWignerMinusDraw",likelihoodBreitWignerMinusDraw,3,4,85,95,0);
   breitwigner4->Write();
   


TCanvas *chMasseCut1 = new TCanvas("c1","c1",1600,1000);
hMasseCut->Draw();
gauss1->Draw("same");
chMasseCut1->SaveAs("CANVASpicture_hMasseCut.jpg");
chMasseCut1->Write();

TCanvas *c3 = new TCanvas("c3","c3",1600,1000);
Likelihood1->Draw("LEGO");
c3->SaveAs("CANVASpictureLikelihood.jpg");
c3->Write();

TCanvas *c31 = new TCanvas("c31","c31",1600,1000);
breitwigner4->Draw("LEGO");
c31->SaveAs("likelihoodBreitWignerMinusDraw.jpg");
c31->Write();

TCanvas *c4 = new TCanvas("c4","c4",1600,1000);
f2->Draw("SURF2");         //  Optionen f�r Draw: LEGO,LEGO2, SURF2
c4->SaveAs("CANVASpictureRastrigin.jpg");
c4->Write();

TCanvas *c6 = new TCanvas("c6","c6",1600,1000);
f3->Draw("SURF2");         //  Optionen f�r Draw: LEGO,LEGO2, SURF2
c6->SaveAs("CANVASpictureRastrigin12.jpg");
c6->Write();

TCanvas *c5 = new TCanvas("c1","c1",1600,1000);
f3->Draw("LEGO");         //  Optionen f�r Draw: LEGO,LEGO2, SURF2
c5->SaveAs("CANVASpictureTF2_2.jpg");
c5->Write();

  
  
}




   // Normalverteilung, normiert auf 1
Double_t gauss(Double_t *x,Double_t *par){
      
      Double_t fitval = 800*0.39894228*(TMath::Exp(-0.5*((x[0]-par[1])*(x[0]-par[1]))/(par[2]*par[2]) ))/par[2];
      return fitval;
   }  
   
   
//relativistische Breit-Wigner-Formel:

Double_t breitwigner(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}

Double_t breitwigner2(const Double_t *x, const Double_t *par)
{
  double f;
  double gamma,M,E,k,g; //1/2pi =0,159154943
  E=x[0];
  gamma=par[1];
  M=par[2];
   //f= gamma/( 0,159154943*(E-M)*(E-M)+ gamma*gamma/4);  (liefert die falsche Funktion)
  //https://en.wikipedia.org/wiki/Relativistic_Breit%E2%80%93Wigner_distribution:   (liefert die richtige BreitWigner Funktion)
  // This equation is written using natural units, ? = c = 1
  //funktioniert:-----------------------------------
  //g=sqrt(M*M*(M*M+gamma*gamma));
  //k=2*sqrt(2)*M*gamma*g/(3.141592654*sqrt(M*M+g));
 // f=300*k/( (E*E -M*M)*(E*E -M*M) +M*M*gamma*gamma);
  //---------------------------------------------------------------------------------------------
  double s;
  s=E*E;
  //f= 300*s/(M*M *((s-M*M)*(s-M*M)+s*s*gamma*gamma/(M*M)  ));
  //f=300*s*gamma*gamma/( (s-M*M)*(s-M*M)+M*M*gamma*gamma);
  //f=1/ ( (E*E -M*M)*(E*E -M*M) +M*M*gamma*gamma);
  //(1/2Pi) gamma/((x-mean)*(x-mean) + gamma*gamma/4)
  //f= 0,159154943 *gamma/((E-M)*(E-M) +gamma*gamma/4);
  
  //funktioniert:------------------------------
  //double gammahalf = gamma/2.0;
  //f=  300*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf));
  //----------------------------------
    double gammahalf = gamma/2.0;
  f=  300*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf));
   
  return f;
}

Double_t breitwignerRoot(const Double_t *x, const Double_t *par)
{
  double f;
  double gamma,M,E,k,g; //1/2pi =0,159154943
  E=x[0];
  gamma=par[1];
  M=par[2];

  //funktioniert:------------------------------
  double gammahalf = gamma/2.0;
  f=  par[0]*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf));
  //----------------------------------

  return f;
}

Double_t breitwignerRootGaus(const Double_t *x, const Double_t *par)
{
  double f;
  double gamma,M,E,k,g; //1/2pi =0,159154943
  E=x[0];
  gamma=par[1];
  M=par[2];
  double sigma=par[3];
  //double mean=par[4];
  //double eta=0.85;   //Chi2 Prob 0.76
 // double eta=0.95;  Prob 0.66
  double eta=0.80;

  //funktioniert:------------------------------func3->SetParameters(788.0,2.49,90.56,3.0585,90.56);
  double gammahalf = gamma/2.0;
  f=  par[0]*  (eta*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf)) +(1-eta)*0.39894228*(TMath::Exp(-0.5*((x[0]-M)*(x[0]-M))/(sigma*sigma) ))/sigma);
  
       
  //----------------------------------

  return f;
}

Double_t breitwignerRoot2(const Double_t *x, const Double_t *par)
{
  double f;
  double gamma,M,E,k,g; //1/2pi =0,159154943
  E=x[0];
  gamma=par[1];
  M=par[2];
   //f= gamma/( 0,159154943*(E-M)*(E-M)+ gamma*gamma/4);  (liefert die falsche Funktion)
  //https://en.wikipedia.org/wiki/Relativistic_Breit%E2%80%93Wigner_distribution:   (liefert die richtige BreitWigner Funktion)
  // This equation is written using natural units, ? = c = 1
  //funktioniert:-----------------------------------
  /*
  g=sqrt(M*M*(M*M+gamma*gamma));
  k=2*sqrt(2)*M*gamma*g/(3.141592654*sqrt(M*M+g));
  f=par[0]*k/( (E*E -M*M)*(E*E -M*M) +M*M*gamma*gamma);
  */
  
  //---------------------------------------------------------------------------------------------
  double s;
  s=E*E;
  //f= 300*s/(M*M *((s-M*M)*(s-M*M)+s*s*gamma*gamma/(M*M)  ));
  //f=300*s*gamma*gamma/( (s-M*M)*(s-M*M)+M*M*gamma*gamma);
  //f=1/ ( (E*E -M*M)*(E*E -M*M) +M*M*gamma*gamma);
  //(1/2Pi) gamma/((x-mean)*(x-mean) + gamma*gamma/4)
  //f= 0,159154943 *gamma/((E-M)*(E-M) +gamma*gamma/4);
  
  //funktioniert:------------------------------
  //double gammahalf = gamma/2.0;
  //f=  300*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf));
  //----------------------------------
   double gammahalf = gamma/2.0;
  f=  par[3]+par[0]*gammahalf/(M_PI * ((E-M)*(E-M) + gammahalf*gammahalf));
  return f;
}

Double_t breitwignerMinuit2(const Double_t *x, const Double_t *par)
{
  double f;
  double gamma,M,E,k,g; //1/2pi =0,159154943
  E=x[0];
  gamma=par[0];
  M=par[1];
   //f= gamma/( 0,159154943*(E-M)*(E-M)+ gamma*gamma/4);  (liefert die falsche Funktion)
  //https://en.wikipedia.org/wiki/Relativistic_Breit%E2%80%93Wigner_distribution:   (liefert die richtige BreitWigner Funktion)
  g=sqrt(M*M*(M*M+gamma*gamma));
  k=2*sqrt(2)*M*gamma*g/(3.141592654*sqrt(M*M+g));
  
  f=300*k/( (E*E -M*M)*(E*E -M*M) +M*M*gamma*gamma);
  
  return f;
}
 
   
Double_t likelihoodBreitWigner(const Double_t *x){
    Double_t Likelihood=0.0;
    double xwert[2];
    

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enth�lt den Wert M f�r Event i
       xwert[0]=xwerte[i];
  Likelihood=Likelihood +  log(breitwigner2(xwert,x));//x sind hier die beiden Fitparameter
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}

Double_t likelihoodBreitWignerMinus(const Double_t *x){
    Double_t Likelihood=0.0;
    double xwert[2];
    

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enth�lt den Wert M f�r Event i
       xwert[0]=xwerte[i];
  Likelihood=Likelihood -  log(breitwignerMinuit2(xwert,x));//x sind hier die beiden Fitparameter
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}

Double_t likelihoodBreitWignerMinusDraw(const Double_t *x,const Double_t *par){
    Double_t Likelihood=0.0;
    double xwert[2];
    

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enth�lt den Wert M f�r Event i
       xwert[0]=xwerte[i];
  Likelihood=Likelihood -  log(breitwignerMinuit2(xwert,x));//x sind hier die beiden Fitparameter
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}
   
   
   
   
// Funktion berechnet die Log-Likelihood-Funktion mit *x==Array der x-Werte;*par==Parameter f�r Optimierung
// Achtung: x ist hier eigentlich der(die) Parameter
Double_t likelihood(Double_t *x,Double_t *par){
    Double_t Likelihood=0.0;
  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

      xwert[0] =xwerte[i]; //xwerte[i] enth�lt den Wert M f�r Event i
 
  Likelihood=Likelihood+  log(gauss(xwert,x));  //  Likelihoodfunktion als Produkt aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
  }
  return Likelihood;
}

Double_t likelihood1d(const Double_t *x,const Double_t *par){
    Double_t Likelihood=0.0;

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enth�lt den Wert M f�r Event i
 
  Likelihood=Likelihood+  log(800*0.39894228*(TMath::Exp(-0.5*((xwerte[i]-x[0])*(xwerte[i]-x[0]))/(x[1]*x[1]) ))/x[1]);
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}

// diese Version der Likelihoodfunktion ist f�r Minuit2, welcher das Minimum sucht, deshalb likelihoodxy = Minus likelihood (der eigene Maximierer sucht ein Maximum)
Double_t likelihoodxy(const Double_t *x){
    Double_t Likelihood=0.0;

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enth�lt den Wert M f�r Event i
 
  Likelihood=Likelihood-  log(800*0.39894228*(TMath::Exp(-0.5*((xwerte[i]-x[0])*(xwerte[i]-x[0]))/(x[1]*x[1]) ))/x[1]);
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}
   

   
   
   
Double_t funktion1(Double_t *x,Double_t *par){
Double_t f;

f=sin(x[0])*sin(x[1])/(x[0]*x[1]);
return f;
}


    // eigene Funktion "fitf" definieren: define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par){
      Double_t arg = 0;
      if (par[2]!=0) arg = (x[0] - par[1])/par[2];
      Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
      return fitval;
   } 
   
   

  
   
   
   
Double_t fitf2(Double_t *x,Double_t *par){
      Double_t fitval = 700*0.39894228*(TMath::Exp(-0.5*((x[0]-par[1])*(x[0]-par[1]))/(par[2]*par[2]) ))/par[2];
      return fitval;
   } 
   
   

double RosenBrock(const double *xx )
{
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

double Rastrigin(const Double_t *x,const Double_t *par){
    double f;
    f= 20.0 + x[0]*x[0] - 10.0*cos(2*3.14159265359*x[0]) + x[1]*x[1] - 10.0*cos(2*3.14159265359*x[1]);
    return f;
}
double Rastrigin12(const Double_t *x,const Double_t *par){
    double f;
    f= 20.0 + x[1]*x[1] - 10.0*cos(2*3.14159265359*x[1]) + x[2]*x[2] - 10.0*cos(2*3.14159265359*x[2]);
    return f;
}
// Version f�r Minuit2; 
double Rastrigin01(const Double_t *x){
    double f;
    f= 20.0 + x[0]*x[0] - 10.0*cos(2*3.14159265359*x[0]) + x[1]*x[1] - 10.0*cos(2*3.14159265359*x[1]);
    return f;
}





int NumericalMinimization()
{
   ROOT::Math::GSLSimAnMinimizer min;
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
 
   ROOT::Math::Functor f(&likelihoodxy,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = {85.0,6.0};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout<<endl;
   cout << "GSLSimAnMinimizer:  Minimum: f(" << xs[0] << "," << xs[1] << "): "   << likelihoodxy(xs) << endl;
   cout << " number of function calls to reach the minimum : "<<min.NCalls()<<endl;
   cout << "minimizer provides error and error matrix: (0=false or 1=true) "<<min.ProvidesError()<<endl;
   cout << " Errors: "<<min.Errors()<<endl;
 
   return 0;
}


int NumericalMinimization2()
{

    // Choose method upon creation between:
   // kConjugateFR, kConjugatePR, kVectorBFGS,
   // kVectorBFGS2, kSteepestDescent
   ROOT::Math::GSLMinimizer min( ROOT::Math::kSteepestDescent );
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
 
   ROOT::Math::Functor f(&likelihoodxy,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = {85.0,6.0};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout<<endl;
   cout << "GSLMinimizer_kSteepestDescent :  Minimum: f(" << xs[0] << "," << xs[1] << "): "  << likelihoodxy(xs) << endl;
   cout << "number of function calls to reach the minimum :"<<min.NCalls()<<endl;
   cout << "minimizer provides error and error matrix: (0=false or 1=true) "<<min.ProvidesError()<<endl;
   cout << "min.Errors: "<<min.Errors()<<endl;
 
   return 0;
}

int NumericalMinimization3()
{

    // Choose method upon creation between:
   // kConjugateFR, kConjugatePR, kVectorBFGS,
   // kVectorBFGS2, kSteepestDescent
   ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugateFR );
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
 
   ROOT::Math::Functor f(&likelihoodxy,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = {85.0,6.0};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout<<endl;
   cout << "GSLMinimizer_kConjugateFR :  Minimum: f(" << xs[0] << "," << xs[1] << "): "  << likelihoodxy(xs) << endl;
   cout << "number of function calls to reach the minimum :"<<min.NCalls()<<endl;
   cout << "minimizer provides error and error matrix: (0=false or 1=true) "<<min.ProvidesError()<<endl;
   cout << "min.Errors: "<<min.Errors()<<endl;
 
   return 0;
}

int NumericalMinimization4()
{

    // Choose method upon creation between:
   // kConjugateFR, kConjugatePR, kVectorBFGS,
   // kVectorBFGS2, kSteepestDescent
   ROOT::Math::GSLMinimizer min( ROOT::Math::kConjugatePR );
 
   min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
 
   ROOT::Math::Functor f(&likelihoodxy,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = {85.0,6.0};
 
   min.SetFunction(f);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout<<endl;
   cout << "GSLMinimizer_kConjugatePR :  Minimum: f(" << xs[0] << "," << xs[1] << "): "  << likelihoodxy(xs) << endl;
   cout << "number of function calls to reach the minimum :"<<min.NCalls()<<endl;
   cout << "minimizer provides error and error matrix: (0=false or 1=true) "<<min.ProvidesError()<<endl;
   cout << "min.Errors: "<<min.Errors()<<endl;
 
   return 0;
}

/*--------------------------------------------------------------------*/
Double_t breitgausfun(Double_t *x, Double_t *par) /*--------------------------------------------------------------------*/
{

//Fit parameters:
//par[0]=Width (scale) Breit-Wigner
//par[1]=Most Probable (MP, location) Breit mean
//par[2]=Total area (integral -inf to inf, normalization constant)
//par[3]=Width (sigma) of convoluted Gaussian function
//
//In the Landau distribution (represented by the CERNLIB approximation), //the maximum is located at x=-0.22278298 with the location parameter=0.
//This shift is corrected within this function, so that the actual
//maximum is identical to the MP parameter.

// Numeric constants
Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
Double_t twoPi = 6.2831853071795;//2Pi

// Control constants
Double_t np = 100.0; // number of convolution steps
Double_t sc = .4; // convolution extends to +-sc Gaussian sigmas

// Variables
Double_t xx;
Double_t fland;
Double_t sum = 0.0;
Double_t xlow,xupp;
Double_t step;
Double_t i;

// Range of convolution integral
xlow = x[0] - sc * par[3];
xupp = x[0] + sc * par[3];

step = (xupp-xlow) / np;

// Convolution integral of Breit and Gaussian by sum
for(i=1.0; i<=np/2; i++) {
xx = xlow + (i-.5) * step;
fland = TMath::BreitWigner(xx,par[1],par[0]);
sum += fland * TMath::Gaus(x[0],xx,par[3]);

xx = xupp - (i-.5) * step;
fland = TMath::BreitWigner(xx,par[1],par[0]);
sum += fland * TMath::Gaus(x[0],xx,par[3]);
}

return (par[2] * step * sum * invsq2pi / par[3]);
}



