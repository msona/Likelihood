/*
Information zu Likelihood.cc:

Das Programm liest Events im Format *.csv ein, schreibt die Daten als tree in eine *.root Datei
und schreibt in eine weitere *.root Datei (tree1.root) Histogramme der einzelnen Variablen wie Masse M,
Transversalimpuls pt usw.

Dann werden Cuts in pt durchgeführt und damit besondere Events wie Z0, Charmonium, Bottonium selektiert. Die selektierten
Events werden in einem Stack-Plot mit dem Background dargestellt (der Stack-Plot wird ebenfalls in tree1.root gespeichert).

Im vorliegenden Code ist die *.csv Datei "dielectron100k.csv" zum Einlesen hardgecoded, diese befindet sich im gleichen Verzeichnis
wie die ausführbare Datei "Likelihood".

In einem weiteren Schritt wurde die Maximum-Likelihood-Methode (ab Zeile 503) implementiert, sowie eine Gittersuchmethode (also ein Minimierer)
zum Maximieren der Log-Likelihoodfunktion. Damit wird dann ein Fit der Gaußverteilung an den Z0-Peak 
durchgeführt. Im Detail wird die Variable "M" als Stichproben-Zufallsvariable verwendet, dabei werden
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


//neue includes für Einlesen von *.csv Dateien:
#include "Riostream.h"
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TSystem.h"   // Compilerfehler ist weg sobald diese Datei included wird!!!
#include "TTree.h"
#include "THStack.h"   // Histogramme aufeinanderlegen (z.b. Background und dann Signal darauf)
#include "TF2.h"


//#include "Math/Minimizer.h"
#include "Math/Functor.h"   // aus C-Funktion einen Functor machen für Root-Minimierer

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



#include "TStopwatch.h"

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


//globale Variablen (damit bei Funktionsaufrufen nicht weitere Parameterübergaben nötig sind)
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
   TFile *f = new TFile("dielectron100k.root","RECREATE");
   TTree *tree = new TTree("ntuple","Daten aus einer *.csv Datei");
   //tree->ReadFile("dielectron100k.csv","Type/C:Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');
   tree->ReadFile("dielectron100k.csv","Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');
   
   f->Write();
  
   double_t px1, py1, pz1,pt1,M, E1,E2,pt2;
   
  TH1F *hM   = new TH1F("hM","M Masse distribution",100,-3,3);
  TH2F *hpxpy = new TH2F("hpxpy","py vs px",30,-3,3,30,-3,3);
  
  
  hpxpy->GetXaxis()->SetTitle(" x-Achse hat den Titel: px1 in GeV");
  hpxpy->GetYaxis()->SetTitle(" y-Achse hat den Titel: py1 in GeV");
  
  
  Int_t nentries = (Int_t)tree->GetEntries(); 
   
  tree->SetBranchAddress("pt1",&pt1);   //Setze die Einzulesende Variable pt1 für tree->GetEntry(i) 
 for (Int_t i = 0; i<nentries; i++) {
 tree->GetEntry(i); // Code läuft durch Compiler mit tree->GetEntry
 //cout << " pt1: " << pt1 << " px1: "<<px1;
 }
   
 tree->SetBranchAddress("px1",&px1);
  tree->SetBranchAddress("py1",&py1);
   tree->SetBranchAddress("pz1",&pz1);
    tree->SetBranchAddress("M",&M);
  for (Int_t i = 0; i<nentries; i++) {
 tree->GetEntry(i); // Code läuft durch Compiler mit tree->GetEntry  ; Die Werte von pt1 stimmen mit den echten aus den Daten überein!
 
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
 
 TFile f1("tree1.root","recreate");
 TTree t1("t1","a simple Tree with simple variables");
 
 TH1F *hHt   = new TH1F("hHt"," Transversale Impulssumme Ht=pt1+pt2",500,0,310);
   hHt->GetXaxis()->SetTitle(" Transversale Impulssumme Ht = pt1+pt2 in GeV");
   hHt->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hHt->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hHt->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
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
   hPt1berechnet->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1berechnet->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1berechnet->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1berechnet->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1berechnet->SetLineColor(4);
   hPt1berechnet->SetFillColor(4);
   
  TH1F *hPt1   = new TH1F("hPt1"," Pt1 aus den Daten",400,0,120);
   hPt1->GetXaxis()->SetTitle(" Pt1 in GeV");
   hPt1->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt1->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1->SetLineColor(4);
   hPt1->SetFillColor(4);
   
     TH1F *hPt2   = new TH1F("hPt2"," Pt2 aus den Daten",400,0,120);
   hPt2->GetXaxis()->SetTitle(" Pt2 in GeV");
   hPt2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt2->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt2->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt2->SetLineColor(4);
   hPt2->SetFillColor(4);
   
  TH1F *hMasse   = new TH1F("hMasse"," Masse M",200,0,115);
   hMasse->GetXaxis()->SetTitle(" M in GeV");
   hMasse->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasse->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasse->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
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
  TH1F *hMasseCut   = new TH1F("hMasseCut"," pp->ee Masse M fuer Cut: pt1>20GeV pt2>20GeV und M>75GeV ",1000,0,130);
   hMasseCut->GetXaxis()->SetTitle(" M in GeV");
   hMasseCut->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCut->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCut->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCut->SetLineColor(kBlue);
   //hMasseCut->SetFillColor(kRed);
   
  TH1F *hPt1Cut   = new TH1F("hPt1Cut"," Pt1 nach Cut: pt1>20GeV und pt2>20GeV ",1000,0,130);
   hPt1Cut->GetXaxis()->SetTitle(" Pt1 in GeV");
   hPt1Cut->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hPt1Cut->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1Cut->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hPt1Cut->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hPt1Cut->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hPt1Cut->SetLineColor(4);
   hPt1Cut->SetFillColor(4);
   
 TH1F *hMasseCutBackground   = new TH1F("hMasseCutBackground"," pp->ee Masse M fuer Cut: Background (kein Z0)  ",1000,0,130);
   hMasseCutBackground->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutBackground->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutBackground->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutBackground->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
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
       
      //if(pt2>20&&pt1>20&&E1>20&&E2>20){    // 6675 Events erfüllen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
        if(pt2>20&&pt1>20&&M>75){      // 6675 Events erfüllen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
       hMasseCut->Fill(M);
       hPt1Cut->Fill(pt1);
      }else{//hMasseCutBackground->Fill(M);}
      }
      
      if( !(pt2>20&&pt1>20&&M>75)){
      hMasseCutBackground->Fill(M);    
      }
   }
   

   hStack->Add(hMasseCutBackground);
   hStack->Add(hMasseCut);
   
   
   //Beginn: Cuts für das J/Psi
   TH1F *hMasseCutJpsi   = new TH1F("hMasseCutJpsi"," pp->ee Masse M fuer Cut: Jpsi ",1000,0,130);
   hMasseCutJpsi->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutJpsi->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutJpsi->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutJpsi->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutJpsi->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCutJpsi->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCutJpsi->SetLineColor(kBlue);
   hMasseCutJpsi->SetFillColor(kBlue);
   
   TH1F *hMasseCutJpsiBackground   = new TH1F("hMasseCutJpsiBackground"," pp->ee Masse M fuer Cut: Jpsi (Background) ",1000,0,130);
   hMasseCutJpsiBackground->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutJpsiBackground->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutJpsiBackground->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutJpsiBackground->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutJpsiBackground->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCutJpsiBackground->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCutJpsiBackground->SetLineColor(kYellow);
   hMasseCutJpsiBackground->SetFillColor(kYellow);
   
   
      for (Int_t i=0; i<nentries; i++) {
       tree->GetEntry(i);
       Ht=pt1+pt2;
      
      //  if(M>2,0&&M<5,16&&pt1>2,9&&pt1<6,2&&pt2>2,9&&pt2<6,2){      // selektiert einfach alle 100000 Events! (Fehler) Cut für das J/Psi Teilchen im Massenspektrum pp->ee
        if(M>2.0&&M<5.16&&pt1>2.9&&pt1<19.72&&pt2>2.9&&pt2<19.72){
       hMasseCutJpsi->Fill(M);
       
      }
      else{//hMasseCutJpsiBackground->Fill(M);
          }
      
      
            if( (!(pt2>20&&pt1>20&&M>80)) && !(M>2.0&&M<5.16&&pt1>2.9&&pt1<19.72&&pt2>2.9&&pt2<19.72)  ){
      hMasseCutJpsiBackground->Fill(M);    
      }

      }
      
   hStackJpsi->Add(hMasseCutJpsiBackground);
   hStackJpsi->Add(hMasseCutJpsi);
   hStackJpsi->Add(hMasseCut);
   
   
   //Ende JPsi Cuts
   
   
      //Beginn: Cuts für das Bottonium
   
   
   THStack *hStackBottonium = new THStack("hStackBottonium"," Background (Magenta) Z0 (Rot)  J/Psi (Blau) und Bottonium (Grün) "); // Histogramme aufeinanderlegen
   
   
   TH1F *hMasseCutBottonium   = new TH1F("hMasseCutBottonium"," pp->ee Masse M fuer Cut: Bottonium ",1000,0,130);
   hMasseCutBottonium->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutBottonium->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutBottonium->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutBottonium->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutBottonium->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCutBottonium->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCutBottonium->SetLineColor(kGreen);
   hMasseCutBottonium->SetFillColor(kGreen);
   
   TH1F *hMasseCutBottoniumBackground   = new TH1F("hMasseCutBottoniumBackground"," pp->ee Masse M fuer Cut: Bottonium (Background) ",1000,0,130);
   hMasseCutBottoniumBackground->GetXaxis()->SetTitle(" M in GeV");
   hMasseCutBottoniumBackground->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCutBottoniumBackground->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutBottoniumBackground->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCutBottoniumBackground->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCutBottoniumBackground->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCutBottoniumBackground->SetLineColor(kMagenta);
   hMasseCutBottoniumBackground->SetFillColor(kMagenta);
   
   
         for (Int_t i=0; i<nentries; i++) {
       tree->GetEntry(i);
       Ht=pt1+pt2;
      
      
        if(Ht>6.6&&Ht<38&&pt1<23.3&&pt2<22.2&&M>7.9&&M<12.07){
       hMasseCutBottonium->Fill(M);
       
      }
      else{//hMasseCutBottoniumBackground->Fill(M);
          }
          
      if( (!(pt2>20&&pt1>20&&M>75)) && !(M>2.0&&M<5.16&&pt1>2.9&&pt1<19.72&&pt2>2.9&&pt2<19.72) && !(Ht>6.6&&Ht<38&&pt1<23.3&&pt2<22.2&&M>7.9&&M<12.07)  ){  
       hMasseCutBottoniumBackground->Fill(M);   
      }
          

      }
   
   hStackBottonium->Add(hMasseCutBottoniumBackground);
   hStackBottonium->Add(hMasseCutBottonium);
   hStackBottonium->Add(hMasseCutJpsi);
   hStackBottonium->Add(hMasseCut);
   
   
   
   hStack->Write();
   hStackJpsi->Write();
   hStackBottonium->Write();
   //hMasseCut->Write();
   hMasseCutBackground->Write();
   hPt1Cut->Write();
   
   hMasseCutJpsi->Write();
    hMasseCutJpsiBackground->Write();
    
    hMasseCutBottonium->Write();
    hMasseCutBottoniumBackground->Write();
t1.Write();



//Ende: Versuche aus pp->ee durch Cuts die Z0 Events herauszufiltern. --------------------------------------------------------
   
   // obiger Code führt zu einem Einlesen der ersten 10 Spalten der csv-Datei-Tabelle!!
    // der Inhalt ist im ROOT-TBrowser darstellbar!
   // Type	Run	Event	E1	px1	py1	pz1	pt1	eta1	phi1	Q1	E2	px2	py2	pz2	pt2	eta2	phi2	Q2	M

  // Die Anzahl der Events (663 Stück) stimmt mit der *.csv Datei überein!
   // Plausibilitätstest: Die Variable M (invariante Masse) hat ihr Maximum bei 90 GeV, der Masse des Z0 aus dem die beiden Leptonen erzeugt wurden   
// Root  Version   5.34/19       9 July 2014   *
// g++ Version: 4:4.9.2-2ubuntu2      gcc Version: 4.9.2-2ubuntu2
  
 

// Beginn: Fitten mit der Maximum-Likelihood-Methode (diese soll selbst implementiert werden)
// Aufstellen der Likelihood-Funktion und Maximieren dieser Funktion (oder -Log-Likelihood und minimieren)
// d.h. der Parameter a der Wahrscheinlichkeitsdichtefunktion wird so gewählt, dass
// die Likelihoodfunktion maximal wird.

  //hMasseCut->Fit("funktion1","R");
 
  //hMasseCut->Fit("gaus");
  


   
        // Create a TF1 object using the function defined above.
      // The last three parameters specify the number of parameters
      // for the function.

//gStyle->SetOptFit(); //Fit-Parameter anzeigen

   TF1 *fit1 = new TF1("fit1",fitf,0,150,3); 
   TF1 *fit2 = new TF1("fit2",fitf2,0,150,3); 
   

   
         // set the parameters to the mean and RMS of the histogram
  fit1->SetParameters(500,hMasseCut->GetMean(),hMasseCut->GetRMS());
  fit2->SetParameters(500,hMasseCut->GetMean(),hMasseCut->GetRMS());
  
       // give the parameters meaningful names
  fit1->SetParNames ("Constant","Mean_value","Sigma");
  fit2->SetParNames ("Constant","Mean_value","Sigma");
  
       // call TH1::Fit with the name of the TF1 object
      //hMasseCut->Fit("fit1");
    // hMasseCut->Fit("gaus","L");
  
  hMasseCut->Write();
  fit1->Write();
  fit2->Write();
  
 
  // par[3] ist ein Array mit den 3 Parametern der Gaussfunktion
  double par [3] = {500.0,hMasseCut->GetMean(),hMasseCut->GetRMS()};
  

// Beginn Code für MaximumLikelihood-Methode

// zum Beschleunigen der Funktionsauswertungen werden die Stichprobenwerte X_i (Variable M) in einem Array xwerte[7000] zwischengespeichert

   AnzahlEvents=0;
  for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
       
      if(pt2>20&&pt1>20&&M>75){
      xwerte[AnzahlEvents] =M;AnzahlEvents++;}     // lese die Daten ein (M ist die Masse; Datensatz ist dielectron100k
  }

//Setze Anfangswerte der Parameter par[1] und par[2]:
double par1=86.0;
double par2=6.0;
par[1]=par1;
par[2]=par2;

//setze die x_i in die Likelihoodfunktion ein und gebe die Funktionswerte mit cout aus:
   //  fitf(xwert,par) ist bereits ein Teil der Likelihoodfunktion für den x_i Wert xwert[0]
   Double_t Likelihood=-1000000.0;
   long counter =0;
   
   Double_t MaxLikelihood1=-1000000.0;
   double parameter1=0.0,parameter2=0.0,parameter0=0.0;
   
// erster Teil der Maximierung (Suche Intervalle nach dem Maximum ab und gebe das größte Intervall an den nächsten Maximierungsschritt weiter) 
// 2 for Schleifen für die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   



// erster Teil der Maximierung (Intervallhalbierungsverfahren/Gittersuchmethode): feines Gitter (damit wird das Finden von Nebenmaxima vermieden) ---------------------------------------------------------
double intervallgroesse=5.0;
    
for(int n=0;n<1;n++){  //Schleife: für jedes n wird ein noch kleineres Intervall gewählt

// 2 for Schleifen für die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   
for (double j2=-intervallgroesse; j2<intervallgroesse; j2=j2+intervallgroesse/2) {     
   if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}   
   par[2] = j2+par2;
for (double j=-intervallgroesse; j<intervallgroesse; j=j+intervallgroesse/2) {  
    if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}
    //Likelihood für verschiedene Parameterwerte auswerten:
    par[1] = j+par1;
    Likelihood=0.0;
    
  // Schleife über alle Events für den Likelihood-Fit
  Likelihood=likelihood(par,par);
  //  Likelihood= -Rastrigin12(par,par);
  
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

//zweiter Teil Maximierung: verkleinere die Gitterseitenlänge um Faktor 1/2 und wähle gröberen Abtastschritt als in erstem Teil
for(int n=0;n<15;n++){  //Schleife: für jedes n wird ein noch kleineres Intervall gewählt

// 2 for Schleifen für die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   
for (double j2=-intervallgroesse; j2<intervallgroesse; j2=j2+intervallgroesse/2) {     
   if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}   
   par[2] = j2+par2;
for (double j=-intervallgroesse; j<intervallgroesse; j=j+intervallgroesse/2) {  
    if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}
    //Likelihood für verschiedene Parameterwerte auswerten:
    par[1] = j+par1;
    Likelihood=0.0;
    
  // Schleife über alle Events für den Likelihood-Fit
  Likelihood=likelihood(par,par);
  //  Likelihood= -Rastrigin12(par,par);
  
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
cout << " MaximumLikelihoodfunktion gefunden für Parameter "<<" p1= "<<parameter1<<" p2= "<<parameter2<<"  Likelihoodfunktion= "<<  MaxLikelihood1<<endl;

cout << " GetMean =  "<< hMasseCut->GetMean()<<" GetRMS = " <<hMasseCut->GetRMS()<< "  Anzahl Events: "<< AnzahlEvents <<endl;
cout <<endl;
cout<<endl;



//Minimierung der gleichen Funktion mittels TMinuit2-------------------------------------------------------------
// Zur Info zum Compilerfehler wegen Minuit2:  zusätzlich zum Minuit2-Includefile muss auch die Bibliothek -lMinuit2 und der Pfad zu dieser Bibliothek im Makefile angegeben werden, sonst kann man nicht
// auf Minuit2 zugreifen und es gibt den Compilerfehler dem Sinn nach "Minuit2 nicht definiert".

ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

//ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kSimplex );

min.SetMaxFunctionCalls(1000000);
   min.SetMaxIterations(100000);
   min.SetTolerance(0.001);
   
ROOT::Math::Functor f5(&likelihoodxy,2); //Minuit2 findet hier das richtige Minimum
//   ROOT::Math::Functor f5(&Rastrigin01,2);  // Minuit2 findet nicht das richtige Minimum
   double step[2] = {0.01,0.01};
   double variable[2] = {80.0,6.0};
 
   min.SetFunction(f5);
 
   // Set the free variables to be minimized!
   min.SetVariable(0,"x",variable[0], step[0]);
   min.SetVariable(1,"y",variable[1], step[1]);
 
   min.Minimize(); 
 
   const double *xs = min.X();
   cout << "   Minuit2:  Minimum: f(" << xs[0] << "," << xs[1] << ")  =   " << likelihoodxy(xs) << endl << endl<<endl;
   
   cout<< "Minuit2 Hesse: min.Hesse(); "<<min.Hesse()<<endl;
   
   cout<< "Minuit2 PrintResults: "<<endl;
   min.PrintResults();   //print result of minimization. Quelle: Minuit2-Class-Reference, Members der Minuit2-Klasse
   
   cout<< "Minuit2 Correlation between Variable 0 und 1:   "<<min.Correlation(0,1)<<endl;
   
   cout<< "Minuit2 min.Errors: "<<min.Errors()<<endl;
   cout<< "Minuit2 min.ProvidesError: "<<min.ProvidesError()<<endl;
   cout<< "Minuit2 min.NCalls:  number of function calls to reach the minimum :"<<min.NCalls()<<endl;
   
   cout<<endl;
   
//     cout << "   Minuit2:  Minimum: f(" << xs[0] << "," << xs[1] << ")  =   " << Rastrigin01(xs) << endl << endl<<endl;   

//Minuit1:

  
TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
//TMinuit minuit(2);

// Beginn-Minimierung per Root-------------------------------------------------------------------

  // Choose method upon creation between:
   // kConjugateFR, kConjugatePR, kVectorBFGS,
   // kVectorBFGS2, kSteepestDescent
 // ROOT::Math::GSLSimAnMinimizer min;




// Ende TMinuit-Minimierung-------------------------------------------------------------------

par[1]=parameter1;
par[2]=parameter2;

double para2 [3];;
para2[0]=parameter1;
para2[1]=parameter2;



// Wiso ist p2=0 ???

//Baue eine TF1 Funktion der PDF gauss() mit oben per MaxLikelihood bestimmten Parametern
par[1]=parameter1;
par[2]=parameter2;


   TF1 *gauss1 = new TF1("gauss1",gauss,0,150,3);
   gauss1->SetParameters(500,parameter1,parameter2);
   gauss1->SetParNames ("Constant","Mean_value","Sigma");
   gauss1->Write();
   
   TF2 *f2 = new TF2("f2",Rastrigin,-5,5,-5,5,2);
   f2->Write();
   
   TF2 *f3 = new TF2("f3",Rastrigin12,-5,5,-5,5,2); // xmin,xmax,ymin,ymax,n_par,n_dimensions
   f3->Write();
   
   TF2 *Likelihood1 = new TF2("likelihood1d",likelihood1d,85,95,4,20,0);
   Likelihood1->Write();


TCanvas *chMasseCut1 = new TCanvas("c1","c1",1600,1000);
hMasseCut->Draw();
gauss1->Draw("same");
chMasseCut1->SaveAs("CANVASpicture_hMasseCut.jpg");
chMasseCut1->Write();

TCanvas *c3 = new TCanvas("c3","c3",1600,1000);
Likelihood1->Draw("LEGO");
c3->SaveAs("CANVASpictureLikelihood.jpg");
c3->Write();

TCanvas *c4 = new TCanvas("c4","c4",1600,1000);
f2->Draw("SURF2");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c4->SaveAs("CANVASpictureRastrigin.jpg");
c4->Write();

TCanvas *c6 = new TCanvas("c6","c6",1600,1000);
f3->Draw("SURF2");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c6->SaveAs("CANVASpictureRastrigin12.jpg");
c6->Write();

TCanvas *c5 = new TCanvas("c1","c1",1600,1000);
f3->Draw("LEGO");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c5->SaveAs("CANVASpictureTF2_2.jpg");
c5->Write();

  
  
}




   // Normalverteilung, normiert auf 1
Double_t gauss(Double_t *x,Double_t *par){
      
      Double_t fitval = 700*0.39894228*(TMath::Exp(-0.5*((x[0]-par[1])*(x[0]-par[1]))/(par[2]*par[2]) ))/par[2];
      return fitval;
   }  
// Funktion berechnet die Log-Likelihood-Funktion mit *x==Array der x-Werte;*par==Parameter für Optimierung
// Achtung: x ist hier eigentlich der(die) Parameter
Double_t likelihood(Double_t *x,Double_t *par){
    Double_t Likelihood=0.0;
  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

      xwert[0] =xwerte[i]; //xwerte[i] enthält den Wert M für Event i
 
  Likelihood=Likelihood+  log(gauss(xwert,x));  //  Likelihoodfunktion als Produkt aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
  }
  return Likelihood;
}

Double_t likelihood1d(const Double_t *x,const Double_t *par){
    Double_t Likelihood=0.0;

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enthält den Wert M für Event i
 
  Likelihood=Likelihood+  log(700*0.39894228*(TMath::Exp(-0.5*((xwerte[i]-x[0])*(xwerte[i]-x[0]))/(x[1]*x[1]) ))/x[1]);
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}

// diese Version der Likelihoodfunktion ist für Minuit2, welcher das Minimum sucht, deshalb likelihoodxy = Minus likelihood (der eigene Maximierer sucht ein Maximum)
Double_t likelihoodxy(const Double_t *x){
    Double_t Likelihood=0.0;

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enthält den Wert M für Event i
 
  Likelihood=Likelihood-  log(700*0.39894228*(TMath::Exp(-0.5*((xwerte[i]-x[0])*(xwerte[i]-x[0]))/(x[1]*x[1]) ))/x[1]);
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
// Version für Minuit2; 
double Rastrigin01(const Double_t *x){
    double f;
    f= 20.0 + x[0]*x[0] - 10.0*cos(2*3.14159265359*x[0]) + x[1]*x[1] - 10.0*cos(2*3.14159265359*x[1]);
    return f;
}
