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
 
 
 // Ende der Main-Funktion bei Zeile 951 
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

// für Breit-Wigner-Verteilung:
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


// für Zufallszahlen
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



Double_t breitwigner2(const Double_t* x,const Double_t* par);
Double_t breitwignerRoot(const Double_t *x, const Double_t *par);
Double_t breitwignerRoot2(const Double_t *x, const Double_t *par);

Double_t breitwignerRootGaus(const Double_t *x, const Double_t *par);
Double_t breitgaus(Double_t *x, Double_t *par);


Double_t likelihoodBreitWigner(const Double_t *x);
Double_t likelihoodBreitWignerMinusDraw(const Double_t *x,const Double_t *par);

// Breit-Wigner-Funktion im Format für Minuit2  (par[0] und par[1] )
Double_t likelihoodBreitWignerMinus(const Double_t *x);
Double_t breitwignerMinuit2(const Double_t *x, const Double_t *par);

//globale Variablen (damit bei Funktionsaufrufen nicht weitere Parameterübergaben nötig sind)
int AnzahlEvents=0;
Double_t xwerte[7000];
Double_t xwert [1] = {0.0};



// Beginn der Main-Funktion  *******************************************************************************************************************************************

//    ******************************************************************************************************************************************************************
int main(int argc, char **argv)
{


//# Einlesen der *.csv Datei dielectron100k.csv mit TTree->Readfile und schreiben in Rootfile f (Ergebnisse/dielectron100k.root) per TFile->Write   ******************************************************
  
   //Einlesen und konvertieren von *.csv Datei nach *.root
    TFile *f = new TFile("Ergebnisse/dielectron100k.root","RECREATE");
   TTree *tree = new TTree("ntuple","Daten aus einer *.csv Datei");
   //tree->ReadFile("dielectron100k.csv","Type/C:Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');  
   //Hinweis: für dielectron100k.csv entfällt "Type/C:" weil es diese Spalte in der Datei nicht gibt; Für die anderen Event-Files muss Type/C: angegeben werden.
   tree->ReadFile("Daten/dielectron100k.csv","Run/D:Event/D:E1/D:px1/D:py1/D:pz1/D:pt1/D:eta1/D:phi1/D:Q1/D:E2/D:px2/D:py2/D:pz2/D:pt2/D:eta2/D:phi2/D:Q2/D:M/D",',');
   f->Write();
   Int_t nentries = (Int_t)tree->GetEntries();
   Int_t AnzahlEreignisse = (Int_t)tree->GetEntries();
// Ende Einlesen der Daten ************************************************************************************************************************************************   

   
   
   
// Erzeugen einer zweiten Datei "tree1.root"; Dort werden die weiter unten erstellten Histogramme gespeichert. Die Datei "dielectron100k.root" enthält lediglich die Daten aus der *.csv-Datei.
  TFile f1("Ergebnisse/tree1.root","recreate");
  TTree t1("t1","a simple Tree with simple variables");
   
   
// Erstellen von Histogrammen aus den Daten, die sich jetzt im Tree "tree" befinden****************************************************************************************   
   double_t px1, py1, pz1,pt1,M, E1,E2,pt2,phi1,phi2,Ht, px2,py2,pz2, eta1,eta2;
   
  TH1F *hM   = new TH1F("hM","M Masse distribution",400,0,130);
  TH2F *hpxpy = new TH2F("hpxpy","py vs px",30,-3,3,30,-3,3);
  hpxpy->GetXaxis()->SetTitle(" x-Achse hat den Titel: px1 in GeV");
  hpxpy->GetYaxis()->SetTitle(" y-Achse hat den Titel: py1 in GeV");
  

  tree->SetBranchAddress("pt1",&pt1);   //Setze die Einzulesende Variable pt1 für tree->GetEntry(i) 
  tree->SetBranchAddress("px1",&px1);
  tree->SetBranchAddress("py1",&py1);
  tree->SetBranchAddress("pz1",&pz1);
  tree->SetBranchAddress("M",&M);
  tree->SetBranchAddress("M",&M);
  tree->SetBranchAddress("E1",&E1);
  tree->SetBranchAddress("E2",&E2);
  tree->SetBranchAddress("phi1",&phi1);
  tree->SetBranchAddress("phi2",&phi2);
  tree->SetBranchAddress("pt1",&pt1);
  tree->SetBranchAddress("pt2",&pt2);
  
    for (Int_t i = 0; i<AnzahlEreignisse; i++) {
 tree->GetEntry(i); // Lesen des Eintrags Nummer i auf dem Tree "tree"
     hM->Fill(M);   // Füllen des Histogramms hM mit der Variable M aus dem Tree tree.
    hpxpy->Fill(px1,py1);
  }  
 hM->Write();
 
 // Kopieren der Variable M aus dem Tree tree nach Tree "t1" mit Variable "M1"  funktioniert, die neue Datei heisst tree1.root
 // Bauen eines neuen Trees t1 mit neuen Variablen, die aus den Variablen des alten Trees aus Zmumu.root zusammengesetzt sind


 // Deklarieren verschiedener Histogramme für die Daten aus tree, welche dann im Tree t1 im TFile f1 "Ergebnisse/tree1.root" gespeichert werden
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
   
  TH1F *hMasse   = new TH1F("hMasse"," Masse M",1000,0,130);
   hMasse->GetXaxis()->SetTitle(" M in GeV");
   hMasse->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasse->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasse->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasse->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasse->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasse->SetLineColor(kBlue);
   hMasse->SetFillColor(kBlue);
   
   TH1F *hMasseinvariant   = new TH1F("hMasseinvariant"," invariante Masse M",1000,0,130);
   hMasseinvariant->GetXaxis()->SetTitle(" M in GeV");
   hMasseinvariant->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseinvariant->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseinvariant->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseinvariant->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseinvariant->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseinvariant->SetLineColor(kBlue);
   hMasseinvariant->SetFillColor(kBlue);
   
   
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
 double Mass;
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
       
       //Mass =sqrt( 2*(E1*E2 - px1*px2 - py1*py2-pz1*pz2 ));         //näherungsweise schon gleich M nur ohne Z-Peak
       Mass =sqrt(2*pt1*pt2* (cosh(eta1-eta2) -cos(phi1-phi2)));   //näherungsweise schon gleich M  nur ohne Z-Peak
        
       hMasseinvariant->Fill(Mass);
       
 t1.Fill();
   }
   
hpxpy->Write();   // Histogramm hpxpy in die Datei tree1.root schreiben!
hM->Write();
hsumE1E2->Write();
hHt->Write();
hPt1berechnet->Write();
hPt1->Write();
hPt2->Write();
hMasse->Write();
hMasseinvariant->Write();
t1.Write();
   // Ende Kopieren der Variable
   

//Beginn: Versuche aus pp->ee durch Cuts die Z0 Events herauszufiltern:
// Die Histogramme hMasseCut und hMasseCut2 enthalten die als Z0-Teilchen selektierten Events; Es sind 2 Histogramme weil an hMasseCut die BreitWignerVerteilung und an hMasseCut2 die Voigtfunktion gefittet wird.
  TH1F *hMasseCut   = new TH1F("hMasseCut","  Masse M fuer Cut: pt1>20GeV pt2>20GeV und 130>M>50GeV ",1000,0,130);
   hMasseCut->GetXaxis()->SetTitle(" M in GeV");
   hMasseCut->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCut->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCut->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCut->SetLineColor(kBlue);
   hMasseCut->SetFillColor(kBlue);
   
   TH1F *hMasseCut2   = new TH1F("hMasseCut2"," pp->ee Masse M fuer Cut: pt1>20GeV pt2>20GeV und 130>M>50GeV ",1000,0,130);
   hMasseCut2->GetXaxis()->SetTitle(" M in GeV");
   hMasseCut2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseCut2->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut2->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseCut2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseCut2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseCut2->SetLineColor(kBlue);
   hMasseCut2->SetFillColor(kBlue);
   
   
   TH1F *hMasseRandom   = new TH1F("hMasseRandom","  Masse M fuer Z0-Peak: gaussverteilte Zufallszahlen ",1000,0,130);
   hMasseRandom->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom->SetLineColor(kBlue);
   hMasseRandom->SetFillColor(kBlue);
   
      TH1F *hMasseRandom2   = new TH1F("hMasseRandom2","  simulierter Z0-Peak: Zufallszahlen, BreitWignerVerteilung ",1000,0,130);
   hMasseRandom2->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom2->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom2->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom2->SetLineColor(kBlue);
   hMasseRandom2->SetFillColor(kBlue);
   
      TH1F *hMasseRandom3   = new TH1F("hMasseRandom3","  Masse M fuer Z0-Peak: Zufallszahlen natuerliche Breite als BreitWigner und Addition der Detektoreffekte als Gauß ",1000,0,130);
   hMasseRandom3->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandom3->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandom3->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom3->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandom3->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandom3->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandom3->SetLineColor(kBlue);
   hMasseRandom3->SetFillColor(kBlue);

// Das Histogramm hMasseRandomVoigt wird mit Zufallszahlen der Voigtfunktion gefüllt (zum Generieren der Zufallszahlen verwendete Funktion: breitgaus)   
    TH1F *hMasseRandomVoigt   = new TH1F("hMasseRandomVoigt","  Masse M fuer Z0-Peak: Zufallszahlen per Voigtfunktion ",1000,80,100);
   hMasseRandomVoigt->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandomVoigt->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandomVoigt->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandomVoigt->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandomVoigt->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandomVoigt->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandomVoigt->SetLineColor(kBlue);
   hMasseRandomVoigt->SetFillColor(kBlue);
   
       TH1F *hMasseRandomVoigt2   = new TH1F("hMasseRandomVoigt2","  Masse M fuer Z0-Peak: Zufallszahlen per Voigtfunktion ",1000,0,130);
   hMasseRandomVoigt2->GetXaxis()->SetTitle(" M in GeV");
   hMasseRandomVoigt2->GetYaxis()->SetTitle(" Anzahl Ereignisse");
   hMasseRandomVoigt2->GetXaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandomVoigt2->GetYaxis()->SetLabelSize(0.02);   // Größe der Zahlen des Histogramms
   hMasseRandomVoigt2->GetXaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen
   hMasseRandomVoigt2->GetYaxis()->CenterTitle();  //Achsentitel in die Mitte der Achse legen   
   hMasseRandomVoigt2->SetLineColor(kBlue);
   hMasseRandomVoigt2->SetFillColor(kBlue);
   

   
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
         
// Ende: Erstellen/Deklarieren von Histogrammen ****************************************************************************************  

// Beginn: Füllen der Histogramme mit Daten aus dem Tree "tree" ****************************************************************************************    
   double_t MasseCut;
    
   THStack *hStack = new THStack("hStack","Stacked 1D histograms: Z0 Background (gelb) und Z0 Events (rot); Cut pt1&pt2>20GeV"); // Histogramme aufeinanderlegen
  THStack *hStackJpsi = new THStack("hStackJpsi","Stacked 1D histograms: Background (gelb) und Z0 Events (Z0 rot) und JPsi blau"); // Histogramme aufeinanderlegen
   
   
   for (Int_t i=0; i<nentries; i++) {
       tree->GetEntry(i);
       
      //if(pt2>20&&pt1>20&&E1>20&&E2>20){    // 6675 Events erfüllen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
        if(pt2>20&&pt1>20&&M>50){      // 6675 Events erfüllen diese Bedingung, Massenpeak ist sehr deutlich bei 90 GeV
       hMasseCut->Fill(M);
       hMasseCut2->Fill(M);
       hPt1Cut->Fill(pt1);
      }else{//hMasseCutBackground->Fill(M);}
      }
      
      if( !(pt2>20&&pt1>20&&M>50)){
      hMasseCutBackground->Fill(M);    
      }
   }
   

   hStack->Add(hMasseCutBackground);
   hStack->Add(hMasseCut);
   
   
   //Beginn: Cuts für das J/Psi
 
   
   //Ende JPsi Cuts
   
   
      //Beginn: Cuts für das Bottonium
   
  
t1.Write();
hMasseCut->Write();


// Ende: Füllen der Histogramme mit Daten aus dem Tree "tree" ***********************************************************************

//Ende: Versuche aus pp->ee durch Cuts die Z0 Events herauszufiltern. --------------------------------------------------------

// Root  Version   5.34/19       9 July 2014   *
// g++ Version: 4:4.9.2-2ubuntu2      gcc Version: 4.9.2-2ubuntu2
  
 

// Beginn: Fitten mit der Maximum-Likelihood-Methode (diese soll selbst implementiert werden)
// Aufstellen der Likelihood-Funktion und Maximieren dieser Funktion (oder -Log-Likelihood und minimieren)
// d.h. der Parameter a der Wahrscheinlichkeitsdichtefunktion wird so gewählt, dass
// die Likelihoodfunktion maximal wird.

 
  // par[3] ist ein Array mit den 3 Parametern der Gaussfunktion
  double par [3] = {500.0,hMasseCut->GetMean(),hMasseCut->GetRMS()};
  

// Beginn Code für MaximumLikelihood-Methode

// zum Beschleunigen der Funktionsauswertungen werden die Stichprobenwerte X_i (Variable M) in einem Array xwerte[7000] zwischengespeichert

   AnzahlEvents=0;
  for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
       
      if(pt2>20&&pt1>20&&M>50){
      xwerte[AnzahlEvents] =M;AnzahlEvents++;}     // lese die Daten ein (M ist die Masse; Datensatz ist dielectron100k
  }

//Setze Anfangswerte der Parameter par[1] und par[2]:
double par1=15.0;
double par2=85.0;
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
double intervallgroesse=14.0;
    
for(int n=0;n<1;n++){  //Schleife: für jedes n wird ein noch kleineres Intervall gewählt

// 2 for Schleifen für die 2 Parameter der Wahrscheinlichkeitsdichtefunktion gauss(x,par)   
for (double j2=-intervallgroesse; j2<intervallgroesse; j2=j2+intervallgroesse/10) {     
   if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}   
   par[2] = j2+par2;
for (double j=-intervallgroesse; j<intervallgroesse; j=j+intervallgroesse/10) {  
    if(Likelihood>MaxLikelihood1){     MaxLikelihood1=Likelihood;parameter1=par[1];parameter2=par[2];}
    //Likelihood für verschiedene Parameterwerte auswerten:
    par[1] = j+par1;
    Likelihood=0.0;
    
  // Schleife über alle Events für den Likelihood-Fit
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

//zweiter Teil Maximierung: verkleinere die Gitterseitenlänge um Faktor 1/2 und wähle gröberen Abtastschritt als in erstem Teil
for(int n=0;n<20;n++){  //Schleife: für jedes n wird ein noch kleineres Intervall gewählt

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
// Verhältnis 6,401188647 


//Minimierung der gleichen Funktion mittels TMinuit2-------------------------------------------------------------
// Zur Info zum Compilerfehler wegen Minuit2:  zusätzlich zum Minuit2-Includefile muss auch die Bibliothek -lMinuit2 und der Pfad zu dieser Bibliothek im Makefile angegeben werden, sonst kann man nicht
// auf Minuit2 zugreifen und es gibt den Compilerfehler dem Sinn nach "Minuit2 nicht definiert".

   
//     cout << "   Minuit2:  Minimum: f(" << xs[0] << "," << xs[1] << ")  =   " << Rastrigin01(xs) << endl << endl<<endl;   

//--------------------------------------------------------------------------

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
cMasseCutBreitWigner->SaveAs("Histogramme/cMasseCutBreitWigner.jpg");
cMasseCutBreitWigner->Write();   

// Fitten Breit-Wigner an Z0-Peak mit Root ##############################################################################################

TF1 *func = new TF1("breitwignerRoot",breitwignerRoot,0, 130,3);
   func->SetParameters(788.1,5.543,90.56);
   func->SetParNames ("Normierungskonstante","Gamma","M");
 

TCanvas *cMasseCutBreitWigner2 = new TCanvas("cMasseCutbreitwignerRoot","c1",1600,1000);
hMasseCut->Draw();
cMasseCutBreitWigner2->SaveAs("Histogramme/cMasseCutbreitwignerRootOffsetPlus3.jpg");
cMasseCutBreitWigner2->Write(); 

func->Write();



TF1 *func2 = new TF1("breitwignerRoot2",breitwignerRoot2,0, 130,3);
   func2->SetParameters(1.0,5.0,95.0,3.0);
   func2->SetParNames ("Normierungskonstante","Gamma","M","Offset");
   
TF1 *func3 = new TF1("breitwignerRootGaus1",breitwignerRootGaus,0, 130,4);
   func3->SetParameters(788.0,2.49,90.56,3.0585);
   func3->SetParNames ("Normierungskonstante","Gamma","M","sigma"); 
func3->Write();   

func->SetNpx(10000);
func2->SetNpx(10000);
func3->SetNpx(10000);


// Ende Fitten Breit-Wigner an Z0-Peak mit Root ##############################################################################################
   
// füllen von hMasseRandom mit Zufallszahlen +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TRandom *r0 = new TRandom();

//hMasseRandom->FillRandom("gaus",10000);

func->SetParameters(788.1,5.543,90.56);
double random0,random1;

// Histogramm mit Breit-Wigner verteilten Zufallszahlen füllen:
  for (Int_t i=0; i<50000; i++) {
  //random0 = r0->Rndm(i);
  random0 = func->GetRandom();
  
  hMasseRandom2->Fill(random0);
  }


//hMasseRandom2->Fit("breitwignerRoot","IL");
//hMasseRandom2->Write();

//func->SetParameters(788.1,2.49,91.18);


  TStopwatch sw;
  gRandom = r0;
  sw.Start();
  
  
// Histogramm mit gaussverteilten Zufallszahlen füllen:  
    for (Int_t i=0; i<5868; i++) {
     // Mean=90.6378 Sigma=4.54337
     random1 = r0->Gaus(90.6378,4.54337);
     hMasseRandom->Fill(random1);
    
  }
  
func->SetParameters(788.1,2.49,90.56);   //Die Funktion "func" ist die BreitWigner-Verteilung
  // Histogramm mit Voigtfunktion- verteilten Zufallszahlen füllen:  
   for (Int_t i=0; i<5868; i++) {
     random0 = func->GetRandom();
     random1 = r0->Gaus(0,3.0585);
     
     hMasseRandom3->Fill(random0+random1);
  }
  


  //Voigtfunktion:
TF1 *voigt1 = new TF1("voigt1",breitgaus, 0, 130 ,4);        //Fitrange einschränken auf 50-130
//voigt1->SetParameters(2.49,90.56,826.1,3.578);
//voigt1->SetParameters(2.49,90.56,826.1,9.0);
voigt1->SetParameters(3.500,90.72,826,3.500);
voigt1->SetParNames ("Gamma","Mz","Normierungskonstante","sigma"); 
voigt1->SetNpx(10000); // setze die Anzahl der Plotpunkte auf 10000 (Standart: 100 Punkte). Damit wird der Funktionsplot genauer.
voigt1->Write();  
  
  func3->SetParameters(788.0,2.49,90.56,3.0585);
       for (Int_t i=0; i<6000; i++) {
     
     random1 = voigt1->GetRandom();
     hMasseRandomVoigt->Fill(random1);
     hMasseRandomVoigt2->Fill(voigt1->GetRandom());
  }
  
  

 TF1 * gaus1 = new TF1("gauss", "[0] / sqrt(2.0 * TMath::Pi()) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", 0, 100);
  
hMasseRandom2->Fit("breitwignerRoot","ILR");
hMasseRandom2->Write();

//gStyle->SetOptFit(1111);  
//voigt1->SetParameters(2.49,90.56,826.1,5.0);
voigt1->SetParameters(3.500,90.72,826,3.500);


hMasseRandom3->Fit("voigt1","ILR","",0,130);
hMasseRandom3->Write();



hMasseRandomVoigt->Fit("voigt1","ILLR","",80,100);
  hMasseRandomVoigt->Write();
  
 //voigt1->SetParLimits(0,2.49,2.49);    // Parameter 0 fixieren auf 2.49
 voigt1->SetParLimits(3,4.5,4.5);    // Parameter 3 (sigma) fixieren auf 12.1
hMasseRandomVoigt2->Fit("voigt1","ILLR","",0,130);
  hMasseRandomVoigt2->Write();  
  
 /* If I remember the syntax (maybe you should look it up) to do a _binned_ likelihood fit you just need to do:
   h->Fit("gaus","L") 

  *
  *
  *
  */ 



//hMasseRandom->Fit("breitwignerRoot","IL");
  hMasseRandom->Fit("voigt1","ILR");

cout << " Fit der Voigtfunktion (also Faltung BreitWigner mit Gauss) an den echten Z0-Peak: (hMasseCut2 mit 100 Bins und hMasseCut mit 1000 Bins): "<<endl;
//hMasseCut2->Fit("breitwignerRootGaus1","IL");



//Fitrange:

//voigt1->SetNumberFitPoints(500);
//voigt1->SetNpx(1000);  // Anzahl der zu plottenden Punkte festlegen (100 ist Standart-Wert)
//voigt1->SetParLimits(0,2.49,2.49);  // Parameter 0 fixieren auf 2.49
//voigt1->SetParLimits(3,-8.382,-8.382);  // Parameter 3 fixieren auf 8.382
hMasseCut->Fit("voigt1","ILLREMW","",0,130);  // nehme Funktionsbereich von voigt1 als Fitrange
hMasseCut->Write();
//cout >> " Anzahl leere Bins: " >> hMasseCut->GetNumberFitPoints()>>endl;


hMasseCut2->Fit("voigt1","ILLR","",0,130);  // setze Fitrange 70-110
hMasseCut2->Write();
//cout >> " Anzahl leere Bins: " >> hMasseCut2->GetNumberFitPoints()>>endl;




TCanvas *crandom = new TCanvas("crandom","c1",1600,1000);
hMasseRandom3->Draw();
crandom->SaveAs("Histogramme/hMasseRandom3_BreitWigner_mit_Gauss.jpg");
  
hMasseRandom->Write();  
hMasseRandom3->Write();  
hMasseRandom2->Write();


Double_t res[1000];


   cout << " Chi-Quadrat-Test fuer Vergleich Gauss-Verteilung mit Z0-Peak: (Mean=90.6378 Sigma=4.54337) "<<endl;
hMasseRandom->Chi2Test(hMasseCut,"UW P ",res);
   cout << " Chi-Quadrat-Test fuer Vergleich Breit-Wigner-Verteilung mit Z0-Peak: (Gamma=5.543 M=90.56) "<<endl;
hMasseRandom2->Chi2Test(hMasseCut,"UW P ",res);
 cout << " Chi-Quadrat-Test fuer Vergleich Breit-Wigner-Verteilung (Gamma=2.49 M=90.56) und Gauss-Addition (Sigma=3.0585,Mean=0.0) mit Z0-Peak: "<<endl;
hMasseRandom3->Chi2Test(hMasseCut,"UW P ",res);


cout << " Chi-Quadrat-Test fuer Vergleich Voigtfunktion mit realem Z0-Peak: "<<endl;
//hMasseRandomVoigt->Chi2Test(hMasseCut,"UW P ",res);
voigt1->SetParameters(4.067,90.29,2513,9.266);
cout << " Chi-Quadrat-Test fuer Vergleich Voigtfunktion mit Voigt-Zufallszahlen-Peak: "<<endl;
//hMasseRandomVoigt->Chi2Test(hMasseRandom3,"UW P ",res);


// Chi-Quadrat-Test zwischen Histogramm und Funktion:
cout << endl<<" Chisqare-Test  Z-Peak mit Voigt : "<<hMasseCut->Chisquare(voigt1)<<endl;
cout << endl<<" Chisqare-Test  Z-Peak mit Breit-Wigner : "<<hMasseCut->Chisquare(func)<<endl;
cout << endl<<" Chisqare-Test  Z-Peak mit Gauss : "<<hMasseCut->Chisquare(gaus1)<<endl;



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
chMasseCut1->SaveAs("Histogramme/CANVASpicture_hMasseCut.jpg");
chMasseCut1->Write();

TCanvas *c3 = new TCanvas("c3","c3",1600,1000);
Likelihood1->Draw("LEGO");
c3->SaveAs("Histogramme/CANVASpictureLikelihood.jpg");
c3->Write();

TCanvas *c31 = new TCanvas("c31","c31",1600,1000);
breitwigner4->Draw("LEGO");
c31->SaveAs("Histogramme/likelihoodBreitWignerMinusDraw.jpg");
c31->Write();

TCanvas *c4 = new TCanvas("c4","c4",1600,1000);
f2->Draw("SURF2");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c4->SaveAs("Histogramme/CANVASpictureRastrigin.jpg");
c4->Write();

TCanvas *c6 = new TCanvas("c6","c6",1600,1000);
f3->Draw("SURF2");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c6->SaveAs("Histogramme/CANVASpictureRastrigin12.jpg");
c6->Write();

TCanvas *c5 = new TCanvas("c1","c1",1600,1000);
f3->Draw("LEGO");         //  Optionen für Draw: LEGO,LEGO2, SURF2
c5->SaveAs("Histogramme/CANVASpictureTF2_2.jpg");
c5->Write();

  
  
}

// Ende der Main-Funktion ##########################################################################################################################################################



// Ende der Main-Funktion ##########################################################################################################################################################



   // Normalverteilung, normiert auf 1
Double_t gauss(Double_t *x,Double_t *par){
      
      Double_t fitval = 800*0.39894228*(TMath::Exp(-0.5*((x[0]-par[1])*(x[0]-par[1]))/(par[2]*par[2]) ))/par[2];
      return fitval;
   }  
   
   
//relativistische Breit-Wigner-Formel:



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

       //xwerte[i] enthält den Wert M für Event i
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

       //xwerte[i] enthält den Wert M für Event i
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

       //xwerte[i] enthält den Wert M für Event i
       xwert[0]=xwerte[i];
  Likelihood=Likelihood -  log(breitwignerMinuit2(xwert,x));//x sind hier die beiden Fitparameter
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
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
 
  Likelihood=Likelihood+  log(800*0.39894228*(TMath::Exp(-0.5*((xwerte[i]-x[0])*(xwerte[i]-x[0]))/(x[1]*x[1]) ))/x[1]);
  //  LogLikelihoodfunktion als Summe aller Einzelfunktionswerte der Zufallsvariablen xwert[0]=M
 
  }
  return Likelihood;
}

// diese Version der Likelihoodfunktion ist für Minuit2, welcher das Minimum sucht, deshalb likelihoodxy = Minus likelihood (der eigene Maximierer sucht ein Maximum)
Double_t likelihoodxy(const Double_t *x){
    Double_t Likelihood=0.0;

  for (Int_t i=0; i<AnzahlEvents; i++) {  //6675 ist die Anzahl an Z0-Events

       //xwerte[i] enthält den Wert M für Event i
 
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
// Version für Minuit2; 
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

/*  Voigtfunktion  (eine Faltung aus Gauss und Breit-Wigner Verteilung)--------------------------------------------------------------------*/
Double_t breitgaus(Double_t *x, Double_t *par) /*--------------------------------------------------------------------*/
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



