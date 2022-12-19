/****************************************************************************************************************************
******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
******    Mike Sas, mike.sas@cern.ch                                                                                    *****
******    Lucas Altenkaemper, lucas.altenkaemper@cern.ch                                                                *****
******    Lucia Leardini, lucia.leardini@cern.ch                                                                        *****
*****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
#include "CocktailHEPDataPbPb.h"
#include <iostream>


void  CocktailProduceCompleteInputFilePbPb( TString enableCentralities2760GeV   = "01111111111011",
                                            TString enableCentralities5TeV      = "000000000000",
                                            TString enableParticles             = "11000111010100011000001000",
                                            TString suffix                      = "eps",
                                            Bool_t produceAllSpectraPlots       = kTRUE,
                                            Int_t withBinShift                  = 1,
                                            Bool_t enableFlow                   = kTRUE
                                            //TString convertTo               = ""                      // "dN/d#it{y}dpT" or "1/2pipT*dN/d#it{y}dpT"
)
{
    //================================================================================================================
    // Enable energies and centralities to be included
    //================================================================================================================
    // This is done by changing the "enableCentralitiesXTeV" string
    // "enableCentralitiesXTeV" has 11 digits, set each digit either 0(exclude) or 1(include)
    // see arrays above for order
    // to exclude an energy, disable all centralities
    // Example: enableCentralities5TeV="0000111100" enables PbPb 5TeV and the centralities 10-20, 0-20, 20-40, 20-50

    Bool_t enable2760GeV    = kTRUE;
    Bool_t enable5TeV       = kFALSE;

    Bool_t includeCentrality2760GeV[nCentralities]  = {kFALSE};
    Bool_t includeCentrality5TeV[nCentralities]     = {kFALSE};

    Int_t tempCounter   = 0;
    TString tempString  = "";
    for (Int_t i=0; i<nCentralities; i++) {
        tempString = enableCentralities2760GeV(i,1);

        if (tempString.CompareTo("1") == 0) includeCentrality2760GeV[i] = kTRUE;
        else tempCounter++;
    }
    if (tempCounter == nCentralities) enable2760GeV = kFALSE;

    tempCounter = 0;
    tempString  = "";
    for (Int_t i=0; i<nCentralities; i++) {
        tempString = enableCentralities5TeV(i,1);

        if (tempString.CompareTo("1") == 0) includeCentrality5TeV[i] = kTRUE;
        else tempCounter++;
    }
    if (tempCounter == nCentralities) enable5TeV = kFALSE;

    if (!enable5TeV && !enable2760GeV) {
        cout << "Warning: No energies and centralities selected, stopping!" << endl;
        return;
    }

    //================================================================================================================
    // Enable particles to be included
    //================================================================================================================
    // This is done by changing the "enableParticles" string
    // "enableParticles" has 22 digits, set each digit either 0(exclude) or 1(include)
    // the particles are sorted as follows:
    // 0)   1)   2)     3)    4)         5)      6)     7)        8)     9)   10)   11)    12)      13)      14)        15)    16)     17)      18)        19)        20)     21)   22)     23)     24)     25)
    // pi0, eta, omega, eta', gamma_dir, pi^+/-, K^+/-, p/anti-p, h^+/-, phi, K^*0, rho^0, rho^+/-, Delta^0, Delta^+/-, K^0_s, Lambda, Sigma^0, Sigma^+/-, Omega^+/-, Xi^+/-, J/psi D^0     D^+     D^*+    D_s^+
    // Example: enableParticles="0000011100000000000000" enables pi^+/-, K^+/- and p/anti-p

    Bool_t includeParticle[nParticles] = {kFALSE};
    tempCounter = 0;
    tempString  = "";
    for (Int_t i=0; i<nParticles; i++) {
        tempString = enableParticles(i,1);

        if (tempString.CompareTo("1") == 0) includeParticle[i] = kTRUE;
        else tempCounter++;
    }

    if (tempCounter == nParticles) {
        cout << "Warning: No particles selected, stopping!" << endl;
        return;
    }


    //================================================================================================================
    //Creating output file structure
    //================================================================================================================

    TFile *outputFile                   = new TFile("CocktailInputPbPb.root","RECREATE");
    TList *lists2760GeV[nCentralities]  = {NULL};
    TList *lists5TeV[nCentralities]     = {NULL};

    if (enable2760GeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality2760GeV[i]) {
                lists2760GeV[i] = new TList();
                lists2760GeV[i]->SetName(Form("%s_%s_%s", fCollSys[2].Data(), fEnergy[1].Data(), fCentrality[i].Data()));
            }
        }
    }

    //================================================================================================================
    // Set input files
    //================================================================================================================

    TFile *fileNeutralMeson2760GeV                   = NULL;
    if (includeParticle[0] || includeParticle[1]) {
        if (enable2760GeV) fileNeutralMeson2760GeV   = new TFile("PbPb/2760_GeV/NeutralMesonInputPbPb2760GeV_2017_07_28.root"); //was 14july
    }

    TFile *fileChargedPionKaonProton2760GeV          = NULL;
    TFile *fileChargedPionKaonProton2760GeVRatios    = NULL;
    if (includeParticle[5] || includeParticle[6] || includeParticle[7]) {
        if (enable2760GeV) {
            fileChargedPionKaonProton2760GeV         = new TFile("PbPb/2760_GeV/PbPb276.fullpT.INEL.20140329.root");
            fileChargedPionKaonProton2760GeVRatios   = new TFile("PbPb/2760_GeV/PbPb276.fullpT.RATIOS.20140329.root");
        }
    }

    TFile *filePhiMeson2760GeV                   = NULL;
    if (includeParticle[9]) {
        if (enable2760GeV) filePhiMeson2760GeV   = new TFile("PbPb/2760_GeV/PhiSpectraPbPb2760GeV.root");
    }

    TFile *fileRho0Meson2760GeV                   = NULL;
    if (includeParticle[11]) {
        if (enable2760GeV) fileRho0Meson2760GeV   = new TFile("PbPb/2760_GeV/rho_spectra_PbPb2760GeV.root");
    }

    TFile *fileK0sLambda2760GeV                            = NULL;
    if (includeParticle[15] || includeParticle[16]) {
        if (enable2760GeV) fileK0sLambda2760GeV            = new TFile("PbPb/2760_GeV/k0s_lambda_final_spectra.root");
    }


    //================================================================================================================
    // creating histos and graphs for 2760GeV
    //================================================================================================================
    TDirectoryFile* tempList = NULL;
    if (enable2760GeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality2760GeV[i]) {

                cout << fEnergy[1].Data() << "_" << fCentrality[i].Data() << endl;

                //================================================================================================================
                // reading and writing pi0 to 2760GeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[0]) {
                    TGraphAsymmErrors* graphNPionComb2760GeVStat  = NULL;
                    TGraphAsymmErrors* graphNPionComb2760GeVSys   = NULL;
                    TGraphAsymmErrors* graphNPionPCM2760GeVStat   = NULL;
                    TGraphAsymmErrors* graphNPionPCM2760GeVSys    = NULL;
                    TGraphAsymmErrors* graphNPionPHOS2760GeVStat  = NULL;
                    TGraphAsymmErrors* graphNPionPHOS2760GeVSys   = NULL;
                    TGraphAsymmErrors* graphNPionEMCal2760GeVStat = NULL;
                    TGraphAsymmErrors* graphNPionEMCal2760GeVSys  = NULL;

                    TH1D* histoNPionComb2760GeVStat  = NULL;
                    TH1D* histoNPionComb2760GeVSys   = NULL;
                    TH1D* histoNPionPCM2760GeVStat   = NULL;
                    TH1D* histoNPionPCM2760GeVSys    = NULL;
                    TH1D* histoNPionPHOS2760GeVStat  = NULL;
                    TH1D* histoNPionPHOS2760GeVSys   = NULL;
                    TH1D* histoNPionEMCal2760GeVStat = NULL;
                    TH1D* histoNPionEMCal2760GeVSys  = NULL;
                    if (i==0) {
                        // MB
                    } else {
                        // centrality classes
                        if (fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphInvYieldPi0CombPbPb2760GeVStatErr_%s", fCentrality[i].Data())) && fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphInvYieldPi0CombPbPb2760GeVSysErr_%s", fCentrality[i].Data()))) {

                            // spectra
                            cout << " - pi^0 spectrum" << endl;
                            if(withBinShift){

                                TString shift = "";
                                if(withBinShift == 2) shift = "_yShifted";
                                else shift = "";
                                graphNPionComb2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0CombPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionComb2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0CombPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionPCM2760GeVStat   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PCMPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionPCM2760GeVSys    = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PCMPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionPHOS2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PHOSPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionPHOS2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PHOSPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionEMCal2760GeVStat = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0EMCalPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphNPionEMCal2760GeVSys  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0EMCalPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));


                                graphNPionComb2760GeVStat                                = ConvertYieldGraph(graphNPionComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                graphNPionComb2760GeVSys                                 = ConvertYieldGraph(graphNPionComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                SetGraphProperties(graphNPionComb2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                SetGraphProperties(graphNPionComb2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                lists2760GeV[i]->Add(graphNPionComb2760GeVStat);
                                lists2760GeV[i]->Add(graphNPionComb2760GeVSys);

                                if(graphNPionPCM2760GeVStat){
                                  graphNPionPCM2760GeVStat     = ConvertYieldGraph(graphNPionPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionPCM2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionPCM2760GeVStat);
                                }
                                if(graphNPionPCM2760GeVSys){
                                  graphNPionPCM2760GeVSys       = ConvertYieldGraph(graphNPionPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionPCM2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionPCM2760GeVSys);
                                }
                                if(graphNPionPHOS2760GeVStat){
                                  graphNPionPHOS2760GeVStat   = ConvertYieldGraph(graphNPionPHOS2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionPHOS2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionPHOS2760GeVStat);
                                }
                                if(graphNPionPHOS2760GeVSys){
                                  graphNPionPHOS2760GeVSys     = ConvertYieldGraph(graphNPionPHOS2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionPHOS2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionPHOS2760GeVSys);
                                }
                                if(graphNPionEMCal2760GeVStat){
                                  graphNPionEMCal2760GeVStat = ConvertYieldGraph(graphNPionEMCal2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionEMCal2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionEMCal2760GeVStat);
                                }
                                if(graphNPionEMCal2760GeVSys){
                                  graphNPionEMCal2760GeVSys = ConvertYieldGraph(graphNPionEMCal2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphNPionEMCal2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphNPionEMCal2760GeVSys);
                                }

                            } else {

                                graphNPionComb2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0CombPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphNPionComb2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0CombPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                                graphNPionPCM2760GeVStat   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PCMPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphNPionPCM2760GeVSys    = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PCMPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                                graphNPionPHOS2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PHOSPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphNPionPHOS2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0PHOSPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                                graphNPionEMCal2760GeVStat = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0EMCalPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphNPionEMCal2760GeVSys  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldPi0EMCalPbPb2760GeVSysErr_%s", fCentrality[i].Data()));

                                histoNPionComb2760GeVStat = GraphToHist_withErrors(graphNPionComb2760GeVStat);
                                histoNPionComb2760GeVSys = GraphToHist_withErrors(graphNPionComb2760GeVSys);
                                if(graphNPionPCM2760GeVStat)histoNPionPCM2760GeVStat = GraphToHist_withErrors(graphNPionPCM2760GeVStat);
                                if(graphNPionPCM2760GeVSys)histoNPionPCM2760GeVSys = GraphToHist_withErrors(graphNPionPCM2760GeVSys);
                                if(graphNPionPHOS2760GeVStat)histoNPionPHOS2760GeVStat = GraphToHist_withErrors(graphNPionPHOS2760GeVStat);
                                if(graphNPionPHOS2760GeVSys)histoNPionPHOS2760GeVSys = GraphToHist_withErrors(graphNPionPHOS2760GeVSys);
                                if(graphNPionEMCal2760GeVStat)histoNPionEMCal2760GeVStat = GraphToHist_withErrors(graphNPionEMCal2760GeVStat);
                                if(graphNPionEMCal2760GeVSys)histoNPionEMCal2760GeVSys = GraphToHist_withErrors(graphNPionEMCal2760GeVSys);

                                histoNPionComb2760GeVStat                                = ConvertYieldHisto(histoNPionComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                histoNPionComb2760GeVSys                                 = ConvertYieldHisto(histoNPionComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                SetHistoProperties(histoNPionComb2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                SetHistoProperties(histoNPionComb2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                lists2760GeV[i]->Add(histoNPionComb2760GeVStat);
                                lists2760GeV[i]->Add(histoNPionComb2760GeVSys);

                                if(graphNPionPCM2760GeVStat){
                                  histoNPionPCM2760GeVStat     = ConvertYieldHisto(histoNPionPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionPCM2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionPCM2760GeVStat);
                                }
                                if(graphNPionPCM2760GeVSys){
                                  histoNPionPCM2760GeVSys       = ConvertYieldHisto(histoNPionPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionPCM2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionPCM2760GeVSys);
                                }
                                if(graphNPionPHOS2760GeVStat){
                                  histoNPionPHOS2760GeVStat   = ConvertYieldHisto(histoNPionPHOS2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionPHOS2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionPHOS2760GeVStat);
                                }
                                if(graphNPionPHOS2760GeVSys){
                                  histoNPionPHOS2760GeVSys     = ConvertYieldHisto(histoNPionPHOS2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionPHOS2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionPHOS2760GeVSys);
                                }
                                if(graphNPionEMCal2760GeVStat){
                                  histoNPionEMCal2760GeVStat = ConvertYieldHisto(histoNPionEMCal2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionEMCal2760GeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionEMCal2760GeVStat);
                                }
                                if(graphNPionEMCal2760GeVSys){
                                  histoNPionEMCal2760GeVSys = ConvertYieldHisto(histoNPionEMCal2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoNPionEMCal2760GeVSys, Form("%s%sSys", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoNPionEMCal2760GeVSys);
                                }
                            }
                        }
                    }
                }
                //================================================================================================================
                // reading and writing eta to 2760GeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[1]) {
                    TGraphAsymmErrors* graphEtaComb2760GeVStat  = NULL;
                    TGraphAsymmErrors* graphEtaComb2760GeVSys   = NULL;
                    TGraphAsymmErrors* graphEtaPCM2760GeVStat   = NULL;
                    TGraphAsymmErrors* graphEtaPCM2760GeVSys    = NULL;
                    TGraphAsymmErrors* graphEtaEMCal2760GeVStat = NULL;
                    TGraphAsymmErrors* graphEtaEMCal2760GeVSys  = NULL;

                    TH1D* histoEtaComb2760GeVStat  = NULL;
                    TH1D* histoEtaComb2760GeVSys   = NULL;
                    TH1D* histoEtaPCM2760GeVStat   = NULL;
                    TH1D* histoEtaPCM2760GeVSys    = NULL;
                    TH1D* histoEtaEMCal2760GeVStat = NULL;
                    TH1D* histoEtaEMCal2760GeVSys  = NULL;

                    TGraphAsymmErrors* graphEtaToNPionComb2760GeVStat  = NULL;
                    TGraphAsymmErrors* graphEtaToNPionComb2760GeVSys   = NULL;
                    TGraphAsymmErrors* graphEtaToNPionPCM2760GeVStat   = NULL;
                    TGraphAsymmErrors* graphEtaToNPionPCM2760GeVSys    = NULL;
                    TGraphAsymmErrors* graphEtaToNPionEMCal2760GeVStat = NULL;
                    TGraphAsymmErrors* graphEtaToNPionEMCal2760GeVSys  = NULL;
                    if (i==0) {
                        // MB
                    } else {
                        // centrality classes
                        if (fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphInvYieldEtaCombPbPb2760GeVStatErr_%s", fCentrality[i].Data())) && fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphInvYieldEtaCombPbPb2760GeVSysErr_%s", fCentrality[i].Data()))) {

                            // spectra
                            cout << " - eta spectrum" << endl;
                            if(withBinShift){

                                TString shift = "";
                                if(withBinShift == 2) shift = "_yShifted";
                                else shift = "";

                                graphEtaComb2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaCombPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphEtaComb2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaCombPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphEtaPCM2760GeVStat   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaPCMPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphEtaPCM2760GeVSys    = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaPCMPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphEtaEMCal2760GeVStat = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaEMCalPbPb2760GeVStatErr%s_%s", shift.Data(), fCentrality[i].Data()));
                                graphEtaEMCal2760GeVSys  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaEMCalPbPb2760GeVSysErr%s_%s", shift.Data(), fCentrality[i].Data()));

                                graphEtaComb2760GeVStat                                = ConvertYieldGraph(graphEtaComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                graphEtaComb2760GeVSys                                 = ConvertYieldGraph(graphEtaComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                SetGraphProperties(graphEtaComb2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                SetGraphProperties(graphEtaComb2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                lists2760GeV[i]->Add(graphEtaComb2760GeVStat);
                                lists2760GeV[i]->Add(graphEtaComb2760GeVSys);

                                if(graphEtaPCM2760GeVStat){
                                  graphEtaPCM2760GeVStat     = ConvertYieldGraph(graphEtaPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphEtaPCM2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphEtaPCM2760GeVStat);
                                }
                                if(graphEtaPCM2760GeVSys){
                                  graphEtaPCM2760GeVSys       = ConvertYieldGraph(graphEtaPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphEtaPCM2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphEtaPCM2760GeVSys);
                                }
                                if(graphEtaEMCal2760GeVStat){
                                  graphEtaEMCal2760GeVStat = ConvertYieldGraph(graphEtaEMCal2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphEtaEMCal2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphEtaEMCal2760GeVStat);
                                }
                                if(graphEtaEMCal2760GeVSys){
                                  graphEtaEMCal2760GeVSys = ConvertYieldGraph(graphEtaEMCal2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetGraphProperties(graphEtaEMCal2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(graphEtaEMCal2760GeVSys);
                                }

                            } else {

                                graphEtaComb2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaCombPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphEtaComb2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaCombPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                                graphEtaPCM2760GeVStat   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaPCMPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphEtaPCM2760GeVSys    = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaPCMPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                                graphEtaEMCal2760GeVStat = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaEMCalPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                                graphEtaEMCal2760GeVSys  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphInvYieldEtaEMCalPbPb2760GeVSysErr_%s", fCentrality[i].Data()));

                                histoEtaComb2760GeVStat = GraphToHist_withErrors(graphEtaComb2760GeVStat);
                                histoEtaComb2760GeVSys = GraphToHist_withErrors(graphEtaComb2760GeVSys);
                                if(graphEtaPCM2760GeVStat)histoEtaPCM2760GeVStat = GraphToHist_withErrors(graphEtaPCM2760GeVStat);
                                if(graphEtaPCM2760GeVSys)histoEtaPCM2760GeVSys = GraphToHist_withErrors(graphEtaPCM2760GeVSys);
                                if(graphEtaEMCal2760GeVStat)histoEtaEMCal2760GeVStat = GraphToHist_withErrors(graphEtaEMCal2760GeVStat);
                                if(graphEtaEMCal2760GeVSys)histoEtaEMCal2760GeVSys = GraphToHist_withErrors(graphEtaEMCal2760GeVSys);

                                histoEtaComb2760GeVStat                                = ConvertYieldHisto(histoEtaComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                histoEtaComb2760GeVSys                                 = ConvertYieldHisto(histoEtaComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                SetHistoProperties(histoEtaComb2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                SetHistoProperties(histoEtaComb2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                lists2760GeV[i]->Add(histoEtaComb2760GeVStat);
                                lists2760GeV[i]->Add(histoEtaComb2760GeVSys);

                                if(graphEtaPCM2760GeVStat){
                                  histoEtaPCM2760GeVStat     = ConvertYieldHisto(histoEtaPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoEtaPCM2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoEtaPCM2760GeVStat);
                                }
                                if(graphEtaPCM2760GeVSys){
                                  histoEtaPCM2760GeVSys       = ConvertYieldHisto(histoEtaPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoEtaPCM2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoEtaPCM2760GeVSys);
                                }
                                if(graphEtaEMCal2760GeVStat){
                                  histoEtaEMCal2760GeVStat = ConvertYieldHisto(histoEtaEMCal2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoEtaEMCal2760GeVStat, Form("%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoEtaEMCal2760GeVStat);
                                }
                                if(graphEtaEMCal2760GeVSys){
                                  histoEtaEMCal2760GeVSys = ConvertYieldHisto(histoEtaEMCal2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                                  SetHistoProperties(histoEtaEMCal2760GeVSys, Form("%s%sSys", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                                  lists2760GeV[i]->Add(histoEtaEMCal2760GeVSys);
                                }
                            }
                        }

                        // ratio
                        if (fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphEtaToPi0CombPbPb2760GeVStatErr_%s", fCentrality[i].Data())) && fileNeutralMeson2760GeV->GetListOfKeys()->Contains(Form("graphEtaToPi0CombPbPb2760GeVSysErr_%s", fCentrality[i].Data()))) {

                            // spectra
                            cout << " - eta / pi0 ratio" << endl;
                            graphEtaToNPionComb2760GeVStat  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0CombPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                            graphEtaToNPionComb2760GeVSys   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0CombPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                            graphEtaToNPionPCM2760GeVStat   = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0PCMPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                            graphEtaToNPionPCM2760GeVSys    = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0PCMPbPb2760GeVSysErr_%s", fCentrality[i].Data()));
                            graphEtaToNPionEMCal2760GeVStat = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0EMCalPbPb2760GeVStatErr_%s", fCentrality[i].Data()));
                            graphEtaToNPionEMCal2760GeVSys  = (TGraphAsymmErrors*)fileNeutralMeson2760GeV->Get(Form("graphEtaToPi0EMCalPbPb2760GeVSysErr_%s", fCentrality[i].Data()));

                            TH1D* histoEtaToNPionComb2760GeVStat = GraphToHist_withErrors(graphEtaToNPionComb2760GeVStat);
                            TH1D* histoEtaToNPionComb2760GeVSys = GraphToHist_withErrors(graphEtaToNPionComb2760GeVSys);
                            SetHistoProperties(histoEtaToNPionComb2760GeVStat, Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#eta/#pi^{0}", "");
                            SetHistoProperties(histoEtaToNPionComb2760GeVSys, Form("%sTo%s%sSys", fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#eta/#pi^{0}", "");
                            lists2760GeV[i]->Add(histoEtaToNPionComb2760GeVStat);
                            lists2760GeV[i]->Add(histoEtaToNPionComb2760GeVSys);

                            TH1D* histoEtaToNPionPCM2760GeVStat;
                            TH1D* histoEtaToNPionPCM2760GeVSys;
                            TH1D* histoEtaToNPionEMCal2760GeVStat;
                            TH1D* histoEtaToNPionEMCal2760GeVSys;

                            if(graphEtaToNPionPCM2760GeVStat){
                              histoEtaToNPionPCM2760GeVStat = GraphToHist_withErrors(graphEtaToNPionPCM2760GeVStat);
                              SetHistoProperties(histoEtaToNPionPCM2760GeVStat, Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                              lists2760GeV[i]->Add(histoEtaToNPionPCM2760GeVStat);
                            }
                            if(graphEtaToNPionPCM2760GeVSys){
                              histoEtaToNPionPCM2760GeVSys = GraphToHist_withErrors(graphEtaToNPionPCM2760GeVSys);
                              SetHistoProperties(histoEtaToNPionPCM2760GeVSys, Form("%sTo%s%sSys", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                              lists2760GeV[i]->Add(histoEtaToNPionPCM2760GeVSys);
                            }
                            if(graphEtaToNPionEMCal2760GeVStat){
                              histoEtaToNPionEMCal2760GeVStat = GraphToHist_withErrors(graphEtaToNPionEMCal2760GeVStat);
                              SetHistoProperties(histoEtaToNPionEMCal2760GeVStat, Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                              lists2760GeV[i]->Add(histoEtaToNPionEMCal2760GeVStat);
                            }
                            if(graphEtaToNPionEMCal2760GeVSys){
                              histoEtaToNPionEMCal2760GeVSys = GraphToHist_withErrors(graphEtaToNPionEMCal2760GeVSys);
                              SetHistoProperties(histoEtaToNPionEMCal2760GeVSys, Form("%sTo%s%sSys", fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                              lists2760GeV[i]->Add(histoEtaToNPionEMCal2760GeVSys);
                            }
                        }
                    }
                }

                //================================================================================================================
                // reading and writing omega to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[2]) {
                    cout << " - omega" << endl;
                }

                //================================================================================================================
                // reading and writing eta! to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[3]) {
                    cout << " - eta'" << endl;
                }

                //================================================================================================================
                // reading and writing gamma_dir to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[4]) {
                    cout << " - gamma_dir" << endl;
                }

                //================================================================================================================
                // reading and writing pi+-, K+- and p/bar{p} to 2760GeV list
                // input is given as fully invariant yield
                //================================================================================================================
                if (includeParticle[5]) {
                    // spectra
                    if (fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hstat_PbPb276_%s_pion_sum", fCentrality[i].Data())) && fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hsys_PbPb276_%s_pion_sum", fCentrality[i].Data()))) {
                        cout << " - pi^+/- spectrum" << endl;

                        TH1D* histoCPionComb2760GeVStat                = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hstat_PbPb276_%s_pion_sum", fCentrality[i].Data()));
                        TH1D* histoCPionComb2760GeVSys                 = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hsys_PbPb276_%s_pion_sum", fCentrality[i].Data()));

                        histoCPionComb2760GeVStat                      = ConvertYieldHisto(histoCPionComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoCPionComb2760GeVSys                       = ConvertYieldHisto(histoCPionComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoCPionComb2760GeVStat->Scale(0.5);
                        histoCPionComb2760GeVSys->Scale(0.5);

                        SetHistoProperties(histoCPionComb2760GeVStat, Form("%sStat", fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoCPionComb2760GeVSys, Form("%sSys", fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(histoCPionComb2760GeVStat);
                        lists2760GeV[i]->Add(histoCPionComb2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2CPion2760GeVStat = ReadFileforv2Input("Chargedpions",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2CPion2760GeVSyst = ReadFileforv2Input("Chargedpions",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2CPion2760GeVStat && graphv2CPion2760GeVSyst){
                            SetGraphProperties(graphv2CPion2760GeVStat,    "v2_CPionStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2CPion2760GeVSyst,    "v2_CPionSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2CPion2760GeVStat);
                            lists2760GeV[i]->Add(graphv2CPion2760GeVSyst);
                          }
                        }
                    }
                }

                if (includeParticle[6]) {
                    // spectra
                    if (fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hstat_PbPb276_%s_kaon_sum", fCentrality[i].Data())) && fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hsys_PbPb276_%s_kaon_sum", fCentrality[i].Data()))) {
                        cout << " - K^+/- spectrum" << endl;

                        TH1D* histoCKaonComb2760GeVStat                = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hstat_PbPb276_%s_kaon_sum", fCentrality[i].Data()));
                        TH1D* histoCKaonComb2760GeVSys                 = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hsys_PbPb276_%s_kaon_sum", fCentrality[i].Data()));

                        histoCKaonComb2760GeVStat                      = ConvertYieldHisto(histoCKaonComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoCKaonComb2760GeVSys                       = ConvertYieldHisto(histoCKaonComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoCKaonComb2760GeVStat->Scale(0.5);
                        histoCKaonComb2760GeVSys->Scale(0.5);

                        SetHistoProperties(histoCKaonComb2760GeVStat, Form("%sStat", fParticle[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoCKaonComb2760GeVSys, Form("%sSys", fParticle[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(histoCKaonComb2760GeVStat);
                        lists2760GeV[i]->Add(histoCKaonComb2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2CKaon2760GeVStat = ReadFileforv2Input("Chargedkaons",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2CKaon2760GeVSyst = ReadFileforv2Input("Chargedkaons",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2CKaon2760GeVStat && graphv2CKaon2760GeVSyst){
                            SetGraphProperties(graphv2CKaon2760GeVStat,    "v2_CKaonStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2CKaon2760GeVSyst,    "v2_CKaonSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2CKaon2760GeVStat);
                            lists2760GeV[i]->Add(graphv2CKaon2760GeVSyst);
                          }
                        }
                    }

                    // ratios
                    if (fileChargedPionKaonProton2760GeVRatios->GetListOfKeys()->Contains(Form("hstat_PbPb276_%s_kaon_to_pion_sum", fCentrality[i].Data())) && fileChargedPionKaonProton2760GeVRatios->GetListOfKeys()->Contains(Form("hsys_PbPb276_%s_kaon_to_pion_sum", fCentrality[i].Data()))) {
                        cout << " - K^+/-/pi^+/-  ratio" << endl;

                        TH1D* histoCKaonToCPionComb2760GeVStat         = (TH1D*)fileChargedPionKaonProton2760GeVRatios->Get(Form("hstat_PbPb276_%s_kaon_to_pion_sum", fCentrality[i].Data()));
                        TH1D* histoCKaonToCPionComb2760GeVSys          = (TH1D*)fileChargedPionKaonProton2760GeVRatios->Get(Form("hsys_PbPb276_%s_kaon_to_pion_sum", fCentrality[i].Data()));

                        SetHistoProperties(histoCKaonToCPionComb2760GeVStat, Form("%sTo%sStat", fParticle[6].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");
                        SetHistoProperties(histoCKaonToCPionComb2760GeVSys, Form("%sTo%sSys", fParticle[6].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");

                        lists2760GeV[i]->Add(histoCKaonToCPionComb2760GeVStat);
                        lists2760GeV[i]->Add(histoCKaonToCPionComb2760GeVSys);
                    }
                }

                if (includeParticle[7]) {
                    // spectra
                    if (fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hstat_PbPb276_%s_proton_sum", fCentrality[i].Data())) && fileChargedPionKaonProton2760GeV->GetListOfKeys()->Contains(Form("hsys_PbPb276_%s_proton_sum", fCentrality[i].Data()))) {
                        cout << " - p/anti-p spectrum" << endl;

                        TH1D* histoProtonComb2760GeVStat               = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hstat_PbPb276_%s_proton_sum", fCentrality[i].Data()));
                        TH1D* histoProtonComb2760GeVSys                = (TH1D*)fileChargedPionKaonProton2760GeV->Get(Form("hsys_PbPb276_%s_proton_sum", fCentrality[i].Data()));

                        histoProtonComb2760GeVStat                     = ConvertYieldHisto(histoProtonComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        histoProtonComb2760GeVSys                      = ConvertYieldHisto(histoProtonComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        // scale by 0.5 to get averaged single particle
                        histoProtonComb2760GeVStat->Scale(0.5);
                        histoProtonComb2760GeVSys->Scale(0.5);

                        SetHistoProperties(histoProtonComb2760GeVStat, Form("%sStat", fParticle[7].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoProtonComb2760GeVSys, Form("%sSys", fParticle[7].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(histoProtonComb2760GeVStat);
                        lists2760GeV[i]->Add(histoProtonComb2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2Proton2760GeVStat = ReadFileforv2Input("ProtonsAntiProtons",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2Proton2760GeVSyst = ReadFileforv2Input("ProtonsAntiProtons",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2Proton2760GeVStat && graphv2Proton2760GeVSyst){
                            SetGraphProperties(graphv2Proton2760GeVStat,    "v2_ProtonStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2Proton2760GeVSyst,    "v2_ProtonSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2Proton2760GeVStat);
                            lists2760GeV[i]->Add(graphv2Proton2760GeVSyst);
                          }
                        }
                    }

                    // ratios
                    if (fileChargedPionKaonProton2760GeVRatios->GetListOfKeys()->Contains(Form("hstat_PbPb276_%s_proton_to_pion_sum", fCentrality[i].Data())) && fileChargedPionKaonProton2760GeVRatios->GetListOfKeys()->Contains(Form("hsys_PbPb276_%s_proton_to_pion_sum", fCentrality[i].Data()))) {
                        cout << " - p/pi^+/- ratio" << endl;

                        TH1D* histoProtonToCPionComb2760GeVStat        = (TH1D*)fileChargedPionKaonProton2760GeVRatios->Get(Form("hstat_PbPb276_%s_proton_to_pion_sum", fCentrality[i].Data()));
                        TH1D* histoProtonToCPionComb2760GeVSys         = (TH1D*)fileChargedPionKaonProton2760GeVRatios->Get(Form("hsys_PbPb276_%s_proton_to_pion_sum", fCentrality[i].Data()));

                        SetHistoProperties(histoProtonToCPionComb2760GeVStat, Form("%sTo%sStat", fParticle[7].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");
                        SetHistoProperties(histoProtonToCPionComb2760GeVSys, Form("%sTo%sSys", fParticle[7].Data(), fParticle[5].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");

                        lists2760GeV[i]->Add(histoProtonToCPionComb2760GeVStat);
                        lists2760GeV[i]->Add(histoProtonToCPionComb2760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing h+- to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[8]) {
                    cout << " - h^+/-" << endl;
                }

                //================================================================================================================
                // reading and writing phi to 2760GeV list
                // input is given as 1/2pi 1/pt dN/dydpt and scaled for clarity by 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-3 for each cent. (skipping 40-50%)
                //================================================================================================================
                if (includeParticle[9]) {
                    // spectra
                    if (filePhiMeson2760GeV->GetListOfKeys()->Contains(Form("gr_%s_stat", fCentralityOpt4[i].Data())) && filePhiMeson2760GeV->GetListOfKeys()->Contains(Form("gr_%s_syst", fCentralityOpt4[i].Data()))) {
                        cout << " - phi spectrum" << endl;

                        TGraphErrors* graphPhi2760GeVStat        = (TGraphErrors*)filePhiMeson2760GeV->Get(Form("gr_%s_stat", fCentralityOpt4[i].Data()));
                        TGraphErrors* graphPhi2760GeVSys         = (TGraphErrors*)filePhiMeson2760GeV->Get(Form("gr_%s_syst", fCentralityOpt4[i].Data()));

                        if(i==1){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e3);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e3);
                        } else if(i==2){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e2);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e2);
                        } else if(i==4){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e1);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e1);
                        } else if(i==12){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e0);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e0);
                        } else if(i==13){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e-1);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e-1);
                        } else if(i==9){
                          graphPhi2760GeVStat = ScaleGraph(graphPhi2760GeVStat, 1./1e-3);
                          graphPhi2760GeVSys = ScaleGraph(graphPhi2760GeVSys, 1./1e-3);
                        }
                        for(Int_t n=0; n<graphPhi2760GeVStat->GetN(); n++){
                            graphPhi2760GeVStat->SetPointError(n,graphPhi2760GeVSys->GetErrorX(n),graphPhi2760GeVStat->GetErrorY(n));
                        }
                        graphPhi2760GeVStat                      = ConvertYieldGraph(graphPhi2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        graphPhi2760GeVSys                       = ConvertYieldGraph(graphPhi2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        TH1D* histPhi2760GeVStat = GraphToHist_withErrors(graphPhi2760GeVStat);
                        TH1D* histPhi2760GeVSys = GraphToHist_withErrors(graphPhi2760GeVSys);
                        SetHistoProperties(histPhi2760GeVStat, Form("%sStat", fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histPhi2760GeVSys, Form("%sSys", fParticle[9].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        lists2760GeV[i]->Add(histPhi2760GeVStat);
                        lists2760GeV[i]->Add(histPhi2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2Phi2760GeVStat = ReadFileforv2Input("Phi",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2Phi2760GeVSyst = ReadFileforv2Input("Phi",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2Phi2760GeVStat && graphv2Phi2760GeVSyst){
                            SetGraphProperties(graphv2Phi2760GeVStat,    "v2_PhiStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2Phi2760GeVSyst,    "v2_PhiSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2Phi2760GeVStat);
                            lists2760GeV[i]->Add(graphv2Phi2760GeVSyst);
                          }
                        }
                    }
                }

                //================================================================================================================
                // reading and writing K^*0 to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[10]) {
                    cout << " - K^*0" << endl;

                    //
                    // should exist
                    //
                }

                //================================================================================================================
                // reading and writing rho^0 to 2760GeV list
                // input is given as 1/2pi 1/pt dN/dydpt
                //================================================================================================================
                if (includeParticle[11]) {
                    // spectra
                    if (fileRho0Meson2760GeV->GetListOfKeys()->Contains(Form("c%s_stat", fCentralityOpt4[i].Data())) && fileRho0Meson2760GeV->GetListOfKeys()->Contains(Form("c%s_sys", fCentralityOpt4[i].Data()))) {
                        cout << " - rho^0 spectrum" << endl;

                        TGraphErrors* graphRho02760GeVStat        = (TGraphErrors*)fileRho0Meson2760GeV->Get(Form("c%s_stat", fCentralityOpt4[i].Data()));
                        TGraphAsymmErrors* graphRho02760GeVSys    = (TGraphAsymmErrors*)fileRho0Meson2760GeV->Get(Form("c%s_sys", fCentralityOpt4[i].Data()));
                        if(i==9){
                          graphRho02760GeVStat->RemovePoint(graphRho02760GeVStat->GetN()-1);
                          graphRho02760GeVSys->RemovePoint(graphRho02760GeVSys->GetN()-1);
                        }

                        for(Int_t n=0; n<graphRho02760GeVStat->GetN(); n++){
                            graphRho02760GeVStat->SetPointError(n,graphRho02760GeVSys->GetErrorX(n),graphRho02760GeVStat->GetErrorY(n));
                        }
                        graphRho02760GeVStat                      = ConvertYieldGraph(graphRho02760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                        graphRho02760GeVSys                       = ConvertYieldGraph(graphRho02760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                        TH1D* histNRho2760GeVStat = GraphToHist_withErrors(graphRho02760GeVStat);
                        TH1D* histNRho2760GeVSys = GraphToHist_withErrors(graphRho02760GeVSys);

                        SetHistoProperties(histNRho2760GeVStat, Form("%sStat", fParticle[11].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histNRho2760GeVSys, Form("%sSys", fParticle[11].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        lists2760GeV[i]->Add(histNRho2760GeVStat);
                        lists2760GeV[i]->Add(histNRho2760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing rho^+/- to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[12]) {
                    cout << " - rho^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing Delta^0 to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[13]) {
                    cout << " - Delta^0" << endl;

                }

                //================================================================================================================
                // reading and writing Delta^+/- to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[14]) {
                    cout << " - Delta^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing K^0_s to 2760GeV list
                // input is given as yield (dN/d#it{y}dpT)
                //================================================================================================================
                if (includeParticle[15]) {
                    // spectra
                    if (fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("statonly_cent%s_K0s", fCentrality[i].Data())) && fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("systonly_cent%s_K0s", fCentrality[i].Data()))) {
                        cout << " - K^0_s spectrum" << endl;

                        TH1D* histoK0s2760GeVStat               = (TH1D*)fileK0sLambda2760GeV->Get(Form("statonly_cent%s_K0s", fCentrality[i].Data()));
                        TH1D* histoK0s2760GeVSys                = (TH1D*)fileK0sLambda2760GeV->Get(Form("systonly_cent%s_K0s", fCentrality[i].Data()));

                        SetHistoProperties(histoK0s2760GeVStat, Form("%sStat", fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoK0s2760GeVSys, Form("%sSys", fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(histoK0s2760GeVStat);
                        lists2760GeV[i]->Add(histoK0s2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2K0s2760GeVStat = ReadFileforv2Input("K0s",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2K0s2760GeVSyst = ReadFileforv2Input("K0s",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2K0s2760GeVStat && graphv2K0s2760GeVSyst){
                            SetGraphProperties(graphv2K0s2760GeVStat,    "v2_K0sStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2K0s2760GeVSyst,    "v2_K0sSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2K0s2760GeVStat);
                            lists2760GeV[i]->Add(graphv2K0s2760GeVSyst);
                          }
                        }
                    }
                }

                //================================================================================================================
                // reading and writing Lambda to 2760GeV list
                // input is given as yield (dN/d#it{y}dpT)
                //================================================================================================================
                if (includeParticle[16]) {
                    // spectra
                    if (fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("statonly_cent%s_Lambda", fCentrality[i].Data())) && fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("systonly_cent%s_Lambda", fCentrality[i].Data())) /*&& fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("fHistPtAntiLambda_%s_OnlyStat", fCentrality[i].Data())) && fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("fHistPtAntiLambda_%s_OnlySyst", fCentrality[i].Data()))*/ ) {
                        cout << " - Lambda spectrum" << endl;

                        TH1D* histoLambda2760GeVStat                   = (TH1D*)fileK0sLambda2760GeV->Get(Form("statonly_cent%s_Lambda", fCentrality[i].Data()));
                        TH1D* histoLambda2760GeVSys                    = (TH1D*)fileK0sLambda2760GeV->Get(Form("systonly_cent%s_Lambda", fCentrality[i].Data()));
//                         TH1D* histoAntiLambda2760GeVStat               = (TH1D*)fileK0sLambda2760GeV->Get(Form("fHistPtAntiLambda_%s_OnlyStat", fCentrality[i].Data()));
//                         TH1D* histoAntiLambda2760GeVSys                = (TH1D*)fileK0sLambda2760GeV->Get(Form("fHistPtAntiLambda_%s_OnlySyst", fCentrality[i].Data()));

//                         histoLambda2760GeVStat->Add(histoAntiLambda2760GeVStat);
//                         histoLambda2760GeVSys->Add(histoAntiLambda2760GeVSys);

//                         histoLambda2760GeVStat->Scale(0.5);
//                         histoLambda2760GeVSys->Scale(0.5);

                        SetHistoProperties(histoLambda2760GeVStat, Form("%sStat", fParticle[16].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetHistoProperties(histoLambda2760GeVSys, Form("%sSys", fParticle[16].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(histoLambda2760GeVStat);
                        lists2760GeV[i]->Add(histoLambda2760GeVSys);

                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2Lambda2760GeVStat = ReadFileforv2Input("LambdaAntiLambda",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2Lambda2760GeVSyst = ReadFileforv2Input("LambdaAntiLambda",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2Lambda2760GeVStat && graphv2Lambda2760GeVSyst){
                            SetGraphProperties(graphv2Lambda2760GeVStat,    "v2_LambdaStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2Lambda2760GeVSyst,    "v2_LambdaSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2Lambda2760GeVStat);
                            lists2760GeV[i]->Add(graphv2Lambda2760GeVSyst);
                          }
                        }
                    }

                    // ratios
                    if (fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("statonly_cent%s_LambdaK0s", fCentrality[i].Data())) && fileK0sLambda2760GeV->GetListOfKeys()->Contains(Form("systonly_cent%s_LambdaK0s", fCentrality[i].Data()))) {
                        cout << " - Lambda/K0s ratio" << endl;

                        TH1D* histoLambdaToKaon2760GeVStat             = (TH1D*)fileK0sLambda2760GeV->Get(Form("statonly_cent%s_LambdaK0s", fCentrality[i].Data()));
                        TH1D* histoLambdaToKaon2760GeVSys              = (TH1D*)fileK0sLambda2760GeV->Get(Form("systonly_cent%s_LambdaK0s", fCentrality[i].Data()));

                        SetHistoProperties(histoLambdaToKaon2760GeVStat, Form("%sTo%sStat", fParticle[16].Data(), fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#Lambda / K^{0}_{s}", "");
                        SetHistoProperties(histoLambdaToKaon2760GeVSys, Form("%sTo%sSys", fParticle[16].Data(), fParticle[15].Data()), "#it{p}_{T} (GeV/#it{c})", "#Lambda / K^{0}_{s}", "");

                        lists2760GeV[i]->Add(histoLambdaToKaon2760GeVStat);
                        lists2760GeV[i]->Add(histoLambdaToKaon2760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing Sigma^0 to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[17]) {
                    cout << " - Sigma^0" << endl;

                }

                //================================================================================================================
                // reading and writing Sigma^+/- to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[18]) {
                    cout << " - Sigma^+/-" << endl;

                }

                //================================================================================================================
                // reading and writing Omega^+/- to 2760GeV list
                // input is given as (dN/d#it{y}dpT)
                //================================================================================================================
                if (includeParticle[19]) {
                    // spectra



                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2Omega2760GeVStat = ReadFileforv2Input("OmegaAntiOmega",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2Omega2760GeVSyst = ReadFileforv2Input("OmegaAntiOmega",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2Omega2760GeVStat && graphv2Omega2760GeVSyst){
                            SetGraphProperties(graphv2Omega2760GeVStat,    "v2_XiStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2Omega2760GeVSyst,    "v2_XiSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2Omega2760GeVStat);
                            lists2760GeV[i]->Add(graphv2Omega2760GeVSyst);
                          }
                        }
                }

                //================================================================================================================
                // reading and writing Xi^+/- to 2760GeV list
                // input is given as yield (dN/d#it{y}dpT)
                //================================================================================================================
                if (includeParticle[20]) {
                    // spectra



                        if (enableFlow) {
                          TGraphAsymmErrors* graphv2Xi2760GeVStat = ReadFileforv2Input("XiAntiXi",Form("%s",fCentrality[i].Data()),kTRUE);
                          TGraphAsymmErrors* graphv2Xi2760GeVSyst = ReadFileforv2Input("XiAntiXi",Form("%s",fCentrality[i].Data()),kFALSE);
                          if(graphv2Xi2760GeVStat && graphv2Xi2760GeVSyst){
                            SetGraphProperties(graphv2Xi2760GeVStat,    "v2_XiStat", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            SetGraphProperties(graphv2Xi2760GeVSyst,    "v2_XiSyst", "#it{p}_{T} (GeV/#it{c})", "v_{2}", "");
                            lists2760GeV[i]->Add(graphv2Xi2760GeVStat);
                            lists2760GeV[i]->Add(graphv2Xi2760GeVSyst);
                          }
                        }
                }

                //================================================================================================================
                // reading and writing J/psi to 2760GeV list
                // input is given as ???
                //================================================================================================================
                if (includeParticle[21]) {
                    cout << " - J/psi" << endl;

                }

                //================================================================================================================
                // reading and writing D0 to 2760GeV list
                // input is given as dN/dy/dpt
                //================================================================================================================
                if (includeParticle[22]) {
                    if ( i == 3){
                        cout << " - D0" << endl;
                        TGraphAsymmErrors* graphD02760GeVStat           = (TGraphAsymmErrors*)GetD0MesonPbPb0010(kTRUE);
                        TGraphAsymmErrors* graphD02760GeVSys            = (TGraphAsymmErrors*)GetD0MesonPbPb0010(kFALSE);
                        SetGraphProperties(graphD02760GeVStat,    Form("%sStat", fParticle[22].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphD02760GeVSys,     Form("%sSys",  fParticle[22].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(graphD02760GeVStat);
                        lists2760GeV[i]->Add(graphD02760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing D+ to 2760GeV list
                // input is given as dN/dy/dpt
                //================================================================================================================
                if (includeParticle[23]) {
                    if ( i == 3){
                        cout << " - D+" << endl;
                        TGraphAsymmErrors* graphDPlus2760GeVStat        = (TGraphAsymmErrors*)GetDPlusMesonPbPb0010(kTRUE);
                        TGraphAsymmErrors* graphDPlus2760GeVSys         = (TGraphAsymmErrors*)GetDPlusMesonPbPb0010(kFALSE);
                        SetGraphProperties(graphDPlus2760GeVStat,    Form("%sStat", fParticle[23].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDPlus2760GeVSys,     Form("%sSys",  fParticle[23].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(graphDPlus2760GeVStat);
                        lists2760GeV[i]->Add(graphDPlus2760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing D*+ to 2760GeV list
                // input is given as dN/dy/dpt
                //================================================================================================================
                if (includeParticle[24]) {
                    if ( i == 3){
                        cout << " - D*+" << endl;
                        TGraphAsymmErrors* graphDStarPlus2760GeVStat    = (TGraphAsymmErrors*)GetDStarPlusMesonPbPb0010(kTRUE);
                        TGraphAsymmErrors* graphDStarPlus2760GeVSys     = (TGraphAsymmErrors*)GetDStarPlusMesonPbPb0010(kFALSE);
                        SetGraphProperties(graphDStarPlus2760GeVStat,    Form("%sStat", fParticle[24].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDStarPlus2760GeVSys,     Form("%sSys",  fParticle[24].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(graphDStarPlus2760GeVStat);
                        lists2760GeV[i]->Add(graphDStarPlus2760GeVSys);
                    }
                }

                //================================================================================================================
                // reading and writing D*+ to 2760GeV list
                // input is given as dN/dy/dpt
                //================================================================================================================
                if (includeParticle[25]) {
                    if ( i == 3){
                        cout << " - D_s^+" << endl;
                        TGraphAsymmErrors* graphDSPlus2760GeVStat       = (TGraphAsymmErrors*)GetDSPlusMesonPbPb0010(kTRUE);
                        TGraphAsymmErrors* graphDSPlus2760GeVSys        = (TGraphAsymmErrors*)GetDSPlusMesonPbPb0010(kFALSE);
                        SetGraphProperties(graphDSPlus2760GeVStat,    Form("%sStat", fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDSPlus2760GeVSys,     Form("%sSys",  fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(graphDSPlus2760GeVStat);
                        lists2760GeV[i]->Add(graphDSPlus2760GeVSys);
                    }
                    if ( i == 7){
                        cout << " - D_s^+" << endl;
                        TGraphAsymmErrors* graphDSPlus2760GeVStat             = (TGraphAsymmErrors*)GetDSPlusMesonPbPb2050(kTRUE);
                        TGraphAsymmErrors* graphDSPlus2760GeVSys              = (TGraphAsymmErrors*)GetDSPlusMesonPbPb2050(kFALSE);
                        SetGraphProperties(graphDSPlus2760GeVStat,    Form("%sStat", fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");
                        SetGraphProperties(graphDSPlus2760GeVSys,     Form("%sSys",  fParticle[25].Data()), "#it{p}_{T} (GeV/#it{c})", "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})", "");

                        lists2760GeV[i]->Add(graphDSPlus2760GeVStat);
                        lists2760GeV[i]->Add(graphDSPlus2760GeVSys);
                    }
                }
            }
        }
    }

    //================================================================================================================
    //Produce particle spectra and ratio for 0-10% centrality
    //================================================================================================================

    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[5], "", "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[6], "", "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[7], "", "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[9], "", "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[15], "", "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[16], "", "");

    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[6], fParticle[9], "", "");

    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[6], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[7], fParticle[5], "");
    ProduceSpectrumInCentralityBin(lists2760GeV, fCentrality[3], fParticle[16], fParticle[15], "");

    //================================================================================================================
    //Produce plots containing all particle spectra if specified
    //================================================================================================================
    if (produceAllSpectraPlots) {
        if (enable2760GeV) {
            for (Int_t i=0; i<nCentralities; i++) {
                if (includeCentrality2760GeV[i]) {
                    cout << "plotting " << lists2760GeV[i]->GetName() << endl;
                    ProduceParticleSpectraPlotFromList(lists2760GeV[i], fCollSys[2], fEnergy[1], fCentrality[i], suffix);
                    ProduceParticleSpectraPlotFromListOnlyFinal(lists2760GeV[i], fCollSys[2], fEnergy[1], fCentrality[i], suffix);
                }
            }
        }
    }

    //================================================================================================================
    //Produce plots containing all particle ratio if specified
    //================================================================================================================
    if (produceAllSpectraPlots) {
        if (enable2760GeV) {
            for (Int_t i=0; i<nCentralities; i++) {
                if (includeCentrality2760GeV[i]) {
                    cout << "plotting " << lists2760GeV[i]->GetName() << endl;
                    ProduceParticleRatioPlotFromList(lists2760GeV[i], fCollSys[2], fEnergy[1], fCentrality[i], suffix);
                    ProduceParticleRatioPlotFromListOnlyFinal(lists2760GeV[i], fCollSys[2], fEnergy[1], fCentrality[i], suffix);
                }
            }
        }
    }

    //================================================================================================================
    //Saving the TLists to the final file
    //================================================================================================================
    cout << "writing lists" << endl;
    outputFile->cd();
    if (enable2760GeV) {
        for (Int_t i=0; i<nCentralities; i++) {
            if (includeCentrality2760GeV[i]) {
                lists2760GeV[i]->Write(Form("%s_%s_%s", fCollSys[2].Data(), fEnergy[1].Data(), fCentrality[i].Data()), TObject::kSingleKey);
            }
        }
    }
    outputFile->Close();

}
