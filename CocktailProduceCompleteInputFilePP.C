/****************************************************************************************************************************
******    Friederike Bock, friederike.bock@cern.ch                                                                      *****
******    Mike Sas, mike.sas@cern.ch                                                                                    *****
******    Lucas Altenkaemper, lucas.altenkaemper@cern.ch                                                                *****
*****************************************************************************************************************************/
#include "CocktailFunctions.h"
#include "CocktailPlotting.h"
//#include "CocktailFitting.h" // already included in cocktail plotting
#include "CocktailHEPDataPP.h"
#include "TDatabasePDG.h"
#include <iostream>

void  CocktailProduceCompleteInputFilePP(   TString enableEnergy        = "110110",
                                            Bool_t  enable_NPion        = 1,
                                            Bool_t  enable_CPion        = 1,
                                            Bool_t  enable_Eta          = 1,
                                            Bool_t  enable_CKaon        = 1,
                                            Bool_t  enable_OmegaMeson   = 1,
                                            Bool_t  enable_Proton       = 1,
                                            Bool_t  enable_Phi          = 1,
                                            Bool_t  enable_K0s          = 1,
                                            Bool_t  enable_Lambda       = 1,
                                            Bool_t  enable_K0star       = 1,
                                            Bool_t  enable_COmega       = 1,
                                            Bool_t  enable_CXi          = 1,
                                            Bool_t  enable_CSigma       = 1,
                                            Bool_t  enable_JPsi         = 1,
                                            Bool_t  enable_Rho0         = 1,
                                            Bool_t  enablePtY           = 1,
                                            TString convertTo           = "dN/dydpT",  // "dN/dydpT" or "1/2pipT*dN/dydpT"
                                            TString suffix              = "eps",
                                            Bool_t parametrizeMC        = kFALSE
                                        ){

    //================================================================================================================
    //Enable energies to be included
    //================================================================================================================

    //This is done by changing the "enableEnergy" string
    //"enableEnergy" has 14 digits, set each digit either 0(exclude) or 1(include)
    //Example: "110100" Enables energies 900GeV,2.76TeV,7TeV

    TString String_900GeV   = enableEnergy(0,1);
    TString String_2760GeV  = enableEnergy(1,1);
    TString String_5TeV     = enableEnergy(2,1);
    TString String_7TeV     = enableEnergy(3,1);
    TString String_8TeV     = enableEnergy(4,1);
    TString String_13TeV    = enableEnergy(5,1);

    Bool_t Include_900GeV   = kFALSE;
    Bool_t Include_2760GeV  = kFALSE;
    Bool_t Include_5TeV     = kFALSE;
    Bool_t Include_7TeV     = kFALSE;
    Bool_t Include_8TeV     = kFALSE;
    Bool_t Include_13TeV    = kFALSE;

    if(String_900GeV.CompareTo("1") == 0)   Include_900GeV  = kTRUE;
    if(String_2760GeV.CompareTo("1") == 0)  Include_2760GeV = kTRUE;
    if(String_5TeV.CompareTo("1") == 0)     Include_5TeV    = kTRUE;
    if(String_7TeV.CompareTo("1") == 0)     Include_7TeV    = kTRUE;
    if(String_8TeV.CompareTo("1") == 0)     Include_8TeV    = kTRUE;
    if(String_13TeV.CompareTo("1") == 0)    Include_13TeV   = kTRUE;

    //================================================================================================================
    //Creating output file structure
    //================================================================================================================

    TFile *output_File      = new TFile("CocktailInputPP.root","RECREATE");
    TList *list_900GeV      = NULL;
    if (Include_900GeV){
        list_900GeV         = new TList();
        list_900GeV->SetName("pp_0.9TeV");
    }
    TList *list_2760GeV     = NULL;
    if (Include_2760GeV){
        list_2760GeV        = new TList();
        list_2760GeV->SetName("pp_2.76TeV");
    }
    TList *list_5TeV        = NULL;
    if (Include_5TeV){
        list_5TeV           = new TList();
        list_5TeV->SetName("pp_5TeV");
    }
    TList *list_7TeV        = NULL;
    if (Include_7TeV){
        list_7TeV           = new TList();
        list_7TeV->SetName("pp_7TeV");
    }
    TList *list_8TeV        = NULL;
    if (Include_8TeV){
        list_8TeV           = new TList();
        list_8TeV->SetName("pp_8TeV");
    }
    TList *list_13TeV       = NULL;
    if (Include_13TeV){
        list_13TeV          = new TList();
        list_13TeV->SetName("pp_13TeV");
    }

    //================================================================================================================
    // Which kind of conversion of the spectra is needed?
    //================================================================================================================
    TString globalOutputQuantity                = "";
    Int_t changeQuantity                        = 0;
    if ( convertTo.CompareTo("")==0 ){
        changeQuantity                          = 0;
    } else if ( convertTo.CompareTo("dN/dydpT")==0 || convertTo.CompareTo("dn/dydpt")==0 ){
        changeQuantity                          = 1;
        globalOutputQuantity                    = "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})";
    } else if ( convertTo.CompareTo("1/2pipT*dN/dydpT")==0 || convertTo.CompareTo("1/2pipt*dn/dydpt")==0 ){
        changeQuantity                          = 2;
        globalOutputQuantity                    = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";
    }

    //================================================================================================================
    // Set correct xSection for pp
    //================================================================================================================
    // pp 0.9 TeV
    Double_t xSection900GeVINEL                 = 52.5*1e9;
    Double_t xSection900GeVV0OR                 = 47.78*1e9;
    Double_t xSection900GeVV0AND                = 40.06*1e9;
    // pp 2.76 TeV
    Double_t xSection2760GeVOR                  = 55.416*1e9;
    Double_t xSection2760GeVV0AND               = 47.73*1e9;
    Double_t xSection2760GeVINEL                = 62.8*1e9;
    // pp 5.023TeV
    Double_t xSection5023GeVV0AND               = 51.2*1e9;    // from https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/mgagliar/2016-Aug-22-paper_draft-vdmNote_5TeV.pdf
    // pp 7 TeV
    Double_t xSection7TeVINEL                   = 73.2*1e9;
    Double_t xSection7TeVOR                     = 62.22*1e9;
    Double_t xSection7TeVV0AND                  = 54.31*1e9;
    // pp 8 TeV
    Double_t xSection8TeVINEL                   = 72.3*1e9;   // from https://aliceinfo.cern.ch/Notes/node/665
    Double_t xSection8TeVV0AND                  = 55.8*1e9;   // from https://aliceinfo.cern.ch/Notes/node/531
    Double_t xSection8TeVT0AND                  = 25.5*1e9;   // from https://aliceinfo.cern.ch/Notes/node/531

    //================================================================================================================
    // creating histos and graphs for 900GeV
    //================================================================================================================
    if (Include_900GeV){
        if (enable_NPion || enable_Eta){
            TString localOutputQuantity                                 = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";
            if (changeQuantity != 0){
                localOutputQuantity                                     = globalOutputQuantity;
            }

            TFile*          fFileNeutralMeson900GeV                     = new TFile("pp/CombinedResultsPaperPP900GeV_2017_04_05.root");
            TDirectoryFile* fFileNeutralMeson900GeVPi0                  = (TDirectoryFile*)fFileNeutralMeson900GeV->Get("Pi0900GeV");
            TDirectoryFile* fFileNeutralMeson900GeVEta                  = (TDirectoryFile*)fFileNeutralMeson900GeV->Get("Eta900GeV");

            //================================================================================================================
            // reading and writing pi0 to 900GeV list
            // input from pi0 is given as invariant cross section shifted in X direction
            //================================================================================================================
            if (enable_NPion){
                TGraphAsymmErrors* graphNPionInvYieldComb900GeVStat             = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0CombStat");
                TGraphAsymmErrors* graphNPionInvYieldComb900GeVSys              = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0CombSys");
                TGraphAsymmErrors* graphNPionInvYieldPCM900GeVStat              = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0PCMStat");
                TGraphAsymmErrors* graphNPionInvYieldPCM900GeVSys               = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0PCMSys");
                TGraphAsymmErrors* graphNPionInvYieldPHOS900GeVStat             = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0PHOSStat");
                TGraphAsymmErrors* graphNPionInvYieldPHOS900GeVSys              = (TGraphAsymmErrors*)fFileNeutralMeson900GeVPi0->Get("graphInvYieldPi0PHOSSys");

                if (changeQuantity == 0){
                    TGraphAsymmErrors* graphNPionXSecComb900GeVStat             = (TGraphAsymmErrors*)graphNPionInvYieldComb900GeVStat->Clone(  "graphInvCrosSecPi0CombStat");
                    TGraphAsymmErrors* graphNPionXSecComb900GeVSys              = (TGraphAsymmErrors*)graphNPionInvYieldComb900GeVSys->Clone(   "graphInvCrosSecPi0CombSys");
                    TGraphAsymmErrors* graphNPionXSecPCM900GeVStat              = (TGraphAsymmErrors*)graphNPionInvYieldPCM900GeVStat->Clone(   "graphInvCrosSecPi0PCMStat");
                    TGraphAsymmErrors* graphNPionXSecPCM900GeVSys               = (TGraphAsymmErrors*)graphNPionInvYieldPCM900GeVSys->Clone(    "graphInvCrosSecPi0PCMSys");
                    TGraphAsymmErrors* graphNPionXSecPHOS900GeVStat             = (TGraphAsymmErrors*)graphNPionInvYieldPHOS900GeVStat->Clone(  "graphInvCrosSecPi0PHOSStat");
                    TGraphAsymmErrors* graphNPionXSecPHOS900GeVSys              = (TGraphAsymmErrors*)graphNPionInvYieldPHOS900GeVSys->Clone(   "graphInvCrosSecPi0PHOSSys");

                    graphNPionXSecComb900GeVStat                                = ScaleGraph(graphNPionXSecComb900GeVStat,  xSection900GeVV0OR);
                    graphNPionXSecComb900GeVSys                                 = ScaleGraph(graphNPionXSecComb900GeVSys,   xSection900GeVV0OR);
                    graphNPionXSecPCM900GeVStat                                 = ScaleGraph(graphNPionXSecPCM900GeVStat,   xSection900GeVV0OR);
                    graphNPionXSecPCM900GeVSys                                  = ScaleGraph(graphNPionXSecPCM900GeVSys,    xSection900GeVV0OR);
                    graphNPionXSecPHOS900GeVStat                                = ScaleGraph(graphNPionXSecPHOS900GeVStat,  xSection900GeVV0OR);
                    graphNPionXSecPHOS900GeVSys                                 = ScaleGraph(graphNPionXSecPHOS900GeVSys,   xSection900GeVV0OR);

                    SetGraphProperties(graphNPionXSecComb900GeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecComb900GeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCM900GeVStat,     Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCM900GeVSys,      Form("%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPHOS900GeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPHOS900GeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    // simply add Xsection to the output file
                    list_900GeV->Add(graphNPionXSecComb900GeVStat);
                    list_900GeV->Add(graphNPionXSecComb900GeVSys);
                    list_900GeV->Add(graphNPionXSecPCM900GeVStat);
                    list_900GeV->Add(graphNPionXSecPCM900GeVSys);
                    list_900GeV->Add(graphNPionXSecPHOS900GeVStat);
                    list_900GeV->Add(graphNPionXSecPHOS900GeVSys);
                } else {
                    // convert to INEL events
                    graphNPionInvYieldComb900GeVStat                            = ScaleGraph(graphNPionInvYieldComb900GeVStat,  xSection900GeVV0OR / xSection900GeVINEL);
                    graphNPionInvYieldComb900GeVSys                             = ScaleGraph(graphNPionInvYieldComb900GeVSys,   xSection900GeVV0OR / xSection900GeVINEL);
                    graphNPionInvYieldPCM900GeVStat                             = ScaleGraph(graphNPionInvYieldPCM900GeVStat,   xSection900GeVV0OR / xSection900GeVINEL);
                    graphNPionInvYieldPCM900GeVSys                              = ScaleGraph(graphNPionInvYieldPCM900GeVSys,    xSection900GeVV0OR / xSection900GeVINEL);
                    graphNPionInvYieldPHOS900GeVStat                            = ScaleGraph(graphNPionInvYieldPHOS900GeVStat,  xSection900GeVV0OR / xSection900GeVINEL);
                    graphNPionInvYieldPHOS900GeVSys                             = ScaleGraph(graphNPionInvYieldPHOS900GeVSys,   xSection900GeVV0OR / xSection900GeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphNPionInvYieldComb900GeVStat                            = ConvertYieldGraph(graphNPionInvYieldComb900GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionInvYieldComb900GeVSys                             = ConvertYieldGraph(graphNPionInvYieldComb900GeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionInvYieldPCM900GeVStat                             = ConvertYieldGraph(graphNPionInvYieldPCM900GeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionInvYieldPCM900GeVSys                              = ConvertYieldGraph(graphNPionInvYieldPCM900GeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionInvYieldPHOS900GeVStat                            = ConvertYieldGraph(graphNPionInvYieldPHOS900GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionInvYieldPHOS900GeVSys                             = ConvertYieldGraph(graphNPionInvYieldPHOS900GeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphNPionInvYieldComb900GeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionInvYieldComb900GeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionInvYieldPCM900GeVStat,     Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionInvYieldPCM900GeVSys,      Form("%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionInvYieldPHOS900GeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionInvYieldPHOS900GeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_900GeV->Add(graphNPionInvYieldComb900GeVStat);
                    list_900GeV->Add(graphNPionInvYieldComb900GeVSys);
                    list_900GeV->Add(graphNPionInvYieldPCM900GeVStat);
                    list_900GeV->Add(graphNPionInvYieldPCM900GeVSys);
                    list_900GeV->Add(graphNPionInvYieldPHOS900GeVStat);
                    list_900GeV->Add(graphNPionInvYieldPHOS900GeVSys);
                }
            }
            if (enable_Eta){
                TGraphAsymmErrors* graphEtaInvYieldPCM900GeVStat                = (TGraphAsymmErrors*)fFileNeutralMeson900GeVEta->Get("graphInvYieldEtaPCMStat");
                TGraphAsymmErrors* graphEtaInvYieldPCM900GeVSys                 = (TGraphAsymmErrors*)fFileNeutralMeson900GeVEta->Get("graphInvYieldEtaPCMSys");

                TGraphAsymmErrors* graphEtaToPi0PCM900GeVStat                   = (TGraphAsymmErrors*)fFileNeutralMeson900GeVEta->Get("graphEtaToPi0PCMStat");
                TGraphAsymmErrors* graphEtaToPi0PCM900GeVSys                    = (TGraphAsymmErrors*)fFileNeutralMeson900GeVEta->Get("graphEtaToPi0PCMSys");
                SetGraphProperties(graphEtaToPi0PCM900GeVStat,   Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCM900GeVSys,    Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                list_900GeV->Add(graphEtaToPi0PCM900GeVStat);
                list_900GeV->Add(graphEtaToPi0PCM900GeVSys);

                if (changeQuantity == 0){
                    TGraphAsymmErrors* graphEtaXSecPCM900GeVStat                = (TGraphAsymmErrors*)graphEtaInvYieldPCM900GeVStat->Clone(  "graphInvCrosSecEtaPCMStat");
                    TGraphAsymmErrors* graphEtaXSecPCM900GeVSys                 = (TGraphAsymmErrors*)graphEtaInvYieldPCM900GeVSys->Clone(   "graphInvCrosSecEtaPCMSys");

                    graphEtaXSecPCM900GeVStat                                   = ScaleGraph(graphEtaXSecPCM900GeVStat, xSection900GeVV0OR);
                    graphEtaXSecPCM900GeVSys                                    = ScaleGraph(graphEtaXSecPCM900GeVSys,  xSection900GeVV0OR);

                    SetGraphProperties(graphEtaXSecPCM900GeVStat,    Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCM900GeVSys,     Form("%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    list_900GeV->Add(graphEtaXSecPCM900GeVStat);
                    list_900GeV->Add(graphEtaXSecPCM900GeVSys);

                } else {
                    // convert to INEL events
                    graphEtaInvYieldPCM900GeVStat                               = ScaleGraph(graphEtaInvYieldPCM900GeVStat,   xSection900GeVV0OR / xSection900GeVINEL);
                    graphEtaInvYieldPCM900GeVSys                                = ScaleGraph(graphEtaInvYieldPCM900GeVSys,    xSection900GeVV0OR / xSection900GeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphEtaInvYieldPCM900GeVStat                               = ConvertYieldGraph(graphEtaInvYieldPCM900GeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaInvYieldPCM900GeVSys                                = ConvertYieldGraph(graphEtaInvYieldPCM900GeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphEtaInvYieldPCM900GeVStat,     Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaInvYieldPCM900GeVSys,      Form("%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_900GeV->Add(graphEtaInvYieldPCM900GeVStat);
                    list_900GeV->Add(graphEtaInvYieldPCM900GeVSys);
                }
            }
        }
        if (enable_CPion && enable_CKaon && enable_Proton && enable_Lambda && enable_K0s){
            //================================================================================================================
            // reading and writing pi+-, K+- and p/bar{p} to 7TeV list
            // input is given as fully invariant yield
            //================================================================================================================
            TString localOutputQuantity                                         = globalOutputQuantity;

            TFile *fileIdentifiedCharged900GeV                                  = new TFile("pp/900GeV_pub_CPionCKaonProtonK0sLambda.root");

            TGraphAsymmErrors* graphCPion900GeVStat                             = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("chargedPionStat");
            TGraphAsymmErrors* graphCPion900GeVSys                              = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("chargedPionSys");
            TGraphAsymmErrors* graphCKaon900GeVStat                             = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("chargedKaonStat");
            TGraphAsymmErrors* graphCKaon900GeVSys                              = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("chargedKaonSys");
            TGraphAsymmErrors* graphProton900GeVStat                            = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("protonStat");
            TGraphAsymmErrors* graphProton900GeVSys                             = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("protonSys");
            TGraphAsymmErrors* graphLambda900GeVStat                            = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("lambdaStat");
            TGraphAsymmErrors* graphLambda900GeVSys                             = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("lambdaSys");
            TGraphAsymmErrors* graphK0s900GeVStat                               = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("k0sStat");
            TGraphAsymmErrors* graphK0s900GeVSys                                = (TGraphAsymmErrors*)fileIdentifiedCharged900GeV->Get("k0sSys");

            graphCPion900GeVStat                                                = ScaleGraph(graphCPion900GeVStat,  0.5);
            graphCPion900GeVSys                                                 = ScaleGraph(graphCPion900GeVSys,   0.5);
            graphCKaon900GeVStat                                                = ScaleGraph(graphCKaon900GeVStat,  0.5);
            graphCKaon900GeVSys                                                 = ScaleGraph(graphCKaon900GeVSys,   0.5);
            graphProton900GeVStat                                               = ScaleGraph(graphProton900GeVStat, 0.5);
            graphProton900GeVSys                                                = ScaleGraph(graphProton900GeVSys,  0.5);
            graphLambda900GeVStat                                               = ScaleGraph(graphLambda900GeVStat, 0.5);
            graphLambda900GeVSys                                                = ScaleGraph(graphLambda900GeVSys,  0.5);

            SetGraphProperties(graphCPion900GeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphCPion900GeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphCKaon900GeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphCKaon900GeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphProton900GeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphProton900GeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphLambda900GeVStat, Form("%s%sStat", fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphLambda900GeVSys,  Form("%s%sSys",  fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphK0s900GeVStat,    Form("%s%sStat", fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetGraphProperties(graphK0s900GeVSys,     Form("%s%sSys",  fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

            list_900GeV->Add(graphCPion900GeVStat);
            list_900GeV->Add(graphCPion900GeVSys);
            list_900GeV->Add(graphCKaon900GeVStat);
            list_900GeV->Add(graphCKaon900GeVSys);
            list_900GeV->Add(graphProton900GeVStat);
            list_900GeV->Add(graphProton900GeVSys);
            list_900GeV->Add(graphLambda900GeVStat);
            list_900GeV->Add(graphLambda900GeVSys);
            list_900GeV->Add(graphK0s900GeVStat);
            list_900GeV->Add(graphK0s900GeVSys);
        }
    }

    //================================================================================================================
    // creating histos and graphs for 2.76TeV
    //================================================================================================================
    if (Include_2760GeV) {

        //================================================================================================================
        // reading and writing pi0/eta to 2.76TeV list
        // input from pi0/eta is given as invariant cross section shifted in X direction
        //================================================================================================================
        if (enable_NPion || enable_Eta) {

            TString localOutputQuantity                                 = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";
            if (changeQuantity != 0){
                localOutputQuantity                                     = globalOutputQuantity;
            }

            TFile* fFileNeutralMeson2760GeV                             = new TFile("pp/CombinedResultsPaperPP2760GeV_2017_07_10_FrediV2Clusterizer.root");

            //================================================================================================================
            // pi0
            //================================================================================================================
            if (enable_NPion){

                TDirectoryFile* fFilePi02760GeV                             = (TDirectoryFile*)fFileNeutralMeson2760GeV->Get("Pi02.76TeV");

                TGraphAsymmErrors* graphNPionXSecComb2760GeVStat            = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr");
                TGraphAsymmErrors* graphNPionXSecComb2760GeVSys             = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
                TGraphAsymmErrors* graphNPionXSecPCM2760GeVStat             = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCM2760GeVSys              = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMCAL2760GeVStat        = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMCAL2760GeVSys         = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr");
                TGraphAsymmErrors* graphNPionXSecEMCAL2760GeVStat           = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
                TGraphAsymmErrors* graphNPionXSecEMCAL2760GeVSys            = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
                TGraphAsymmErrors* graphNPionXSecEMCALMerged2760GeVStat     = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0EMCALMerged2760GeVStatErr");
                TGraphAsymmErrors* graphNPionXSecEMCALMerged2760GeVSys      = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0EMCALMerged2760GeVSysErr");
                TGraphAsymmErrors* graphNPionXSecPHOS2760GeVStat            = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPHOS2760GeVSys             = (TGraphAsymmErrors*)fFilePi02760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");

                if (changeQuantity == 0){
                    // simply add Xsection to the output file
                    list_2760GeV->Add(graphNPionXSecComb2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecComb2760GeVSys);
                    list_2760GeV->Add(graphNPionXSecPCM2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecPCM2760GeVSys);
                    list_2760GeV->Add(graphNPionXSecPCMEMCAL2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecPCMEMCAL2760GeVSys);
                    list_2760GeV->Add(graphNPionXSecEMCAL2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecEMCAL2760GeVSys);
                    list_2760GeV->Add(graphNPionXSecEMCALMerged2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecEMCALMerged2760GeVSys);
                    list_2760GeV->Add(graphNPionXSecPHOS2760GeVStat);
                    list_2760GeV->Add(graphNPionXSecPHOS2760GeVSys);
                } else {
                    // convert XSection to invariant yield
                    TGraphAsymmErrors* graphNPionComb2760GeVStat            = ScaleGraph(graphNPionXSecComb2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionComb2760GeVSys             = ScaleGraph(graphNPionXSecComb2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPCM2760GeVStat             = ScaleGraph(graphNPionXSecPCM2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPCM2760GeVSys              = ScaleGraph(graphNPionXSecPCM2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPCMEMCAL2760GeVStat        = ScaleGraph(graphNPionXSecPCMEMCAL2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPCMEMCAL2760GeVSys         = ScaleGraph(graphNPionXSecPCMEMCAL2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionEMCAL2760GeVStat           = ScaleGraph(graphNPionXSecEMCAL2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionEMCAL2760GeVSys            = ScaleGraph(graphNPionXSecEMCAL2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionEMCALMerged2760GeVStat     = ScaleGraph(graphNPionXSecEMCALMerged2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionEMCALMerged2760GeVSys      = ScaleGraph(graphNPionXSecEMCALMerged2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPHOS2760GeVStat            = ScaleGraph(graphNPionXSecPHOS2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphNPionPHOS2760GeVSys             = ScaleGraph(graphNPionXSecPHOS2760GeVSys,1./xSection2760GeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphNPionComb2760GeVStat                               = ConvertYieldGraph(graphNPionComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionComb2760GeVSys                                = ConvertYieldGraph(graphNPionComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCM2760GeVStat                                = ConvertYieldGraph(graphNPionPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCM2760GeVSys                                 = ConvertYieldGraph(graphNPionPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCMEMCAL2760GeVStat                           = ConvertYieldGraph(graphNPionPCMEMCAL2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCMEMCAL2760GeVSys                            = ConvertYieldGraph(graphNPionPCMEMCAL2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCAL2760GeVStat                              = ConvertYieldGraph(graphNPionEMCAL2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCAL2760GeVSys                               = ConvertYieldGraph(graphNPionEMCAL2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCALMerged2760GeVStat                        = ConvertYieldGraph(graphNPionEMCALMerged2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCALMerged2760GeVSys                         = ConvertYieldGraph(graphNPionEMCALMerged2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPHOS2760GeVStat                               = ConvertYieldGraph(graphNPionPHOS2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPHOS2760GeVSys                                = ConvertYieldGraph(graphNPionPHOS2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphNPionComb2760GeVStat,           Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionComb2760GeVSys,            Form("%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCM2760GeVStat,            Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCM2760GeVSys,             Form("%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCMEMCAL2760GeVStat,       Form("%s%sStat", fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCMEMCAL2760GeVSys,        Form("%s%sSys",  fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCAL2760GeVStat,          Form("%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCAL2760GeVSys,           Form("%s%sSys",  fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCALMerged2760GeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[7].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCALMerged2760GeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[7].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPHOS2760GeVStat,           Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPHOS2760GeVSys,            Form("%s%sSys",  fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_2760GeV->Add(graphNPionComb2760GeVStat);
                    list_2760GeV->Add(graphNPionComb2760GeVSys);
                    list_2760GeV->Add(graphNPionPCM2760GeVStat);
                    list_2760GeV->Add(graphNPionPCM2760GeVSys);
                    list_2760GeV->Add(graphNPionPCMEMCAL2760GeVStat);
                    list_2760GeV->Add(graphNPionPCMEMCAL2760GeVSys);
                    list_2760GeV->Add(graphNPionEMCAL2760GeVStat);
                    list_2760GeV->Add(graphNPionEMCAL2760GeVSys);
                    list_2760GeV->Add(graphNPionEMCALMerged2760GeVStat);
                    list_2760GeV->Add(graphNPionEMCALMerged2760GeVSys);
                    list_2760GeV->Add(graphNPionPHOS2760GeVStat);
                    list_2760GeV->Add(graphNPionPHOS2760GeVSys);
                }
            }

            //================================================================================================================
            // eta
            //================================================================================================================
            if (enable_Eta){

                TDirectoryFile* fFileEta2760GeV                             = (TDirectoryFile*)fFileNeutralMeson2760GeV->Get("Eta2.76TeV");

                // eta to pi0
                TH1D*               histoEtaToPi0PCM2760GeVStat             = (TH1D*)fFileEta2760GeV->Get("histoRatioEtaToPi0PCM2760GeVStatErr");
                TGraphAsymmErrors*  graphEtaToPi0PCM2760GeVSys              = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphRatioEtaToPi0PCM2760GeVSysErr");
                TH1D*               histoEtaToPi0PCMEMCAL2760GeVStat        = (TH1D*)fFileEta2760GeV->Get("histoRatioEtaToPi0PCMEMCAL2760GeVStatErr");
                TGraphAsymmErrors*  graphEtaToPi0PCMEMCAL2760GeVSys         = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphRatioEtaToPi0PCMEMCAL2760GeVSysErr");
                TH1D*               histoEtaToPi0EMCAL2760GeVStat           = (TH1D*)fFileEta2760GeV->Get("histoRatioEtaToPi0EMCAL2760GeVStatErr");
                TGraphAsymmErrors*  graphEtaToPi0EMCAL2760GeVSys            = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphRatioEtaToPi0EMCAL2760GeVSysErr");
                TGraphAsymmErrors*  graphEtaToPi0Comb2760GeVStat            = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphRatioEtaToPi0Comb2760GeVStatErr");
                TGraphAsymmErrors*  graphEtaToPi0Comb2760GeVSys             = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphRatioEtaToPi0Comb2760GeVSysErr");

                SetGraphProperties(graphEtaToPi0Comb2760GeVStat,            Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0Comb2760GeVSys,             Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetHistoProperties(histoEtaToPi0PCM2760GeVStat,             Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCM2760GeVSys,              Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetHistoProperties(histoEtaToPi0PCMEMCAL2760GeVStat,        Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMEMCAL2760GeVSys,         Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetHistoProperties(histoEtaToPi0EMCAL2760GeVStat,           Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0EMCAL2760GeVSys,            Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");

                list_2760GeV->Add(graphEtaToPi0Comb2760GeVStat);
                list_2760GeV->Add(graphEtaToPi0Comb2760GeVSys);
                list_2760GeV->Add(histoEtaToPi0PCM2760GeVStat);
                list_2760GeV->Add(graphEtaToPi0PCM2760GeVSys);
                list_2760GeV->Add(histoEtaToPi0PCMEMCAL2760GeVStat);
                list_2760GeV->Add(graphEtaToPi0PCMEMCAL2760GeVSys);
                list_2760GeV->Add(histoEtaToPi0EMCAL2760GeVStat);
                list_2760GeV->Add(graphEtaToPi0EMCAL2760GeVSys);

                // eta
                TGraphAsymmErrors* graphEtaXSecComb2760GeVStat              = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaComb2760GeVAStatErr");
                TGraphAsymmErrors* graphEtaXSecComb2760GeVSys               = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaComb2760GeVASysErr");
                TGraphAsymmErrors* graphEtaXSecPCM2760GeVStat               = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCM2760GeVSys                = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVSysErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMCAL2760GeVStat          = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaPCMEMCAL2760GeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMCAL2760GeVSys           = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaPCMEMCAL2760GeVSysErr");
                TGraphAsymmErrors* graphEtaXSecEMCAL2760GeVStat             = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
                TGraphAsymmErrors* graphEtaXSecEMCAL2760GeVSys              = (TGraphAsymmErrors*)fFileEta2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr");

                if (changeQuantity == 0){
                    // simply add Xsection to the output file
                    list_2760GeV->Add(graphEtaXSecComb2760GeVStat);
                    list_2760GeV->Add(graphEtaXSecComb2760GeVSys);
                    list_2760GeV->Add(graphEtaXSecPCM2760GeVStat);
                    list_2760GeV->Add(graphEtaXSecPCM2760GeVSys);
                    list_2760GeV->Add(graphEtaXSecPCMEMCAL2760GeVStat);
                    list_2760GeV->Add(graphEtaXSecPCMEMCAL2760GeVSys);
                    list_2760GeV->Add(graphEtaXSecEMCAL2760GeVStat);
                    list_2760GeV->Add(graphEtaXSecEMCAL2760GeVSys);
                } else {
                    // convert XSection to invariant yield
                    TGraphAsymmErrors* graphEtaComb2760GeVStat              = ScaleGraph(graphEtaXSecComb2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaComb2760GeVSys               = ScaleGraph(graphEtaXSecComb2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaPCM2760GeVStat               = ScaleGraph(graphEtaXSecPCM2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaPCM2760GeVSys                = ScaleGraph(graphEtaXSecPCM2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaPCMEMCAL2760GeVStat          = ScaleGraph(graphEtaXSecPCMEMCAL2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaPCMEMCAL2760GeVSys           = ScaleGraph(graphEtaXSecPCMEMCAL2760GeVSys,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaEMCAL2760GeVStat             = ScaleGraph(graphEtaXSecEMCAL2760GeVStat,1./xSection2760GeVINEL);
                    TGraphAsymmErrors* graphEtaEMCAL2760GeVSys              = ScaleGraph(graphEtaXSecEMCAL2760GeVSys,1./xSection2760GeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphEtaComb2760GeVStat                                 = ConvertYieldGraph(graphEtaComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaComb2760GeVSys                                  = ConvertYieldGraph(graphEtaComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCM2760GeVStat                                  = ConvertYieldGraph(graphEtaPCM2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCM2760GeVSys                                   = ConvertYieldGraph(graphEtaPCM2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCMEMCAL2760GeVStat                             = ConvertYieldGraph(graphEtaPCMEMCAL2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCMEMCAL2760GeVSys                              = ConvertYieldGraph(graphEtaPCMEMCAL2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaEMCAL2760GeVStat                                = ConvertYieldGraph(graphEtaEMCAL2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaEMCAL2760GeVSys                                 = ConvertYieldGraph(graphEtaEMCAL2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphEtaComb2760GeVStat,             Form("%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaComb2760GeVSys,              Form("%s%sSys",  fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCM2760GeVStat,              Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCM2760GeVSys,               Form("%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCMEMCAL2760GeVStat,         Form("%s%sStat", fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCMEMCAL2760GeVSys,          Form("%s%sSys",  fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaEMCAL2760GeVStat,            Form("%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaEMCAL2760GeVSys,             Form("%s%sSys",  fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_2760GeV->Add(graphEtaComb2760GeVStat);
                    list_2760GeV->Add(graphEtaComb2760GeVSys);
                    list_2760GeV->Add(graphEtaPCM2760GeVStat);
                    list_2760GeV->Add(graphEtaPCM2760GeVSys);
                    list_2760GeV->Add(graphEtaPCMEMCAL2760GeVStat);
                    list_2760GeV->Add(graphEtaPCMEMCAL2760GeVSys);
                    list_2760GeV->Add(graphEtaEMCAL2760GeVStat);
                    list_2760GeV->Add(graphEtaEMCAL2760GeVSys);
                }
            }
        }

        if (enable_CPion && enable_CKaon && enable_Proton ){
            //================================================================================================================
            // reading and writing pi+-, K+- and p/bar{p} to 2.76TeV list
            // input is given as fully invariant yield
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile *fileIdentifiedCharged2760GeV                     = new TFile("pp/pp276.fullpT.INEL.20140504.root");

            TH1D* histoCPionComb2760GeVStat                         = (TH1D*)fileIdentifiedCharged2760GeV->Get("hstat_pp276_pion_sum");
            TH1D* histoCPionComb2760GeVSys                          = (TH1D*)fileIdentifiedCharged2760GeV->Get("hsys_pp276_pion_sum");
            TH1D* histoCKaonComb2760GeVStat                         = (TH1D*)fileIdentifiedCharged2760GeV->Get("hstat_pp276_kaon_sum");
            TH1D* histoCKaonComb2760GeVSys                          = (TH1D*)fileIdentifiedCharged2760GeV->Get("hsys_pp276_kaon_sum");
            TH1D* histoProtonComb2760GeVStat                        = (TH1D*)fileIdentifiedCharged2760GeV->Get("hstat_pp276_proton_sum");
            TH1D* histoProtonComb2760GeVSys                         = (TH1D*)fileIdentifiedCharged2760GeV->Get("hsys_pp276_proton_sum");

            // scale by 0.5 to get averaged single particle
            histoCPionComb2760GeVStat->Scale(0.5);
            histoCPionComb2760GeVSys->Scale(0.5);
            histoCKaonComb2760GeVStat->Scale(0.5);
            histoCKaonComb2760GeVSys->Scale(0.5);
            histoProtonComb2760GeVStat->Scale(0.5);
            histoProtonComb2760GeVSys->Scale(0.5);

            if (changeQuantity == 0 || changeQuantity == 2){
                SetHistoProperties(histoCPionComb2760GeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCPionComb2760GeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCKaonComb2760GeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCKaonComb2760GeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoProtonComb2760GeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoProtonComb2760GeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_2760GeV->Add(histoCPionComb2760GeVStat);
                list_2760GeV->Add(histoCPionComb2760GeVSys);
                list_2760GeV->Add(histoCKaonComb2760GeVStat);
                list_2760GeV->Add(histoCKaonComb2760GeVSys);
                list_2760GeV->Add(histoProtonComb2760GeVStat);
                list_2760GeV->Add(histoProtonComb2760GeVSys);
            } else {
                histoCPionComb2760GeVStat                              = ConvertYieldHisto(histoCPionComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                histoCPionComb2760GeVSys                               = ConvertYieldHisto(histoCPionComb2760GeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                histoCKaonComb2760GeVStat                              = ConvertYieldHisto(histoCKaonComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                histoCKaonComb2760GeVSys                               = ConvertYieldHisto(histoCKaonComb2760GeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                histoProtonComb2760GeVStat                             = ConvertYieldHisto(histoProtonComb2760GeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                histoProtonComb2760GeVSys                              = ConvertYieldHisto(histoProtonComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                SetHistoProperties(histoCPionComb2760GeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCPionComb2760GeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCKaonComb2760GeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCKaonComb2760GeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoProtonComb2760GeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoProtonComb2760GeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_2760GeV->Add(histoCPionComb2760GeVStat);
                list_2760GeV->Add(histoCPionComb2760GeVSys);
                list_2760GeV->Add(histoCKaonComb2760GeVStat);
                list_2760GeV->Add(histoCKaonComb2760GeVSys);
                list_2760GeV->Add(histoProtonComb2760GeVStat);
                list_2760GeV->Add(histoProtonComb2760GeVSys);
            }
        }

        if (enable_K0s) {
            //================================================================================================================
            // reading and writing K0s to 2.76TeV list
            // input from K0s is given as 1/N_inel dN^2/dydpt
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileK0s2760GeV                                     = new TFile("pp/pp_spectra_18_02_16_extrapol_etaErr_2760GeV.root");
            TH1F* histK0s2760GeVStat                                   = (TH1F*)fFileK0s2760GeV->Get("K0s_corr_staterr_centpp_extrapolated");
            TH1F* histK0s2760GeVSys                                    = (TH1F*)fFileK0s2760GeV->Get("K0s_corr_systerr_centpp_extrapolated");

            SetHistoProperties(histK0s2760GeVStat,  Form("%s%sStat", fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histK0s2760GeVSys,   Form("%s%sSys",  fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_2760GeV->Add(histK0s2760GeVStat);
            list_2760GeV->Add(histK0s2760GeVSys);
        }

        if (enable_Lambda) {
            //================================================================================================================
            // reading and writing Lambda to 2.76TeV list
            // input from Lambda is given as 1/N_inel dN^2/dydpt
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileLambda2760GeV                                  = new TFile("pp/pp_spectra_18_02_16_extrapol_etaErr_2760GeV.root");
            TH1F* histLambda2760GeVStat                                = (TH1F*)fFileLambda2760GeV->Get("L_corr_staterr_centpp_extrapolated");
            histLambda2760GeVStat->SetName("histLambda2760GeVStat");
            histLambda2760GeVStat->Sumw2();
            TH1F* histLambda2760GeVSys                                 = (TH1F*)fFileLambda2760GeV->Get("L_corr_systerr_centpp_extrapolated");
            histLambda2760GeVSys->SetName("histLambda2760GeVSys");
            histLambda2760GeVSys->Sumw2();

            TH1F* histAntiLambda2760GeVStat                            = (TH1F*)fFileLambda2760GeV->Get("AL_corr_staterr_centpp_extrapolated");
            histAntiLambda2760GeVStat->SetName("histAntiLambda2760GeVStat");
            histAntiLambda2760GeVStat->Sumw2();
            TH1F* histAntiLambda2760GeVSys                             = (TH1F*)fFileLambda2760GeV->Get("AL_corr_systerr_centpp_extrapolated");
            histAntiLambda2760GeVSys->SetName("histAntiLambda2760GeVSys");
            histAntiLambda2760GeVSys->Sumw2();

            TH1F* histLambdaComb2760GeVStat                            = (TH1F*)histLambda2760GeVStat->Clone("histLambdaComb2760GeVStat");
            histLambdaComb2760GeVStat->Add(histAntiLambda2760GeVStat);
            histLambdaComb2760GeVStat->Scale(0.5);
            TH1F* histLambdaComb2760GeVSys                             = (TH1F*)histLambda2760GeVSys->Clone("histLambdaComb2760GeVSys");
            histLambdaComb2760GeVSys->Add(histAntiLambda2760GeVSys);
            histLambdaComb2760GeVSys->Scale(0.5);

            SetHistoProperties(histLambdaComb2760GeVStat,  Form("%s%sStat", fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histLambdaComb2760GeVSys,   Form("%s%sSys",  fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_2760GeV->Add(histLambdaComb2760GeVStat);
            list_2760GeV->Add(histLambdaComb2760GeVSys);
        }

        if (enable_Rho0){
            //================================================================================================================
            // reading and writing rho0 to 2.76TeV list
            // input from rho0 is given as invariant x-section
            //================================================================================================================

            TFile* fileRho02760GeV                                  = new TFile("pp/pp2760GeV_rho_spectra_prelimViktor20170617.root");
            TGraphErrors* graphNRhoXSec2760GeVStat                  = (TGraphErrors*)fileRho02760GeV->Get("pp_stat");
            TGraphAsymmErrors* graphNRhoXSec2760GeVSys              = (TGraphAsymmErrors*)fileRho02760GeV->Get("pp_sys");

            for (Int_t i = 0; i < graphNRhoXSec2760GeVStat->GetN(); i++){
                graphNRhoXSec2760GeVStat->SetPointError(i, graphNRhoXSec2760GeVSys->GetEXlow()[i], graphNRhoXSec2760GeVStat->GetEY()[i] );
            }

            if (changeQuantity == 0){
                // simply add Xsection to the output file
                list_2760GeV->Add(graphNRhoXSec2760GeVStat);
                list_2760GeV->Add(graphNRhoXSec2760GeVSys);
            } else {
                // convert XSection to invariant yield
                TGraphErrors* graphNRho2760GeVStat                  = ScaleGraph(graphNRhoXSec2760GeVStat,1./xSection2760GeVINEL);
                TGraphAsymmErrors* graphNRho2760GeVSys              = ScaleGraph(graphNRhoXSec2760GeVSys,1./xSection2760GeVINEL);

                // convert inv yield to 1/N d^2N/dydpT
                graphNRho2760GeVStat                                = ConvertYieldGraph(graphNRho2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphNRho2760GeVSys                                 = ConvertYieldGraph(graphNRho2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                SetGraphProperties(graphNRho2760GeVStat,            Form("%s%sStat", fParticle[11].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphNRho2760GeVSys,             Form("%s%sSys",  fParticle[11].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_2760GeV->Add(graphNRho2760GeVStat);
                list_2760GeV->Add(graphNRho2760GeVSys);
            }
        }

        if (enable_Phi) {
            //================================================================================================================
            // reading and writing Phi to 2.76TeV list
            // input from K0s is given as invariant yield (arxiv1702.00555)
            //================================================================================================================
            TString phiHEPDataFile                                  = "pp/HEPdata/pp2760GeV_Phi_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphPhiYieldComb2760GeVStat         = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphPhiYieldComb2760GeVSys          = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
            graphPhiYieldComb2760GeVStat                            = ConvertYieldGraph(graphPhiYieldComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            graphPhiYieldComb2760GeVSys                             = ConvertYieldGraph(graphPhiYieldComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphPhiXSecComb2760GeVStat      = ScaleGraph(graphPhiYieldComb2760GeVStat,xSection2760GeVINEL);
                TGraphAsymmErrors* graphPhiXSecComb2760GeVSys       = ScaleGraph(graphPhiYieldComb2760GeVSys,xSection2760GeVINEL);

                graphPhiXSecComb2760GeVStat                         = ConvertYieldGraph(graphPhiXSecComb2760GeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphPhiXSecComb2760GeVSys                          = ConvertYieldGraph(graphPhiXSecComb2760GeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphPhiXSecComb2760GeVStat,  Form("Xsec_%s%sStat", fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphPhiXSecComb2760GeVSys,   Form("Xsec_%s%sSys",  fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_2760GeV->Add(graphPhiXSecComb2760GeVStat);
                list_2760GeV->Add(graphPhiXSecComb2760GeVSys);
            } else {
                SetGraphProperties(graphPhiYieldComb2760GeVStat,    Form("%s%sStat", fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphPhiYieldComb2760GeVSys,     Form("%s%sSys",  fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_2760GeV->Add(graphPhiYieldComb2760GeVStat);
                list_2760GeV->Add(graphPhiYieldComb2760GeVSys);
            }

            TString phiToPiHEPDataFile                              = "pp/HEPdata/pp2760GeV_PhiToPi_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphPhiToPiComb2760GeVStat          = ParseHEPData(phiToPiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphPhiToPiComb2760GeVSys           = ParseHEPData(phiToPiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphPhiToPiComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[9].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/#pi", "");
            SetGraphProperties(graphPhiToPiComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[9].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/#pi", "");
            list_2760GeV->Add(graphPhiToPiComb2760GeVStat);
            list_2760GeV->Add(graphPhiToPiComb2760GeVSys);

            TString phiToKHEPDataFile                               = "pp/HEPdata/pp2760GeV_PhiToK_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphPhiToKComb2760GeVStat           = ParseHEPData(phiToKHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphPhiToKComb2760GeVSys            = ParseHEPData(phiToKHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphPhiToKComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[9].Data(), fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/K", "");
            SetGraphProperties(graphPhiToKComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[9].Data(), fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/K", "");
            list_2760GeV->Add(graphPhiToKComb2760GeVStat);
            list_2760GeV->Add(graphPhiToKComb2760GeVSys);

            TString pToPhiHEPDataFile                               = "pp/HEPdata/pp2760GeV_PToPhi_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphpToPhiComb2760GeVStat           = ParseHEPData(pToPhiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphpToPhiComb2760GeVSys            = ParseHEPData(pToPhiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphpToPhiComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[7].Data(), fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "p/#phi", "");
            SetGraphProperties(graphpToPhiComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[7].Data(), fParticle[9].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "p/#phi", "");
            list_2760GeV->Add(graphpToPhiComb2760GeVStat);
            list_2760GeV->Add(graphpToPhiComb2760GeVSys);
        }

        if (enable_K0star) {
            //================================================================================================================
            // reading and writing K0star to 2.76TeV list
            // input from K0s is given as invariant yield (arxiv1702.00555)
            //================================================================================================================
            TString k0starHEPDataFile                                  = "pp/HEPdata/pp2760GeV_K0Star_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphK0StarYieldComb2760GeVStat         = ParseHEPData(k0starHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphK0StarYieldComb2760GeVSys          = ParseHEPData(k0starHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);
            graphK0StarYieldComb2760GeVStat                            = ConvertYieldGraph(graphK0StarYieldComb2760GeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            graphK0StarYieldComb2760GeVSys                             = ConvertYieldGraph(graphK0StarYieldComb2760GeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphK0StarXSecComb2760GeVStat      = ScaleGraph(graphK0StarYieldComb2760GeVStat,xSection2760GeVINEL);
                TGraphAsymmErrors* graphK0StarXSecComb2760GeVSys       = ScaleGraph(graphK0StarYieldComb2760GeVSys,xSection2760GeVINEL);

                graphK0StarXSecComb2760GeVStat                         = ConvertYieldGraph(graphK0StarXSecComb2760GeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphK0StarXSecComb2760GeVSys                          = ConvertYieldGraph(graphK0StarXSecComb2760GeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphK0StarXSecComb2760GeVStat,  Form("Xsec_%s%sStat", fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphK0StarXSecComb2760GeVSys,   Form("Xsec_%s%sSys",  fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_2760GeV->Add(graphK0StarXSecComb2760GeVStat);
                list_2760GeV->Add(graphK0StarXSecComb2760GeVSys);
            } else {
                SetGraphProperties(graphK0StarYieldComb2760GeVStat,    Form("%s%sStat", fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphK0StarYieldComb2760GeVSys,     Form("%s%sSys",  fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_2760GeV->Add(graphK0StarYieldComb2760GeVStat);
                list_2760GeV->Add(graphK0StarYieldComb2760GeVSys);
            }
            TString k0StarToPiHEPDataFile                              = "pp/HEPdata/pp2760GeV_K0StarToPi_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphK0StarToPiComb2760GeVStat          = ParseHEPData(k0StarToPiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphK0StarToPiComb2760GeVSys           = ParseHEPData(k0StarToPiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphK0StarToPiComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[10].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K^{0*}/#pi", "");
            SetGraphProperties(graphK0StarToPiComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[10].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K^{0*}/#pi", "");
            list_2760GeV->Add(graphK0StarToPiComb2760GeVStat);
            list_2760GeV->Add(graphK0StarToPiComb2760GeVSys);

            TString k0StarToKHEPDataFile                               = "pp/HEPdata/pp2760GeV_K0StarToK_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphK0StarToKComb2760GeVStat           = ParseHEPData(k0StarToKHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphK0StarToKComb2760GeVSys            = ParseHEPData(k0StarToKHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphK0StarToKComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[9].Data(), fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K^{0*}/K", "");
            SetGraphProperties(graphK0StarToKComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[9].Data(), fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K^{0*}/K", "");
            list_2760GeV->Add(graphK0StarToKComb2760GeVStat);
            list_2760GeV->Add(graphK0StarToKComb2760GeVSys);

            TString pToK0StarHEPDataFile                               = "pp/HEPdata/pp2760GeV_PToK0Star_arxiv1702.00555.csv";

            TGraphAsymmErrors* graphpToK0StarComb2760GeVStat           = ParseHEPData(pToK0StarHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE, kFALSE);
            TGraphAsymmErrors* graphpToK0StarComb2760GeVSys            = ParseHEPData(pToK0StarHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            SetGraphProperties(graphpToK0StarComb2760GeVStat,   Form("%sTo%s%sStat", fParticle[7].Data(), fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "p/K^{0*}", "");
            SetGraphProperties(graphpToK0StarComb2760GeVSys,    Form("%sTo%s%sSys",  fParticle[7].Data(), fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "p/K^{0*}", "");
            list_2760GeV->Add(graphpToK0StarComb2760GeVStat);
            list_2760GeV->Add(graphpToK0StarComb2760GeVSys);
        }

        if (enable_CXi && enable_COmega) {
            //================================================================================================================
            // reading and writing Xi+ + Xi- to 2.76TeV list
            // input from Xis is given as yield https://aliceinfo.cern.ch/ArtSubmission/node/2211
            //================================================================================================================
            TString localOutputQuantity                                = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileCXi2760GeV                                     = new TFile("pp/XiAndOmega_fittedgfcorrectedspectrum_woSDD_weightedonZdistrib.root");
            TH1F* histCXi2760GeVStat                                   = (TH1F*)fFileCXi2760GeV->Get("hptcorr_norm_4");
            histCXi2760GeVStat->Scale(0.5);
            TH1F* histCXi2760GeVSys                                    = (TH1F*)histCXi2760GeVStat->Clone("hptcorr_norm_4_sys");
            for (Int_t i = 1; i< histCXi2760GeVSys->GetNbinsX()+1; i++){
                histCXi2760GeVSys->SetBinError(i, histCXi2760GeVSys->GetBinContent(i)*0.05);
            }

            SetHistoProperties(histCXi2760GeVStat,  Form("%s%sStat", fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histCXi2760GeVSys,   Form("%s%sSys",  fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_2760GeV->Add(histCXi2760GeVStat);
            list_2760GeV->Add(histCXi2760GeVSys);

            //================================================================================================================
            // reading and writing Omega+ + Omega- to 2.76TeV list
            // input from Xis is given as yield https://aliceinfo.cern.ch/ArtSubmission/node/2211
            //================================================================================================================
            TH1F* histCOmega2760GeVStat                                = (TH1F*)fFileCXi2760GeV->Get("hptcorr_norm_5");
            histCOmega2760GeVStat->Scale(0.5);
            TH1F* histCOmega2760GeVSys                                 = (TH1F*)histCOmega2760GeVStat->Clone("hptcorr_norm_5_sys");
            for (Int_t i = 1; i< histCOmega2760GeVSys->GetNbinsX()+1; i++){
                histCOmega2760GeVSys->SetBinError(i, histCOmega2760GeVSys->GetBinContent(i)*0.09);
            }

            SetHistoProperties(histCOmega2760GeVStat,  Form("%s%sStat", fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histCOmega2760GeVSys,   Form("%s%sSys",  fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_2760GeV->Add(histCOmega2760GeVStat);
            list_2760GeV->Add(histCOmega2760GeVSys);
        }

        if (enablePtY) {
            //================================================================================================================
            // read in pythia pt y distributions for particles (i.e. PWGGA/GammaConv/AliAnalysisTaskGammaPureMC output)
            //================================================================================================================
            TFile*          fFilePtY2760GeV                             = new TFile("pp/pp2760GeV_PtYDistributions_Pythia8.root");
            TDirectoryFile* fDirPtY2760GeV                              = (TDirectoryFile*)fFilePtY2760GeV->Get("GammaPureMC");
            TList*          fListPtY2760GeV                             = (TList*)fDirPtY2760GeV->Get("GammaPureMC");

            Int_t       nBinsRebinX                                     = 5;
            Int_t       nBinsRebinY                                     = 4;
            Double_t    relStatErrThresh                                = 0.01;

            TH2F* histPtYPi0                                            = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Pi0");
            histPtYPi0->Sumw2();
            histPtYPi0->SetName("111_pt_y");
            NormalizePtYHistogram(histPtYPi0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYEta                                            = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Eta");
            histPtYEta->Sumw2();
            histPtYEta->SetName("221_pt_y");
            NormalizePtYHistogram(histPtYEta, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYEtaPrim                                        = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_EtaPrim");
            histPtYEtaPrim->Sumw2();
            histPtYEtaPrim->SetName("331_pt_y");
            NormalizePtYHistogram(histPtYEtaPrim, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYOmega                                          = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Omega");
            histPtYOmega->Sumw2();
            histPtYOmega->SetName("223_pt_y");
            NormalizePtYHistogram(histPtYOmega, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYRho0                                           = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Rho0");
            histPtYRho0->Sumw2();
            histPtYRho0->SetName("113_pt_y");
            NormalizePtYHistogram(histPtYRho0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYRhoPl                                          = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_RhoPl");
            histPtYRhoPl->Sumw2();
            histPtYRhoPl->SetName("213_pt_y");
            NormalizePtYHistogram(histPtYRhoPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYRhoMi                                          = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_RhoMi");
            histPtYRhoMi->Sumw2();
            histPtYRhoMi->SetName("-213_pt_y");
            NormalizePtYHistogram(histPtYRhoMi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYPhi                                            = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Phi");
            histPtYPhi->Sumw2();
            histPtYPhi->SetName("333_pt_y");
            NormalizePtYHistogram(histPtYPhi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYJPsi                                           = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_JPsi");
            histPtYJPsi->Sumw2();
            histPtYJPsi->SetName("443_pt_y");
            NormalizePtYHistogram(histPtYJPsi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYSigma0                                         = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Sigma0");
            histPtYSigma0->Sumw2();
            histPtYSigma0->SetName("3212_pt_y");
            NormalizePtYHistogram(histPtYSigma0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYK0s                                            = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_K0s");
            histPtYK0s->Sumw2();
            histPtYK0s->SetName("310_pt_y");
            NormalizePtYHistogram(histPtYK0s, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYK0l                                            = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_K0l");
            histPtYK0l->Sumw2();
            histPtYK0l->SetName("130_pt_y");
            NormalizePtYHistogram(histPtYK0l, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYK0star                                         = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_K0star");
            histPtYK0star->Sumw2();
            histPtYK0star->SetName("313_pt_y");
            NormalizePtYHistogram(histPtYK0star, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYDeltaPlPl                                      = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_DeltaPlPl");
            histPtYDeltaPlPl->Sumw2();
            histPtYDeltaPlPl->SetName("2224_pt_y");
            NormalizePtYHistogram(histPtYDeltaPlPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYDeltaPl                                        = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_DeltaPl");
            histPtYDeltaPl->Sumw2();
            histPtYDeltaPl->SetName("2214_pt_y");
            NormalizePtYHistogram(histPtYDeltaPl, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYDeltaMi                                        = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_DeltaMi");
            histPtYDeltaMi->Sumw2();
            histPtYDeltaMi->SetName("1114_pt_y");
            NormalizePtYHistogram(histPtYDeltaMi, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYDelta0                                         = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Delta0");
            histPtYDelta0->Sumw2();
            histPtYDelta0->SetName("2114_pt_y");
            NormalizePtYHistogram(histPtYDelta0, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            TH2F* histPtYLambda                                         = (TH2F*)fListPtY2760GeV->FindObject("Pt_Y_Lambda");
            histPtYLambda->Sumw2();
            histPtYLambda->SetName("3122_pt_y");
            NormalizePtYHistogram(histPtYLambda, nBinsRebinX, nBinsRebinY, relStatErrThresh);

            list_2760GeV->Add(histPtYPi0);
            list_2760GeV->Add(histPtYEta);
            list_2760GeV->Add(histPtYEtaPrim);
            list_2760GeV->Add(histPtYOmega);
            list_2760GeV->Add(histPtYRho0);
            list_2760GeV->Add(histPtYRhoPl);
            list_2760GeV->Add(histPtYRhoMi);
            list_2760GeV->Add(histPtYPhi);
            list_2760GeV->Add(histPtYJPsi);
            list_2760GeV->Add(histPtYSigma0);
            list_2760GeV->Add(histPtYK0s);
            list_2760GeV->Add(histPtYK0l);
            list_2760GeV->Add(histPtYK0star);
            list_2760GeV->Add(histPtYDeltaPlPl);
            list_2760GeV->Add(histPtYDeltaPl);
            list_2760GeV->Add(histPtYDeltaMi);
            list_2760GeV->Add(histPtYDelta0);
            list_2760GeV->Add(histPtYLambda);
        }
    }

    //================================================================================================================
    // creating histos and graphs for 7TeV
    //================================================================================================================
    if (Include_7TeV){

        //================================================================================================================
        // reading and writing pi0 and eta to 7TeV list
        // input from pi0 and eta is given as invariant cross section shifted in X direction
        //================================================================================================================
        if (enable_NPion || enable_Eta){
            TString localOutputQuantity                                 = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

            TFile*          fFileNeutralMeson7TeV                       = new TFile("pp/CombinedResultsPaperPP7TeV_2017_07_14.root");
            TDirectoryFile* fFileNeutralMeson7TeVPi0                    = (TDirectoryFile*)fFileNeutralMeson7TeV->Get("Pi07TeV");
            TDirectoryFile* fFileNeutralMeson7TeVEta                    = (TDirectoryFile*)fFileNeutralMeson7TeV->Get("Eta7TeV");

            TFile*          fFileNeutralMeson7TeVPass2                  = new TFile("pp/data_GammaConversionResultsFullCorrectionNoBinShifting_PCM_020712.root");
            TDirectoryFile* fFileNeutralMeson7TeVPass2Pi0               = (TDirectoryFile*)fFileNeutralMeson7TeVPass2->Get("Pi07TeV");
            TDirectoryFile* fFileNeutralMeson7TeVPass2Eta               = (TDirectoryFile*)fFileNeutralMeson7TeVPass2->Get("Eta7TeV");

            if (enable_NPion){
                TGraphAsymmErrors* graphNPionXSecComb7TeVStat           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0Comb7TeVAStatErr");
                TGraphAsymmErrors* graphNPionXSecComb7TeVSys            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0Comb7TeVASysErr");
                TGraphAsymmErrors* graphNPionXSecPCM7TeVStat            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCM7TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCM7TeVSys             = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCM7TeVSysErr");
                TGraphAsymmErrors* graphNPionXSecPHOS7TeVStat           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PHOS7TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPHOS7TeVSys            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PHOS7TeVSysErr");
                TGraphAsymmErrors* graphNPionXSecEMCal7TeVStat          = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0EMCAL7TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecEMCal7TeVSys           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0EMCAL7TeVSysErr");
                //TGraphAsymmErrors* graphNPionXSecPCMPHOS7TeVStat        = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCMPHOS7TeVStatErr");
                //TGraphAsymmErrors* graphNPionXSecPCMPHOS7TeVSys         = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCMPHOS7teVSysErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMC7TeVStat         = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCMEMCAL7TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMC7TeVSys          = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPi0->Get("graphInvCrossSectionPi0PCMEMCAL7TeVSysErr");

                // pass2
                TH1D*               histoNPionXSecPCMPass27TeVStat      = (TH1D*)fFileNeutralMeson7TeVPass2Pi0->Get("InvCrossSectionPi0");
                TGraphAsymmErrors*  graphNPionXSecPCMPass27TeVStat      = HistToGraph(histoNPionXSecPCMPass27TeVStat);
                graphNPionXSecPCMPass27TeVStat->RemovePoint(0);
                TGraphAsymmErrors*  graphNPionXSecPCMPass27TeVSys       = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPass2Pi0->Get("InvCrossSectionPi0Sys");

                if (changeQuantity == 0) {
                    // simply add Xsection to the output file
                    SetGraphProperties(graphNPionXSecComb7TeVStat,      Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecComb7TeVSys,       Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCM7TeVStat,       Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCM7TeVSys,        Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPHOS7TeVStat,      Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPHOS7TeVSys,       Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecEMCal7TeVStat,     Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecEMCal7TeVSys,      Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    //SetGraphProperties(graphNPionXSecPCMPHOS7TeVStat,   Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[5].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    //SetGraphProperties(graphNPionXSecPCMPHOS7TeVSys,    Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[5].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMEMC7TeVStat,    Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMEMC7TeVSys,     Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMPass27TeVStat,  Form("Xsec_%s%sStat", fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMPass27TeVSys,   Form("Xsec_%s%sSys",  fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    list_7TeV->Add(graphNPionXSecComb7TeVStat);
                    list_7TeV->Add(graphNPionXSecComb7TeVSys);
                    list_7TeV->Add(graphNPionXSecPCM7TeVStat);
                    list_7TeV->Add(graphNPionXSecPCM7TeVSys);
                    list_7TeV->Add(graphNPionXSecPHOS7TeVStat);
                    list_7TeV->Add(graphNPionXSecPHOS7TeVSys);
                    list_7TeV->Add(graphNPionXSecEMCal7TeVStat);
                    list_7TeV->Add(graphNPionXSecEMCal7TeVSys);
                    //list_7TeV->Add(graphNPionXSecPCMPHOS7TeVStat);
                    //list_7TeV->Add(graphNPionXSecPCMPHOS7TeVSys);
                    list_7TeV->Add(graphNPionXSecPCMEMC7TeVStat);
                    list_7TeV->Add(graphNPionXSecPCMEMC7TeVSys);
                    list_7TeV->Add(graphNPionXSecPCMPass27TeVStat);
                    list_7TeV->Add(graphNPionXSecPCMPass27TeVSys);
                } else {
                    TGraphAsymmErrors* graphNPionYieldComb7TeVStat      = (TGraphAsymmErrors*)graphNPionXSecComb7TeVStat->Clone(   "graphYieldPi0CombStat");
                    TGraphAsymmErrors* graphNPionYieldComb7TeVSys       = (TGraphAsymmErrors*)graphNPionXSecComb7TeVSys->Clone(    "graphYieldPi0CombSys");
                    TGraphAsymmErrors* graphNPionYieldPCM7TeVStat       = (TGraphAsymmErrors*)graphNPionXSecPCM7TeVStat->Clone(    "graphYieldPi0PCMStat");
                    TGraphAsymmErrors* graphNPionYieldPCM7TeVSys        = (TGraphAsymmErrors*)graphNPionXSecPCM7TeVSys->Clone(     "graphYieldPi0PCMSys");
                    TGraphAsymmErrors* graphNPionYieldPHOS7TeVStat      = (TGraphAsymmErrors*)graphNPionXSecPHOS7TeVStat->Clone(   "graphYieldPi0PHOSStat");
                    TGraphAsymmErrors* graphNPionYieldPHOS7TeVSys       = (TGraphAsymmErrors*)graphNPionXSecPHOS7TeVSys->Clone(    "graphYieldPi0PHOSSys");
                    TGraphAsymmErrors* graphNPionYieldEMCal7TeVStat     = (TGraphAsymmErrors*)graphNPionXSecEMCal7TeVStat->Clone(  "graphYieldPi0EMCalStat");
                    TGraphAsymmErrors* graphNPionYieldEMCal7TeVSys      = (TGraphAsymmErrors*)graphNPionXSecEMCal7TeVSys->Clone(   "graphYieldPi0EMCalSys");
                    //TGraphAsymmErrors* graphNPionYieldPCMPHOS7TeVStat   = (TGraphAsymmErrors*)graphNPionXSecPCMPHOS7TeVStat->Clone("graphYieldPi0PCM-PHOSStat");
                    //TGraphAsymmErrors* graphNPionYieldPCMPHOS7TeVSys    = (TGraphAsymmErrors*)graphNPionXSecPCMPHOS7TeVSys->Clone( "graphYieldPi0PCM-PHOSSys");
                    TGraphAsymmErrors* graphNPionYieldPCMEMC7TeVStat    = (TGraphAsymmErrors*)graphNPionXSecPCMEMC7TeVStat->Clone( "graphYieldPi0PCMEMCStat");
                    TGraphAsymmErrors* graphNPionYieldPCMEMC7TeVSys     = (TGraphAsymmErrors*)graphNPionXSecPCMEMC7TeVSys->Clone(  "graphYieldPi0PCMEMCSys");
                    TGraphAsymmErrors* graphNPionYieldPCMPass27TeVStat  = (TGraphAsymmErrors*)graphNPionXSecPCMPass27TeVStat->Clone( "graphYieldPi0PCMPass2Stat");
                    TGraphAsymmErrors* graphNPionYieldPCMPass27TeVSys   = (TGraphAsymmErrors*)graphNPionXSecPCMPass27TeVSys->Clone(  "graphYieldPi0PCMPass2Sys");


                    // convert XSection to invariant yield
                    graphNPionYieldComb7TeVStat                         = ScaleGraph(graphNPionYieldComb7TeVStat,   1./xSection7TeVINEL);
                    graphNPionYieldComb7TeVSys                          = ScaleGraph(graphNPionYieldComb7TeVSys,    1./xSection7TeVINEL);
                    graphNPionYieldPCM7TeVStat                          = ScaleGraph(graphNPionYieldPCM7TeVStat,    1./xSection7TeVINEL);
                    graphNPionYieldPCM7TeVSys                           = ScaleGraph(graphNPionYieldPCM7TeVSys,     1./xSection7TeVINEL);
                    graphNPionYieldPHOS7TeVStat                         = ScaleGraph(graphNPionYieldPHOS7TeVStat,   1./xSection7TeVINEL);
                    graphNPionYieldPHOS7TeVSys                          = ScaleGraph(graphNPionYieldPHOS7TeVSys,    1./xSection7TeVINEL);
                    graphNPionYieldEMCal7TeVStat                        = ScaleGraph(graphNPionYieldEMCal7TeVStat,  1./xSection7TeVINEL);
                    graphNPionYieldEMCal7TeVSys                         = ScaleGraph(graphNPionYieldEMCal7TeVSys,   1./xSection7TeVINEL);
                    //graphNPionYieldPCMPHOS7TeVStat                      = ScaleGraph(graphNPionYieldPCMPHOS7TeVStat,1./xSection7TeVINEL);
                    //graphNPionYieldPCMPHOS7TeVSys                       = ScaleGraph(graphNPionYieldPCMPHOS7TeVSys, 1./xSection7TeVINEL);
                    graphNPionYieldPCMEMC7TeVStat                       = ScaleGraph(graphNPionYieldPCMEMC7TeVStat, 1./xSection7TeVINEL);
                    graphNPionYieldPCMEMC7TeVSys                        = ScaleGraph(graphNPionYieldPCMEMC7TeVSys,  1./xSection7TeVINEL);
                    graphNPionYieldPCMPass27TeVStat                     = ScaleGraph(graphNPionYieldPCMPass27TeVStat, 1./xSection7TeVINEL);
                    graphNPionYieldPCMPass27TeVSys                      = ScaleGraph(graphNPionYieldPCMPass27TeVSys,  1./xSection7TeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphNPionYieldComb7TeVStat                         = ConvertYieldGraph(graphNPionYieldComb7TeVStat,    kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldComb7TeVSys                          = ConvertYieldGraph(graphNPionYieldComb7TeVSys,     kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCM7TeVStat                          = ConvertYieldGraph(graphNPionYieldPCM7TeVStat,     kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCM7TeVSys                           = ConvertYieldGraph(graphNPionYieldPCM7TeVSys,      kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPHOS7TeVStat                         = ConvertYieldGraph(graphNPionYieldPHOS7TeVStat,    kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPHOS7TeVSys                          = ConvertYieldGraph(graphNPionYieldPHOS7TeVSys,     kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldEMCal7TeVStat                        = ConvertYieldGraph(graphNPionYieldEMCal7TeVStat,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldEMCal7TeVSys                         = ConvertYieldGraph(graphNPionYieldEMCal7TeVSys,    kFALSE, kFALSE, kTRUE, kTRUE);
                    //graphNPionYieldPCMPHOS7TeVStat                      = ConvertYieldGraph(graphNPionYieldPCMPHOS7TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                    //graphNPionYieldPCMPHOS7TeVSys                       = ConvertYieldGraph(graphNPionYieldPCMPHOS7TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCMEMC7TeVStat                       = ConvertYieldGraph(graphNPionYieldPCMEMC7TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCMEMC7TeVSys                        = ConvertYieldGraph(graphNPionYieldPCMEMC7TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCMPass27TeVStat                     = ConvertYieldGraph(graphNPionYieldPCMPass27TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionYieldPCMPass27TeVSys                      = ConvertYieldGraph(graphNPionYieldPCMPass27TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphNPionYieldComb7TeVStat,     Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldComb7TeVSys,      Form("%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCM7TeVStat,      Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCM7TeVSys,       Form("%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPHOS7TeVStat,     Form("%s%sStat", fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPHOS7TeVSys,      Form("%s%sSys",  fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldEMCal7TeVStat,    Form("%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldEMCal7TeVSys,     Form("%s%sSys",  fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    //SetGraphProperties(graphNPionYieldPCMPHOS7TeVStat,  Form("%s%sStat", fParticle[0].Data(), fMethod[5].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    //SetGraphProperties(graphNPionYieldPCMPHOS7TeVSys,   Form("%s%sSys",  fParticle[0].Data(), fMethod[5].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCMEMC7TeVStat,   Form("%s%sStat", fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCMEMC7TeVSys,    Form("%s%sSys",  fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCMPass27TeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionYieldPCMPass27TeVSys,  Form("%s%sSys",  fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_7TeV->Add(graphNPionYieldComb7TeVStat);
                    list_7TeV->Add(graphNPionYieldComb7TeVSys);
                    list_7TeV->Add(graphNPionYieldPCM7TeVStat);
                    list_7TeV->Add(graphNPionYieldPCM7TeVSys);
                    list_7TeV->Add(graphNPionYieldPHOS7TeVStat);
                    list_7TeV->Add(graphNPionYieldPHOS7TeVSys);
                    list_7TeV->Add(graphNPionYieldEMCal7TeVStat);
                    list_7TeV->Add(graphNPionYieldEMCal7TeVSys);
                    //list_7TeV->Add(graphNPionYieldPCMPHOS7TeVStat);
                    //list_7TeV->Add(graphNPionYieldPCMPHOS7TeVSys);
                    list_7TeV->Add(graphNPionYieldPCMEMC7TeVStat);
                    list_7TeV->Add(graphNPionYieldPCMEMC7TeVSys);
                    list_7TeV->Add(graphNPionYieldPCMPass27TeVStat);
                    list_7TeV->Add(graphNPionYieldPCMPass27TeVSys);
                }
            }
            if (enable_Eta){
                TGraphAsymmErrors* graphEtaXSecComb7TeVStat             = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaComb7TeVAStatErr");
                TGraphAsymmErrors* graphEtaXSecComb7TeVSys              = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaComb7TeVASysErr");
                TGraphAsymmErrors* graphEtaXSecPCM7TeVStat              = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaPCM7TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCM7TeVSys               = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaPCM7TeVSysErr");
                TGraphAsymmErrors* graphEtaXSecEMCal7TeVStat            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaEMCAL7TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecEMCal7TeVSys             = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaEMCAL7TeVSysErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMC7TeVStat           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaPCMEMCAL7TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMC7TeVSys            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphInvCrossSectionEtaPCMEMCAL7TeVSysErr");

                // pass2
                TH1D*               histoEtaXSecPCMPass27TeVStat      = (TH1D*)fFileNeutralMeson7TeVPass2Eta->Get("InvCrossSectionEta");
                TGraphAsymmErrors*  graphEtaXSecPCMPass27TeVStat      = HistToGraph(histoEtaXSecPCMPass27TeVStat);
                graphEtaXSecPCMPass27TeVStat->RemovePoint(0);
                TGraphAsymmErrors*  graphEtaXSecPCMPass27TeVSys       = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPass2Eta->Get("InvCrossSectionEtaSys");

                // remove first 5 points from eta stat EMCal (are zero)
                graphEtaXSecEMCal7TeVStat->RemovePoint(0);
                graphEtaXSecEMCal7TeVStat->RemovePoint(0);
                graphEtaXSecEMCal7TeVStat->RemovePoint(0);
                graphEtaXSecEMCal7TeVStat->RemovePoint(0);
                graphEtaXSecEMCal7TeVStat->RemovePoint(0);
//                cout << "graphInvCrossSectionEtaEMCAL7TeVStatErr points:" << endl;
//                for (Int_t i=0; i<graphEtaXSecEMCal7TeVStat->GetN(); i++) {
//                    if (i < graphEtaXSecEMCal7TeVSys->GetN())
//                        cout << i << ": stat = " << graphEtaXSecEMCal7TeVStat->GetX()[i] << "\t sys = " << graphEtaXSecEMCal7TeVSys->GetX()[i] << endl;
//                    else
//                        cout << i << ": stat = " << graphEtaXSecEMCal7TeVStat->GetX()[i] << endl;
//                }

                TGraphAsymmErrors* graphEtaToPi0Comb7TeVStat            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0Comb7TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0Comb7TeVSys             = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0Comb7TeVSysErr");
                TGraphAsymmErrors* graphEtaToPi0PCM7TeVStat             = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0PCM7TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0PCM7TeVSys              = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0PCM7TeVSysErr");
                TGraphAsymmErrors* graphEtaToPi0EMCal7TeVStat           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0EMCAL7TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0EMCal7TeVSys            = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0EMCAL7TeVSysErr");
                TGraphAsymmErrors* graphEtaToPi0PCMEMC7TeVStat          = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0PCMEMCAL7TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0PCMEMC7TeVSys           = (TGraphAsymmErrors*)fFileNeutralMeson7TeVEta->Get("graphRatioEtaToPi0PCMEMCAL7TeVSysErr");

                // pass2
                //TH1D*               histoEtaToPi0PCMPass27TeVStat       = (TH1D*)fFileNeutralMeson7TeVPass2Eta->Get("EtatoPi0RatioConversionBinShifted");
                TH1D*               histoEtaToPi0PCMPass27TeVStat       = (TH1D*)fFileNeutralMeson7TeVPass2Eta->Get("EtatoPi0RatioConversion");
                TGraphAsymmErrors*  graphEtaToPi0PCMPass27TeVStat       = HistToGraph(histoEtaToPi0PCMPass27TeVStat);
                graphEtaToPi0PCMPass27TeVStat->RemovePoint(0);
                TGraphAsymmErrors*  graphEtaToPi0PCMPass27TeVSys        = (TGraphAsymmErrors*)fFileNeutralMeson7TeVPass2Eta->Get("EtatoPi0RatioConversionBinShiftedSys");

                SetGraphProperties(graphEtaToPi0Comb7TeVStat,       Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0Comb7TeVSys,        Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCM7TeVStat,        Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCM7TeVSys,         Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0EMCal7TeVStat,      Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0EMCal7TeVSys,       Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMEMC7TeVStat,     Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMEMC7TeVSys,      Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMPass27TeVStat,   Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMPass27TeVSys,    Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");

                list_7TeV->Add(graphEtaToPi0Comb7TeVStat);
                list_7TeV->Add(graphEtaToPi0Comb7TeVSys);
                list_7TeV->Add(graphEtaToPi0PCM7TeVStat);
                list_7TeV->Add(graphEtaToPi0PCM7TeVSys);
                list_7TeV->Add(graphEtaToPi0EMCal7TeVStat);
                list_7TeV->Add(graphEtaToPi0EMCal7TeVSys);
                list_7TeV->Add(graphEtaToPi0PCMEMC7TeVStat);
                list_7TeV->Add(graphEtaToPi0PCMEMC7TeVSys);
                list_7TeV->Add(graphEtaToPi0PCMPass27TeVStat);
                list_7TeV->Add(graphEtaToPi0PCMPass27TeVSys);

                if (changeQuantity == 0){
                    // simply add Xsection to the output file
                    SetGraphProperties(graphEtaXSecComb7TeVStat,        Form("Xsec_%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecComb7TeVSys,         Form("Xsec_%s%sSys",  fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCM7TeVStat,         Form("Xsec_%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCM7TeVSys,          Form("Xsec_%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecEMCal7TeVStat,       Form("Xsec_%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecEMCal7TeVSys,        Form("Xsec_%s%sSys",  fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMEMC7TeVStat,      Form("Xsec_%s%sStat", fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMEMC7TeVSys,       Form("Xsec_%s%sSys",  fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMPass27TeVStat,    Form("Xsec_%s%sStat", fParticle[1].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMPass27TeVSys,     Form("Xsec_%s%sSys",  fParticle[1].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    list_7TeV->Add(graphEtaXSecComb7TeVStat);
                    list_7TeV->Add(graphEtaXSecComb7TeVStat);
                    list_7TeV->Add(graphEtaXSecPCM7TeVStat);
                    list_7TeV->Add(graphEtaXSecPCM7TeVSys);
                    list_7TeV->Add(graphEtaXSecEMCal7TeVStat);
                    list_7TeV->Add(graphEtaXSecEMCal7TeVSys);
                    list_7TeV->Add(graphEtaXSecPCMEMC7TeVStat);
                    list_7TeV->Add(graphEtaXSecPCMEMC7TeVSys);
                    list_7TeV->Add(graphEtaXSecPCMPass27TeVStat);
                    list_7TeV->Add(graphEtaXSecPCMPass27TeVSys);
                } else {
                    TGraphAsymmErrors* graphEtaYieldComb7TeVStat        = (TGraphAsymmErrors*)graphEtaXSecComb7TeVStat->Clone(      "graphYieldEtaCombStat");
                    TGraphAsymmErrors* graphEtaYieldComb7TeVSys         = (TGraphAsymmErrors*)graphEtaXSecComb7TeVSys->Clone(       "graphYieldEtaCombSys");
                    TGraphAsymmErrors* graphEtaYieldPCM7TeVStat         = (TGraphAsymmErrors*)graphEtaXSecPCM7TeVStat->Clone(       "graphYieldEtaPCMStat");
                    TGraphAsymmErrors* graphEtaYieldPCM7TeVSys          = (TGraphAsymmErrors*)graphEtaXSecPCM7TeVSys->Clone(        "graphYieldEtaPCMSys");
                    TGraphAsymmErrors* graphEtaYieldEMCal7TeVStat       = (TGraphAsymmErrors*)graphEtaXSecEMCal7TeVStat->Clone(     "graphYieldEtaEMCalStat");
                    TGraphAsymmErrors* graphEtaYieldEMCal7TeVSys        = (TGraphAsymmErrors*)graphEtaXSecEMCal7TeVSys->Clone(      "graphYieldEtaEMCalSys");
                    TGraphAsymmErrors* graphEtaYieldPCMEMC7TeVStat      = (TGraphAsymmErrors*)graphEtaXSecPCMEMC7TeVStat->Clone(    "graphYieldEtaPCMEMCStat");
                    TGraphAsymmErrors* graphEtaYieldPCMEMC7TeVSys       = (TGraphAsymmErrors*)graphEtaXSecPCMEMC7TeVSys->Clone(     "graphYieldEtaPCMEMCSys");
                    TGraphAsymmErrors* graphEtaYieldPCMPass27TeVStat    = (TGraphAsymmErrors*)graphEtaXSecPCMPass27TeVStat->Clone(  "graphYieldEtaPCMPass2Stat");
                    TGraphAsymmErrors* graphEtaYieldPCMPass27TeVSys     = (TGraphAsymmErrors*)graphEtaXSecPCMPass27TeVSys->Clone(   "graphYieldEtaPCMPass2Sys");

                    // convert XSection to invariant yield
                    graphEtaYieldComb7TeVStat                           = ScaleGraph(graphEtaYieldComb7TeVStat,     1./xSection7TeVINEL);
                    graphEtaYieldComb7TeVSys                            = ScaleGraph(graphEtaYieldComb7TeVSys,      1./xSection7TeVINEL);
                    graphEtaYieldPCM7TeVStat                            = ScaleGraph(graphEtaYieldPCM7TeVStat,      1./xSection7TeVINEL);
                    graphEtaYieldPCM7TeVSys                             = ScaleGraph(graphEtaYieldPCM7TeVSys,       1./xSection7TeVINEL);
                    graphEtaYieldEMCal7TeVStat                          = ScaleGraph(graphEtaYieldEMCal7TeVStat,    1./xSection7TeVINEL);
                    graphEtaYieldEMCal7TeVSys                           = ScaleGraph(graphEtaYieldEMCal7TeVSys,     1./xSection7TeVINEL);
                    graphEtaYieldPCMEMC7TeVStat                         = ScaleGraph(graphEtaYieldPCMEMC7TeVStat,   1./xSection7TeVINEL);
                    graphEtaYieldPCMEMC7TeVSys                          = ScaleGraph(graphEtaYieldPCMEMC7TeVSys,    1./xSection7TeVINEL);
                    graphEtaYieldPCMPass27TeVStat                       = ScaleGraph(graphEtaYieldPCMPass27TeVStat, 1./xSection7TeVINEL);
                    graphEtaYieldPCMPass27TeVSys                        = ScaleGraph(graphEtaYieldPCMPass27TeVSys,  1./xSection7TeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphEtaYieldComb7TeVStat                           = ConvertYieldGraph(graphEtaYieldComb7TeVStat,      kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldComb7TeVSys                            = ConvertYieldGraph(graphEtaYieldComb7TeVSys,       kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCM7TeVStat                            = ConvertYieldGraph(graphEtaYieldPCM7TeVStat,       kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCM7TeVSys                             = ConvertYieldGraph(graphEtaYieldPCM7TeVSys,        kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldEMCal7TeVStat                          = ConvertYieldGraph(graphEtaYieldEMCal7TeVStat,     kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldEMCal7TeVSys                           = ConvertYieldGraph(graphEtaYieldEMCal7TeVSys,      kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCMEMC7TeVStat                         = ConvertYieldGraph(graphEtaYieldPCMEMC7TeVStat,    kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCMEMC7TeVSys                          = ConvertYieldGraph(graphEtaYieldPCMEMC7TeVSys,     kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCMPass27TeVStat                       = ConvertYieldGraph(graphEtaYieldPCMPass27TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaYieldPCMPass27TeVSys                        = ConvertYieldGraph(graphEtaYieldPCMPass27TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphEtaYieldComb7TeVStat,       Form("%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldComb7TeVSys,        Form("%s%sSys",  fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCM7TeVStat,        Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCM7TeVSys,         Form("%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldEMCal7TeVStat,      Form("%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldEMCal7TeVSys,       Form("%s%sSys",  fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCMEMC7TeVStat,     Form("%s%sStat", fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCMEMC7TeVSys,      Form("%s%sSys",  fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCMPass27TeVStat,   Form("%s%sStat", fParticle[1].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaYieldPCMPass27TeVSys,    Form("%s%sSys",  fParticle[1].Data(), fMethod[11].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_7TeV->Add(graphEtaYieldComb7TeVStat);
                    list_7TeV->Add(graphEtaYieldComb7TeVSys);
                    list_7TeV->Add(graphEtaYieldPCM7TeVStat);
                    list_7TeV->Add(graphEtaYieldPCM7TeVSys);
                    list_7TeV->Add(graphEtaYieldEMCal7TeVStat);
                    list_7TeV->Add(graphEtaYieldEMCal7TeVSys);
                    list_7TeV->Add(graphEtaYieldPCMEMC7TeVStat);
                    list_7TeV->Add(graphEtaYieldPCMEMC7TeVSys);
                    list_7TeV->Add(graphEtaYieldPCMPass27TeVStat);
                    list_7TeV->Add(graphEtaYieldPCMPass27TeVSys);
                }
            }
        }

        if (enable_CPion && enable_CKaon && enable_Proton ){
            //================================================================================================================
            // reading and writing pi+-, K+- and p/bar{p} to 7TeV list
            // input is given as fully invariant yield
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile *fileIdentifiedCharged7TeV                        = new TFile("pp/Pi_K_P_7TeV_INELSpectra_Paper2016_20150803.root");

            TH1D* histoCPionComb7TeVStat                            = (TH1D*)fileIdentifiedCharged7TeV->Get("hstat_pp7_pion_sum");
            TH1D* histoCPionComb7TeVSys                             = (TH1D*)fileIdentifiedCharged7TeV->Get("hsys_pp7_pion_sum");
            TH1D* histoCKaonComb7TeVStat                            = (TH1D*)fileIdentifiedCharged7TeV->Get("hstat_pp7_kaon_sum");
            TH1D* histoCKaonComb7TeVSys                             = (TH1D*)fileIdentifiedCharged7TeV->Get("hsys_pp7_kaon_sum");
            TH1D* histoProtonComb7TeVStat                           = (TH1D*)fileIdentifiedCharged7TeV->Get("hstat_pp7_proton_sum");
            TH1D* histoProtonComb7TeVSys                            = (TH1D*)fileIdentifiedCharged7TeV->Get("hsys_pp7_proton_sum");

            // scale by 0.5 to get averaged single particle
            histoCPionComb7TeVStat->Scale(  0.5);
            histoCPionComb7TeVSys->Scale(   0.5);
            histoCKaonComb7TeVStat->Scale(  0.5);
            histoCKaonComb7TeVSys->Scale(   0.5);
            histoProtonComb7TeVStat->Scale( 0.5);
            histoProtonComb7TeVSys->Scale(  0.5);

            if (changeQuantity == 0 || changeQuantity == 2){
                SetHistoProperties(histoCPionComb7TeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCPionComb7TeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCKaonComb7TeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoCKaonComb7TeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoProtonComb7TeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetHistoProperties(histoProtonComb7TeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(histoCPionComb7TeVStat);
                list_7TeV->Add(histoCPionComb7TeVSys);
                list_7TeV->Add(histoCKaonComb7TeVStat);
                list_7TeV->Add(histoCKaonComb7TeVSys);
                list_7TeV->Add(histoProtonComb7TeVStat);
                list_7TeV->Add(histoProtonComb7TeVSys);
            } else {
                histoCPionComb7TeVStat                              = ConvertYieldHisto(histoCPionComb7TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                histoCPionComb7TeVSys                               = ConvertYieldHisto(histoCPionComb7TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                histoCKaonComb7TeVStat                              = ConvertYieldHisto(histoCKaonComb7TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                histoCKaonComb7TeVSys                               = ConvertYieldHisto(histoCKaonComb7TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                histoProtonComb7TeVStat                             = ConvertYieldHisto(histoProtonComb7TeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                histoProtonComb7TeVSys                              = ConvertYieldHisto(histoProtonComb7TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                SetHistoProperties(histoCPionComb7TeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCPionComb7TeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCKaonComb7TeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoCKaonComb7TeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoProtonComb7TeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetHistoProperties(histoProtonComb7TeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(histoCPionComb7TeVStat);
                list_7TeV->Add(histoCPionComb7TeVSys);
                list_7TeV->Add(histoCKaonComb7TeVStat);
                list_7TeV->Add(histoCKaonComb7TeVSys);
                list_7TeV->Add(histoProtonComb7TeVStat);
                list_7TeV->Add(histoProtonComb7TeVSys);
            }

            TFile *fileIdentifiedChargedRatios7TeV                  = new TFile("pp/Pi_K_P_7TeV_mb_RATIOS_Paper2016_20150803.root");

            TH1D* histoCKaonToCPionComb7TeVStat                     = (TH1D*)fileIdentifiedChargedRatios7TeV->Get("hstat_pp7_kaon_to_pion_sum");
            TH1D* histoCKaonToCPionComb7TeVSys                      = (TH1D*)fileIdentifiedChargedRatios7TeV->Get("hsys_pp7_kaon_to_pion_sum");
            TH1D* histoProtonToCKaonComb7TeVStat                    = (TH1D*)fileIdentifiedChargedRatios7TeV->Get("hstat_pp7_proton_to_pion_sum");
            TH1D* histoProtonToCKaonComb7TeVSys                     = (TH1D*)fileIdentifiedChargedRatios7TeV->Get("hsys_pp7_proton_to_pion_sum");

            SetHistoProperties(histoCKaonToCPionComb7TeVStat,   Form("%sTo%s%sStat", fParticle[6].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");
            SetHistoProperties(histoCKaonToCPionComb7TeVSys,    Form("%sTo%s%sSys",  fParticle[6].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "(K^{+} + K^{-}) / (#pi^{+} + #pi^{-})", "");
            SetHistoProperties(histoProtonToCKaonComb7TeVStat,  Form("%sTo%s%sStat", fParticle[7].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");
            SetHistoProperties(histoProtonToCKaonComb7TeVSys,   Form("%sTo%s%sSys",  fParticle[7].Data(), fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "(p + #bar{p}) / (#pi^{+} + #pi^{-})", "");

            list_7TeV->Add(histoCKaonToCPionComb7TeVStat);
            list_7TeV->Add(histoCKaonToCPionComb7TeVSys);
            list_7TeV->Add(histoProtonToCKaonComb7TeVStat);
            list_7TeV->Add(histoProtonToCKaonComb7TeVSys);
        }

        if (enable_OmegaMeson ){
            //================================================================================================================
            // reading and writing omega to 7TeV list
            // input is given as invariant cross section
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileOmegaMesonPHOS7TeV                          = new TFile("pp/PHOS_pp_omega_7TeV_07082012.root");

            TGraphAsymmErrors* graphOmegaXSecPHOS7TeVStat           = (TGraphAsymmErrors*)fFileOmegaMesonPHOS7TeV->Get("graphOmegaStat");
            TGraphAsymmErrors* graphOmegaXSecPHOS7TeVSys            = (TGraphAsymmErrors*)fFileOmegaMesonPHOS7TeV->Get("graphOmegaSyst");

            TGraphAsymmErrors* graphOmegaToPi0Comb7TeVStatTemp      = (TGraphAsymmErrors*)fFileOmegaMesonPHOS7TeV->Get("omega_to_pi0_stat_err");
            Int_t nPoints                                           = graphOmegaToPi0Comb7TeVStatTemp->GetN();
            Double_t* xVal                                          = graphOmegaToPi0Comb7TeVStatTemp->GetX();
            Double_t* xErrUp                                        = graphOmegaXSecPHOS7TeVStat->GetEXhigh();
            Double_t* xErrDown                                      = graphOmegaXSecPHOS7TeVStat->GetEXlow();
            Double_t* yVal                                          = graphOmegaToPi0Comb7TeVStatTemp->GetY();
            Double_t* yErr                                          = graphOmegaToPi0Comb7TeVStatTemp->GetEY();
            TGraphAsymmErrors* graphOmegaToPi0Comb7TeVStat          = new TGraphAsymmErrors(nPoints, xVal, yVal, xErrDown, xErrUp, yErr, yErr);
            TGraphAsymmErrors* graphOmegaToPi0Comb7TeVSys           = (TGraphAsymmErrors*)graphOmegaToPi0Comb7TeVStat->Clone("clone");

            SetGraphProperties(graphOmegaToPi0Comb7TeVStat,   Form("%sTo%s%sStat", fParticle[2].Data(), fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#omega/#pi^{0}", "");
            SetGraphProperties(graphOmegaToPi0Comb7TeVSys,    Form("%sTo%s%sSys",  fParticle[2].Data(), fParticle[0].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", "#omega/#pi^{0}", "");
            list_7TeV->Add(graphOmegaToPi0Comb7TeVStat);
            list_7TeV->Add(graphOmegaToPi0Comb7TeVSys);
            SetGraphProperties(graphOmegaToPi0Comb7TeVStat,   Form("%sTo%s%sStat", fParticle[2].Data(), fParticle[0].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#omega/#pi^{0}", "");
            SetGraphProperties(graphOmegaToPi0Comb7TeVSys,    Form("%sTo%s%sSys",  fParticle[2].Data(), fParticle[0].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "#omega/#pi^{0}", "");
            list_7TeV->Add(graphOmegaToPi0Comb7TeVStat);
            list_7TeV->Add(graphOmegaToPi0Comb7TeVSys);

            if (changeQuantity == 0){
                // simply add Xsection to the output file
                SetGraphProperties(graphOmegaXSecPHOS7TeVStat,    Form("Xsec_%s%sStat", fParticle[2].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphOmegaXSecPHOS7TeVSys,     Form("Xsec_%s%sSys",  fParticle[2].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphOmegaXSecPHOS7TeVStat);
                list_7TeV->Add(graphOmegaXSecPHOS7TeVSys);

                SetGraphProperties(graphOmegaXSecPHOS7TeVStat,    Form("Xsec_%s%sStat", fParticle[2].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphOmegaXSecPHOS7TeVSys,     Form("Xsec_%s%sSys",  fParticle[2].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphOmegaXSecPHOS7TeVStat);
                list_7TeV->Add(graphOmegaXSecPHOS7TeVSys);
            } else {
                // convert XSection to invariant yield
                TGraphAsymmErrors* graphOmegaPHOS7TeVStat           = ScaleGraph(graphOmegaXSecPHOS7TeVStat,1./xSection7TeVINEL);
                TGraphAsymmErrors* graphOmegaPHOS7TeVSys            = ScaleGraph(graphOmegaXSecPHOS7TeVSys,1./xSection7TeVINEL);

                // convert inv yield to 1/N d^2N/dydpT
                graphOmegaPHOS7TeVStat                              = ConvertYieldGraph(graphOmegaPHOS7TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphOmegaPHOS7TeVSys                               = ConvertYieldGraph(graphOmegaPHOS7TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);

                SetGraphProperties(graphOmegaPHOS7TeVStat,    Form("%s%sStat", fParticle[2].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphOmegaPHOS7TeVSys,     Form("%s%sSys",  fParticle[2].Data(), fMethod[3].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphOmegaPHOS7TeVStat);
                list_7TeV->Add(graphOmegaPHOS7TeVSys);
                SetGraphProperties(graphOmegaPHOS7TeVStat,    Form("%s%sStat", fParticle[2].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphOmegaPHOS7TeVSys,     Form("%s%sSys",  fParticle[2].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphOmegaPHOS7TeVStat);
                list_7TeV->Add(graphOmegaPHOS7TeVSys);
            }
        }

        if (enable_Phi) {
            //================================================================================================================
            // reading and writing phi to 7TeV list
            // input from phi is given as yield per inel. event
            //================================================================================================================
            TString phiHEPDataFile                                  = "pp/HEPdata/pp7TeV_phi_arxiv1208.5717.txt";

            TGraphAsymmErrors* graphPhiYieldComb7TeVStat            = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphPhiYieldComb7TeVSys             = ParseHEPData(phiHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphPhiXSecComb7TeVStat         = ScaleGraph(graphPhiYieldComb7TeVStat,xSection7TeVINEL);
                TGraphAsymmErrors* graphPhiXSecComb7TeVSys          = ScaleGraph(graphPhiYieldComb7TeVSys,xSection7TeVINEL);

                graphPhiXSecComb7TeVStat                            = ConvertYieldGraph(graphPhiXSecComb7TeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphPhiXSecComb7TeVSys                             = ConvertYieldGraph(graphPhiXSecComb7TeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphPhiXSecComb7TeVStat,  Form("Xsec_%s%sStat", fParticle[9].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphPhiXSecComb7TeVSys,   Form("Xsec_%s%sSys",  fParticle[9].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphPhiXSecComb7TeVStat);
                list_7TeV->Add(graphPhiXSecComb7TeVSys);
            } else {
                SetGraphProperties(graphPhiYieldComb7TeVStat,    Form("%s%sStat", fParticle[9].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphPhiYieldComb7TeVSys,     Form("%s%sSys",  fParticle[9].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphPhiYieldComb7TeVStat);
                list_7TeV->Add(graphPhiYieldComb7TeVSys);

                // get phi to pi0 ratio from fit to pi0
                if ( list_7TeV->FindObject(Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data())) ) {
                    TGraphAsymmErrors* graphNPionPCM7TeVStat       = (TGraphAsymmErrors*)list_7TeV->FindObject(Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()));

                    // fit pi0 spectrum
                    Double_t mass                                   = TDatabasePDG::Instance()->GetParticle(111)->Mass();
                    TF1* graphNPionPCM7TeVStatFit                   = new TF1("graphNPionPCM7TeVStatFit", Form("x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])"), 0.0, 25.0);
                    graphNPionPCM7TeVStatFit->SetParameters(12.8, 0.5, 2.2e-10, 0.8, 6.2);
                    graphNPionPCM7TeVStatFit->SetParLimits(4, 6.05, 6.5);
                    graphNPionPCM7TeVStat->Fit(graphNPionPCM7TeVStatFit,"QNRME+","",0.4,16.0);
                    graphNPionPCM7TeVStatFit->SetName(Form("%s_tempFit", graphNPionPCM7TeVStat->GetName()));

                    // get ratio and add to list
                    TGraphAsymmErrors* graphPhiToPi0Comb7TeVStat    = (TGraphAsymmErrors*)CalculateParticleRatioWithFit(graphPhiYieldComb7TeVStat, graphNPionPCM7TeVStatFit);
                    TGraphAsymmErrors* graphPhiToPi0Comb7TeVSys     = (TGraphAsymmErrors*)CalculateParticleRatioWithFit(graphPhiYieldComb7TeVSys,  graphNPionPCM7TeVStatFit);

                    SetGraphProperties(graphPhiToPi0Comb7TeVStat,   Form("%sTo%s%sStat", fParticle[9].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/#pi^{0}", "");
                    SetGraphProperties(graphPhiToPi0Comb7TeVSys,    Form("%sTo%s%sSys",  fParticle[9].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#phi/#pi^{0}", "");

                    list_7TeV->Add(graphNPionPCM7TeVStatFit);
                    list_7TeV->Add(graphPhiToPi0Comb7TeVStat);
                    list_7TeV->Add(graphPhiToPi0Comb7TeVSys);
                }
            }
        }

        if (enable_K0s) {
            //================================================================================================================
            // reading and writing K0s to 7TeV list
            // input from K0s is given as 1/N_inel dN^2/dydpt
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileK0s7TeV                                     = new TFile("pp/K0s-pp7TeV-Preliminary.root");
            TH1F* histK0s7TeVStat                                   = (TH1F*)fFileK0s7TeV->Get("fHistPtK0ShortStatOnly");
            TH1F* histK0s7TeVSys                                    = (TH1F*)fFileK0s7TeV->Get("fHistPtK0ShortStatAndSystExceptNormalization");
            Double_t tempErrTot                                     = 0.;
            Double_t tempErrStat                                    = 0.;
            Double_t tempErrSys                                     = 0.;
            for (Int_t i=1; i<=histK0s7TeVSys->GetNbinsX(); i++) {
                tempErrTot                                          = histK0s7TeVSys->GetBinError(i);
                tempErrStat                                         = histK0s7TeVStat->GetBinError(i);
                tempErrSys                                          = TMath::Sqrt(tempErrTot*tempErrTot - tempErrStat*tempErrStat);
                histK0s7TeVSys->SetBinError(i, tempErrSys);
            }

            SetHistoProperties(histK0s7TeVStat,  Form("%s%sStat", fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histK0s7TeVSys,   Form("%s%sSys",  fParticle[15].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_7TeV->Add(histK0s7TeVStat);
            list_7TeV->Add(histK0s7TeVSys);
        }

        if (enable_Lambda) {
            //================================================================================================================
            // reading and writing Lambda to 7TeV list
            // input from Lambda is given as 1/N_inel dN^2/dydpt
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile* fFileLambda7TeV                                  = new TFile("pp/Lambda-pp7TeV-Preliminary.root");
            TH1F* histLambda7TeVStat                                = (TH1F*)fFileLambda7TeV->Get("fHistPtLambdaStatOnly");
            TH1F* histLambda7TeVSys                                 = (TH1F*)fFileLambda7TeV->Get("fHistPtLambdaStatAndSystExceptNormalization");
            histLambda7TeVStat->Sumw2();
            histLambda7TeVSys->Sumw2();
            Double_t tempErrTot                                     = 0.;
            Double_t tempErrStat                                    = 0.;
            Double_t tempErrSys                                     = 0.;
            for (Int_t i=1; i<=histLambda7TeVSys->GetNbinsX(); i++) {
                tempErrTot                                          = histLambda7TeVSys->GetBinError(i);
                tempErrStat                                         = histLambda7TeVStat->GetBinError(i);
                tempErrSys                                          = TMath::Sqrt(tempErrTot*tempErrTot - tempErrStat*tempErrStat);
                histLambda7TeVSys->SetBinError(i, tempErrSys);
            }

            TFile* fFileAntiLambda7TeV                              = new TFile("pp/AntiLambda-pp7TeV-Preliminary.root");
            TH1F* histAntiLambda7TeVStat                            = (TH1F*)fFileAntiLambda7TeV->Get("fHistPtAntiLambdaStatOnly");
            TH1F* histAntiLambda7TeVSys                             = (TH1F*)fFileAntiLambda7TeV->Get("fHistPtAntiLambdaStatAndSystExceptNormalization");
            histAntiLambda7TeVStat->Sumw2();
            histAntiLambda7TeVSys->Sumw2();
            tempErrTot                                              = 0.;
            tempErrStat                                             = 0.;
            tempErrSys                                              = 0.;
            for (Int_t i=1; i<=histAntiLambda7TeVSys->GetNbinsX(); i++) {
                tempErrTot                                          = histAntiLambda7TeVSys->GetBinError(i);
                tempErrStat                                         = histAntiLambda7TeVStat->GetBinError(i);
                tempErrSys                                          = TMath::Sqrt(tempErrTot*tempErrTot - tempErrStat*tempErrStat);
                histAntiLambda7TeVSys->SetBinError(i, tempErrSys);
            }

            TH1F* histLambdaComb7TeVStat                            = (TH1F*)histLambda7TeVStat->Clone("histLambdaComb7TeVStat");
            histLambdaComb7TeVStat->Add(histAntiLambda7TeVStat);
            histLambdaComb7TeVStat->Scale(0.5);
            TH1F* histLambdaComb7TeVSys                             = (TH1F*)histLambda7TeVSys->Clone("histLambdaComb7TeVSys");
            histLambdaComb7TeVSys->Add(histAntiLambda7TeVSys);
            histLambdaComb7TeVSys->Scale(0.5);

            SetHistoProperties(histLambdaComb7TeVStat,  Form("%s%sStat", fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histLambdaComb7TeVSys,   Form("%s%sSys",  fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_7TeV->Add(histLambdaComb7TeVStat);
            list_7TeV->Add(histLambdaComb7TeVSys);
        }

        if (enable_K0star) {
            //================================================================================================================
            // reading and writing K0* to 7TeV list
            // input from K0* is given as yield per inel. event
            //================================================================================================================
            TString K0starHEPDataFile                               = "pp/HEPdata/pp7TeV_K0star_arxiv1208.5717.txt";

            TGraphAsymmErrors* graphK0starYieldComb7TeVStat         = ParseHEPData(K0starHEPDataFile, 8, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphK0starYieldComb7TeVSys          = ParseHEPData(K0starHEPDataFile, 8, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphK0starXSecComb7TeVStat      = ScaleGraph(graphK0starYieldComb7TeVStat,xSection7TeVINEL);
                TGraphAsymmErrors* graphK0starXSecComb7TeVSys       = ScaleGraph(graphK0starYieldComb7TeVSys,xSection7TeVINEL);

                graphK0starXSecComb7TeVStat                         = ConvertYieldGraph(graphK0starXSecComb7TeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphK0starXSecComb7TeVSys                          = ConvertYieldGraph(graphK0starXSecComb7TeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphK0starXSecComb7TeVStat,  Form("Xsec_%s%sStat", fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphK0starXSecComb7TeVSys,   Form("Xsec_%s%sSys",  fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphK0starXSecComb7TeVStat);
                list_7TeV->Add(graphK0starXSecComb7TeVSys);
            } else {
                SetGraphProperties(graphK0starYieldComb7TeVStat,    Form("%s%sStat", fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphK0starYieldComb7TeVSys,     Form("%s%sSys",  fParticle[10].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphK0starYieldComb7TeVStat);
                list_7TeV->Add(graphK0starYieldComb7TeVSys);

                // get K0star to pi0 ratio from fit to pi0
                if ( list_7TeV->FindObject(Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data())) ) {
                    TGraphAsymmErrors* graphNPionPCM7TeVStat       = (TGraphAsymmErrors*)list_7TeV->FindObject(Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()));

                    // fit pi0 spectrum
                    Double_t mass                                   = TDatabasePDG::Instance()->GetParticle(111)->Mass();
                    TF1* graphNPionPCM7TeVStatFit                   = new TF1("graphNPionPCM7TeVStatFit", Form("x*[0]*pow(exp(-[1]*x-abs([2])*x*x)+x/[3],-[4])"), 0.0, 25.0);
                    graphNPionPCM7TeVStatFit->SetParameters(12.8, 0.5, 2.2e-10, 0.8, 6.2);
                    graphNPionPCM7TeVStatFit->SetParLimits(4, 6.05, 6.5);
                    graphNPionPCM7TeVStat->Fit(graphNPionPCM7TeVStatFit,"QNRME+","",0.4,16.0);
                    graphNPionPCM7TeVStatFit->SetName(Form("%s_tempFit", graphNPionPCM7TeVStat->GetName()));

                    // get ratio and add to list
                    TGraphAsymmErrors* graphK0starToPi0Comb7TeVStat = (TGraphAsymmErrors*)CalculateParticleRatioWithFit(graphK0starYieldComb7TeVStat, graphNPionPCM7TeVStatFit);
                    TGraphAsymmErrors* graphK0starToPi0Comb7TeVSys  = (TGraphAsymmErrors*)CalculateParticleRatioWithFit(graphK0starYieldComb7TeVSys,  graphNPionPCM7TeVStatFit);

                    SetGraphProperties(graphK0starToPi0Comb7TeVStat,   Form("%sTo%s%sStat", fParticle[10].Data(), fParticle[0].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K*^{0}/#pi^{0}", "");
                    SetGraphProperties(graphK0starToPi0Comb7TeVSys,    Form("%sTo%s%sSys",  fParticle[10].Data(), fParticle[0].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", "K*^{0}/#pi^{0}", "");

                    //list_7TeV->Add(graphNPionPCM7TeVStatFit);
                    list_7TeV->Add(graphK0starToPi0Comb7TeVStat);
                    list_7TeV->Add(graphK0starToPi0Comb7TeVSys);
                }
            }
        }

        if (enable_COmega) {
            //================================================================================================================
            // reading and writing Omega+- to 7TeV list
            // input from Omega+- is given as yield per inel. event
            //================================================================================================================
            TString chargedOmegaHEPDataFile                         = "pp/HEPdata/pp7TeV_Omega+-_arxiv1204.0282.txt";

            TGraphAsymmErrors* graphOmegaMinusYield7TeVStat         = ParseHEPData(chargedOmegaHEPDataFile, 13, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphOmegaMinusYield7TeVSys          = ParseHEPData(chargedOmegaHEPDataFile, 13, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            TGraphAsymmErrors* graphOmegaPlusYield7TeVStat          = ParseHEPData(chargedOmegaHEPDataFile, 13, 0, 1, 2, 8, 9, 10, kFALSE, kTRUE);
            TGraphAsymmErrors* graphOmegaPlusYield7TeVSys           = ParseHEPData(chargedOmegaHEPDataFile, 13, 0, 1, 2, 8, 11, 12, kFALSE, kTRUE);

            TGraphAsymmErrors* graphChargedOmegaYield7TeVStat       = (TGraphAsymmErrors*)graphOmegaMinusYield7TeVStat->Clone("graphChargedOmegaYield7TeVStat");
            graphChargedOmegaYield7TeVStat                          = Add2TGraphErrorsSameBinning(graphChargedOmegaYield7TeVStat, graphOmegaPlusYield7TeVStat);
            graphChargedOmegaYield7TeVStat                          = ScaleGraph(graphChargedOmegaYield7TeVStat, 0.5);

            TGraphAsymmErrors* graphChargedOmegaYield7TeVSys        = (TGraphAsymmErrors*)graphOmegaMinusYield7TeVSys->Clone("graphChargedOmegaYield7TeVSys");
            graphChargedOmegaYield7TeVSys                           = Add2TGraphErrorsSameBinning(graphChargedOmegaYield7TeVSys, graphOmegaPlusYield7TeVSys);
            graphChargedOmegaYield7TeVSys                           = ScaleGraph(graphChargedOmegaYield7TeVSys, 0.5);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphChargedOmegaXSec7TeVStat    = ScaleGraph(graphChargedOmegaYield7TeVStat,xSection7TeVINEL);
                TGraphAsymmErrors* graphChargedOmegaXSec7TeVSys     = ScaleGraph(graphChargedOmegaYield7TeVSys,xSection7TeVINEL);

                graphChargedOmegaXSec7TeVStat                       = ConvertYieldGraph(graphChargedOmegaXSec7TeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphChargedOmegaXSec7TeVSys                        = ConvertYieldGraph(graphChargedOmegaXSec7TeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphChargedOmegaXSec7TeVStat,  Form("Xsec_%s%sStat", fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphChargedOmegaXSec7TeVSys,   Form("Xsec_%s%sSys",  fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphChargedOmegaXSec7TeVStat);
                list_7TeV->Add(graphChargedOmegaXSec7TeVSys);
            } else {
                SetGraphProperties(graphChargedOmegaYield7TeVStat,    Form("%s%sStat", fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphChargedOmegaYield7TeVSys,     Form("%s%sSys",  fParticle[19].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphChargedOmegaYield7TeVStat);
                list_7TeV->Add(graphChargedOmegaYield7TeVSys);
            }
        }

        if (enable_CXi) {
            //================================================================================================================
            // reading and writing Xi+- to 7TeV list
            // input from Xi+- is given as yield per inel. event
            //================================================================================================================
            TString chargedXiHEPDataFile                            = "pp/HEPdata/pp7TeV_Xi+-_arxiv1406.3206v2.txt";

            TGraphAsymmErrors* graphXiMinusYield7TeVStat            = ParseHEPData(chargedXiHEPDataFile, 13, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphXiMinusYield7TeVSys             = ParseHEPData(chargedXiHEPDataFile, 13, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            TGraphAsymmErrors* graphXiPlusYield7TeVStat             = ParseHEPData(chargedXiHEPDataFile, 13, 0, 1, 2, 8, 9, 10, kFALSE, kTRUE);
            TGraphAsymmErrors* graphXiPlusYield7TeVSys              = ParseHEPData(chargedXiHEPDataFile, 13, 0, 1, 2, 8, 11, 12, kFALSE, kTRUE);

            TGraphAsymmErrors* graphChargedXiYield7TeVStat          = (TGraphAsymmErrors*)graphXiMinusYield7TeVStat->Clone("graphChargedXiYield7TeVStat");
            graphChargedXiYield7TeVStat                             = Add2TGraphErrorsSameBinning(graphChargedXiYield7TeVStat, graphXiPlusYield7TeVStat);
            graphChargedXiYield7TeVStat                             = ScaleGraph(graphChargedXiYield7TeVStat, 0.5);

            TGraphAsymmErrors* graphChargedXiYield7TeVSys           = (TGraphAsymmErrors*)graphXiMinusYield7TeVSys->Clone("graphChargedXiYield7TeVSys");
            graphChargedXiYield7TeVSys                              = Add2TGraphErrorsSameBinning(graphChargedXiYield7TeVSys, graphXiPlusYield7TeVSys);
            graphChargedXiYield7TeVSys                              = ScaleGraph(graphChargedXiYield7TeVSys, 0.5);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphChargedXiXSec7TeVStat       = ScaleGraph(graphChargedXiYield7TeVStat,xSection7TeVINEL);
                TGraphAsymmErrors* graphChargedXiXSec7TeVSys        = ScaleGraph(graphChargedXiYield7TeVSys,xSection7TeVINEL);

                graphChargedXiXSec7TeVStat                          = ConvertYieldGraph(graphChargedXiXSec7TeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphChargedXiXSec7TeVSys                           = ConvertYieldGraph(graphChargedXiXSec7TeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphChargedXiXSec7TeVStat,  Form("Xsec_%s%sStat", fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphChargedXiXSec7TeVSys,   Form("Xsec_%s%sSys",  fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphChargedXiXSec7TeVStat);
                list_7TeV->Add(graphChargedXiXSec7TeVSys);
            } else {
                SetGraphProperties(graphChargedXiYield7TeVStat,    Form("%s%sStat", fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphChargedXiYield7TeVSys,     Form("%s%sSys",  fParticle[20].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphChargedXiYield7TeVStat);
                list_7TeV->Add(graphChargedXiYield7TeVSys);
            }
        }

        if (enable_CSigma) {
            //================================================================================================================
            // reading and writing Sigma(1385)+- to 7TeV list
            // input from Sigma(1385)+- is given as yield per inel. event
            //================================================================================================================
            TString chargedSigmaHEPDataFile                         = "pp/HEPdata/pp7TeV_Sigma+-_arxiv1406.3206v2.txt";

            TGraphAsymmErrors* graphSigmaPlusYield7TeVStat          = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphSigmaPlusYield7TeVSys           = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 3, 6, 7, kFALSE, kTRUE);

            TGraphAsymmErrors* graphSigmaMinusYield7TeVStat         = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 8, 9, 10, kFALSE, kTRUE);
            TGraphAsymmErrors* graphSigmaMinusYield7TeVSys          = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 8, 11, 12, kFALSE, kTRUE);

            TGraphAsymmErrors* graphSigmaBarMinusYield7TeVStat      = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 13, 14, 15, kFALSE, kTRUE);
            TGraphAsymmErrors* graphSigmaBarMinusYield7TeVSys       = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 13, 16, 17, kFALSE, kTRUE);

            TGraphAsymmErrors* graphSigmaBarPlusYield7TeVStat       = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 18, 19, 20, kFALSE, kTRUE);
            TGraphAsymmErrors* graphSigmaBarPlusYield7TeVSys        = ParseHEPData(chargedSigmaHEPDataFile, 23, 0, 1, 2, 18, 21, 22, kFALSE, kTRUE);

            // charged sigma
            TGraphAsymmErrors* graphChargedSigmaYield7TeVStat       = (TGraphAsymmErrors*)graphSigmaPlusYield7TeVStat->Clone("graphChargedSigmaYield7TeVStat");
            graphChargedSigmaYield7TeVStat                          = Add2TGraphErrorsSameBinning(graphChargedSigmaYield7TeVStat, graphSigmaMinusYield7TeVStat);
            graphChargedSigmaYield7TeVStat                          = ScaleGraph(graphChargedSigmaYield7TeVStat, 0.5);

            TGraphAsymmErrors* graphChargedSigmaYield7TeVSys        = (TGraphAsymmErrors*)graphSigmaPlusYield7TeVSys->Clone("graphChargedSigmaYield7TeVSys");
            graphChargedSigmaYield7TeVSys                           = Add2TGraphErrorsSameBinning(graphChargedSigmaYield7TeVSys, graphSigmaMinusYield7TeVSys);
            graphChargedSigmaYield7TeVSys                           = ScaleGraph(graphChargedSigmaYield7TeVSys, 0.5);

            // charged sigma bar
            TGraphAsymmErrors* graphChargedSigmaBarYield7TeVStat    = (TGraphAsymmErrors*)graphSigmaBarPlusYield7TeVStat->Clone("graphChargedSigmaBarYield7TeVStat");
            graphChargedSigmaBarYield7TeVStat                       = Add2TGraphErrorsSameBinning(graphChargedSigmaBarYield7TeVStat, graphSigmaBarMinusYield7TeVStat);
            graphChargedSigmaBarYield7TeVStat                       = ScaleGraph(graphChargedSigmaBarYield7TeVStat, 0.5);

            TGraphAsymmErrors* graphChargedSigmaBarYield7TeVSys     = (TGraphAsymmErrors*)graphSigmaBarPlusYield7TeVSys->Clone("graphChargedSigmaBarYield7TeVSys");
            graphChargedSigmaBarYield7TeVSys                        = Add2TGraphErrorsSameBinning(graphChargedSigmaBarYield7TeVSys, graphSigmaBarMinusYield7TeVSys);
            graphChargedSigmaBarYield7TeVSys                        = ScaleGraph(graphChargedSigmaBarYield7TeVSys, 0.5);

            // mean
            TGraphAsymmErrors* graphChargedSigmaCombYield7TeVStat   = (TGraphAsymmErrors*)graphChargedSigmaYield7TeVStat->Clone("graphChargedSigmaCombYield7TeVStat");
            graphChargedSigmaCombYield7TeVStat                      = Add2TGraphErrorsSameBinning(graphChargedSigmaCombYield7TeVStat, graphChargedSigmaBarYield7TeVStat);
            graphChargedSigmaCombYield7TeVStat                      = ScaleGraph(graphChargedSigmaCombYield7TeVStat, 0.5);

            TGraphAsymmErrors* graphChargedSigmaCombYield7TeVSys    = (TGraphAsymmErrors*)graphChargedSigmaYield7TeVSys->Clone("graphChargedSigmaCombYield7TeVSys");
            graphChargedSigmaCombYield7TeVSys                       = Add2TGraphErrorsSameBinning(graphChargedSigmaCombYield7TeVSys, graphChargedSigmaBarYield7TeVSys);
            graphChargedSigmaCombYield7TeVSys                       = ScaleGraph(graphChargedSigmaCombYield7TeVSys, 0.5);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                // calculate xSec
                TGraphAsymmErrors* graphChargedSigmaCombXSec7TeVStat    = ScaleGraph(graphChargedSigmaCombYield7TeVStat,xSection7TeVINEL);
                TGraphAsymmErrors* graphChargedSigmaCombXSec7TeVSys     = ScaleGraph(graphChargedSigmaCombYield7TeVSys,xSection7TeVINEL);

                graphChargedSigmaCombXSec7TeVStat                   = ConvertYieldGraph(graphChargedSigmaCombXSec7TeVStat, kTRUE, kTRUE, kFALSE, kFALSE);
                graphChargedSigmaCombXSec7TeVSys                    = ConvertYieldGraph(graphChargedSigmaCombXSec7TeVSys,  kTRUE, kTRUE, kFALSE, kFALSE);

                SetGraphProperties(graphChargedSigmaCombXSec7TeVStat,  Form("Xsec_%s%sStat", fParticle[18].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphChargedSigmaCombXSec7TeVSys,   Form("Xsec_%s%sSys",  fParticle[18].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphChargedSigmaCombXSec7TeVStat);
                list_7TeV->Add(graphChargedSigmaCombXSec7TeVSys);
            } else {
                SetGraphProperties(graphChargedSigmaCombYield7TeVStat,    Form("%s%sStat", fParticle[18].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphChargedSigmaCombYield7TeVSys,     Form("%s%sSys",  fParticle[18].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphChargedSigmaCombYield7TeVStat);
                list_7TeV->Add(graphChargedSigmaCombYield7TeVSys);
            }
        }

        if (enable_JPsi) {
            //================================================================================================================
            // reading and writing J/Psi to 7TeV list
            // input from J/Psi is given as inv. cross section per inel. event
            //================================================================================================================
            TString JPsiHEPDataFile                                 = "pp/HEPdata/pp7TeV_JPsi_y0.9_arxiv1105.0380.txt";

            TGraphAsymmErrors* graphJPsiXSec7TeVStat                = ParseHEPData(JPsiHEPDataFile, 14, 0, 1, 2, 3, 4, 5, kFALSE, kTRUE);
            TGraphAsymmErrors* graphJPsiXSec7TeVSys                 = ParseHEPData(JPsiHEPDataFile, 14, 0, 1, 2, 3, 8, 9, kFALSE, kTRUE);

            if (changeQuantity == 0){
                TString localOutputQuantity                         = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

                SetGraphProperties(graphJPsiXSec7TeVStat,  Form("Xsec_%s%sStat", fParticle[21].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphJPsiXSec7TeVSys,   Form("Xsec_%s%sSys",  fParticle[21].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_7TeV->Add(graphJPsiXSec7TeVStat);
                list_7TeV->Add(graphJPsiXSec7TeVSys);
            } else {
                TGraphAsymmErrors* graphJPsiYield7TeVStat           = ScaleGraph(graphJPsiXSec7TeVStat, 1/xSection7TeVINEL);
                TGraphAsymmErrors* graphJPsiYield7TeVSys            = ScaleGraph(graphJPsiXSec7TeVSys,  1/xSection7TeVINEL);

                graphJPsiYield7TeVStat                              = ConvertYieldGraph(graphJPsiYield7TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphJPsiYield7TeVSys                               = ConvertYieldGraph(graphJPsiYield7TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);

                SetGraphProperties(graphJPsiYield7TeVStat,    Form("%s%sStat", fParticle[21].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphJPsiYield7TeVSys,     Form("%s%sSys",  fParticle[21].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_7TeV->Add(graphJPsiYield7TeVStat);
                list_7TeV->Add(graphJPsiYield7TeVSys);
            }
        }
    }

    //================================================================================================================
    // creating histos and graphs for 7TeV
    //================================================================================================================
    if (Include_8TeV){

        //================================================================================================================
        // reading and writing pi0 and eta to 8TeV list
        // input from pi0 and eta is given as inv. cross section or inv. yield (no bin shift available)
        //================================================================================================================
        if (!parametrizeMC && (enable_NPion || enable_Eta)){
            TString localOutputQuantity                                 = "#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3})";

            TFile* fFileNeutralMeson8TeVComb                            = new TFile("pp/CombinedResultsPaperPP8TeV_2017_07_12.root");
            TDirectoryFile* fFileNeutralMeson8TeVCombPi0                = (TDirectoryFile*)fFileNeutralMeson8TeVComb->Get("Pi08TeV");
            TDirectoryFile* fFileNeutralMeson8TeVCombEta                = (TDirectoryFile*)fFileNeutralMeson8TeVComb->Get("Eta8TeV");

            if (enable_NPion){
                TGraphAsymmErrors* graphNPionXSecPCM8TeVStat            = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0PCM8TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCM8TeVSys             = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0PCM8TeVSysErr");
                TGraphAsymmErrors* graphNPionXSecEMCAL8TeVStat          = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0EMCAL8TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecEMCAL8TeVSys           = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0EMCAL8TeVSysErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMCAL8TeVStat       = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0PCMEMCAL8TeVStatErr");
                TGraphAsymmErrors* graphNPionXSecPCMEMCAL8TeVSys        = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0PCMEMCAL8TeVSysErr");
                TGraphAsymmErrors* graphNPionXSecComb8TeVStat           = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0Comb8TeVAStatErr");
                TGraphAsymmErrors* graphNPionXSecComb8TeVSys            = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombPi0->Get("graphInvCrossSectionPi0Comb8TeVASysErr");

                if (changeQuantity == 0){
                    // simply add Xsection to the output file
                    SetGraphProperties(graphNPionXSecPCM8TeVStat,     Form("Xsec_%s%sStat",    fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCM8TeVSys,      Form("Xsec_%s%sSys",     fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecEMCAL8TeVStat,   Form("Xsec_%s%sStat",    fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecEMCAL8TeVSys,    Form("Xsec_%s%sSys",     fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMEMCAL8TeVStat,Form("Xsec_%s%sStat",    fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecPCMEMCAL8TeVSys, Form("Xsec_%s%sSys",     fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecComb8TeVStat,    Form("Xsec_%s%sStat",    fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphNPionXSecComb8TeVSys,     Form("Xsec_%s%sSys",     fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    list_8TeV->Add(graphNPionXSecPCM8TeVStat);
                    list_8TeV->Add(graphNPionXSecPCM8TeVSys);
                    list_8TeV->Add(graphNPionXSecEMCAL8TeVStat);
                    list_8TeV->Add(graphNPionXSecEMCAL8TeVSys);
                    list_8TeV->Add(graphNPionXSecPCMEMCAL8TeVStat);
                    list_8TeV->Add(graphNPionXSecPCMEMCAL8TeVSys);
                    list_8TeV->Add(graphNPionXSecComb8TeVStat);
                    list_8TeV->Add(graphNPionXSecComb8TeVSys);
                } else {
                    // convert XSection to invariant yield
                    TGraphAsymmErrors* graphNPionPCM8TeVStat        = ScaleGraph(graphNPionXSecPCM8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionPCM8TeVSys         = ScaleGraph(graphNPionXSecPCM8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionEMCAL8TeVStat      = ScaleGraph(graphNPionXSecEMCAL8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionEMCAL8TeVSys       = ScaleGraph(graphNPionXSecEMCAL8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionPCMEMCAL8TeVStat   = ScaleGraph(graphNPionXSecPCMEMCAL8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionPCMEMCAL8TeVSys    = ScaleGraph(graphNPionXSecPCMEMCAL8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionComb8TeVStat       = ScaleGraph(graphNPionXSecComb8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphNPionComb8TeVSys        = ScaleGraph(graphNPionXSecComb8TeVSys,1./xSection8TeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphNPionPCM8TeVStat                           = ConvertYieldGraph(graphNPionPCM8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCM8TeVSys                            = ConvertYieldGraph(graphNPionPCM8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCAL8TeVStat                         = ConvertYieldGraph(graphNPionEMCAL8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionEMCAL8TeVSys                          = ConvertYieldGraph(graphNPionEMCAL8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCMEMCAL8TeVStat                      = ConvertYieldGraph(graphNPionPCMEMCAL8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionPCMEMCAL8TeVSys                       = ConvertYieldGraph(graphNPionPCMEMCAL8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionComb8TeVStat                          = ConvertYieldGraph(graphNPionComb8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphNPionComb8TeVSys                           = ConvertYieldGraph(graphNPionComb8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphNPionPCM8TeVStat,       Form("%s%sStat", fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCM8TeVSys,        Form("%s%sSys",  fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCAL8TeVStat,     Form("%s%sStat", fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionEMCAL8TeVSys,      Form("%s%sSys",  fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCMEMCAL8TeVStat,  Form("%s%sStat", fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionPCMEMCAL8TeVSys,   Form("%s%sSys",  fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionComb8TeVStat,      Form("%s%sStat", fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphNPionComb8TeVSys,       Form("%s%sSys",  fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_8TeV->Add(graphNPionPCM8TeVStat);
                    list_8TeV->Add(graphNPionPCM8TeVSys);
                    list_8TeV->Add(graphNPionEMCAL8TeVStat);
                    list_8TeV->Add(graphNPionEMCAL8TeVSys);
                    list_8TeV->Add(graphNPionPCMEMCAL8TeVStat);
                    list_8TeV->Add(graphNPionPCMEMCAL8TeVSys);
                    list_8TeV->Add(graphNPionComb8TeVStat);
                    list_8TeV->Add(graphNPionComb8TeVSys);
                }
            }
            if (enable_Eta){
                TGraphAsymmErrors* graphEtaXSecPCM8TeVStat        = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaPCM8TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCM8TeVSys         = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaPCM8TeVSysErr");
                TGraphAsymmErrors* graphEtaXSecEMCAL8TeVStat      = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaEMCAL8TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecEMCAL8TeVSys       = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaEMCAL8TeVSysErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMCAL8TeVStat   = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaPCMEMCAL8TeVStatErr");
                TGraphAsymmErrors* graphEtaXSecPCMEMCAL8TeVSys    = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaPCMEMCAL8TeVSysErr");
                TGraphAsymmErrors* graphEtaXSecComb8TeVStat       = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaComb8TeVAStatErr");
                TGraphAsymmErrors* graphEtaXSecComb8TeVSys        = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphInvCrossSectionEtaComb8TeVASysErr");

                TH1F* graphEtaToPi0PCM8TeVStat      			 = (TH1F*)fFileNeutralMeson8TeVCombEta->Get("histoRatioEtaToPi0PCM8TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0PCM8TeVSys       = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphRatioEtaToPi0PCM8TeVSysErr");
                TH1F* graphEtaToPi0EMCAL8TeVStat    		     = (TH1F*)fFileNeutralMeson8TeVCombEta->Get("histoRatioEtaToPi0EMCAL8TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0EMCAL8TeVSys     = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphRatioEtaToPi0EMCAL8TeVSysErr");
                TH1F* graphEtaToPi0PCMEMCAL8TeVStat 			 = (TH1F*)fFileNeutralMeson8TeVCombEta->Get("histoRatioEtaToPi0PCMEMCAL8TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0PCMEMCAL8TeVSys  = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphRatioEtaToPi0PCMEMCAL8TeVSysErr");
                TGraphAsymmErrors* graphEtaToPi0Comb8TeVStat     = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphRatioEtaToPi0Comb8TeVStatErr");
                TGraphAsymmErrors* graphEtaToPi0Comb8TeVSys      = (TGraphAsymmErrors*)fFileNeutralMeson8TeVCombEta->Get("graphRatioEtaToPi0Comb8TeVSysErr");

                SetHistoProperties(graphEtaToPi0PCM8TeVStat,     Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCM8TeVSys,      Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetHistoProperties(graphEtaToPi0EMCAL8TeVStat,   Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0EMCAL8TeVSys,    Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetHistoProperties(graphEtaToPi0PCMEMCAL8TeVStat,Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0PCMEMCAL8TeVSys, Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0Comb8TeVStat,    Form("%sTo%s%sStat", fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                SetGraphProperties(graphEtaToPi0Comb8TeVSys,     Form("%sTo%s%sSys",  fParticle[1].Data(), fParticle[0].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", "#eta/#pi^{0}", "");
                list_8TeV->Add(graphEtaToPi0PCM8TeVStat);
                list_8TeV->Add(graphEtaToPi0PCM8TeVSys);
                list_8TeV->Add(graphEtaToPi0EMCAL8TeVStat);
                list_8TeV->Add(graphEtaToPi0EMCAL8TeVSys);
                list_8TeV->Add(graphEtaToPi0PCMEMCAL8TeVStat);
                list_8TeV->Add(graphEtaToPi0PCMEMCAL8TeVSys);
                list_8TeV->Add(graphEtaToPi0Comb8TeVStat);
                list_8TeV->Add(graphEtaToPi0Comb8TeVSys);

                if (changeQuantity == 0){
                    // simply add Xsection to the output file
                    SetGraphProperties(graphEtaXSecPCM8TeVStat,     Form("Xsec_%s%sStat",    fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCM8TeVSys,      Form("Xsec_%s%sSys",     fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecEMCAL8TeVStat,   Form("Xsec_%s%sStat",    fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecEMCAL8TeVSys,    Form("Xsec_%s%sSys",     fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMEMCAL8TeVStat,Form("Xsec_%s%sStat",    fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecPCMEMCAL8TeVSys, Form("Xsec_%s%sSys",     fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecComb8TeVStat,    Form("Xsec_%s%sStat",    fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                    SetGraphProperties(graphEtaXSecComb8TeVSys,     Form("Xsec_%s%sSys",     fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                    list_8TeV->Add(graphEtaXSecPCM8TeVStat);
                    list_8TeV->Add(graphEtaXSecPCM8TeVSys);
                    list_8TeV->Add(graphEtaXSecEMCAL8TeVStat);
                    list_8TeV->Add(graphEtaXSecEMCAL8TeVSys);
                    list_8TeV->Add(graphEtaXSecPCMEMCAL8TeVStat);
                    list_8TeV->Add(graphEtaXSecPCMEMCAL8TeVSys);
                    list_8TeV->Add(graphEtaXSecComb8TeVStat);
                    list_8TeV->Add(graphEtaXSecComb8TeVSys);
                } else {
                    // convert XSection to invariant yield
                    TGraphAsymmErrors* graphEtaPCM8TeVStat        = ScaleGraph(graphEtaXSecPCM8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaPCM8TeVSys         = ScaleGraph(graphEtaXSecPCM8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaEMCAL8TeVStat      = ScaleGraph(graphEtaXSecEMCAL8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaEMCAL8TeVSys       = ScaleGraph(graphEtaXSecEMCAL8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaPCMEMCAL8TeVStat   = ScaleGraph(graphEtaXSecPCMEMCAL8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaPCMEMCAL8TeVSys    = ScaleGraph(graphEtaXSecPCMEMCAL8TeVSys,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaComb8TeVStat       = ScaleGraph(graphEtaXSecComb8TeVStat,1./xSection8TeVINEL);
                    TGraphAsymmErrors* graphEtaComb8TeVSys        = ScaleGraph(graphEtaXSecComb8TeVSys,1./xSection8TeVINEL);

                    // convert inv yield to 1/N d^2N/dydpT
                    graphEtaPCM8TeVStat                           = ConvertYieldGraph(graphEtaPCM8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCM8TeVSys                            = ConvertYieldGraph(graphEtaPCM8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaEMCAL8TeVStat                         = ConvertYieldGraph(graphEtaEMCAL8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaEMCAL8TeVSys                          = ConvertYieldGraph(graphEtaEMCAL8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCMEMCAL8TeVStat                      = ConvertYieldGraph(graphEtaPCMEMCAL8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaPCMEMCAL8TeVSys                       = ConvertYieldGraph(graphEtaPCMEMCAL8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaComb8TeVStat                          = ConvertYieldGraph(graphEtaComb8TeVStat,  kFALSE, kFALSE, kTRUE, kTRUE);
                    graphEtaComb8TeVSys                           = ConvertYieldGraph(graphEtaComb8TeVSys,   kFALSE, kFALSE, kTRUE, kTRUE);

                    SetGraphProperties(graphEtaPCM8TeVStat,       Form("%s%sStat", fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCM8TeVSys,        Form("%s%sSys",  fParticle[1].Data(), fMethod[2].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaEMCAL8TeVStat,     Form("%s%sStat", fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaEMCAL8TeVSys,      Form("%s%sSys",  fParticle[1].Data(), fMethod[4].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCMEMCAL8TeVStat,  Form("%s%sStat", fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaPCMEMCAL8TeVSys,   Form("%s%sSys",  fParticle[1].Data(), fMethod[6].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaComb8TeVStat,      Form("%s%sStat", fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                    SetGraphProperties(graphEtaComb8TeVSys,       Form("%s%sSys",  fParticle[1].Data(), fMethod[1].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                    list_8TeV->Add(graphEtaPCM8TeVStat);
                    list_8TeV->Add(graphEtaPCM8TeVSys);
                    list_8TeV->Add(graphEtaEMCAL8TeVStat);
                    list_8TeV->Add(graphEtaEMCAL8TeVSys);
                    list_8TeV->Add(graphEtaPCMEMCAL8TeVStat);
                    list_8TeV->Add(graphEtaPCMEMCAL8TeVSys);
                    list_8TeV->Add(graphEtaComb8TeVStat);
                    list_8TeV->Add(graphEtaComb8TeVSys);
                }
            }
        }
        if (!parametrizeMC && enable_CPion && enable_CKaon && enable_Proton && enable_Lambda){
            //================================================================================================================
            // reading and writing pi+-, K+- and p/bar{p} to 8TeV list
            // input is given as fully invariant yield
            //================================================================================================================
            TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

            TFile *fileExtrapolated8TeV                        		= new TFile("pp/8TeV_extrapolated_spectra.root");

            TGraphErrors* graphCPionInter8TeVStat              = (TGraphErrors*)fileExtrapolated8TeV->Get("graphStatErr_Pion_8TeV");
            TGraphErrors* graphCKaonInter8TeVStat              = (TGraphErrors*)fileExtrapolated8TeV->Get("graphStatErr_Kaon_8TeV");
            TGraphErrors* graphProtonInter8TeVStat             = (TGraphErrors*)fileExtrapolated8TeV->Get("graphStatErr_Proton_8TeV");
            TGraphErrors* graphLambdaInter8TeVStat             = (TGraphErrors*)fileExtrapolated8TeV->Get("graphStatErr_Lambda_8TeV");

            TGraphErrors* graphCPionInter8TeVSys               = (TGraphErrors*)fileExtrapolated8TeV->Get("graphSysErr_Pion_8TeV");
            TGraphErrors* graphCKaonInter8TeVSys               = (TGraphErrors*)fileExtrapolated8TeV->Get("graphSysErr_Kaon_8TeV");
            TGraphErrors* graphProtonInter8TeVSys              = (TGraphErrors*)fileExtrapolated8TeV->Get("graphSysErr_Proton_8TeV");
            TGraphErrors* graphLambdaInter8TeVSys              = (TGraphErrors*)fileExtrapolated8TeV->Get("graphSysErr_Lambda_8TeV");


            if (changeQuantity == 0 || changeQuantity == 2){
                SetGraphProperties(graphCPionInter8TeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphCPionInter8TeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphCKaonInter8TeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphCKaonInter8TeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphProtonInter8TeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphProtonInter8TeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphLambdaInter8TeVStat, Form("%s%sStat", fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
                SetGraphProperties(graphLambdaInter8TeVSys,  Form("%s%sSys",  fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

                list_8TeV->Add(graphCPionInter8TeVStat);
                list_8TeV->Add(graphCPionInter8TeVSys);
                list_8TeV->Add(graphCKaonInter8TeVStat);
                list_8TeV->Add(graphCKaonInter8TeVSys);
                list_8TeV->Add(graphProtonInter8TeVStat);
                list_8TeV->Add(graphProtonInter8TeVSys);
                list_8TeV->Add(graphLambdaInter8TeVStat);
                list_8TeV->Add(graphLambdaInter8TeVSys);
            } else { //TGraphErrors* inputGraph,                                  Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt
                graphCPionInter8TeVStat                              = ConvertYieldGraph(graphCPionInter8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphCPionInter8TeVSys                               = ConvertYieldGraph(graphCPionInter8TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                graphCKaonInter8TeVStat                              = ConvertYieldGraph(graphCKaonInter8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphCKaonInter8TeVSys                               = ConvertYieldGraph(graphCKaonInter8TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                graphProtonInter8TeVStat                             = ConvertYieldGraph(graphProtonInter8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
                graphProtonInter8TeVSys                              = ConvertYieldGraph(graphProtonInter8TeVSys,  kFALSE, kFALSE, kTRUE, kTRUE);
                graphLambdaInter8TeVStat                             = ConvertYieldGraph(graphLambdaInter8TeVStat,kFALSE, kFALSE, kTRUE, kTRUE);
                graphLambdaInter8TeVSys                              = ConvertYieldGraph(graphLambdaInter8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

                SetGraphProperties(graphCPionInter8TeVStat,  Form("%s%sStat", fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphCPionInter8TeVSys,   Form("%s%sSys",  fParticle[5].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphCKaonInter8TeVStat,  Form("%s%sStat", fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphCKaonInter8TeVSys,   Form("%s%sSys",  fParticle[6].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphProtonInter8TeVStat, Form("%s%sStat", fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphProtonInter8TeVSys,  Form("%s%sSys",  fParticle[7].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphLambdaInter8TeVStat, Form("%s%sStat", fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
                SetGraphProperties(graphLambdaInter8TeVSys,  Form("%s%sSys",  fParticle[16].Data(), fMethod[0].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

                list_8TeV->Add(graphCPionInter8TeVStat);
                list_8TeV->Add(graphCPionInter8TeVSys);
                list_8TeV->Add(graphCKaonInter8TeVStat);
                list_8TeV->Add(graphCKaonInter8TeVSys);
                list_8TeV->Add(graphProtonInter8TeVStat);
                list_8TeV->Add(graphProtonInter8TeVSys);
                list_8TeV->Add(graphLambdaInter8TeVStat);
                list_8TeV->Add(graphLambdaInter8TeVSys);
            }
        }

        //================================================================================================================
        // reading and writing MC spectra of pi0, eta, K0s, K0l, lambda to 8TeV list
        // input given as inv. yield
        //================================================================================================================
        if(parametrizeMC && enable_NPion && enable_Eta && enable_K0s && enable_CKaon && enable_Lambda ){

          TString localOutputQuantity                             = "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-2})";

          TFile* fFileMCParam8TeV                                 = new TFile("pp/MCInputCompilationLHC15h_pp8TeV_2.root");

          TH1F* histoNPionMC8TeVStat                 = (TH1F*)fFileMCParam8TeV->Get("MC_Pi0_Pt_Rebinned");
          TH1F* histoNPionMC8TeVSys                  = (TH1F*)fFileMCParam8TeV->Get("MC_Pi0_Pt_Rebinned");
          for(Int_t i=1; i<=histoNPionMC8TeVStat->GetNbinsX(); i++) histoNPionMC8TeVSys->SetBinError(i,histoNPionMC8TeVStat->GetBinContent(i)*0.05);
          TH1F* histoEtaMC8TeVStat                   = (TH1F*)fFileMCParam8TeV->Get("MC_Eta_Pt_Rebinned");
          TH1F* histoEtaMC8TeVSys                    = (TH1F*)fFileMCParam8TeV->Get("MC_Eta_Pt_Rebinned");
          for(Int_t i=1; i<=histoEtaMC8TeVStat->GetNbinsX(); i++) histoEtaMC8TeVSys->SetBinError(i,histoEtaMC8TeVStat->GetBinContent(i)*0.05);
          TH1F* histoK0sMC8TeVStat                   = (TH1F*)fFileMCParam8TeV->Get("MC_K0s_Pt_Rebinned");
          TH1F* histoK0sMC8TeVSys                    = (TH1F*)fFileMCParam8TeV->Get("MC_K0s_Pt_Rebinned");
          for(Int_t i=1; i<=histoK0sMC8TeVStat->GetNbinsX(); i++) histoK0sMC8TeVSys->SetBinError(i,histoK0sMC8TeVStat->GetBinContent(i)*0.05);
          TH1F* histoK0lMC8TeVStat                   = (TH1F*)fFileMCParam8TeV->Get("MC_K0l_Pt_Rebinned");
          TH1F* histoK0lMC8TeVSys                    = (TH1F*)fFileMCParam8TeV->Get("MC_K0l_Pt_Rebinned");
          for(Int_t i=1; i<=histoK0lMC8TeVStat->GetNbinsX(); i++) histoK0lMC8TeVSys->SetBinError(i,histoK0lMC8TeVStat->GetBinContent(i)*0.05);
          TH1F* histoLambdaMC8TeVStat                = (TH1F*)fFileMCParam8TeV->Get("MC_Lambda_Pt_Rebinned");
          TH1F* histoLambdaMC8TeVSys                 = (TH1F*)fFileMCParam8TeV->Get("MC_Lambda_Pt_Rebinned");
          for(Int_t i=1; i<=histoLambdaMC8TeVStat->GetNbinsX(); i++) histoLambdaMC8TeVSys->SetBinError(i,histoLambdaMC8TeVStat->GetBinContent(i)*0.05);

          if (changeQuantity == 0 || changeQuantity == 2){
            SetHistoProperties(histoNPionMC8TeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoNPionMC8TeVSys,  Form("%s%sSys", fParticle[0].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoEtaMC8TeVStat,   Form("%s%sStat", fParticle[1].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoEtaMC8TeVSys,    Form("%s%sSys", fParticle[1].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoK0sMC8TeVStat,   Form("%s%sStat", fParticle[15].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoK0sMC8TeVSys,    Form("%s%sSys", fParticle[15].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoK0lMC8TeVStat,   Form("%s%sStat", fParticle[6].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoK0lMC8TeVSys,    Form("%s%sSys", fParticle[6].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoLambdaMC8TeVStat,Form("%s%sStat", fParticle[16].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");
            SetHistoProperties(histoLambdaMC8TeVSys ,Form("%s%sSys", fParticle[16].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", localOutputQuantity, "");

            list_8TeV->Add(histoNPionMC8TeVStat);
            list_8TeV->Add(histoNPionMC8TeVSys);
            list_8TeV->Add(histoEtaMC8TeVStat);
            list_8TeV->Add(histoEtaMC8TeVSys);
            list_8TeV->Add(histoK0sMC8TeVStat);
            list_8TeV->Add(histoK0sMC8TeVSys);
            list_8TeV->Add(histoK0lMC8TeVStat);
            list_8TeV->Add(histoK0lMC8TeVSys);
            list_8TeV->Add(histoLambdaMC8TeVStat);
            list_8TeV->Add(histoLambdaMC8TeVSys);
          }else{
            histoNPionMC8TeVStat                               = ConvertYieldHisto(histoNPionMC8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            histoNPionMC8TeVSys                                = ConvertYieldHisto(histoNPionMC8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
            histoEtaMC8TeVStat                                 = ConvertYieldHisto(histoEtaMC8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            histoEtaMC8TeVSys                                  = ConvertYieldHisto(histoEtaMC8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
            histoK0sMC8TeVStat                                 = ConvertYieldHisto(histoK0sMC8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            histoK0sMC8TeVSys                                  = ConvertYieldHisto(histoK0sMC8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
            histoK0lMC8TeVStat                                 = ConvertYieldHisto(histoK0lMC8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            histoK0lMC8TeVSys                                  = ConvertYieldHisto(histoK0lMC8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);
            histoLambdaMC8TeVStat                              = ConvertYieldHisto(histoLambdaMC8TeVStat, kFALSE, kFALSE, kTRUE, kTRUE);
            histoLambdaMC8TeVSys                               = ConvertYieldHisto(histoLambdaMC8TeVSys, kFALSE, kFALSE, kTRUE, kTRUE);

            SetHistoProperties(histoNPionMC8TeVStat, Form("%s%sStat", fParticle[0].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoNPionMC8TeVSys,  Form("%s%sSys", fParticle[0].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoEtaMC8TeVStat,   Form("%s%sStat", fParticle[1].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoEtaMC8TeVSys,    Form("%s%sSys", fParticle[1].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoK0sMC8TeVStat,   Form("%s%sStat", fParticle[15].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoK0sMC8TeVSys,    Form("%s%sSys", fParticle[15].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoK0lMC8TeVStat,   Form("%s%sStat", fParticle[6].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoK0lMC8TeVSys,    Form("%s%sSys", fParticle[6].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoLambdaMC8TeVStat,Form("%s%sStat", fParticle[16].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");
            SetHistoProperties(histoLambdaMC8TeVSys ,Form("%s%sSys", fParticle[16].Data(), fMethod[12].Data()), "#it{p}_{T} (GeV/#it{c})", globalOutputQuantity, "");

            list_8TeV->Add(histoNPionMC8TeVStat);
            list_8TeV->Add(histoNPionMC8TeVSys);
            list_8TeV->Add(histoEtaMC8TeVStat);
            list_8TeV->Add(histoEtaMC8TeVSys);
            list_8TeV->Add(histoK0sMC8TeVStat);
            list_8TeV->Add(histoK0sMC8TeVSys);
            list_8TeV->Add(histoK0lMC8TeVStat);
            list_8TeV->Add(histoK0lMC8TeVSys);
            list_8TeV->Add(histoLambdaMC8TeVStat);
            list_8TeV->Add(histoLambdaMC8TeVSys);
          }
        }
    }

    //================================================================================================================
    //Produce plots containing all particle spectra if specified
    //================================================================================================================
    if (Include_900GeV)     ProduceParticleSpectraPlotFromList(list_900GeV,  fCollSys[0], fEnergy[0], "", suffix);
    if (Include_2760GeV)    ProduceParticleSpectraPlotFromList(list_2760GeV, fCollSys[0], fEnergy[1], "", suffix);
    if (Include_5TeV)       ProduceParticleSpectraPlotFromList(list_5TeV,    fCollSys[0], fEnergy[2], "", suffix);
    if (Include_7TeV)       ProduceParticleSpectraPlotFromList(list_7TeV,    fCollSys[0], fEnergy[3], "", suffix);
    if (Include_8TeV)       ProduceParticleSpectraPlotFromList(list_8TeV,    fCollSys[0], fEnergy[4], "", suffix);
    if (Include_13TeV)      ProduceParticleSpectraPlotFromList(list_13TeV,   fCollSys[0], fEnergy[5], "", suffix);


    //================================================================================================================
    //Produce plots containing all particle spectra if specified
    //================================================================================================================
    if (Include_900GeV)     ProduceParticleSpectraPlotFromListOnlyFinal(list_900GeV,  fCollSys[0], fEnergy[0], "", suffix);
    if (Include_2760GeV)    ProduceParticleSpectraPlotFromListOnlyFinal(list_2760GeV, fCollSys[0], fEnergy[1], "", suffix);
    if (Include_5TeV)       ProduceParticleSpectraPlotFromListOnlyFinal(list_5TeV,    fCollSys[0], fEnergy[2], "", suffix);
    if (Include_7TeV)       ProduceParticleSpectraPlotFromListOnlyFinal(list_7TeV,    fCollSys[0], fEnergy[3], "", suffix);
    if (Include_8TeV)       ProduceParticleSpectraPlotFromListOnlyFinal(list_8TeV,    fCollSys[0], fEnergy[4], "", suffix);
    if (Include_13TeV)      ProduceParticleSpectraPlotFromListOnlyFinal(list_13TeV,   fCollSys[0], fEnergy[5], "", suffix);

    //================================================================================================================
    //Produce plots containing all particle spectra if specified
    //================================================================================================================
    if (Include_900GeV)     ProduceParticleRatioPlotFromListOnlyFinal(list_900GeV,  fCollSys[0], fEnergy[0], "", suffix);
    if (Include_2760GeV)    ProduceParticleRatioPlotFromListOnlyFinal(list_2760GeV, fCollSys[0], fEnergy[1], "", suffix);
    if (Include_5TeV)       ProduceParticleRatioPlotFromListOnlyFinal(list_5TeV,    fCollSys[0], fEnergy[2], "", suffix);
    if (Include_7TeV)       ProduceParticleRatioPlotFromListOnlyFinal(list_7TeV,    fCollSys[0], fEnergy[3], "", suffix, kTRUE);
    if (Include_8TeV)       ProduceParticleRatioPlotFromListOnlyFinal(list_8TeV,    fCollSys[0], fEnergy[4], "", suffix);
    if (Include_13TeV)      ProduceParticleRatioPlotFromListOnlyFinal(list_13TeV,   fCollSys[0], fEnergy[5], "", suffix);

    //================================================================================================================
    //Produce plots containing all particle ratios if specified
    //================================================================================================================
    if (Include_900GeV)     ProduceParticleRatioPlotFromList(list_900GeV,  fCollSys[0], fEnergy[0], "", suffix);
    if (Include_2760GeV)    ProduceParticleRatioPlotFromList(list_2760GeV, fCollSys[0], fEnergy[1], "", suffix);
    if (Include_5TeV)       ProduceParticleRatioPlotFromList(list_5TeV,    fCollSys[0], fEnergy[2], "", suffix);
    if (Include_7TeV)       ProduceParticleRatioPlotFromList(list_7TeV,    fCollSys[0], fEnergy[3], "", suffix);
    if (Include_8TeV)       ProduceParticleRatioPlotFromList(list_8TeV,    fCollSys[0], fEnergy[4], "", suffix);
    if (Include_13TeV)      ProduceParticleRatioPlotFromList(list_13TeV,   fCollSys[0], fEnergy[5], "", suffix);

    output_File->cd();
    //================================================================================================================
    //Saving the TList to the final file
    //================================================================================================================
    if(Include_900GeV  && list_900GeV->GetEntries())  list_900GeV->Write( "pp_0.9TeV",  TObject::kSingleKey);
    if(Include_2760GeV && list_2760GeV->GetEntries()) list_2760GeV->Write("pp_2.76TeV", TObject::kSingleKey);
    if(Include_5TeV    && list_5TeV->GetEntries())    list_5TeV->Write(   "pp_5TeV",    TObject::kSingleKey);
    if(Include_7TeV    && list_7TeV->GetEntries())    list_7TeV->Write(   "pp_7TeV",    TObject::kSingleKey);
    if(Include_8TeV    && list_8TeV->GetEntries())    list_8TeV->Write(   "pp_8TeV",    TObject::kSingleKey);
    if(Include_13TeV   && list_13TeV->GetEntries())   list_13TeV->Write(  "pp_13TeV",   TObject::kSingleKey);

    output_File->Close();
}


