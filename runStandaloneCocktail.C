// runStandaloneCocktail.C
// =====================
// This macro runs cocktail generation with AliMCGen framework
//
// orginal author: M. Verweij
// modified by: F.Bock, L.Altenkaemper, N.Schmidt

#include <ctime>
#include "TGrid.h"

const Bool_t saveManager = kFALSE;

void runStandaloneCocktail( Long64_t nEvents        = 10,
                            Int_t collisionsSystem  = 0,            // 0: pp 7TeV, 1: pPb 5TeV
                            Int_t motherSelect      = 67108863,     // includes all particles in generator (see end of macro for further examples)
                            Int_t decayMode         = 0             // 0: kAll, 1: kGammaEM, 2: kElectronEM, 3: kDiElectronEM
){

    AliAnalysisManager* mgr         = new AliAnalysisManager("MCGenEMCocktail");
    AliDummyHandler*    dumH        = new AliDummyHandler();

    // create blank ESD event
    AliESDEvent *esdE               = new AliESDEvent();
    esdE->CreateStdContent();
    AliESDVertex *vtx               = new AliESDVertex(0.,0.,100);
    vtx->SetName("VertexTracks");
    vtx->SetTitle("VertexTracks");
    esdE->SetPrimaryVertexTracks(vtx);
    if(esdE->GetPrimaryVertex()) Printf("vtx set");
    dumH->SetEvent(esdE);
    mgr->SetInputEventHandler(dumH);

    // greate MC input Handler
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    mgr->SetMCtruthEventHandler(mcInputHandler);

    // cocktail generator
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Cocktail/macros/AddMCEMCocktailV2.C");
    AliGenerator* gener             = 0x0;
    if (collisionsSystem == 0) {
        gener = AddMCEMCocktailV2(200,0,decayMode,motherSelect,"parametrizations/pp.root","7TeV_Comb",100,0.,50,2000,0,1,0,0);
    } else if (collisionsSystem == 1) {
        gener = AddMCEMCocktailV2(300,0,decayMode,motherSelect,"parametrizations/pPb.root","5TeV_MB_PCM",1000,0.,50,2000,0,1,0,0);
    } else {
        cout << "collision system not recognized" << endl;
        return;
    }

    mcInputHandler->SetGenerator(gener);
    mcInputHandler->SetSeedMode(2);

    // gamma cocktail consumer task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaCocktailMC.C");
    AliAnalysisTask *taskA = AddTask_GammaCocktailMC(kFALSE, "80");

    // hadronic cocktail consumer task
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_HadronicCocktailMC.C");
    AliAnalysisTask *taskB = AddTask_HadronicCocktailMC(0, kFALSE, "80"); // pi0
    AliAnalysisTask *taskC = AddTask_HadronicCocktailMC(1, kFALSE, "80"); // eta
    AliAnalysisTask *taskC = AddTask_HadronicCocktailMC(2, kFALSE, "80"); // pi+-

    // give some indication everything is running fine
    mgr->SetUseProgressBar(1, 10);

    // Run analysis
    if(saveManager) {
        TFile *fileMgr = new TFile("AnalysisManager.root","RECREATE");
        mgr->Write();

        fileMgr->Write();
        fileMgr->Close();

        runUseSavedManager("AnalysisManager.root",nEvents);
    } else {
        mgr->InitAnalysis();
        mgr->PrintStatus();
        mgr->EventLoop(nEvents);
    }

}

void runUseSavedManager(TString fileMgr = "AnalysisManager.root", Int_t nEvents = 1000, Bool_t useGrid = kFALSE) {

    Printf("run saved");
    LoadLibs();

    TFile *f = new TFile(fileMgr.Data());
    AliAnalysisManager *mgr = dynamic_cast<AliAnalysisManager*>f->Get("MCGenThermalModel");
    if(!mgr) {
        Printf("Did not find manager");
        return;
    }
    mgr->SetDebugLevel(11);
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->EventLoop(nEvents);
}

//================= examples for inclusion of mother particles ================================
// (62591)_10 = (1111010001111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, J/psi, Delta+, Delta0, Rho+, Rho-, K0*

// (29887)_10 = (0111010010111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Sigma0, Delta+, Delta0, Rho+, Rho-, K0*

// (262143)_10 = (111111111111111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Jpsi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0*, k0l, lambda

// (229375)_10 = (110111111111111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, Jpsi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0l, lambda

// (262079)_10 = (111111111110111111)_2
// includes: pi0, eta, rho0, omega, etaprime, phi, sigma0, k0s, delta++, delta+, delta-, delta0, rho+, rho-, k0*, k0l, lambda

// (131072)_10 = (100000000000000000)_2
// includes: lambda

// (65536)_10 = (010000000000000000)_2
// includes: k0l

// (256)_10 = (000000000100000000)_2
// includes: k0s

// (257)_10 = (000000000100000001)_2
// includes: k0s, pi0

// (196864)_10 = (110000000100000000)_2
// includes: k0s, k0l, lambda

// (196875)_10 = (110000000100001011)_2
// includes: pi0, eta, omega, k0s, k0l, lambda

// (67108863)_10 = (11111111111111111111111111)_2
// includes: pi0, eta, rhoß, imega, eta', phi, JPsi, Sigma0, K0s, Delta++, Delta+, Delta-, Delta0, rho+, rho-, K0*, K0l, Lambda, K+, K-, Omega+, Omega-, Xi+, Xi-, Sigma+, Sigma-
//=============================================================================================

