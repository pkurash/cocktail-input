if [ $1 = "produce" ]; then
    root -l -b -q -x 'CocktailProduceCompleteInputFilePP.C+("100110",1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,"dN/dydpT","pdf")'
fi
if [ $1 = "produceMCParam" ]; then
    root -l -b -q -x 'CocktailProduceCompleteInputFilePP.C+("000010",1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,"dN/dydpT","eps",kTRUE)'
fi
if [ $1 = "parameter900" ]; then
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_0.9TeV","pdf","parametrizationSettings/pp_0.9TeV_standard.dat","","PCM")'
fi
if [ $1 = "parameter7" ]; then
    #root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","","PCM")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","PCMnoOmega")'
    #root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","PCMEtaToPi0")'
    #root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","","EMCal")'
    #root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","parametrizationSettings/pp_7TeV_ratio_standard.dat","EMCalEtaToPi0")'
    #root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_7TeV","pdf","parametrizationSettings/pp_7TeV_standard.dat","","PCMEMCal")'
fi
if [ $1 = "parameter8" ]; then
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","","Comb")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","CombEtaToPi0")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","","PCM")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","PCMEtaToPi0")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","","EMCAL")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","EMCALEtaToPi0")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","","PCMEMCAL")'
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","pdf","parametrizationSettings/pp_8TeV_standard.dat","parametrizationSettings/pp_8TeV_ratio_standard.dat","PCMEMCALEtaToPi0")'
fi

if [ $1 = "parameter8MC" ]; then
    root -x -l -b -q 'CocktailInputParametrization.C+("CocktailInputPP.root","pp_8TeV","eps","parametrizationSettings/pp_8TeV_standardMC.dat","","MCParam")'
fi

if [ $1 = "parameterIntYield7" ]; then
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","pdf","integratedyieldParametrizations/pp_7TeV_Comb_Fit_oHagPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","pdf","integratedyieldParametrizations/pp_7TeV_Comb_Fit_oHagPt_2.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","pdf","integratedyieldParametrizations/pp_7TeV_Comb_Fit_tmpt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","7TeV","pdf","integratedyieldParametrizations/pp_7TeV_Comb_Fit_tcmpt.dat")'
fi
if [ $1 = "parameterIntYield8" ]; then
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","pdf","integratedyieldParametrizations/pp_8TeV_Comb_Fit_oHagPt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","pdf","integratedyieldParametrizations/pp_8TeV_Comb_Fit_tmpt.dat")'
    root -x -l -b -q 'CalculateIntegratedYield.C+("CocktailInputPP.root","8TeV","pdf","integratedyieldParametrizations/pp_8TeV_Comb_Fit_tcmpt.dat")'
fi

if [ $1 = "plotIntYield7" ]; then
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_7TeV.root","7TeV","","pdf","paramSettingsIntegYield/pp_7TeV_YieldComp.dat")'
fi
if [ $1 = "plotIntYield8" ]; then
    root -x -l -b -q 'PlotIntegratedYieldsAndMeanPt.C+("CocktailInputPP_IntegYield_8TeV.root","8TeV","","pdf","paramSettingsIntegYield/pp_8TeV_YieldComp.dat")'
fi
