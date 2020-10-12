#include "interface/weightDefinitions.hpp"
#include "interface/functions.hpp"
RNode weightDefinitions::run(RNode d)
{

    auto definePUweights = [](float weightUp, float weightDown) {
        ROOT::VecOps::RVec<float> PUVars;
        PUVars.emplace_back(weightUp);
        PUVars.emplace_back(weightDown);
        return PUVars;
    };

    auto defineWHSF = [this](float pt, float eta, float charge) {
        int binReco = _Reco->FindBin(eta, pt);
        int binTrigger = _TriggerPlus->FindBin(eta, pt);
        if (charge > 0)
        {
            return _Reco->GetBinContent(binReco) * _TriggerPlus->GetBinContent(binTrigger);
        }
        else
        {
            return _Reco->GetBinContent(binReco) * _TriggerMinus->GetBinContent(binTrigger);
        }
    };

    auto h = new TH1F("h", "h", 48., -2.4, 2.4);
    // Define SF: WHSF = Trigger * ISO * ID
    auto defineWHSFVars = [this, h](float pt, float eta, float charge) {
        ROOT::VecOps::RVec<float> WHSF1;
        ROOT::VecOps::RVec<float> WHSF2;
        ROOT::VecOps::RVec<float> WHSF3;
        WHSF1.resize(96); //48 eta bins x 2 charges
        WHSF2.resize(96); //48 eta bins x 2 charges
        WHSF3.resize(96); //48 eta bins x 2 charges

        int binReco = _Reco->FindBin(eta, pt);
        int binTrigger = _TriggerPlus->FindBin(eta, pt);
        int binSyst = _TriggerPlusSyst0->FindBin(eta, pt);
        int binEta = h->FindBin(eta);
        float nomSF = 0.;
        charge > 0 ? nomSF = _Reco->GetBinContent(binReco) * _TriggerPlus->GetBinContent(binTrigger) : _Reco->GetBinContent(binReco) * _TriggerMinus->GetBinContent(binTrigger);

        for (int i = 0; i < 48; i++)
        {
            //condition ? result_if_true : result_if_false
            WHSF1[i] = nomSF;
            WHSF2[i] = nomSF;
            WHSF3[i] = nomSF;

            if (binEta == i + 1) // apply variation only if in the matching SF
            {
                if (i % 2 == 0)
                { // if even, up variation
                    if (charge > 0)
                    {
                        WHSF1[i] *= (1 + TMath::Sqrt(2) * _TriggerPlusSyst0->GetBinContent(binSyst)); //up
                        WHSF2[i] *= (1 + TMath::Sqrt(2) * _TriggerPlusSyst1->GetBinContent(binSyst)); //up
                        WHSF3[i] *= (1 + TMath::Sqrt(2) * _TriggerPlusSyst2->GetBinContent(binSyst)); //up
                    }
                    else
                    {
                        WHSF1[i + 48] *= (1 + TMath::Sqrt(2) * _TriggerMinusSyst0->GetBinContent(binSyst)); //up
                        WHSF2[i + 48] *= (1 + TMath::Sqrt(2) * _TriggerMinusSyst1->GetBinContent(binSyst)); //up
                        WHSF3[i + 48] *= (1 + TMath::Sqrt(2) * _TriggerMinusSyst2->GetBinContent(binSyst)); //up
                    }
                }
                else
                {
                    if (charge > 0)
                    {
                        WHSF1[i] *= (1 - TMath::Sqrt(2) * _TriggerPlusSyst0->GetBinContent(binSyst)); //down
                        WHSF2[i] *= (1 - TMath::Sqrt(2) * _TriggerPlusSyst1->GetBinContent(binSyst)); //down
                        WHSF3[i] *= (1 - TMath::Sqrt(2) * _TriggerPlusSyst2->GetBinContent(binSyst)); //down
                    }
                    else
                    {
                        WHSF1[i + 48] *= (1 - TMath::Sqrt(2) * _TriggerMinusSyst0->GetBinContent(binSyst)); //down
                        WHSF2[i + 48] *= (1 - TMath::Sqrt(2) * _TriggerMinusSyst1->GetBinContent(binSyst)); //down
                        WHSF3[i + 48] *= (1 - TMath::Sqrt(2) * _TriggerMinusSyst2->GetBinContent(binSyst)); //down
                    }
                }
            }
        }

        std::vector<float> WHSF;
        WHSF.insert(WHSF.end(), WHSF1.begin(), WHSF1.end());
        WHSF.insert(WHSF.end(), WHSF2.begin(), WHSF2.end());
        WHSF.insert(WHSF.end(), WHSF3.begin(), WHSF3.end());

        return WHSF;
    };

    auto defineWHSFVarsFlat = [this](float pt, float eta, float charge) {
        ROOT::VecOps::RVec<float> WHSF;
        WHSF.resize(6);
        int binReco = _Reco->FindBin(eta, pt);
        int binTrigger = _TriggerPlus->FindBin(eta, pt);

        float nomSF = 0.;
        charge > 0 ? nomSF = _Reco->GetBinContent(binReco) * _TriggerPlus->GetBinContent(binTrigger) : _Reco->GetBinContent(binReco) * _TriggerMinus->GetBinContent(binTrigger);

        for (int i = 0; i < 6; i++)
            WHSF[i] = nomSF;

        float flatVar = 0;
        if (fabs(eta) < 1)
        {
            flatVar = 0.002;
            WHSF[0] *= (1 + flatVar);
            WHSF[1] *= (1 - flatVar);
        }
        else if (fabs(eta) < 1.5)
        {
            flatVar = 0.004;
            WHSF[2] *= (1 + flatVar);
            WHSF[3] *= (1 - flatVar);
        }
        else
        {
            flatVar = 0.014;
            WHSF[4] *= (1 + flatVar);
            WHSF[5] *= (1 - flatVar);
        }

        return WHSF;
    };

    auto d1 = d.Define("WHSF", defineWHSF, {"Mu1_pt", "Mu1_eta", "Mu1_charge"})
    .Define("WHSFVarsUncorr", defineWHSFVars, {"Mu1_pt", "Mu1_eta", "Mu1_charge"})
    .Define("WHSFVarsFlat", defineWHSFVarsFlat, {"Mu1_pt", "Mu1_eta", "Mu1_charge"})
    .Define("puWeightVars", definePUweights, {"puWeightUp", "puWeightDown"})
    .Alias("PrefireWeightUp", "PrefireWeight_Up")
    .Alias("PrefireWeightDown", "PrefireWeight_Down")
    .Define("PrefireWeightVars", definePUweights, {"PrefireWeightUp", "PrefireWeightDown"}); //same function can be used since only 2 vars like PU

    return d1;
}
