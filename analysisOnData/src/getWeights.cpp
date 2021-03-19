#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/getWeights.hpp"


RNode getWeights::run(RNode d){

    auto getNorm = [](float Wpt, float Wrap, const ROOT::VecOps::RVec<float> &AngCoeff, ROOT::VecOps::RVec<float> harmonicsVec, float totMap){

        float norm = harmonicsVec[8]*totMap;
        
        harmonicsVec[0]/=2.;
        harmonicsVec[1]/=(2.*TMath::Sqrt(2));
        harmonicsVec[2]/=4.;
        harmonicsVec[3]/=(4.*TMath::Sqrt(2));
        harmonicsVec[4]/=2.;
        harmonicsVec[5]/=2.;
        harmonicsVec[6]/=(2.*TMath::Sqrt(2));
        harmonicsVec[7]/=(4.*TMath::Sqrt(2));

        for(unsigned int i=0; i<harmonicsVec.size()-1; i++){ //loop over angular coefficients

             norm+=(AngCoeff[i]*harmonicsVec[i]*totMap); //sum Ai*Pi

        }

        float fact = 3./(16.*TMath::Pi());
        norm*=fact;

        return norm;

    };

    auto getWeights = [](float norm, const ROOT::VecOps::RVec<float>& harmonicsVec){

        return (harmonicsVec/norm);

    };
    
    
    auto getNormVars = [this](float Wpt, float Wrap, const ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> &AngCoeff, ROOT::VecOps::RVec<float> harmonicsVec, ROOT::VecOps::RVec<float> totMap){
        
        harmonicsVec[0]/=2.;
        harmonicsVec[1]/=(2.*TMath::Sqrt(2));
        harmonicsVec[2]/=4.;
        harmonicsVec[3]/=(4.*TMath::Sqrt(2));
        harmonicsVec[4]/=2.;
        harmonicsVec[5]/=2.;
        harmonicsVec[6]/=(2.*TMath::Sqrt(2));
        harmonicsVec[7]/=(4.*TMath::Sqrt(2));
        float fact = 3./(16.*TMath::Pi());
        
         ROOT::VecOps::RVec<float> norm;

        for(long unsigned int v=0; v!=_syst_name.size(); v++) {
            float norm_var = harmonicsVec[8]*totMap[v];
        
            for(unsigned int i=0; i<harmonicsVec.size()-1; i++){ //loop over angular coefficients
                norm_var+=(AngCoeff[v][i]*harmonicsVec[i]*totMap[v]); //sum Ai*Pi
            }
            norm_var*=fact;
            norm.emplace_back(norm_var);
        }        
        return norm;

    };
    
    auto getWeightsVars = [this](ROOT::VecOps::RVec<float> &norm, const ROOT::VecOps::RVec<float>& harmonicsVec){
        
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> harmonicsWeights;
        
        for(long unsigned int v=0; v!=_syst_name.size(); v++) {
            harmonicsWeights.emplace_back(harmonicsVec/norm[v]);
        }
        
        return harmonicsWeights;

    };
    
    
    if(_syst_kind==""){
        auto d1 = d.Define("norm", getNorm, {"Wpt_preFSR", "Wrap_preFSR_abs", "AngCoeffVec", "harmonicsVec", "totMap"})
                  .Define("harmonicsWeights", getWeights, {"norm", "harmonicsVec"});
        return d1;
    }
    else {
        auto d1 = d.Define("norm", getNormVars, {"Wpt_preFSR", "Wrap_preFSR_abs", "AngCoeffVec", "harmonicsVec", "totMap"})
                  .Define("harmonicsWeights", getWeightsVars, {"norm", "harmonicsVec"});
        return d1;
    }
    
}
