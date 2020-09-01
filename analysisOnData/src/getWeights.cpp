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
    auto d1 = d.Define("norm", getNorm, {"GenV_preFSR_qt", "GenV_preFSR_yabs", "AngCoeffVec", "harmonicsVec", "totMap"})
                  .Define("harmonicsWeights", getWeights, {"norm", "harmonicsVec"});

    return d1;
    
}