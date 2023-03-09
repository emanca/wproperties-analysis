#include "interface/getWeights.hpp"

RNode getWeights::run(RNode d)
{

    auto getNorm = [](float Wpt, float Wrap, const ROOT::VecOps::RVec<float> &AngCoeff, float P0, float P1, float P2, float P3, float P4, float P5, float P6, float P7, float PUL, float totMap) {
        
        float norm = 0.;
        norm += P0 * AngCoeff[0] * totMap;
        norm += P1 * AngCoeff[1] * totMap;
        norm += P2 * AngCoeff[2] * totMap;
        norm += P3 * AngCoeff[3] * totMap;
        norm += P4 * AngCoeff[4] * totMap;
        norm += P5 * AngCoeff[5] * totMap;
        norm += P6 * AngCoeff[6] * totMap;
        norm += P7 * AngCoeff[7] * totMap;
        norm += PUL * totMap;

        float fact = 3. / (16. * TMath::Pi());
        norm *= fact;
        return norm;
    };

    auto d1 = d.Define("norm", getNorm, {"Wpt_preFSR", "Wrap_preFSR_abs", "AngCoeffVec", "P0", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "PUL", "totMap"})
                  .Define("weightL", "float(3. / (16. * TMath::Pi()) * totMap * P0 * AngCoeffVec[0]/norm)")
                  .Define("weightI", "float(3. / (16. * TMath::Pi()) * totMap * P1 * AngCoeffVec[1]/norm)")
                  .Define("weightT", "float(3. / (16. * TMath::Pi()) * totMap * P2 * AngCoeffVec[2]/norm)")
                  .Define("weightA", "float(3. / (16. * TMath::Pi()) * totMap * P3 * AngCoeffVec[3]/norm)")
                  .Define("weightP", "float(3. / (16. * TMath::Pi()) * totMap * P4 * AngCoeffVec[4]/norm)")
                  .Define("weightUL", "float(3. / (16. * TMath::Pi()) * totMap * PUL /norm)");

    return d1;
}
