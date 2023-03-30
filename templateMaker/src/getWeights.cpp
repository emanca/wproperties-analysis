#include "interface/getWeights.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

RNode getWeights::run(RNode d)
{

    auto getNorm = [](float Wpt, float Wrap, const ROOT::VecOps::RVec<float> &AngCoeff, const ROOT::VecOps::RVec<float> &harmVec)
    {
        float norm = 0.;
        float totMap = AngCoeff[8];
        norm += harmVec[0] * AngCoeff[0] * totMap;
        norm += harmVec[1] * AngCoeff[1] * totMap;
        norm += harmVec[2] * AngCoeff[2] * totMap;
        norm += harmVec[3] * AngCoeff[3] * totMap;
        norm += harmVec[4] * AngCoeff[4] * totMap;
        norm += harmVec[5] * AngCoeff[5] * totMap;
        norm += harmVec[6] * AngCoeff[6] * totMap;
        norm += harmVec[7] * AngCoeff[7] * totMap;
        norm += harmVec[8] * totMap;

        float fact = 3. / (16. * TMath::Pi());
        norm *= fact;
        return norm;
    };

    auto d1 = d.Define("norm", getNorm, {"Vpt_preFSR", "Vrap_preFSR_abs", "helVec", "harmonicsVec"})
                  .Define("weightL", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[0] * helVec[0]/norm)")
                  .Define("weightI", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[1] * helVec[1]/norm)")
                  .Define("weightT", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[2] * helVec[2]/norm)")
                  .Define("weightA", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[3] * helVec[3]/norm)")
                  .Define("weightP", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[4] * helVec[4]/norm)")
                  .Define("weightUL", "float(3. / (16. * TMath::Pi()) * helVec[8] * harmonicsVec[8]/norm)");

    auto getTensor = [](float weightL, float weightI, float weightT, float weightA, float weightP, float weightUL, double or_weight)
    {
        Eigen::TensorFixedSize<double, Eigen::Sizes<6>> helTensor;
        ROOT::VecOps::RVec<double> weight_array = {weightL, weightI, weightT, weightA, weightP, weightUL};
        auto w = weight_array * or_weight;
        std::copy(std::begin(w), std::end(w), helTensor.data());

        return helTensor;
    };

    auto d2 = d1.Define("helWeightTensor", getTensor, {"weightL", "weightI", "weightT", "weightA", "weightP", "weightUL", "nominal_weight"});
    return d2;
}
