#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/TH2weightsHelper.hpp"
#include "interface/TH1weightsHelper.hpp"
#include "interface/AngCoeff.hpp"

std::vector<std::string> AngCoeff::stringMultiplication(const std::vector<std::string> &v1, const std::vector<std::string> &v2)
{

    std::vector<std::string> products;

    if (v1.size() == 0)
        return v2;

    else
    {

        products.reserve(v1.size() * v2.size());
        for (auto e1 : v1)
        {
            for (auto e2 : v2)
            {
                products.push_back(e2 + e1);
            }
        }

        return products;
    }
}

RNode AngCoeff::run(RNode d)
{
    // for (float i=0.; i<=200.1; i=i+0.1){ _ptArr.push_back(i); }
    // int _nBinsPt = _ptArr.size()-1;
    // for (float i=0.; i<=6.1; i=i+0.1){ _yArr.push_back(i); }
    // int _nBinsY = _yArr.size()-1;
    
    std::vector<std::string> total = stringMultiplication(_syst_name, _coeff);

    auto vecMultiplication = [](const ROOT::VecOps::RVec<float> &v1, const ROOT::VecOps::RVec<float> &v2) {
        ROOT::VecOps::RVec<float> products;

        products.reserve(v1.size() * v2.size());
        for (auto e1 : v1)
            for (auto e2 : v2)
                products.push_back(e1 * e2);

        return products;
    };
    // systematic variations
    if (_syst_name.size() > 0)
    {
        auto d1 = d.Filter(_filter).Define(Form("%sharmonicsVec", _syst_weight.c_str()), vecMultiplication, {_syst_weight, "harmonicsVec"})
                   .Define(Form("%sharmonicsVecSq", _syst_weight.c_str()), vecMultiplication, {_syst_weight, "harmonicsVecSq"});

        TH2weightsHelper helper(std::string("harmonics"), std::string("harmonics"), _nBinsY, _yArr, _nBinsPt, _ptArr, total);
        auto helXsecs = d1.Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), { "Wrap_preFSR_abs", "Wpt_preFSR", "weight", Form("%sharmonicsVec", _syst_weight.c_str()) });
        _h2Group.push_back(helXsecs);

        TH2weightsHelper helperSq(std::string("harmonicsSq"), std::string("harmonicsSq"), _nBinsY, _yArr, _nBinsPt, _ptArr, total);
        auto helXsecsSq = d1.Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperSq), {"Wrap_preFSR_abs", "Wpt_preFSR", "weight", Form("%sharmonicsVecSq", _syst_weight.c_str())});
        _h2Group.push_back(helXsecsSq);

        TH2weightsHelper mapTothelper(std::string("mapTot"), std::string("mapTot"), _nBinsY, _yArr, _nBinsPt, _ptArr, _syst_name);
        auto mapTot = d1.Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(mapTothelper), {"Wrap_preFSR_abs", "Wpt_preFSR", "weight", _syst_weight });
        _h2Group.push_back(mapTot);

        TH1weightsHelper helperPt(std::string("harmonicsPt"), std::string("harmonicsPt"), _nBinsPt, _ptArr, total);
        auto helXsecsPt = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperPt), {"Wpt_preFSR", "weight", Form("%sharmonicsVec", _syst_weight.c_str())});
        _h1Group.push_back(helXsecsPt);

        TH1weightsHelper helperY(std::string("harmonicsY"), std::string("harmonicsY"), _nBinsY, _yArr, total);
        auto helXsecsY = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperY), {"Wrap_preFSR_abs", "weight", Form("%sharmonicsVec", _syst_weight.c_str())});
        _h1Group.push_back(helXsecsY);

        TH1weightsHelper helperMapPt(std::string("Pt"), std::string("Pt"), _nBinsPt, _ptArr, _syst_name);
        auto Pt = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperMapPt), {"Wpt_preFSR", "weight", _syst_weight});
        _h1Group.push_back(Pt);

        TH1weightsHelper helperMapY(std::string("Y"), std::string("Y"), _nBinsY, _yArr, _syst_name);
        auto Y = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperMapY), {"Wrap_preFSR_abs", "weight", _syst_weight});
        _h1Group.push_back(Y);
        
        // TH1weightsHelper helperMass(std::string("Wmass_preFSR"), std::string("Wmass_preFSR"), _Marr.size() - 1, _Marr, _syst_name);
        // auto Mass = d1.Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperMass), {"Wmass_preFSR", "weight", _syst_weight});
        // _h1Group.push_back(Mass);

        return d1;
    }
    // nominal
    else
    {
        auto d1 = d.Filter(_filter);
        TH2weightsHelper helper(std::string("harmonics"), std::string("harmonics"), _nBinsY, _yArr, _nBinsPt, _ptArr, total);
        auto helXsecs = d1.Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helper), {"Wrap_preFSR_abs", "Wpt_preFSR", "weight", Form("%sharmonicsVec", _syst_weight.c_str())});
        _h2Group.push_back(helXsecs);

        TH2weightsHelper helperSq(std::string("harmonicsSq"), std::string("harmonicsSq"), _nBinsY, _yArr, _nBinsPt, _ptArr, total);
        auto helXsecsSq = d1.Book<float, float, float, ROOT::VecOps::RVec<float>>(std::move(helperSq), {"Wrap_preFSR_abs", "Wpt_preFSR", "weight", Form("%sharmonicsVecSq", _syst_weight.c_str())});
        _h2Group.push_back(helXsecsSq);

        auto mapTot = d1.Histo2D(TH2D("mapTot", "mapTot", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "weight");
        _h2List.push_back(mapTot);

        TH1weightsHelper helperPt(std::string("harmonicsPt"), std::string("harmonicsPt"), _nBinsPt, _ptArr, total);
        auto helXsecsPt = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperPt), { "Wpt_preFSR", "weight", Form("%sharmonicsVec", _syst_weight.c_str()) });
        _h1Group.push_back(helXsecsPt);

        TH1weightsHelper helperY(std::string("harmonicsY"), std::string("harmonicsY"), _nBinsY, _yArr, total);
        auto helXsecsY = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Book<float, float, ROOT::VecOps::RVec<float>>(std::move(helperY), { "Wrap_preFSR_abs", "weight", Form("%sharmonicsVec", _syst_weight.c_str()) });
        _h1Group.push_back(helXsecsY);

        auto Pt = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Histo1D(TH1D("Pt", "Pt", _nBinsPt, _ptArr.data()), "Wpt_preFSR", "weight");
        _h1List.push_back(Pt);

        auto Y = d1.Filter("Wpt_preFSR<60. && Wrap_preFSR_abs<2.4").Histo1D(TH1D("Y", "Y", _nBinsY, _yArr.data()), "Wrap_preFSR_abs", "weight");
        _h1List.push_back(Y);
        
        // auto Mass = d1.Histo1D(TH1D("Wmass_preFSR", "Wmass_preFSR",_Marr.size() - 1, _Marr.data()), "Wmass_preFSR", "weight");
        // _h1List.push_back(Mass);
        
        // auto YqTcT = d1.Histo3D(TH3D("YqTcT", "YqTcT", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data(),  _cosThetaArr.size() -1, _cosThetaArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR","CStheta_preFSR","weight");
        // _h3List.push_back(YqTcT);    
                
        // auto norm = d1.Histo2D(TH2D("norm", "norm", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR","norm");
        // _h2List.push_back(norm);    
        return d1;
    }
}

// void AngCoeff::setAxisarrays()
// {
// //   for (int i=0; i<1001;i++) _Marr[i] = 75.+i*(85.-75.)/1000.;
//   for (int i=0; i<101;i++) _cosThetaArr[i] = -1.+i*(1.+1.)/100.;
// }