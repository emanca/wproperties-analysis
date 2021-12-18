#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "interface/accMap.hpp"

RNode accMap::run(RNode d){
    
    // for (float i=0.; i<=200.1; i=i+0.1){ _ptArr.push_back(i); }
    // int _nBinsPt = _ptArr.size()-1;
    // for (float i=0.; i<=6.1; i=i+0.1){ _yArr.push_back(i); }
    // int _nBinsY = _yArr.size()-1;
    
    auto d1 = d.Filter(_filter);

    auto mapTot = d1.Histo2D(TH2D("mapTot", "mapTot", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "weight");
    auto mapAccEta = d1.Filter("fabs(Mueta_preFSR)<2.4").Histo2D(TH2D("mapAccEta", "mapAccEta", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "weight");
    auto mapAcc = d1.Filter("fabs(Mueta_preFSR)<2.4 && Mupt_preFSR>25. && Mupt_preFSR<55.").Histo2D(TH2D("mapAcc", "mapAcc", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "weight");
    // auto sumw = d1.Define("genratio", "Generator_weight/genEventSumw").Define("genratioRew", "genratio*weightPt*weightY").Histo2D(TH2D("sumw", "sumw", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "genratioRew");
    auto sumw = d1.Define("genratio", "Generator_weight/genEventSumw").Define("genratioRew", "genratio*weightPt").Histo2D(TH2D("sumw", "sumw", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "genratioRew");

    _h2List.push_back(mapTot);
    _h2List.push_back(mapAccEta);
    _h2List.push_back(mapAcc);
    _h2List.push_back(sumw);
    
        
    // auto mapTotLowPtCut = d1.Filter("Mupt_preFSR>25.").Histo2D(TH2D("mapTotLowPtCut", "mapTotLowPtCut", _nBinsY, _yArr.data(), _nBinsPt, _ptArr.data()), "Wrap_preFSR_abs", "Wpt_preFSR", "weight");
    // auto mapQtPt = d1.Filter("fabs(Mueta_preFSR)<2.4 && Mupt_preFSR>25. && Mupt_preFSR<55.").Histo2D(TH2D("mapQtPt", "mapQtPt", _nBinsPt, _ptArr.data(), _nBinsPtMu, _ptMuArr.data()), "Wpt_preFSR", "Mupt_preFSR", "weight");
    // auto mapTotYEta = d1.Define("absEtaMu","fabs(Mueta_preFSR)").Filter("fabs(Mueta_preFSR)<2.4 && Mupt_preFSR>25. && Mupt_preFSR<55.").Histo2D(TH2D("mapTotYEta", "mapTotYEta", _nBinsY, _yArr.data(), _nBinsEtaMu, _etaMuArr.data()), "Wrap_preFSR_abs", "absEtaMu", "weight");
    // _h2List.push_back(mapTotLowPtCut);
    // _h2List.push_back(mapQtPt);
    // _h2List.push_back(mapTotYEta);
    
    
    return d;

}