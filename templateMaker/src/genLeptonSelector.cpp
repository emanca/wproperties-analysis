#include "interface/genLeptonSelector.hpp"
#include "interface/functions.hpp"
#include <algorithm>
#include <utility>

RNode genLeptonSelector::run(RNode d)
{
  //Logic taken from https://github.com/WMass/nanoAOD-tools/blob/master/python/postprocessing/wmass/genLepSelection.py
  //have to think of refinement
  auto getGenLeptonIdx = [](const RVec<int> &pdgId, const RVec<int> &status, const RVec<int> &genPartIdxMother, const RVec<int> &statusFlag, const RVec<float> &pt)
  {
    RVec<std::pair<int, float>> status746;
    RVec<std::pair<int, float>> other;

    for (unsigned int i = 0; i < pdgId.size(); i++)
    {
      if (genPartIdxMother[i] < 0)
        continue;
      int mompdgId = pdgId[genPartIdxMother[i]];
      if (abs(pdgId[i]) < 11 || abs(pdgId[i]) > 16)
        continue;
      if (abs(mompdgId) == 23 || abs(mompdgId) == 24)
      {                       //
        if (status[i] == 746) //status 746 is pre photos FSR in new powheg samples. so if they exist we want those.
          status746.emplace_back(std::make_pair(i, pt[i]));
        else
          other.emplace_back(std::make_pair(i, pt[i]));
      }
      else if (status[i] == 23) //in madgraph5_aMC@NLO there are some events without a W/Z, those have status 23 leptons (tested on >1M events)
        other.emplace_back(std::make_pair(i, pt[i]));
    } //loop over genP

    if (other.size() > 2)
    {
      RVec<std::pair<int, float>> otherTmp;
      for (auto &gp : other)
      {
        if ((statusFlag[gp.first] >> 8) & 1)
          otherTmp.emplace_back(gp);
      }
      other = otherTmp;
    }
    //if 746 size is size2, then we return it by default
    std::pair<int, int> prefsrlepidx(-1, -1);
    if (status746.size() == 2)
    {
      prefsrlepidx.first = status746[0].second > status746[1].second ? status746[0].first : status746[1].first;
      prefsrlepidx.second = status746[0].second > status746[1].second ? status746[1].first : status746[0].first;
    }
    else if (status746.size() == 1 || other.size() == 1)
    {
      prefsrlepidx.first = status746[0].second > other[0].second ? status746[0].first : other[0].first;
      prefsrlepidx.second = status746[0].second > other[0].second ? other[0].first : status746[0].first;
    }
    else if (status746.size() == 0 and other.size() == 2)
    {
      prefsrlepidx.first = other[0].second > other[1].second ? other[0].first : other[1].first;
      prefsrlepidx.second = other[0].second > other[1].second ? other[1].first : other[0].first;
    }

    //swap indices to save +ve lepton as first index
    std::vector<int> genLep;
    genLep.push_back(pdgId[prefsrlepidx.first] < 0 ? prefsrlepidx.second : prefsrlepidx.first);
    genLep.push_back(pdgId[prefsrlepidx.first] < 0 ? prefsrlepidx.first : prefsrlepidx.second);

    //return prefsrlepidx;
    return genLep;
  }; //function

  auto getVtype = [](const int idx1, const int idx2, const RVec<int> &pdg)
  {
    int vtype = -1;
    if (idx1 != -1 && idx2 == -1) //only 1 lepton
      vtype = pdg[idx1];
    else if (idx1 == -1 && idx2 == -1)
    { //https://github.com/WMass/nanoAOD-tools/blob/master/python/postprocessing/wmass/genLepSelection.py#L78-L86
      // int vcharge = pdg[idx1] % 2 ? -1 * pdg[idx1] / std::abs(pdg[idx1]) : -1 * pdg[idx2] / std::abs(pdg[idx2]);
      // vtype = vcharge * int((abs(pdg[idx1]) + abs(pdg[idx1])) / 2.);
      vtype = -99;
    }
    else if (idx1 != -1 && idx2 != -1)
    {
      vtype = pdg[idx1] % 2 ? -1 * pdg[idx1] : -1 * pdg[idx2];
    }
      // std::cout << vtype<< std::endl;
    return vtype;
  };

  auto d1 = d.Define("SelectedGenPartIdxs", getGenLeptonIdx, {"GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother", "GenPart_statusFlags", "GenPart_pt"})
                .Define("GenPart_preFSRLepIdx1", [](const std::vector<int> &gpIdxs)
                        { return gpIdxs[0]; },
                        {"SelectedGenPartIdxs"})
                .Define("GenPart_preFSRLepIdx2", [](const std::vector<int> &gpIdxs)
                        { return gpIdxs[1]; },
                        {"SelectedGenPartIdxs"})
                .Define("genVtype", getVtype, {"GenPart_preFSRLepIdx1", "GenPart_preFSRLepIdx2", "GenPart_pdgId"});

  return d1;
}
