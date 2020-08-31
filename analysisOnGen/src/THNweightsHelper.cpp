#include "THn.h"
#include "interface/THNweightsHelper.hpp"
   /// This constructor takes all the parameters necessary to build the THnTs. In addition, it requires the names of
   /// the columns which will be used.
using THn_t = THnT<float>;
THNweightsHelper::THNweightsHelper(std::string name, std::string title,
                    std::array<int, 4> nbins, std::array<double, 4> xmins,
                    std::array<double, 4> xmax,
                    std::vector<std::string> weightNames)
   {
      _name = name;
      _nbins = nbins;
      _xmins = xmins;
      _xmax = xmax;
      _weightNames = weightNames;

      const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
      for (auto slot : ROOT::TSeqU(nSlots)) {
	fHistos.emplace_back(std::make_shared<std::vector<std::unique_ptr<THn_t>>>());
         (void)slot;

         std::vector<std::unique_ptr<THn_t>> &histos = *fHistos[slot];
         auto n_histos = _weightNames.size();

         std::string slotnum = "";
         slotnum = slot > 0 ? std::to_string(slot) : "";

         for (unsigned int i = 0; i < n_histos; ++i){
	   histos.push_back(std::make_unique<THn_t>(std::string(_name + _weightNames[i] + slotnum).c_str(),
				    std::string(_name + _weightNames[i] + slotnum).c_str(),
				    4,
				    _nbins.data(),
				    _xmins.data(),
				    _xmax.data()));
	   //histos.back().SetDirectory(nullptr);
        }

      }
   }

std::shared_ptr<std::vector<std::unique_ptr<THn_t>>> THNweightsHelper::GetResultPtr() const { return fHistos[0]; }
   void THNweightsHelper::Initialize() {}
   void THNweightsHelper::InitTask(TTreeReader *, unsigned int) {}
   /// This is a method executed at every entry

void THNweightsHelper::Exec(unsigned int slot, const float &var1, const float &var2, const float &var3, const float &var4, const float &weight, const ROOT::VecOps::RVec<float> &weights)
   {
      auto &histos = *fHistos[slot];
      const auto n_histos = histos.size();
      //std::array<double, sizeof...(ColumnTypes)> valuesArr{static_cast<double>(values)...};
      std::vector<double> values;
      values.emplace_back(var1);
      values.emplace_back(var2);
      values.emplace_back(var3);
      values.emplace_back(var4);
      for (unsigned int i = 0; i < n_histos; ++i)
	histos[i]->Fill(values.data(), weight * weights[i]);
}
   /// This method is called at the end of the event loop. It is used to merge all the internal THnTs which
   /// were used in each of the data processing slots.
   void THNweightsHelper::Finalize()
   {
      auto &res_vec = *fHistos[0];
      for (auto slot : ROOT::TSeqU(1, fHistos.size())) {
         auto& histo_vec = *fHistos[slot];
         for (auto i : ROOT::TSeqU(0, res_vec.size()))
           res_vec[i]->Add(&(*histo_vec[i]));
      }
   }
   std::string THNweightsHelper::GetActionName(){
      return "THNweightsHelper";
   }
