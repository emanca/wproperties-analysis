#include "interface/customizeforUL.hpp"
#include "interface/functions.hpp"
#include "ROOT/RDFHelpers.hxx"

RNode customizeforUL::run(RNode d)
{
  //for both data/mc
  auto d1 = d.Alias("MET_pt_nom", "MET_pt")
                .Alias("MET_phi_nom", "MET_phi");

  //for all MC
  /*if (isMC_)
    {
        d1 = d1.Alias("MET_pt_jesTotalUp", "MET_pt_jesTotalUp")
                 .Alias("MET_phi_jesTotalUp", "MET_phi_jesTotalUp")
                 .Alias("MET_pt_jesTotalDown", "MET_pt_jesTotalDown")
                 .Alias("MET_phi_jesTotalDown", "MET_phi_jesTotalDown");
    }*/
  //For only W&Z MC
  if (isWorZMC_)
  { //basic gen
    //Have to define the post nano sequence
    d1 = d1.Alias("Wrap_preFSR", "Vrap_preFSR")
             .Alias("Wpt_preFSR", "Vpt_preFSR")
             .Alias("Wmass_preFSR", "Vmass_preFSR")
             .Alias("GenPart_preFSRMuonIdx", "GenPart_preFSRLepIdx1");

    auto getSameVec = [](ROOT::VecOps::RVec<float> red) {
      return red;
    };

    /*
      Old Nano defn
      Float_t LHE scale variation weights (w_var / w_nominal); 
      [0] is muR=0.50000E+00 muF=0.50000E+00 ; 
      [1] is muR=0.50000E+00 muF=0.10000E+01 ; 
      [2] is muR=0.50000E+00 muF=0.20000E+01 ; 
      [3] is muR=0.10000E+01 muF=0.50000E+00 ; 
      [4] is muR=0.10000E+01 muF=0.10000E+01 ; 
      [5] is muR=0.10000E+01 muF=0.20000E+01 ; 
      [6] is muR=0.20000E+01 muF=0.50000E+00 ; 
      [7] is muR=0.20000E+01 muF=0.10000E+01 ; 
      [8] is muR=0.20000E+01 muF=0.20000E+01 *
      @ForSRC:Once we switch to new MC fully, we create the needed vector of 6 directly
    */
    // trigger matching not in legacy rereco
    //Will be defined later
    //d1 = d1.Define("Mu1_hasTriggerMatch", getIntFromIdx, {"Muon_hasTriggerMatch", "Idx_mu1"});

    /*no longer needed
        d1 = d1.Define("LHEScaleWeight", ROOT::Internal::RDF::PassAsVec<9, float>(getSameVec),
                       {"scaleWeightMuR05MuF05",
                        "scaleWeightMuR05MuF1",
                        "scaleWeightMuR05MuF2",
                        "scaleWeightMuR1MuF05",
                        "scaleWeightMuR1MuF1",
                        "scaleWeightMuR1MuF2",
                        "scaleWeightMuR2MuF05",
                        "scaleWeightMuR2MuF1",
                        "scaleWeightMuR2MuF2"});

        //pdf weights//
        d1 = d1.Define("nLHEPdfWeight", []() { unsigned int np = 102; return np; })
                 .Define("LHEPdfWeight", ROOT::Internal::RDF::PassAsVec<102, float>(getSameVec),
                         {
                             "pdfWeightNNPDF1",
                             "pdfWeightNNPDF2",
                             "pdfWeightNNPDF3",
                             "pdfWeightNNPDF4",
                             "pdfWeightNNPDF5",
                             "pdfWeightNNPDF6",
                             "pdfWeightNNPDF7",
                             "pdfWeightNNPDF8",
                             "pdfWeightNNPDF9",
                             "pdfWeightNNPDF10",
                             "pdfWeightNNPDF11",
                             "pdfWeightNNPDF12",
                             "pdfWeightNNPDF13",
                             "pdfWeightNNPDF14",
                             "pdfWeightNNPDF15",
                             "pdfWeightNNPDF16",
                             "pdfWeightNNPDF17",
                             "pdfWeightNNPDF18",
                             "pdfWeightNNPDF19",
                             "pdfWeightNNPDF20",
                             "pdfWeightNNPDF21",
                             "pdfWeightNNPDF22",
                             "pdfWeightNNPDF23",
                             "pdfWeightNNPDF24",
                             "pdfWeightNNPDF25",
                             "pdfWeightNNPDF26",
                             "pdfWeightNNPDF27",
                             "pdfWeightNNPDF28",
                             "pdfWeightNNPDF29",
                             "pdfWeightNNPDF30",
                             "pdfWeightNNPDF31",
                             "pdfWeightNNPDF32",
                             "pdfWeightNNPDF33",
                             "pdfWeightNNPDF34",
                             "pdfWeightNNPDF35",
                             "pdfWeightNNPDF36",
                             "pdfWeightNNPDF37",
                             "pdfWeightNNPDF38",
                             "pdfWeightNNPDF39",
                             "pdfWeightNNPDF40",
                             "pdfWeightNNPDF41",
                             "pdfWeightNNPDF42",
                             "pdfWeightNNPDF43",
                             "pdfWeightNNPDF44",
                             "pdfWeightNNPDF45",
                             "pdfWeightNNPDF46",
                             "pdfWeightNNPDF47",
                             "pdfWeightNNPDF48",
                             "pdfWeightNNPDF49",
                             "pdfWeightNNPDF50",
                             "pdfWeightNNPDF51",
                             "pdfWeightNNPDF52",
                             "pdfWeightNNPDF53",
                             "pdfWeightNNPDF54",
                             "pdfWeightNNPDF55",
                             "pdfWeightNNPDF56",
                             "pdfWeightNNPDF57",
                             "pdfWeightNNPDF58",
                             "pdfWeightNNPDF59",
                             "pdfWeightNNPDF60",
                             "pdfWeightNNPDF61",
                             "pdfWeightNNPDF62",
                             "pdfWeightNNPDF63",
                             "pdfWeightNNPDF64",
                             "pdfWeightNNPDF65",
                             "pdfWeightNNPDF66",
                             "pdfWeightNNPDF67",
                             "pdfWeightNNPDF68",
                             "pdfWeightNNPDF69",
                             "pdfWeightNNPDF70",
                             "pdfWeightNNPDF71",
                             "pdfWeightNNPDF72",
                             "pdfWeightNNPDF73",
                             "pdfWeightNNPDF74",
                             "pdfWeightNNPDF75",
                             "pdfWeightNNPDF76",
                             "pdfWeightNNPDF77",
                             "pdfWeightNNPDF78",
                             "pdfWeightNNPDF79",
                             "pdfWeightNNPDF80",
                             "pdfWeightNNPDF81",
                             "pdfWeightNNPDF82",
                             "pdfWeightNNPDF83",
                             "pdfWeightNNPDF84",
                             "pdfWeightNNPDF85",
                             "pdfWeightNNPDF86",
                             "pdfWeightNNPDF87",
                             "pdfWeightNNPDF88",
                             "pdfWeightNNPDF89",
                             "pdfWeightNNPDF90",
                             "pdfWeightNNPDF91",
                             "pdfWeightNNPDF92",
                             "pdfWeightNNPDF93",
                             "pdfWeightNNPDF94",
                             "pdfWeightNNPDF95",
                             "pdfWeightNNPDF96",
                             "pdfWeightNNPDF97",
                             "pdfWeightNNPDF98",
                             "pdfWeightNNPDF99",
                             "pdfWeightNNPDF100",
                             "pdfWeightNNPDF101",
                             "pdfWeightNNPDF102",
                         });

        // add variations for PDF
        std::vector<std::string> PDFVars;
        for (int i = 1; i < 103; i++)
            PDFVars.push_back("pdf_" + std::to_string(i));
        customizeforUL::vary("pdfWeightNNPDF0", PDFVars);
      */
  }

  return d1;
}
