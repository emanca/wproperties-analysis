#ifndef MODULE_H
#define MODULE_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "boost/rank_mod.hpp"
#include <boost/histogram.hpp>
#include <memory>
#include <tuple>

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using boost_histogram = boost::histogram::histogram<std::vector<boost::histogram::axis::variable<>>, boost::histogram::storage_adaptor<std::vector<boost::histogram::accumulators::weighted_sum<>, std::allocator<boost::histogram::accumulators::weighted_sum<>>>>>;

enum HistoCategory
{
  Nominal = 0,
  Corrected,
  JME,
  Weights
};

class Module
{

private:
public:
  virtual ~Module(){};
  virtual RNode run(RNode d) = 0;

  std::vector<ROOT::RDF::RResultPtr<TH1D>> _h1List;
  std::vector<ROOT::RDF::RResultPtr<TH2D>> _h2List;
  std::vector<ROOT::RDF::RResultPtr<TH3D>> _h3List;

  // groups of histos
  std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D *>>> _h1Group;
  std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D *>>> _h2Group;
  std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D *>>> _h3Group;

  std::vector<ROOT::RDF::RResultPtr<std::map<std::string, boost_histogram>>> _hNGroup;

  std::vector<ROOT::RDF::RResultPtr<TH1D>> getTH1();
  std::vector<ROOT::RDF::RResultPtr<TH2D>> getTH2();
  std::vector<ROOT::RDF::RResultPtr<TH3D>> getTH3();

  std::vector<ROOT::RDF::RResultPtr<std::vector<TH1D *>>> getGroupTH1();
  std::vector<ROOT::RDF::RResultPtr<std::vector<TH2D *>>> getGroupTH2();
  std::vector<ROOT::RDF::RResultPtr<std::vector<TH3D *>>> getGroupTH3();
  std::vector<ROOT::RDF::RResultPtr<std::map<std::string, boost_histogram>>> getGroupTHN();

  void reset();
};

#endif
