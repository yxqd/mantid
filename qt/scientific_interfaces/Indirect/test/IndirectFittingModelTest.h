#ifndef MANTID_INDIRECTFITTINGMODELTEST_H_
#define MANTID_INDIRECTFITTINGMODELTEST_H_

#include <cxxtest/TestSuite.h>

#include "IndirectFitData.h"
#include "IndirectFittingModel.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"

#include <iostream>

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace MantidQt::CustomInterfaces::IDA;

namespace {

/// Simple class to set up the ADS with the configuration required
struct SetUpADSWithWorkspace {

  SetUpADSWithWorkspace(std::string const &inputWSName,
                        Workspace2D_sptr const &workspace) {
    AnalysisDataService::Instance().addOrReplace(inputWSName, workspace);
  };

  ~SetUpADSWithWorkspace() { AnalysisDataService::Instance().clear(); };
};

/// This is used to compare Spectra which is implemented as a boost::variant
struct AreSpectraEqual : public boost::static_visitor<bool> {

  template <typename T, typename U>
  bool operator()(const T &, const U &) const {
    return false; // cannot compare different types
  }

  template <typename T> bool operator()(const T &lhs, const T &rhs) const {
    return lhs == rhs;
  }
};

class DummyModel : public IndirectFittingModel {
public:
  ~DummyModel(){};

private:
  std::string sequentialFitOutputName() const override { return ""; };
  std::string simultaneousFitOutputName() const override { return ""; };
  std::string singleFitOutputName(std::size_t index,
                                  std::size_t spectrum) const override {
    return "";
  };
};

std::unique_ptr<DummyModel> getModel() {
  return std::make_unique<DummyModel>();
}

std::unique_ptr<DummyModel>
getModelWithWorkspace(std::string const &workspaceName,
                      std::size_t const &numberOfSpectra,
                      std::size_t const &numberOfBins = 3) {
  auto model = getModel();
  auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(
      numberOfSpectra, numberOfBins);
  SetUpADSWithWorkspace ads(workspaceName, workspace);
  model->addWorkspace(workspaceName);
  return model;
}

// std::vector<std::unique_ptr<IndirectFitData>>
// getModelWithFittingData(std::string const &workspaceName,
//                        std::size_t const &numberOfSpectra,
//                        std::size_t const &numberOfBins) {
//  auto model = getModel();
//  auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(
//      numberOfSpectra, numberOfBins);
//  SetUpADSWithWorkspace ads(workspaceName, workspace);
//  //Spectra const spec = std::make_pair(0u, workspace->getNumberHistograms() -
//  1);
//  //IndirectFitData fitData(workspace, spec);
//
//  //std::vector<std::unique_ptr<IndirectFitData>> data{
//  std::make_unique<IndirectFitData>(fitData) };
//
//  //PrivateFittingData privateFitData({
//  std::make_unique<IndirectFitData>(fitData) });
//
//  return data;
//}

} // namespace

class IndirectFittingModelTest : public CxxTest::TestSuite {
public:
  static IndirectFittingModelTest *createSuite() {
    return new IndirectFittingModelTest();
  }

  static void destroySuite(IndirectFittingModelTest *suite) { delete suite; }

  void test_model_is_instantiated_correctly() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT(model->getWorkspace(0));
    TS_ASSERT_EQUALS(model->numberOfWorkspaces(), 1);
    TS_ASSERT_EQUALS(model->getNumberOfSpectra(0), 3);
  }

  void test_that_a_workspace_is_stored_correctly_in_the_ADS() {
    auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads("WorkspaceName", workspace);

    TS_ASSERT(AnalysisDataService::Instance().doesExist("WorkspaceName"));
    MatrixWorkspace_sptr storedWorkspace =
        boost::dynamic_pointer_cast<MatrixWorkspace>(
            AnalysisDataService::Instance().retrieve("WorkspaceName"));
    TS_ASSERT_EQUALS(storedWorkspace->getNumberHistograms(), 3);
  }

  void
  test_that_addWorkspace_will_add_a_workspace_to_the_fittingData_using_the_workspace_name() {
    auto model = getModel();
    auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads("WorkspaceName", workspace);

    model->addWorkspace("WorkspaceName");

    TS_ASSERT_EQUALS(model->getWorkspace(0), workspace);
  }

  void
  test_that_addWorkspace_throws_when_provided_a_workspace_name_and_an_empty_spectraString() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);

    std::string const spectraString("");

    TS_ASSERT_THROWS(model->addWorkspace("WorkspaceName", spectraString),
                     std::runtime_error);
  }

  void
  test_that_addWorkspace_combines_an_input_workspace_with_a_workspace_that_already_exists_if_the_workspaces_have_the_same_name() {
    auto model = getModelWithWorkspace("WorkspaceSameName", 3);
    auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads("WorkspaceSameName", workspace);

    model->addWorkspace("WorkspaceSameName");

    TS_ASSERT(model->getWorkspace(0));
    TS_ASSERT(!model->getWorkspace(1));
  }

  void
  test_that_addWorkspace_does_not_combine_an_input_workspace_with_a_workspace_that_already_exists_if_the_workspaces_are_differently_named() {
    auto model = getModel();
    auto const workspace1 = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    auto const workspace2 = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads1("WorkspaceName1", workspace1);
    SetUpADSWithWorkspace ads2("WorkspaceName2", workspace2);

    model->addWorkspace("WorkspaceName1");
    model->addWorkspace("WorkspaceName2");

    TS_ASSERT_EQUALS(model->getWorkspace(0), workspace1);
    TS_ASSERT_EQUALS(model->getWorkspace(1), workspace2);
  }

  void
  test_that_getWorkspace_returns_a_nullptr_when_getWorkspace_is_provided_an_out_of_range_index() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT_EQUALS(model->getWorkspace(1), nullptr);
  }

  void
  test_that_getSpectra_returns_a_correct_spectra_when_the_index_provided_is_valid() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);

    Spectra const inputSpectra = DiscontinuousSpectra<std::size_t>("0-1");
    model->setSpectra(inputSpectra, 0);

    TS_ASSERT(boost::apply_visitor(AreSpectraEqual(), model->getSpectra(0),
                                   inputSpectra));
  }

  void
  test_that_getSpectra_returns_an_empty_DiscontinuousSpectra_when_provided_an_out_of_range_index() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);

    Spectra const spec(DiscontinuousSpectra<std::size_t>(""));

    TS_ASSERT(
        boost::apply_visitor(AreSpectraEqual(), model->getSpectra(3), spec));
  }

  void
  test_that_getFittingRange_returns_correct_range_when_provided_a_valid_index_and_spectrum() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setStartX(1.2, 0, 0);
    model->setEndX(5.6, 0, 0);

    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).first, 1.2);
    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).second, 5.6);
  }

  void
  test_that_getFittingRange_returns_empty_range_when_provided_an_out_of_range_dataIndex() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setStartX(1.2, 0, 0);
    model->setEndX(5.6, 0, 0);

    TS_ASSERT_EQUALS(model->getFittingRange(1, 0).first, 0.0);
    TS_ASSERT_EQUALS(model->getFittingRange(1, 0).second, 0.0);
  }

  void
  test_that_getFittingRange_returns_empty_range_when_there_are_zero_spectra() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setStartX(1.2, 0, 0);
    model->setEndX(5.6, 0, 0);
    DiscontinuousSpectra<std::size_t> const emptySpec("");
    model->setSpectra(emptySpec, 0);

    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).first, 0.0);
    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).second, 0.0);
  }

  void
  test_that_getExcludeRegion_returns_correct_range_when_provided_a_valid_index_and_spectrum() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setExcludeRegion("0,1,3,4", 0, 0);

    TS_ASSERT_EQUALS(model->getExcludeRegion(0, 0), "0.0,1.0,3.0,4.0");
  }

  void
  test_that_getExcludeRegion_returns_empty_range_when_provided_an_out_of_range_dataIndex() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setExcludeRegion("0,1,3,4", 0, 0);

    TS_ASSERT_EQUALS(model->getExcludeRegion(1, 0), "");
  }

  void
  test_that_getExcludeRegion_returns_empty_range_when_there_are_zero_spectra() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setExcludeRegion("0,1,3,4", 0, 0);
    DiscontinuousSpectra<std::size_t> const emptySpec("");
    model->setSpectra(emptySpec, 0);

    TS_ASSERT_EQUALS(model->getExcludeRegion(1, 0), "");
  }

  void
  test_that_getExcludeRegion_returns_a_region_where_each_range_is_in_order_after_setExcludeRegion_is_given_an_unordered_region_string() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    model->setExcludeRegion("0,1,6,4", 0, 0);

    TS_ASSERT_EQUALS(model->getExcludeRegion(0, 0), "0.0,1.0,4.0,6.0");
  }

  void
  test_that_createDisplayName_returns_valid_string_when_provided_an_in_range_dataIndex() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);

    std::string const formatString = "%1%_s%2%_Gaussian";
    std::string const rangeDelimiter = "_to_";

    TS_ASSERT_EQUALS(model->createOutputName(formatString, rangeDelimiter, 0),
                     "WorkspaceName_s0_Gaussian_Result");
  }

  void
  test_that_createDisplayName_returns_string_with_red_removed_from_the_workspace_name() {
    auto model = getModelWithWorkspace("Workspace_3456_red", 1);

    std::string const formatString = "%1%_s%2%_Gaussian";
    std::string const rangeDelimiter = "_to_";

    TS_ASSERT_EQUALS(model->createOutputName(formatString, rangeDelimiter, 0),
                     "Workspace_3456_s0_Gaussian_Result");
  }

  void
  test_that_createDisplayName_returns_correct_name_when_provided_a_valid_rangeDelimiter_and_formatString() {
    auto model = getModelWithWorkspace("Workspace_3456_red", 1);

    std::vector<std::string> const formatStrings{
        "%1%_s%2%_Gaussian", "%1%_f%2%,s%2%_MSD", "%1%_s%2%_TeixeiraWater"};
    std::string const rangeDelimiter = "_to_";

    TS_ASSERT_EQUALS(
        model->createOutputName(formatStrings[0], rangeDelimiter, 0),
        "Workspace_3456_s0_Gaussian_Result");
    TS_ASSERT_EQUALS(
        model->createOutputName(formatStrings[1], rangeDelimiter, 0),
        "Workspace_3456_f0+s0_MSD_Result");
    TS_ASSERT_EQUALS(
        model->createOutputName(formatStrings[2], rangeDelimiter, 0),
        "Workspace_3456_s0_TeixeiraWater_Result");
  }

  void
  test_that_isMultiFit_returns_true_when_there_are_more_than_one_workspaces_stored_in_the_model() {
    auto model = getModelWithWorkspace("Workspace1", 1);
    auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads("Workspace2", workspace);

    model->addWorkspace("Workspace2");

    TS_ASSERT(model->isMultiFit());
  }

  void
  test_that_isMultiFit_returns_false_when_there_is_one_workspace_stored_in_the_model() {
    auto model = getModelWithWorkspace("Workspace1", 1);
    TS_ASSERT(!model->isMultiFit());
  }

  void
  test_that_isPreviouslyFit_returns_false_if_there_is_no_previous_fit_output_data() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);
    TS_ASSERT(!model->isPreviouslyFit(0, 0));
  }

  void
  test_that_isPreviouslyFit_returns_false_if_the_dataIndex_is_out_of_range() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);
    TS_ASSERT(!model->isPreviouslyFit(4, 0));
  }

  void
  test_that_isPreviouslyFit_returns_true_if_the_spectrum_has_been_fitted_previously() {
    // NOT FINISHED - how to set active function?

    // auto model = getModelWithWorkspace("WorkspaceName", 5);
    // MatrixWorkspace_sptr work = model->getWorkspace(0);

    // std::string const fun =
    //    "name=FlatBackground,A0=0;(composite=Convolution,FixResolution=true,"
    //    "NumDeriv=true;name=Resolution,Workspace=WorkspaceName,WorkspaceIndex="
    //    "0,X=(),Y=();(name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0.0175))";

    // IAlgorithm_sptr convolutionFit =
    //    AlgorithmManager::Instance().create("ConvolutionFitSequential");
    // convolutionFit->initialize();
    // convolutionFit->setProperty("InputWorkspace	", work);
    // convolutionFit->setProperty("StartX", -0.55);
    // convolutionFit->setProperty("EndX", 0.55);
    // convolutionFit->setProperty("Function", fun);
    // convolutionFit->setProperty("OutputWorkspace", "outputws");
    // convolutionFit->execute();

    // Spectra const inputSpectra = DiscontinuousSpectra<std::size_t>("0-1");
    // std::string fun =
    //    "name=FlatBackground,A0=0;(composite=Convolution,FixResolution=true,"
    //    "NumDeriv=true;name=Resolution,Workspace=WorkspaceName,WorkspaceIndex="
    //    "0,X=(),Y=();(name=Lorentzian,Amplitude=1,PeakCentre=0,FWHM=0.0175))";
    // IFunction_sptr fitFunction =
    //    FunctionFactory::Instance().createInitialized(fun);
    // model->setStartX(0.0, 0, 0);
    // model->setEndX(5.3, 0, 0);
    // model->setSpectra(inputSpectra, 0);
    // model->setFitFunction(fitFunction);
    // auto fitAlgorithm = model->getFittingAlgorithm();

    // model->addOutput(fitAlgorithm);

    // TS_ASSERT(!model->isPreviouslyFit(4, 0));
  }

  void test_that_hasZeroSpectra_returns_true_if_workspace_has_zero_spectra() {
    auto model = getModel();
    auto const workspace = boost::make_shared<Workspace2D>();
    SetUpADSWithWorkspace ads("WorkspaceEmpty", workspace);

    model->addWorkspace("WorkspaceEmpty");

    TS_ASSERT(model->hasZeroSpectra(0));
  }

  void
  test_that_hasZeroSpectra_returns_true_if_the_dataIndex_provided_is_out_of_range() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);
    TS_ASSERT(model->hasZeroSpectra(1));
  }

  void
  test_that_hasZeroSpectra_returns_false_if_workspace_contains_one_or_more_spectra() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);
    TS_ASSERT(!model->hasZeroSpectra(0));
  }

  void
  test_that_isInvalidFunction_returns_a_message_when_no_activeFunction_exists() {
    auto model = getModelWithWorkspace("WorkspaceName", 1);
    TS_ASSERT(model->isInvalidFunction());
  }

  void
  test_that_isInvalidFunction_returns_a_message_when_the_activeFunction_contains_zero_parameters_or_functions() {
    // NOT FINISHED - how to set active function?
  }

  void test_isInvalidFunction_returns_none_if_the_activeFunction_is_valid() {
    // NOT FINISHED - how to set active function?
  }

  void
  test_that_numberOfWorkspace_returns_the_number_of_workspace_stored_by_model() {
    auto model = getModelWithWorkspace("Workspace1", 1);
    auto const workspace = WorkspaceCreationHelper::create2DWorkspace123(3, 3);
    SetUpADSWithWorkspace ads("Workspace2", workspace);

    model->addWorkspace("Workspace2");

    TS_ASSERT_EQUALS(model->numberOfWorkspaces(), 2);
  }

  void
  test_that_getNumberOfSpectra_returns_zero_if_dataIndex_is_out_of_range() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT_EQUALS(model->getNumberOfSpectra(1), 0);
  }

  void
  test_that_getNumberOfSpectra_returns_the_number_of_spectra_stored_in_the_workspace_given() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT_EQUALS(model->getNumberOfSpectra(0), 3);
  }

  void
  test_that_getFitParameterNames_returns_an_empty_vector_if_the_fitOutput_is_empty() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT_EQUALS(model->getFitParameterNames(), std::vector<std::string>());
  }

  void
  test_that_getFitParameterNames_returns_a_vector_of_fit_parameters_if_the_fitOutput_contains_parameters() {
    // NOT FINISHED - how to get fitOutput to store output?
  }

  void test_getFittingFunction_returns_null_if_there_is_no_fitting_function() {
    auto model = getModelWithWorkspace("WorkspaceName", 3);
    TS_ASSERT_EQUALS(model->getFittingFunction(), nullptr);
  }

  void
  test_that_getFittingFunction_returns_the_fitting_function_which_is_set() {
    // NOT FINISHED - how to set fitting function?
  }

  void test_that_PrivateFittingData_is_set_correctly() {
    // NOT FINISHED - how to create and set PrivateFittingData?
  }

  void
  test_that_setFittingData_will_set_the_fittingData_to_the_data_provided() {
    // NOT FINISHED - how to create and set PrivateFittingData, and then set
    // this.
  }

  void
  test_that_setSpectra_will_set_the_spectra_to_the_provided_inputSpectra() {
    auto model = getModelWithWorkspace("WorkspaceName", 10);

    Spectra const inputSpectra = DiscontinuousSpectra<std::size_t>("2,4,6-8");
    model->setSpectra(inputSpectra, 0);

    TS_ASSERT(boost::apply_visitor(AreSpectraEqual(), model->getSpectra(0),
                                   inputSpectra));
  }

  void
  test_that_setSpectra_will_set_the_spectra_when_provided_a_spectra_pair() {
    auto model = getModelWithWorkspace("WorkspaceName", 10);

    Spectra const inputSpectra = std::make_pair(0u, 5u);
    model->setSpectra(inputSpectra, 0);

    TS_ASSERT(boost::apply_visitor(AreSpectraEqual(), model->getSpectra(0),
                                   inputSpectra));
  }

  void
  test_that_setSpectra_does_not_throw_when_provided_an_out_of_range_dataIndex() {
    auto model = getModelWithWorkspace("WorkspaceName", 5);
    TS_ASSERT_THROWS_NOTHING(model->getSpectra(1));
  }

  void
  test_that_setStartX_will_set_the_startX_at_the_first_dataIndex_when_the_fit_is_sequential() {
    auto model = getModelWithWorkspace("WorkspaceName", 5);

    model->setStartX(4.0, 3, 0);

    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).first, 4.0);
  }

  void
  test_that_setEndX_will_set_the_endX_at_the_first_dataIndex_when_the_fit_is_sequential() {
    auto model = getModelWithWorkspace("WorkspaceName", 5);

    model->setEndX(4.0, 3, 0);

    TS_ASSERT_EQUALS(model->getFittingRange(0, 0).second, 4.0);
  }

  void
  test_that_setExcludeRegion_set_the_excludeRegion_at_the_first_dataIndex_when_the_fit_is_sequential() {
    auto model = getModelWithWorkspace("WorkspaceName", 5);

    model->setExcludeRegion("0,1,3,4", 3, 0);

    TS_ASSERT_EQUALS(model->getExcludeRegion(0, 0), "0.0,1.0,3.0,4.0");
  }
};

#endif
