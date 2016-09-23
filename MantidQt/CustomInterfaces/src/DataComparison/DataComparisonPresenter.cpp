#include "MantidQtCustomInterfaces/DataComparison/DataComparisonPresenter.h"
#include "MantidQtCustomInterfaces/DataComparison/DataComparisonHelper.h"
#include "MantidQtCustomInterfaces/DataComparison/IDataComparisonView.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/MatrixWorkspace.h"
#include <algorithm>

using namespace Mantid::API;
using namespace MantidQt::CustomInterfaces;

/** Constructor
*
* @param view :: a pointer to the interface for the view we are managing
*/
DataComparisonPresenter::DataComparisonPresenter(IDataComparisonView *view)
    : WorkspaceObserver(), m_view(view), m_diffWsName() {

  observeAfterReplace();
  observeRename();
  observePreDelete();
}

/** Destructor
*
*/
DataComparisonPresenter::~DataComparisonPresenter() {}

/** Deals with notifications triggered by the view
*
* @param notification :: a flag indicating the type of notification
*/
void DataComparisonPresenter::notify(
    IDataComparisonPresenter::Notification notification) {

  switch (notification) {

  case IDataComparisonPresenter::AddWorkspace:
    addWorkspace();
    break;
  case IDataComparisonPresenter::PlotWorkspaces:
    plotWorkspaces();
    break;
  case IDataComparisonPresenter::PlotDiffWorkspaces:
    plotDiffWorkspace();
    break;
  case IDataComparisonPresenter::RemoveAllWorkspaces:
    removeAllWorkspaces();
    break;
  case IDataComparisonPresenter::RemoveSelectedWorkspaces:
    removeSelectedWorkspaces();
    break;
  case IDataComparisonPresenter::RemoveDiffWorkspace:
    removeDiffWorkspace();
    break;
  case IDataComparisonPresenter::WorkspaceIndexChanged:
    plotWorkspaces();
    break;
  case IDataComparisonPresenter::ColorChanged:
    plotWorkspaces();
    break;
  }
}

/** Add data to the interface
*
*/
void DataComparisonPresenter::addWorkspace() {

  // First get from the view the selected data name the user wants to load
  auto dataName = m_view->getSelectedWorkspaceName();

  // Check if the ws exists in the ADS
  if (!AnalysisDataService::Instance().doesExist(dataName)) {
    return;
  }

  // Get the workspace
  Workspace_const_sptr ws =
      AnalysisDataService::Instance().retrieveWS<Workspace>(dataName);
  WorkspaceGroup_const_sptr wsGroup =
      boost::dynamic_pointer_cast<const WorkspaceGroup>(ws);

  // Tell the view to block signals emitted by the table
  m_view->blockTableSignals(true);

  // If this is a WorkspaceGroup then add all items
  if (wsGroup != NULL) {
    size_t numWs = wsGroup->size();
    for (size_t wsIdx = 0; wsIdx < numWs; wsIdx++) {
      addWorkspace(wsGroup->getItem(wsIdx));
    }
  }
  // Otherwise just add the single workspace
  else {
    addWorkspace(ws);
  }

  // Tell the view to unblock signals
  m_view->blockTableSignals(false);

  // Replot the workspaces
  plotWorkspaces();
}

/** Add a workspace to the data table
*
*/
void DataComparisonPresenter::addWorkspace(Workspace_const_sptr ws) {

  // First check that the workspace is the correct type
  // If it is not, tell the view to print an error message

  MatrixWorkspace_const_sptr matrixWs =
      boost::dynamic_pointer_cast<const MatrixWorkspace>(ws);
  if (!matrixWs) {
    m_view->printError("Workspace is of incorrect type");
    return;
  }

  // Then check that the workspace does not already exist in the data table

  if (m_view->containsWorkspace(matrixWs->name())) {
    m_view->printInformation("Workspace already shown in comparison");
    return;
  }

  // Update the table
  m_view->addWorkspace(matrixWs->name(), getInitialColourIndex());
  // Plot the data
  m_view->plotCurve(matrixWs->name(),
                    DataComparisonHelper::curveDataFromWs(matrixWs, 0));
}

/** Plot workspaces loaded into the interface
*
*/
void DataComparisonPresenter::plotWorkspaces() {

  int wsIndex = m_view->getSelectedWorkspaceIndex();
  auto wsColors = m_view->getWorkspaceColors();
  auto wsNames = m_view->getWorkspaceNames();

  for (int ws = 0; ws < static_cast<int>(wsNames.size()); ws++) {

    std::string wsName = wsNames[ws];
    std::string wsColor = wsColors[ws];

    // Check if ws exists in the ADS
    if (!AnalysisDataService::Instance().doesExist(wsName)) {
      continue;
    }

    // Retrieve the workspace
    MatrixWorkspace_const_sptr workspace =
        AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(wsName);
    // and the number of histograms
    int numSpec = static_cast<int>(workspace->getNumberHistograms());

    if (wsIndex >= numSpec) {
      m_view->printError("Workspace index for workspace " + wsName +
                         " is out of range");
      m_view->detachCurve(wsName);
      continue;
    }

    // Detach the old curve from the plot if it exists
    m_view->detachCurve(wsName);
    // Plot the data
    m_view->plotCurve(wsName,
                      DataComparisonHelper::curveDataFromWs(workspace, wsIndex),
                      wsColor);
  }
}

/** Plot diff workspace
*
*/
void DataComparisonPresenter::plotDiffWorkspace() {

  // Detach old curve
  m_view->detachCurve("__diff");

  auto wsNames = m_view->getSelectedWorkspaceNames();
  auto wsIndex = m_view->getSelectedWorkspaceIndex();

  // Print error if there are not two workspaces
  if (wsNames.size() != 2) {
    m_view->printError(
        "Need to have exactly 2 workspaces selected for diff (have " +
        std::to_string(wsNames.size()) + ")");
    return;
  }

  if (!AnalysisDataService::Instance().doesExist(wsNames[0])) {
    m_view->printError("Workspace " + wsNames[0] +
                       " does not exist in the ADS");
  }
  if (!AnalysisDataService::Instance().doesExist(wsNames[1])) {
    m_view->printError("Workspace " + wsNames[1] +
                       " does not exist in the ADS");
  }

  // Get pointers to the workspaces to be diffed
  MatrixWorkspace_sptr ws1 =
      AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(wsNames[0]);
  MatrixWorkspace_sptr ws2 =
      AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(wsNames[1]);

  // Extract the current spectrum for both workspaces
  IAlgorithm_sptr extractWs1Alg =
      AlgorithmManager::Instance().create("ExtractSingleSpectrum");
  extractWs1Alg->setChild(true);
  extractWs1Alg->initialize();
  extractWs1Alg->setProperty("InputWorkspace", ws1);
  extractWs1Alg->setProperty("OutputWorkspace", "__ws1_spec");
  extractWs1Alg->setProperty("WorkspaceIndex", wsIndex);
  extractWs1Alg->execute();
  MatrixWorkspace_sptr ws1SpecWs =
      extractWs1Alg->getProperty("OutputWorkspace");

  IAlgorithm_sptr extractWs2Alg =
      AlgorithmManager::Instance().create("ExtractSingleSpectrum");
  extractWs2Alg->setChild(true);
  extractWs2Alg->initialize();
  extractWs2Alg->setProperty("InputWorkspace", ws2);
  extractWs2Alg->setProperty("OutputWorkspace", "__ws2_spec");
  extractWs2Alg->setProperty("WorkspaceIndex", wsIndex);
  extractWs2Alg->execute();
  MatrixWorkspace_sptr ws2SpecWs =
      extractWs2Alg->getProperty("OutputWorkspace");

  // Rebin the second workspace to the first
  // (needed for identical binning for Minus algorithm)
  IAlgorithm_sptr rebinAlg =
      AlgorithmManager::Instance().create("RebinToWorkspace");
  rebinAlg->setChild(true);
  rebinAlg->initialize();
  rebinAlg->setProperty("WorkspaceToRebin", ws2SpecWs);
  rebinAlg->setProperty("WorkspaceToMatch", ws1SpecWs);
  rebinAlg->setProperty("OutputWorkspace", "__ws2_spec_rebin");
  rebinAlg->execute();
  MatrixWorkspace_sptr rebinnedWs2SpecWs =
      rebinAlg->getProperty("OutputWorkspace");

  // Subtract the two extracted spectra
  IAlgorithm_sptr minusAlg = AlgorithmManager::Instance().create("Minus");
  minusAlg->setChild(true);
  minusAlg->initialize();
  minusAlg->setProperty("LHSWorkspace", ws1SpecWs);
  minusAlg->setProperty("RHSWorkspace", rebinnedWs2SpecWs);
  minusAlg->setProperty("OutputWorkspace", "__diff");
  minusAlg->execute();
  MatrixWorkspace_sptr diffWorkspace = minusAlg->getProperty("OutputWorkspace");

  // Add diff workspace to plot window
  m_diffWsName = wsNames[0] + wsNames[1];
  m_view->plotCurve(m_diffWsName,
                    DataComparisonHelper::curveDataFromWs(diffWorkspace, 0));
}

/** Remove all workspaces from the interface
*
*/
void DataComparisonPresenter::removeAllWorkspaces() {

  auto wsNames = m_view->getWorkspaceNames();
  for (const auto &name : wsNames) {
    m_view->removeWorkspace(name);
    m_view->detachCurve(name);
  }
  m_view->detachCurve(m_diffWsName);
  m_diffWsName.clear();
}

/** Remove selected workspaces from the interface
*
*/
void DataComparisonPresenter::removeSelectedWorkspaces() {

  auto wsNames = m_view->getSelectedWorkspaceNames();
  for (const auto &name : wsNames) {
    m_view->removeWorkspace(name);
    m_view->detachCurve(name);
  }

  if (!m_diffWsName.empty()) {
    for (const auto &name : wsNames) {
      if (m_diffWsName.find(name) != std::string::npos) {
        m_view->detachCurve(m_diffWsName);
        m_diffWsName.clear();
        break;
      }
    }
  }
}

/** Remove diff workspace from the interface
*
*/
void DataComparisonPresenter::removeDiffWorkspace() {

  if (m_diffWsName.empty())
    return;

  m_view->detachCurve(m_diffWsName);
  m_diffWsName.clear();
}

/**
* Gets a colour as an index for the combo box for a new workspace.
* Looks for the lowest unused index, if all colours are used then returns 0.
*
* @return An index to set for the conbo box
*/
int DataComparisonPresenter::getInitialColourIndex() {

  auto availableColors = m_view->getAvailableColors();
  auto currentColors = m_view->getWorkspaceColors();

  // Just use the first colour if this is the first row
  if (availableColors.size() <= 1)
    return 0;

  for (size_t i = 0; i < availableColors.size(); i++) {
    if (std::find(currentColors.begin(), currentColors.end(),
                  availableColors[i]) == currentColors.end()) {
      return static_cast<int>(i);
    }
  }

  return 0;
}

/**
* Handles removing a workspace when it is deleted from ADS.
*
* @param wsName Name of the workspace being deleted
* @param ws Pointer to the workspace
*/
void DataComparisonPresenter::preDeleteHandle(
    const std::string &wsName,
    const boost::shared_ptr<Mantid::API::Workspace> ws) {
  UNUSED_ARG(ws);

  m_view->removeWorkspace(wsName);
  m_view->detachCurve(wsName);

  if (!m_diffWsName.empty()) {
    if (m_diffWsName.find(wsName) != std::string::npos) {
      m_view->detachCurve(m_diffWsName);
      m_diffWsName.clear();
    }
  }
}

/**
* Handle a workspace being renamed.
*
* @param oldName Old name for the workspace
* @param newName New name for the workspace
*/
void DataComparisonPresenter::renameHandle(const std::string &oldName,
                                           const std::string &newName) {

  m_view->detachCurve(oldName);
  m_view->removeWorkspace(oldName);

  // Update diff workspace

  m_view->addWorkspace(newName, getInitialColourIndex());
}

/**
* Handle replotting after a workspace has been changed.
*
* @param wsName Name of changed workspace
* @param ws Pointer to changed workspace
*/
void DataComparisonPresenter::afterReplaceHandle(
    const std::string &wsName,
    const boost::shared_ptr<Mantid::API::Workspace> ws) {
  UNUSED_ARG(wsName);
  UNUSED_ARG(ws);

  MatrixWorkspace_const_sptr matrixWs =
      boost::dynamic_pointer_cast<const MatrixWorkspace>(ws);
  if (!matrixWs)
    return;

  // Update the plot
  auto index = m_view->getSelectedWorkspaceIndex();
  auto colors = m_view->getAvailableColors();

  m_view->detachCurve(wsName);
  m_view->plotCurve(wsName,
                    DataComparisonHelper::curveDataFromWs(matrixWs, index),
                    colors[getInitialColourIndex()]);
}
