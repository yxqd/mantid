//----------------------
// Includes
//----------------------
#include "MantidQtCustomInterfaces/DataComparison/DataComparisonView.h"
#include "MantidQtCustomInterfaces/DataComparison/DataComparisonPresenter.h"

#include "MantidAPI/MatrixWorkspace.h"
#include "MantidQtAPI/QwtWorkspaceSpectrumData.h"

namespace {
Mantid::Kernel::Logger g_log("DataComparisonView");
}

// Add this class to the list of specialised dialogs in this namespace
namespace MantidQt {
namespace CustomInterfaces {
DECLARE_SUBWINDOW(DataComparisonView)
}
}

using namespace MantidQt::CustomInterfaces;
using namespace Mantid::API;

//----------------------
// Public member functions
//----------------------
/// Constructor
DataComparisonView::DataComparisonView(QWidget *parent)
    : UserSubWindow(parent), m_plot(new QwtPlot(parent)),
      m_zoomTool(NULL), m_panTool(NULL), m_magnifyTool(NULL) {

  // Create the presenter
  m_presenter.reset(new DataComparisonPresenter(this));
}

/// Set up the dialog layout
void DataComparisonView::initLayout() {
  m_uiForm.setupUi(this);

  m_zoomTool =
      new QwtPlotZoomer(QwtPlot::xBottom, QwtPlot::yLeft,
                        QwtPicker::DragSelection | QwtPicker::CornerToCorner,
                        QwtPicker::AlwaysOff, m_plot->canvas());
  m_zoomTool->setEnabled(false);

  m_panTool = new QwtPlotPanner(m_plot->canvas());
  m_panTool->setEnabled(false);

  m_magnifyTool = new QwtPlotMagnifier(m_plot->canvas());
  m_magnifyTool->setEnabled(false);

  // Add the plot to the UI
  m_plot->setCanvasBackground(Qt::white);
  m_uiForm.loPlot->addWidget(m_plot);

  // Connect push buttons
  connect(m_uiForm.pbAddData, SIGNAL(clicked()), this, SLOT(addDataClicked()));

  connect(m_uiForm.pbRemoveSelectedData, SIGNAL(clicked()), this,
          SLOT(removeSelectedDataClicked()));
  connect(m_uiForm.pbRemoveAllData, SIGNAL(clicked()), this,
          SLOT(removeAllDataClicked()));

  connect(m_uiForm.pbDiffSelected, SIGNAL(clicked()), this,
          SLOT(diffSelectedClicked()));
  connect(m_uiForm.pbClearDiff, SIGNAL(clicked()), this, SLOT(clearDiffClicked()));

  connect(m_uiForm.pbPan, SIGNAL(toggled(bool)), this, SLOT(togglePan(bool)));
  connect(m_uiForm.pbZoom, SIGNAL(toggled(bool)), this, SLOT(toggleZoom(bool)));
  connect(m_uiForm.pbResetView, SIGNAL(clicked()), this, SLOT(resetView()));

  // Replot spectra when the workspace index is changed
  connect(m_uiForm.sbSpectrum, SIGNAL(valueChanged(int)), this,
          SLOT(workspaceIndexChanged()));

  // Add headers to data table
  QStringList headerLabels;
  headerLabels << "Colour"
               << "Workspace";
  m_uiForm.twCurrentData->setColumnCount(headerLabels.size());
  m_uiForm.twCurrentData->setHorizontalHeaderLabels(headerLabels);

  // Select entire rows when a cell is selected
  m_uiForm.twCurrentData->setSelectionBehavior(QAbstractItemView::SelectRows);

  // Stretch last column
  m_uiForm.twCurrentData->horizontalHeader()->setStretchLastSection(true);
}

/** Slot triggered when 'Add Data' is clicked
 *
 */
void DataComparisonView::addDataClicked() {

  m_presenter->notify(IDataComparisonPresenter::AddWorkspace);
}

/** Slot triggered when 'Remove Selected Data' is clicked
*
*/
void DataComparisonView::removeSelectedDataClicked() {

  m_presenter->notify(IDataComparisonPresenter::RemoveSelectedWorkspaces);
}

/** Slot triggered when 'Remove All Data' is clicked
*
*/
void DataComparisonView::removeAllDataClicked() {

  m_presenter->notify(IDataComparisonPresenter::RemoveAllWorkspaces);
}

/** Slot triggered when 'Diff Selected' is clicked
*
*/
void DataComparisonView::diffSelectedClicked() {

	m_presenter->notify(IDataComparisonPresenter::PlotDiffWorkspaces);
}

/** Slot triggered when 'Clear Diff' is clicked
*
*/
void DataComparisonView::clearDiffClicked() {

	m_presenter->notify(IDataComparisonPresenter::RemoveDiffWorkspace);
}

/** Slot triggered when a different color was selected
*
*/
void DataComparisonView::colorChanged() {

	m_presenter->notify(IDataComparisonPresenter::ColorChanged);
}

/** Slot triggered when global workspace index changed
*
*/
void DataComparisonView::workspaceIndexChanged() {

	m_presenter->notify(IDataComparisonPresenter::WorkspaceIndexChanged);

	bool maintainZoom = m_uiForm.cbMaintainZoom->isChecked();
	if (!maintainZoom)
		resetView();
}

/**
 * Adds a MatrixWorkspace by name to the data table.
 *
 * @param wsName :: the name of the workspace to add.
 */
void DataComparisonView::addWorkspace(const std::string &wsName,
                                      int colorIndex) {

  // Append a new row to the data table
  int currentRows = m_uiForm.twCurrentData->rowCount();
  m_uiForm.twCurrentData->insertRow(currentRows);

  // Insert the colour selector
  QComboBox *colourCombo = new QComboBox();
  auto colors = getAvailableColors();
  for (const auto &color : colors)
    colourCombo->addItem(QString::fromStdString(color));
  colourCombo->setCurrentIndex(colorIndex);
  m_uiForm.twCurrentData->setCellWidget(currentRows, COLOUR, colourCombo);

  // Update plots when colour changed
  connect(colourCombo, SIGNAL(currentIndexChanged(int)), this,
          SLOT(colorChanged()));

  // Insert the workspace name
  QTableWidgetItem *wsNameItem = new QTableWidgetItem(tr(wsName.c_str()));
  wsNameItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
  m_uiForm.twCurrentData->setItem(currentRows, WORKSPACE_NAME, wsNameItem);
}

/**
 * Determines if a given workspace is currently shown in the UI.
 *
 * @param wsName :: the name of the workspace to check
 */
bool DataComparisonView::containsWorkspace(const std::string &wsName) const {
  QString testWsName = QString::fromStdString(wsName);

  int numRows = m_uiForm.twCurrentData->rowCount();
  for (int row = 0; row < numRows; row++) {
    QString workspaceName =
        m_uiForm.twCurrentData->item(row, WORKSPACE_NAME)->text();
    if (workspaceName == testWsName)
      return true;
  }

  return false;
}

/** Return available colors
*
* @return :: available colors
*/
std::vector<std::string> DataComparisonView::getAvailableColors() const {

  std::vector<std::string> colors(16);
  colors[0] = "Black";
  colors[1] = "Red";
  colors[2] = "Green";
  colors[3] = "Blue";
  colors[4] = "Cyan";
  colors[5] = "Magenta";
  colors[6] = "Yellow";
  colors[7] = "Light Gray";
  colors[8] = "Gray";
  colors[9] = "Dark Red";
  colors[10] = "Dark Green";
  colors[11] = "Dark Blue";
  colors[12] = "Dark Cyan";
  colors[13] = "Dark Magenta";
  colors[14] = "Dark Yellow";
  colors[15] = "Dark Gray";
  return colors;
}

/** Remove a workspace from the data table
*
* @param wsName :: the workspace name to remove from the table
*/
void DataComparisonView::removeWorkspace(const std::string &wsName) {

  int numRows = m_uiForm.twCurrentData->rowCount();
  for (int row = 0; row < numRows; row++) {
    if (m_uiForm.twCurrentData->item(row, WORKSPACE_NAME)
            ->text()
            .toStdString() == wsName) {
      m_uiForm.twCurrentData->removeRow(row);
      break;
    }
  }
}

/** Return selected workspace index
*
* @return :: global workspace index
*/
int DataComparisonView::getSelectedWorkspaceIndex() const {

  return m_uiForm.sbSpectrum->value();
}

/** Return workspace names as a vector of strings
*
* @return :: workspace names
*/
std::vector<std::string> DataComparisonView::getWorkspaceNames() const {

  std::vector<std::string> wsNames;

  int numRows = m_uiForm.twCurrentData->rowCount();
  for (int row = 0; row < numRows; row++) {
    // Get workspace
    wsNames.push_back(m_uiForm.twCurrentData->item(row, WORKSPACE_NAME)
                          ->text()
                          .toStdString());
  }

  return wsNames;
}

/** Return workspace colors as a vector of strings
*
* @return :: workspace colors
*/
std::vector<std::string> DataComparisonView::getWorkspaceColors() const {

  std::vector<std::string> wsColors;

  int numRows = m_uiForm.twCurrentData->rowCount();

  for (int row = 0; row < numRows; row++) {
    QComboBox *color = dynamic_cast<QComboBox *>(
        m_uiForm.twCurrentData->cellWidget(row, COLOUR));
    wsColors.push_back(color->currentText().toStdString());
  }
  return wsColors;
}

/**
 * Return names of the workspaces currently selected in the data table
 *O
 * @return :: workspace names
 */
std::vector<std::string> DataComparisonView::getSelectedWorkspaceNames() const {

  std::vector<std::string> wsNames;

  auto selectedItems = m_uiForm.twCurrentData->selectedItems();
  std::set<int> selectedRows;

  // Generate a vector of selected workspaces
  for (auto it = selectedItems.begin(); it != selectedItems.end(); ++it) {

    if (!selectedRows.count((*it)->row())) {
      wsNames.push_back(
          m_uiForm.twCurrentData->item((*it)->row(), WORKSPACE_NAME)
              ->text()
              .toStdString());
	  selectedRows.insert((*it)->row());
    }
  }

  return wsNames;
}

/**
 * Toggles the pan plot tool.
 *
 * @param enabled If the tool should be enabled
 */
void DataComparisonView::togglePan(bool enabled) {
  // First disbale the zoom tool
  if (enabled && m_uiForm.pbZoom->isChecked())
    m_uiForm.pbZoom->setChecked(false);

  g_log.debug() << "Pan tool enabled: " << enabled << '\n';

  m_panTool->setEnabled(enabled);
  m_magnifyTool->setEnabled(enabled);
}

/**
 * Toggles the zoom plot tool.
 *
 * @param enabled If the tool should be enabled
 */
void DataComparisonView::toggleZoom(bool enabled) {
  // First disbale the pan tool
  if (enabled && m_uiForm.pbPan->isChecked())
    m_uiForm.pbPan->setChecked(false);

  g_log.debug() << "Zoom tool enabled: " << enabled << '\n';

  m_zoomTool->setEnabled(enabled);
  m_magnifyTool->setEnabled(enabled);
}

/**
 * Rests the zoom level to fit all curves on the plot.
 */
void DataComparisonView::resetView() {
  g_log.debug("Reset plot view");

  // Auto scale the axis
  m_plot->setAxisAutoScale(QwtPlot::xBottom);
  m_plot->setAxisAutoScale(QwtPlot::yLeft);

  // Set this as the default zoom level
  m_zoomTool->setZoomBase(true);
}

/**
* Print error message
*
* @param messge :: the message to be printed
*/
void DataComparisonView::printError(const std::string &message) {

	g_log.error() << message << "\n";
}

/**
* Print information message
*
* @param messge :: the message to be printed
*/
void DataComparisonView::printInformation(const std::string &message) {

	g_log.information() << message << "\n";
}

/**
* Print debug message
*
* @param messge :: the message to be printed
*/
void DataComparisonView::printDebug(const std::string &message) {

	g_log.debug() << message << "\n";
}

/**
* Return data name currently selected in the DataSelector widget
*
* @return :: the data name
*/
std::string DataComparisonView::getSelectedWorkspaceName() const {

  return m_uiForm.dsData->getCurrentDataName().toStdString();
}

/** Block/unblock signals emitted by the table
*
* @param block :: true if signals must be blocked. False otherwise
*/
void DataComparisonView::blockTableSignals(bool block) {

  m_uiForm.twCurrentData->blockSignals(block);
}

/** Detach a workspace from plot widget
*
* @param wsName :: the name of the workspace to detach from plot
*/
void DataComparisonView::detachCurve(const std::string &wsName) {

  auto name = QString::fromStdString(wsName);
  if (m_curves.contains(name)) {
	  m_curves[name]->attach(NULL);
	  m_curves.remove(name);
  }
  m_plot->replot();
}

/** Plot a curve
*
* @param wsName :: the name of the workspace corresponding to the data
* @param index :: the workspace index to plot
* @param color :: the color (green if empty string)
*/
void DataComparisonView::plotCurve(const std::string &wsName,
                                   const QwtArrayData &curve,
                                   const std::string &color) {

  std::string useColor = color.empty() ? "Green" : color;

  // Create a new curve and attach it to the plot
  auto plotCurve = boost::make_shared<QwtPlotCurve>();
  plotCurve->setData(curve);
  plotCurve->setPen(QColor(QString::fromStdString(useColor)));
  plotCurve->attach(m_plot);
  m_curves[QString::fromStdString(wsName)] = plotCurve;
  m_plot->replot();
}
