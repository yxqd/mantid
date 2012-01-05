//----------------------
// Includes
//----------------------

#include "MantidQtCustomInterfaces/CreateMDWorkspace.h"
#include "MantidQtCustomInterfaces/WorkspaceMemento.h"
#include "MantidQtCustomInterfaces/WorkspaceInADS.h"
#include "MantidQtCustomInterfaces/WorkspaceOnDisk.h"
#include "MantidQtCustomInterfaces/WorkspaceMemento.h"

#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/AnalysisDataService.h"


// Suppress a warning coming out of code that isn't ours
#if defined(__INTEL_COMPILER)
  #pragma warning disable 1125
#elif defined(__GNUC__)
  #if (__GNUC__ >= 4 && __GNUC_MINOR__ >= 6 )
    #pragma GCC diagnostic push
  #endif
  #pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
#include "qttreepropertybrowser.h"
#include "qtpropertymanager.h"
#include "qteditorfactory.h"
#include "DoubleEditorFactory.h"
#if defined(__INTEL_COMPILER)
  #pragma warning enable 1125
#elif defined(__GNUC__)
  #if (__GNUC__ >= 4 && __GNUC_MINOR__ >= 6 )
    #pragma GCC diagnostic pop
  #endif
#endif

#include <QtCheckBoxFactory>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QRadioButton>
#include <strstream>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

using namespace Mantid::API;

namespace MantidQt
{
namespace CustomInterfaces
{

//Add this class to the list of specialised dialogs in this namespace
DECLARE_SUBWINDOW(CreateMDWorkspace); 

  /**
  Helper type to perform comparisons between WorkspaceMementos
  */
  class IdComparitor : public std::unary_function<WorkspaceMemento_sptr,bool>
  {
  private:
    WorkspaceMemento_sptr m_benchmark;
  public:
    IdComparitor(WorkspaceMemento_sptr benchmark) : m_benchmark(benchmark){}
    bool operator()(WorkspaceMemento_sptr a) const
    {
      std::vector<std::string> strs;
      std::string id =  m_benchmark->getId();
      boost::split(strs, id, boost::is_any_of("/,\\"));

      std::stringstream streamPattern;
      streamPattern << "(" << strs.back() << ")$";
      boost::regex pattern(streamPattern.str(), boost::regex_constants::icase); 

      return boost::regex_search(a->getId(), pattern);

    }
  };

/*
Constructor taking a WorkspaceMementoCollection, which acts as the model.
*/
CreateMDWorkspace::CreateMDWorkspace(QWidget *) 
{
  //Generate memento view model.
  m_model = new QtWorkspaceMementoModel(m_data);
}

/*
Initalize the layout.
*/
void CreateMDWorkspace::initLayout()
{
  m_uiForm.setupUi(this);
  connect(m_uiForm.btn_create, SIGNAL(clicked()), this, SLOT(createMDWorkspaceClicked()));
  connect(m_uiForm.btn_add_workspace, SIGNAL(clicked()), this, SLOT(addWorkspaceClicked()));
  connect(m_uiForm.btn_add_file, SIGNAL(clicked()), this, SLOT(addFileClicked()));
  connect(m_uiForm.btn_remove_workspace, SIGNAL(clicked()), this, SLOT(removeSelectedClicked()));
  connect(m_uiForm.btn_set_ub_matrix, SIGNAL(clicked()), this, SLOT(setUBMatrixClicked()));
  connect(m_uiForm.btn_find_ub_matrix, SIGNAL(clicked()), this, SLOT(findUBMatrixClicked()));
  connect(m_uiForm.btn_create, SIGNAL(clicked()), this, SLOT(createMDWorkspaceClicked()));
  connect(m_uiForm.btn_set_goniometer, SIGNAL(clicked()), this, SLOT(setGoniometerClicked()));
  //Set MVC Model
  m_uiForm.tableView->setModel(m_model);
}

/*
Event handler for finding the UBMatrix as part of a peak finding action.
*/
void CreateMDWorkspace::findUBMatrixClicked()
{
  QString command, args, result;
  QStringList bMatrixArgs;

  try
  {
    WorkspaceMemento_sptr memento = getFirstSelected();
    memento->fetchIt();

    // Find the peaks workspace in detector space
    command = "from mantidsimple import *\n"
      "import sys\n"
      "try:\n"
      "    FindSXPeaksDialog(InputWorkspace='%1', OutputWorkspace='%1_peaks')\n"
      "    print 'SUCCESS'\n"
      "except:\n"
      "    print 'FAIL'";
    args = QString(memento->getId().c_str());
    command = command.arg(args);
    result = runPythonCode(command).trimmed();
    if(result == "FAIL")
    {
      runConfirmation("Aborted during PeakFinding.");
      return;
    }

    // Calculate the u matrix and then copy the ub matrix result from the peaksworkspace to the matrix workspace
    command = "try:\n"
      "    alg = CalculateUMatrixDialog(PeaksWorkspace='%1_peaks')\n"
      "    a = alg.getProperty('a').value\n"
      "    b = alg.getProperty('b')\n"
      "    c = alg.getProperty('c')\n"
      "    alpha = alg.getProperty('alpha')\n"
      "    beta = alg.getProperty('beta')\n"
      "    gamma = alg.getProperty('gamma')\n"
      "    CopySample(InputWorkspace='%1_peaks',OutputWorkspace='%1',CopyName='0',CopyMaterial='0',CopyEnvironment='0',CopyShape='0',CopyLattice='1')\n"
      "    print '%(a)s, %(b)s, %(c)s, %(alpha)s, %(beta)s, %(gamma)s' % {'a': a, 'b' : b, 'c' : c, 'alpha' : alpha, 'beta' : beta, 'gamma' : gamma}\n"
      "except:\n"
      "    print 'FAIL'";

    command = command.arg(args);
    result = runPythonCode(command).trimmed();
    if(result == "FAIL")
    {
      runConfirmation("Aborted during calculating and copying the UB matrix.");
      return;
    }
    else
    {
      bMatrixArgs = result.split(',');
    }

    // Index the peaks on the workspace
    if(m_uiForm.ckIndexSXPeaks->isChecked())
    {
      command = "IndexSXPeaksDialog(PeaksWorkspace='%1_peaks', a='%2', b='%3', c='%4', alpha='%5', beta='%6', gamma='%7')";
      args = args + "," + result;
    }
    else
    {
      command = "IndexPeaksDialog(PeaksWorkspace='%1_peaks')";
    }

    command = command.arg(args);
    // Run peak indexing
    runPythonCode(command).trimmed();
  }
  catch(std::invalid_argument& ex)
  {
    runConfirmation(ex.what());
  }
}

/*
Getter for the first selected memento.
@return the first selected memento
@throw invalid argument if nothing is selected
*/
WorkspaceMemento_sptr CreateMDWorkspace::getFirstSelected()
{
  QTableView* view = m_uiForm.tableView;
  QModelIndexList indexes = view->selectionModel()->selection().indexes();
  if(indexes.size() > 0)
  {
    int index = indexes.front().row();
    return m_data[index];
  }
  else
  {
    throw std::invalid_argument("Nothing selected");
  }
}

/*
Event handler for setting the UB Matrix
*/
void CreateMDWorkspace::setUBMatrixClicked()
{
  try
  {
    WorkspaceMemento_sptr memento = getFirstSelected();
    Mantid::API::MatrixWorkspace_sptr ws = boost::dynamic_pointer_cast<MatrixWorkspace>( memento->fetchIt() );
    QString id = QString(memento->getId().c_str());

    QString pyInput =
      "from mantidsimple import *\n"
      "import sys\n"
      "wsName='%1'\n"
      "SetUBDialog(Workspace=wsName)\n"
      "msg = 'ws is: ' + wsName\n"
      "mtd.sendLogMessage(msg)\n"
      "ws = mtd[wsName]\n"
      "lattice = ws.getSample().getOrientedLattice()\n"
      "ub = lattice.getUB()\n"
      "print '%(u00)d, %(u01)d, %(u02)d, %(u10)d, %(u11)d, %(u12)d, %(u20)d, %(u21)d, %(u22)d' "
      "% {'u00': ub[0][0], 'u01' : ub[0][1], 'u02' : ub[0][2], 'u10': ub[1][0], 'u11' : ub[1][1], 'u12' : ub[1][2], 'u20' : ub[2][0], 'u21' : ub[2][1], 'u22' : ub[2][2]}\n";

    pyInput = pyInput.arg(id);
    QString pyOutput = runPythonCode(pyInput).trimmed();

    QStringList ub = pyOutput.split(',');
    if (ub.size() == 9 )
    {
      memento->setUB(ub[0].toDouble(), ub[1].toDouble(), ub[2].toDouble(), ub[3].toDouble(), ub[4].toDouble(), ub[5].toDouble(), ub[6].toDouble(),  ub[7].toDouble(),  ub[8].toDouble());
      memento->cleanUp(); 
      m_model->update();
    }
  }
  catch(std::invalid_argument& ex)
  {
    runConfirmation(ex.what());
  }
}

/*
Add a workspace from the ADS
*/
void CreateMDWorkspace::addWorkspaceClicked()
{
  std::string name = m_uiForm.workspaceSelector->currentText().toStdString();
  if(!name.empty())
  {
    WorkspaceMemento_sptr candidate(new WorkspaceInADS(name));
    addUniqueMemento(candidate);
  }
}

/**
Adds a memento to the existing data, if it's id has not already been used in the data list.
@param candidate : candidate memento to add.
*/
void CreateMDWorkspace::addUniqueMemento(WorkspaceMemento_sptr candidate)
{
  IdComparitor comparitor(candidate);
  WorkspaceMementoCollection::iterator pos = std::find_if(m_data.begin(), m_data.end(), comparitor);
  if(pos == m_data.end())
  {
    m_data.push_back(candidate);
    m_model->update();
  }
  else
  {
    std::string msg = "Already have a workspace by that name loaded. Cannot add another.";
    runConfirmation(msg);
  }
}

/*
Add a raw file from the ADS
*/
void CreateMDWorkspace::addFileClicked()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
    "/home",
    tr("Raw Files (*.raw)"));

  std::string name = fileName.toStdString();
  if(!name.empty())
  {
    try
    {
      WorkspaceMemento_sptr candidate(new WorkspaceOnDisk(name));
      addUniqueMemento(candidate);
    }
    catch(std::invalid_argument& arg)
    {
      this->runConfirmation(arg.what());
    }
  }
}

/**
Handler for setting the goniometer.
*/
void CreateMDWorkspace::setGoniometerClicked()
{
  try
  {
    WorkspaceMemento_sptr memento = getFirstSelected();
    if(memento->locationType() != WorkspaceInADS::locType())
    {
      runConfirmation("Currently, Goniometer settings may only be applied to Workspace in memory");
      return;
    }
    Mantid::API::MatrixWorkspace_sptr ws = boost::dynamic_pointer_cast<MatrixWorkspace>( memento->fetchIt() );
    QString id = QString(memento->getId().c_str());

    QString pyInput =
      "from mantidsimple import *\n"
      "import sys\n"
      "try:\n"
      "    wsName='%1'\n"
      "    SetGoniometerDialog(Workspace=wsName)\n"
      "    print 'SUCCESS'\n"
      "except:\n"
      "    print 'FAIL'\n";

    pyInput = pyInput.arg(id);
    QString pyOutput = runPythonCode(pyInput).trimmed();

  }
  catch(std::invalid_argument& ex)
  {
    runConfirmation(ex.what());
  }
}


/*
Remove any selected workspace mementos
*/
void CreateMDWorkspace::removeSelectedClicked()
{
  QTableView* view = m_uiForm.tableView;
  QModelIndexList indexes = view->selectionModel()->selection().indexes();

  //Copy the original collection
  WorkspaceMementoCollection temp(m_data.size());
  std::copy(m_data.begin(), m_data.end(), temp.begin());

  for (int i = 0; i < indexes.count(); i++)
  {
    QModelIndex index = indexes.at(i);
    int row = index.row();
    //Copy and swap trick.
    WorkspaceMemento_sptr dead = m_data[row];
    WorkspaceMementoCollection::iterator pos = std::find(temp.begin(), temp.end(), dead);
    if(pos != temp.end())
    {
      temp.erase(pos);
    }
  }
  //Assign using the copied/trimmed temp collection.
  m_data = temp;
  //Update the model.
  m_model->update();
}

void CreateMDWorkspace::initLocalPython()
{
}

/*
Run a generic confirmation dialog.
@param message : The message to display.
*/
int CreateMDWorkspace::runConfirmation(const std::string& message)
{
  QMessageBox msgBox;
  msgBox.setText(message.c_str());
  msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
  msgBox.setDefaultButton(QMessageBox::Cancel);
  return msgBox.exec();
}

void CreateMDWorkspace::createMDWorkspaceClicked()
{
  //Launch dialog
}


/// Destructor
CreateMDWorkspace::~CreateMDWorkspace()
{
}

} //namespace CustomInterfaces
} //namespace MantidQt
