#include "MantidQtCustomInterfaces/IndirectCalibration.h"

#include "MantidKernel/Logger.h"

#include <QFileInfo>

namespace
{
  Mantid::Kernel::Logger g_log("IndirectCalibration");
}

namespace MantidQt
{
namespace CustomInterfaces
{
  //----------------------------------------------------------------------------------------------
  /** Constructor
   */
  IndirectCalibration::IndirectCalibration(Ui::IndirectDataReduction& uiForm, QWidget * parent) :
      IndirectDataReductionTab(uiForm, parent)
  {
    m_propTrees["CalPropTree"] = new QtTreePropertyBrowser();
    m_uiForm.cal_treeCal->addWidget(m_propTrees["CalPropTree"]);

    DoubleEditorFactory *doubleEditorFactory = new DoubleEditorFactory();
    m_propTrees["CalPropTree"]->setFactoryForManager(m_dblManager, doubleEditorFactory);

    m_properties["CalPeakMin"] = m_dblManager->addProperty("Peak Min");
    m_properties["CalPeakMax"] = m_dblManager->addProperty("Peak Max");
    m_properties["CalBackMin"] = m_dblManager->addProperty("Back Min");
    m_properties["CalBackMax"] = m_dblManager->addProperty("Back Max");

    m_propTrees["CalPropTree"]->addProperty(m_properties["CalPeakMin"]);
    m_propTrees["CalPropTree"]->addProperty(m_properties["CalPeakMax"]);
    m_propTrees["CalPropTree"]->addProperty(m_properties["CalBackMin"]);
    m_propTrees["CalPropTree"]->addProperty(m_properties["CalBackMax"]);

    m_plots["CalPlot"] = new QwtPlot(m_parentWidget);
    m_plots["CalPlot"]->setAxisFont(QwtPlot::xBottom, parent->font());
    m_plots["CalPlot"]->setAxisFont(QwtPlot::yLeft, parent->font());
    m_uiForm.cal_plotCal->addWidget(m_plots["CalPlot"]);
    m_plots["CalPlot"]->setCanvasBackground(Qt::white);

    m_rangeSelectors["CalPeak"] = new MantidWidgets::RangeSelector(m_plots["CalPlot"]);
    connect(m_rangeSelectors["CalPeak"], SIGNAL(minValueChanged(double)), this, SLOT(calMinChanged(double)));
    connect(m_rangeSelectors["CalPeak"], SIGNAL(maxValueChanged(double)), this, SLOT(calMaxChanged(double)));
    m_rangeSelectors["CalBackground"] = new MantidWidgets::RangeSelector(m_plots["CalPlot"]);
    m_rangeSelectors["CalBackground"]->setColour(Qt::darkGreen); // dark green to signify background range
    connect(m_rangeSelectors["CalBackground"], SIGNAL(minValueChanged(double)), this, SLOT(calMinChanged(double)));
    connect(m_rangeSelectors["CalBackground"], SIGNAL(maxValueChanged(double)), this, SLOT(calMaxChanged(double)));

    // Res
    m_propTrees["ResPropTree"] = new QtTreePropertyBrowser();
    m_uiForm.cal_treeRes->addWidget(m_propTrees["ResPropTree"]);

    m_propTrees["ResPropTree"]->setFactoryForManager(m_dblManager, doubleEditorFactory);

    // Res - Spectra Selection
    m_properties["ResSpecMin"] = m_dblManager->addProperty("Spectra Min");
    m_properties["ResSpecMax"] = m_dblManager->addProperty("Spectra Max");
    m_propTrees["ResPropTree"]->addProperty(m_properties["ResSpecMin"]);
    m_dblManager->setDecimals(m_properties["ResSpecMin"], 0);
    m_propTrees["ResPropTree"]->addProperty(m_properties["ResSpecMax"]);
    m_dblManager->setDecimals(m_properties["ResSpecMax"], 0);

    // Res - Background Properties
    QtProperty* resBG = m_grpManager->addProperty("Background");
    m_properties["ResStart"] = m_dblManager->addProperty("Start");
    m_properties["ResEnd"] = m_dblManager->addProperty("End");
    resBG->addSubProperty(m_properties["ResStart"]);
    resBG->addSubProperty(m_properties["ResEnd"]);
    m_propTrees["ResPropTree"]->addProperty(resBG);

    // Res - rebinning
    const int NUM_DECIMALS = 3;
    QtProperty* resRB = m_grpManager->addProperty("Rebinning");

    m_properties["ResELow"] = m_dblManager->addProperty("Low");
    m_dblManager->setDecimals(m_properties["ResELow"], NUM_DECIMALS);
    m_dblManager->setValue(m_properties["ResELow"], -0.2);

    m_properties["ResEWidth"] = m_dblManager->addProperty("Width");
    m_dblManager->setDecimals(m_properties["ResEWidth"], NUM_DECIMALS);
    m_dblManager->setValue(m_properties["ResEWidth"], 0.002);
    m_dblManager->setMinimum(m_properties["ResEWidth"], 0.001);

    m_properties["ResEHigh"] = m_dblManager->addProperty("High");
    m_dblManager->setDecimals(m_properties["ResEHigh"], NUM_DECIMALS);
    m_dblManager->setValue(m_properties["ResEHigh"], 0.2);

    resRB->addSubProperty(m_properties["ResELow"]);
    resRB->addSubProperty(m_properties["ResEWidth"]);
    resRB->addSubProperty(m_properties["ResEHigh"]);

    m_propTrees["ResPropTree"]->addProperty(resRB);

    m_plots["ResPlot"] = new QwtPlot(m_parentWidget);
    m_plots["ResPlot"]->setAxisFont(QwtPlot::xBottom, parent->font());
    m_plots["ResPlot"]->setAxisFont(QwtPlot::yLeft, parent->font());
    m_uiForm.cal_plotRes->addWidget(m_plots["ResPlot"]);
    m_plots["ResPlot"]->setCanvasBackground(Qt::white);

    // Create ResR2 first so ResR1 is drawn above it.
    m_rangeSelectors["ResBackground"] = new MantidWidgets::RangeSelector(m_plots["ResPlot"], 
        MantidQt::MantidWidgets::RangeSelector::XMINMAX, true, false);
    m_rangeSelectors["ResBackground"]->setColour(Qt::darkGreen);
    m_rangeSelectors["ResPeak"] = new MantidWidgets::RangeSelector(m_plots["ResPlot"], MantidQt::MantidWidgets::RangeSelector::XMINMAX, true, true);
      
    m_uiForm.cal_leIntensityScaleMultiplier->setValidator(m_valDbl);
    m_uiForm.cal_leResScale->setValidator(m_valDbl);

    m_uiForm.cal_valIntensityScaleMultiplier->setText(" ");

    connect(m_rangeSelectors["ResPeak"], SIGNAL(minValueChanged(double)), this, SLOT(calMinChanged(double)));
    connect(m_rangeSelectors["ResPeak"], SIGNAL(maxValueChanged(double)), this, SLOT(calMaxChanged(double)));
    connect(m_rangeSelectors["ResBackground"], SIGNAL(minValueChanged(double)), this, SLOT(calMinChanged(double)));
    connect(m_rangeSelectors["ResBackground"], SIGNAL(maxValueChanged(double)), this, SLOT(calMaxChanged(double)));
    connect(m_rangeSelectors["ResPeak"], SIGNAL(rangeChanged(double, double)), m_rangeSelectors["ResBackground"], SLOT(setRange(double, double)));
    connect(m_dblManager, SIGNAL(valueChanged(QtProperty*, double)), this, SLOT(calUpdateRS(QtProperty*, double)));
    connect(m_dblManager, SIGNAL(valueChanged(QtProperty*, double)), this, SLOT(calUpdateRS(QtProperty*, double)));
    connect(m_dblManager, SIGNAL(valueChanged(QtProperty*, double)), this, SLOT(calUpdateRS(QtProperty*, double)));

    connect(m_uiForm.cal_leRunNo, SIGNAL(filesFound()), this, SLOT(calPlotRaw()));
    connect(m_uiForm.cal_pbPlot, SIGNAL(clicked()), this, SLOT(calPlotRaw()));
    connect(m_uiForm.cal_ckRES, SIGNAL(toggled(bool)), this, SLOT(resCheck(bool)));
    connect(m_uiForm.cal_ckRES, SIGNAL(toggled(bool)), m_uiForm.cal_ckResScale, SLOT(setEnabled(bool)));
    connect(m_uiForm.cal_ckResScale, SIGNAL(toggled(bool)), m_uiForm.cal_leResScale, SLOT(setEnabled(bool)));
    connect(m_uiForm.cal_ckIntensityScaleMultiplier, SIGNAL(toggled(bool)), this, SLOT(intensityScaleMultiplierCheck(bool)));
    connect(m_uiForm.cal_leIntensityScaleMultiplier, SIGNAL(textChanged(const QString &)), this, SLOT(calibValidateIntensity(const QString &)));
  }
    
  //----------------------------------------------------------------------------------------------
  /** Destructor
   */
  IndirectCalibration::~IndirectCalibration()
  {
  }
  
  void IndirectCalibration::setup()
  {
  }

  void IndirectCalibration::run()
  {
    QString file = m_uiForm.cal_leRunNo->getFirstFilename();
    QString filenames = "[r'"+m_uiForm.cal_leRunNo->getFilenames().join("', r'")+"']";

    QString reducer = "from mantid.simpleapi import SaveNexus\n"
      "from inelastic_indirect_reduction_steps import CreateCalibrationWorkspace\n"
      "calib = CreateCalibrationWorkspace()\n"
      "calib.set_files(" + filenames + ")\n"
      "calib.set_detector_range(" + m_uiForm.leSpectraMin->text() + "-1, " + m_uiForm.leSpectraMax->text() + "-1)\n"
      "calib.set_parameters(" + m_properties["CalBackMin"]->valueText() + "," 
      + m_properties["CalBackMax"]->valueText() + ","
      + m_properties["CalPeakMin"]->valueText() + ","
      + m_properties["CalPeakMax"]->valueText() + ")\n"
      "calib.set_analyser('" + m_uiForm.cbAnalyser->currentText() + "')\n"
      "calib.set_reflection('" + m_uiForm.cbReflection->currentText() + "')\n";

    //Scale values by arbitrary scalar if requested
    if(m_uiForm.cal_ckIntensityScaleMultiplier->isChecked())
    {
      QString scale = m_uiForm.cal_leIntensityScaleMultiplier->text(); 
      if(scale.isEmpty())
      {
          scale = "1.0";
        }
        reducer += "calib.set_intensity_scale("+scale+")\n";
      }

      reducer += "calib.execute(None, None)\n"
        "result = calib.result_workspace()\n"
        "print result\n";

      if( m_uiForm.cal_ckSave->isChecked() )
      {
        reducer +=
          "SaveNexus(InputWorkspace=result, Filename=result+'.nxs')\n";
      }

      if ( m_uiForm.cal_ckPlotResult->isChecked() )
      {
        reducer += "from mantidplot import plotTimeBin\n"
          "plotTimeBin(result, 0)\n";
      }

      QString pyOutput = m_pythonRunner.runPythonCode(reducer).trimmed();

      if ( pyOutput == "" )
      {
        emit showMessageBox("An error occurred creating the calib file.\n");
      }
      else
      {
        if ( m_uiForm.cal_ckRES->isChecked() )
        {
          createRESfile(filenames);
        }
        m_uiForm.ind_calibFile->setFileTextWithSearch(pyOutput + ".nxs");
        m_uiForm.ckUseCalib->setChecked(true);
      }
  }

  bool IndirectCalibration::validate()
  {
    MantidQt::CustomInterfaces::UserInputValidator uiv;

    uiv.checkMWRunFilesIsValid("Run", m_uiForm.cal_leRunNo);

    auto peakRange = std::make_pair(m_dblManager->value(m_properties["CalPeakMin"]), m_dblManager->value(m_properties["CalPeakMax"]));
    auto backRange = std::make_pair(m_dblManager->value(m_properties["CalBackMin"]), m_dblManager->value(m_properties["CalBackMax"]));

    uiv.checkValidRange("Peak Range", peakRange);
    uiv.checkValidRange("Back Range", backRange);
    uiv.checkRangesDontOverlap(peakRange, backRange);

    if ( m_uiForm.cal_ckRES->isChecked() )
    {
      auto backgroundRange = std::make_pair(m_dblManager->value(m_properties["ResStart"]), m_dblManager->value(m_properties["ResEnd"]));
      uiv.checkValidRange("Background", backgroundRange);

      double eLow   = m_dblManager->value(m_properties["ResELow"]);
      double eHigh  = m_dblManager->value(m_properties["ResEHigh"]);
      double eWidth = m_dblManager->value(m_properties["ResEWidth"]);

      uiv.checkBins(eLow, eWidth, eHigh);
    }

    if( m_uiForm.cal_ckIntensityScaleMultiplier->isChecked()
        && m_uiForm.cal_leIntensityScaleMultiplier->text().isEmpty() )
    {
      uiv.addErrorMessage("You must enter a scale for the calibration file");
    }

    if( m_uiForm.cal_ckResScale->isChecked() && m_uiForm.cal_leResScale->text().isEmpty() )
    {
      uiv.addErrorMessage("You must enter a scale for the resolution file");
    }

    QString error = uiv.generateErrorMessage();

    if(error != "")
      g_log.warning(error.toStdString());

    return (error == "");
  }

  void IndirectCalibration::calPlotRaw()
  {
    QString filename = m_uiForm.cal_leRunNo->getFirstFilename();

    if ( filename.isEmpty() )
    {
      return;
    }

    QFileInfo fi(filename);
    QString wsname = fi.baseName();
    QString pyInput = "Load(Filename=r'" + filename + "', OutputWorkspace='" + wsname + "', SpectrumMin="
      + m_uiForm.leSpectraMin->text() + ", SpectrumMax="
      + m_uiForm.leSpectraMax->text() + ")\n";

    pyInput = "try:\n  " +
      pyInput +
      "except ValueError as ve:" +
      "  print str(ve)";

    QString pyOutput = m_pythonRunner.runPythonCode(pyInput);

    if(!pyOutput.isEmpty())
    {
      emit showMessageBox("Unable to load file.  Error: \n\n" + pyOutput + "\nCheck whether your file exists and matches the selected instrument in the Energy Transfer tab.");
      return;
    }

    plotMiniPlot(wsname, 0, "CalPlot", "CalCurve");

    //Also replot the energy
    calPlotEnergy();
  }

  void IndirectCalibration::calPlotEnergy()
  {
    if ( ! m_uiForm.cal_leRunNo->isValid() )
    {
      emit showMessageBox("Run number not valid.");
      return;
    }

    QString files = "[r'" + m_uiForm.cal_leRunNo->getFilenames().join("', r'") + "']";
    QString pyInput =
      "from IndirectEnergyConversion import resolution\n"
      "iconOpt = { 'first': " +QString::number(m_dblManager->value(m_properties["ResSpecMin"]))+
      ", 'last': " +QString::number(m_dblManager->value(m_properties["ResSpecMax"]))+ "}\n"
      "instrument = '" + m_uiForm.cbInst->currentText() + "'\n"
      "analyser = '" + m_uiForm.cbAnalyser->currentText() + "'\n"
      "reflection = '" + m_uiForm.cbReflection->currentText() + "'\n"
      "files = " + files + "\n"
      "outWS = resolution(files, iconOpt, '', '', instrument, analyser, reflection, Res=False)\n"
      "print outWS\n";
    QString pyOutput = m_pythonRunner.runPythonCode(pyInput).trimmed();

    //Something went wrong in the Python
    if(pyOutput == "None")
    {
      emit showMessageBox("Failed to convert to energy. See log for details.");
      return;
    }

    plotMiniPlot(pyOutput, 0, "ResPlot", "ResCurve");;
  }

  void IndirectCalibration::calSetDefaultResolution(Mantid::API::MatrixWorkspace_const_sptr ws)
  {
    auto inst = ws->getInstrument();
    auto analyser = inst->getStringParameter("analyser");

    if(analyser.size() > 0)
    {
      auto comp = inst->getComponentByName(analyser[0]);
      auto params = comp->getNumberParameter("resolution", true);

      //set the default instrument resolution
      if(params.size() > 0)
      {
        double res = params[0];
        m_dblManager->setValue(m_properties["ResELow"], -res*10);
        m_dblManager->setValue(m_properties["ResEHigh"], res*10);

        m_dblManager->setValue(m_properties["ResStart"], -res*9);
        m_dblManager->setValue(m_properties["ResEnd"], -res*8);
      }
    }
  }

  void IndirectCalibration::calMinChanged(double val)
  {
    MantidWidgets::RangeSelector* from = qobject_cast<MantidWidgets::RangeSelector*>(sender());
    if ( from == m_rangeSelectors["CalPeak"] )
    {
      m_dblManager->setValue(m_properties["CalPeakMin"], val);
    }
    else if ( from == m_rangeSelectors["CalBackground"] )
    {
      m_dblManager->setValue(m_properties["CalBackMin"], val);
    }
    else if ( from == m_rangeSelectors["ResPeak"] )
    {
      m_dblManager->setValue(m_properties["ResELow"], val);
    }
    else if ( from == m_rangeSelectors["ResBackground"] )
    {
      m_dblManager->setValue(m_properties["ResStart"], val);
    }
  }

  void IndirectCalibration::calMaxChanged(double val)
  {
    MantidWidgets::RangeSelector* from = qobject_cast<MantidWidgets::RangeSelector*>(sender());
    if ( from == m_rangeSelectors["CalPeak"] )
    {
      m_dblManager->setValue(m_properties["CalPeakMax"], val);
    }
    else if ( from == m_rangeSelectors["CalBackground"] )
    {
      m_dblManager->setValue(m_properties["CalBackMax"], val);
    }
    else if ( from == m_rangeSelectors["ResPeak"] )
    {
      m_dblManager->setValue(m_properties["ResEHigh"], val);
    }
    else if ( from == m_rangeSelectors["ResBackground"] )
    {
      m_dblManager->setValue(m_properties["ResEnd"], val);
    }
  }

  void IndirectCalibration::calUpdateRS(QtProperty* prop, double val)
  {
    if ( prop == m_properties["CalPeakMin"] ) m_rangeSelectors["CalPeak"]->setMinimum(val);
    else if ( prop == m_properties["CalPeakMax"] ) m_rangeSelectors["CalPeak"]->setMaximum(val);
    else if ( prop == m_properties["CalBackMin"] ) m_rangeSelectors["CalBackground"]->setMinimum(val);
    else if ( prop == m_properties["CalBackMax"] ) m_rangeSelectors["CalBackground"]->setMaximum(val);
    else if ( prop == m_properties["ResStart"] ) m_rangeSelectors["ResBackground"]->setMinimum(val);
    else if ( prop == m_properties["ResEnd"] ) m_rangeSelectors["ResBackground"]->setMaximum(val);
    else if ( prop == m_properties["ResELow"] ) m_rangeSelectors["ResPeak"]->setMinimum(val);
    else if ( prop == m_properties["ResEHigh"] ) m_rangeSelectors["ResPeak"]->setMaximum(val);
  }

  /**
  * This function enables/disables the display of the options involved in creating the RES file.
  * @param state :: whether checkbox is checked or unchecked
  */
  void IndirectCalibration::resCheck(bool state)
  {
    m_rangeSelectors["ResPeak"]->setVisible(state);
    m_rangeSelectors["ResBackground"]->setVisible(state);
  }

  /**
   * This function is called after calib has run and creates a RES file for use in later analysis (Fury,etc)
   * @param file :: the input file (WBV run.raw)
   */
  void IndirectCalibration::createRESfile(const QString& file)
  {
    QString scaleFactor("1.0");
    if(m_uiForm.cal_ckResScale->isChecked())
    {
      if(!m_uiForm.cal_leResScale->text().isEmpty())
      {
        scaleFactor = m_uiForm.cal_leResScale->text();
      }
    }

    QString pyInput =
      "from IndirectEnergyConversion import resolution\n"
      "iconOpt = { 'first': " +QString::number(m_dblManager->value(m_properties["ResSpecMin"]))+
      ", 'last': " +QString::number(m_dblManager->value(m_properties["ResSpecMax"]))+"}\n"

      "instrument = '" + m_uiForm.cbInst->currentText() + "'\n"
      "analyser = '" + m_uiForm.cbAnalyser->currentText() + "'\n"
      "reflection = '" + m_uiForm.cbReflection->currentText() + "'\n";

    if ( m_uiForm.cal_ckPlotResult->isChecked() ) { pyInput +=	"plot = True\n"; }
    else { pyInput += "plot = False\n"; }

    if ( m_uiForm.cal_ckVerbose->isChecked() ) { pyInput +=  "verbose = True\n"; }
    else { pyInput += "verbose = False\n"; }

    if ( m_uiForm.cal_ckSave->isChecked() ) { pyInput +=  "save = True\n"; }
    else { pyInput += "save = False\n"; }

    QString rebinParam = QString::number(m_dblManager->value(m_properties["ResELow"])) + "," +
      QString::number(m_dblManager->value(m_properties["ResEWidth"])) + "," +
      QString::number(m_dblManager->value(m_properties["ResEHigh"]));

    QString background = "[ " +QString::number(m_dblManager->value(m_properties["ResStart"]))+ ", " +QString::number(m_dblManager->value(m_properties["ResEnd"]))+"]";

    QString scaled = m_uiForm.cal_ckIntensityScaleMultiplier->isChecked() ? "True" : "False";
    pyInput +=
      "background = " + background + "\n"
      "rebinParam = '" + rebinParam + "'\n"
      "file = " + file + "\n"
      "ws = resolution(file, iconOpt, rebinParam, background, instrument, analyser, reflection, Verbose=verbose, Plot=plot, Save=save, factor="+scaleFactor+")\n"
      "scaled = "+ scaled +"\n"
      "scaleFactor = "+m_uiForm.cal_leIntensityScaleMultiplier->text()+"\n"
      "backStart = "+QString::number(m_dblManager->value(m_properties["CalBackMin"]))+"\n"
      "backEnd = "+QString::number(m_dblManager->value(m_properties["CalBackMax"]))+"\n"
      "rebinLow = "+QString::number(m_dblManager->value(m_properties["ResELow"]))+"\n"
      "rebinWidth = "+QString::number(m_dblManager->value(m_properties["ResEWidth"]))+"\n"
      "rebinHigh = "+QString::number(m_dblManager->value(m_properties["ResEHigh"]))+"\n"
      "AddSampleLog(Workspace=ws, LogName='scale', LogType='String', LogText=str(scaled))\n"
      "if scaled:"
      "  AddSampleLog(Workspace=ws, LogName='scale_factor', LogType='Number', LogText=str(scaleFactor))\n"
      "AddSampleLog(Workspace=ws, LogName='back_start', LogType='Number', LogText=str(backStart))\n"
      "AddSampleLog(Workspace=ws, LogName='back_end', LogType='Number', LogText=str(backEnd))\n"
      "AddSampleLog(Workspace=ws, LogName='rebin_low', LogType='Number', LogText=str(rebinLow))\n"
      "AddSampleLog(Workspace=ws, LogName='rebin_width', LogType='Number', LogText=str(rebinWidth))\n"
      "AddSampleLog(Workspace=ws, LogName='rebin_high', LogType='Number', LogText=str(rebinHigh))\n";

    QString pyOutput = m_pythonRunner.runPythonCode(pyInput).trimmed();

    if ( pyOutput != "" )
    {
      emit showMessageBox("Unable to create RES file: \n" + pyOutput);
    }
  }

  void IndirectCalibration::intensityScaleMultiplierCheck(bool state)
  {
    m_uiForm.cal_leIntensityScaleMultiplier->setEnabled(state);
  }

  void IndirectCalibration::calibValidateIntensity(const QString & text)
  {
    if(!text.isEmpty())
    {
      m_uiForm.cal_valIntensityScaleMultiplier->setText(" ");
    }
    else
    {
      m_uiForm.cal_valIntensityScaleMultiplier->setText("*");
    }
  }

  /**
   * Update values for plot bounds when the instument/analyser/reflection is changed in
   * the Convert to Energy tab
   *
   * @param values :: Map of new plot data
   */
  void IndirectCalibration::newPlotValues(QMap<QString, double> &values)
  {
    m_dblManager->setValue(m_properties["ResSpecMin"], values["SpecMin"]);
    m_dblManager->setValue(m_properties["ResSpecMax"], values["SpecMax"]);
    m_dblManager->setValue(m_properties["CalPeakMin"], values["PeakMin"]);
    m_dblManager->setValue(m_properties["CalPeakMax"], values["PeakMax"]);
    m_dblManager->setValue(m_properties["CalBackMin"], values["BackMin"]);
    m_dblManager->setValue(m_properties["CalBackMax"], values["BackMax"]);
  }

} // namespace CustomInterfaces
} // namespace Mantid
