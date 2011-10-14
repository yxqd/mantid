#ifndef SLICEVIEWER_H
#define SLICEVIEWER_H

#include "DimensionSliceWidget.h"
#include "MantidAPI/IMDWorkspace.h"
#include "QwtRasterDataMD.h"
#include "ui_SliceViewer.h"
#include <QtGui/qdialog.h>
#include <QtGui/QWidget>
#include <qwt_plot_spectrogram.h>
#include <qwt_plot.h>
#include <qwt_raster_data.h>
#include <qwt_scale_widget.h>
#include <vector>
#include <qwt_color_map.h>
#include <QtCore/QtCore>
#include "DllOption.h"

class EXPORT_OPT_MANTIDQT_SLICEVIEWER SliceViewer : public QWidget
{
  Q_OBJECT

public:
  SliceViewer(QWidget *parent = 0);
  ~SliceViewer();

  void setWorkspace(Mantid::API::IMDWorkspace_sptr ws);
  void showControls(bool visible);
  void zoomBy(double factor);

public slots:
  void changedShownDim(int index, int dim, int oldDim);
  void updateDisplaySlot(int index, double value);
  void resetZoom();
  void showInfoAt(double, double);
  void colorRangeFullSlot();
  void colorRangeSliceSlot();
  void zoomInSlot();
  void zoomOutSlot();

private:
  void initMenus();
  void initZoomer();

  void updateDisplay(bool resetAxes = false);
  void updateDimensionSliceWidgets();
  void resetAxis(int axis, Mantid::Geometry::IMDDimension_const_sptr dim);

  void findRangeFull();
  void findRangeSlice();


private:
  /// Auto-generated UI controls.
  Ui::SliceViewerClass ui;

  /// Set to true once the first workspace has been loaded in it
  bool m_firstWorkspaceOpen;

  /// Main plot object
  QwtPlot * m_plot;

  /// Spectrogram plot
  QwtPlotSpectrogram * m_spect;

  /// Layout containing the spectrogram
  QHBoxLayout * m_spectLayout;

  /// Color map in use
  QwtLinearColorMap m_colorMap;

  /// Color bar indicating the color scale
  QwtScaleWidget * m_colorBar;

  /// Vector of the widgets for slicing dimensions
  QVector<DimensionSliceWidget *> m_dimWidgets;

  /// Data presenter
  QwtRasterDataMD * m_data;

  /// Workspace being shown
  Mantid::API::IMDWorkspace_sptr m_ws;

  /// The X and Y dimensions being plotted
  Mantid::Geometry::IMDDimension_const_sptr m_X;
  Mantid::Geometry::IMDDimension_const_sptr m_Y;
  size_t m_dimX;
  size_t m_dimY;

  /// The range of values to fit in the color map.
  QwtDoubleInterval m_colorRange;

  /// The calculated range of values in the FULL data set
  QwtDoubleInterval m_colorRangeFull;

  /// The calculated range of values ONLY in the currently viewed part of the slice
  QwtDoubleInterval m_colorRangeSlice;

  /// Use the log of the value for the color scale
  bool m_logColor;

  /// Menus
  QMenu * m_menuColorOptions;
  QMenu * m_menuView;

};

#endif // SLICEVIEWER_H
