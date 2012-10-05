
#include <iostream>
#include  "MantidQtImageViewer/ImageView.h"
#include  "MantidQtImageViewer/ColorMaps.h"

#include "ui_ImageView.h"
#include "MantidQtImageViewer/IVConnections.h"
#include "MantidQtImageViewer/ImageDisplay.h"
#include "MantidQtImageViewer/SliderHandler.h"
#include "MantidQtImageViewer/RangeHandler.h"

namespace MantidQt
{
namespace ImageView
{


/**
 *  Construct an ImageView to display data from the specified data source.
 *  The specified ImageDataSource must be constructed elsewhere and passed
 *  into this ImageView constructor.  Most other components of the ImageView
 *  are managed by this class.  That is the graphs, image display and other
 *  parts of the ImageView are constructed here and are deleted when the
 *  ImageView destructor is called.
 *
 *  @param data_source  The source of the data that will be displayed. 
 */
ImageView::ImageView( ImageDataSource* data_source )
{
  Ui_ImageViewer* ui = new Ui_ImageViewer();
  saved_ui          = ui; 

  QMainWindow* window = this;

  ui->setupUi( window );
  window->resize( 1050, 800 );
  window->show();
  window->setAttribute(Qt::WA_DeleteOnClose);  // We just need to close the
                                               // window to trigger the 
                                               // destructor and clean up

  SliderHandler* slider_handler = new SliderHandler( ui );
  saved_slider_handler = slider_handler;

  RangeHandler* range_handler = new RangeHandler( ui );
  saved_range_handler = range_handler;

  h_graph = new GraphDisplay( ui->h_graphPlot, ui->h_graph_table, false );
  v_graph = new GraphDisplay( ui->v_graphPlot, ui->v_graph_table, true );

  ImageDisplay* image_display = new ImageDisplay( ui->imagePlot,
                                                  slider_handler,
                                                  range_handler,
                                                  h_graph, v_graph,
                                                  ui->image_table );
  saved_image_display = image_display;

  IVConnections* iv_connections = new IVConnections( ui, this, 
                                                     image_display, 
                                                     h_graph, v_graph );
  saved_iv_connections = iv_connections;

  image_display->SetDataSource( data_source );
}


ImageView::~ImageView()
{
//  std::cout << "ImageView destructor called" << std::endl;

  delete  h_graph;
  delete  v_graph;

  ImageDisplay* image_display = static_cast<ImageDisplay*>(saved_image_display);
  delete  image_display;

  SliderHandler* slider_handler = 
                             static_cast<SliderHandler*>(saved_slider_handler);
  delete  slider_handler;

  RangeHandler* range_handler = 
                             static_cast<RangeHandler*>(saved_range_handler);
  delete  range_handler;

  IVConnections* iv_connections = 
                             static_cast<IVConnections*>(saved_iv_connections);
  delete  iv_connections;

  Ui_ImageViewer* ui = static_cast<Ui_ImageViewer*>(saved_ui);
  delete  ui;
}


} // namespace MantidQt 
} // namespace ImageView 

