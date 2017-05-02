class PlotOptionsModel(object):
        
    def __init__(self):
    
        self._plot_options = []
        
    def getListPlotOptions(self, index):
        """
        index: current tab index
        """
        plot_options = self._plot_options
        if index == 0:
            plot_options = ['Spectra','Contour', 'Elwin', 'MSDFit']
        elif index == 1:
            plot_options = ['Spectra', 'Contour', 'Moments']
        elif index == 2:
            plot_options = ['Spectra']
        elif index == 3:
            plot_options = []
  
        self._plot_options = plot_options
        
        return self._plot_options
