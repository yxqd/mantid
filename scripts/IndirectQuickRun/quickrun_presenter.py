
class QuickRunPresenter(object):

    def __init__(self, view, PlotOptionsModel):
        self._view = view
        self._plot_option_model = PlotOptionsModel

        # set up slots
        view.plotClicked.connect(self.onPlotClicked)
        view.tabChanged.connect(self.onTabChanged)

    def onPlotClicked(self):
        view = self._view

        plot_option = view.getPlotOption()
        print(plot_option)

        return None

    def onTabChanged(self, index):
        view = self._view
        model = self._plot_option_model
        
        plot_options = model.getListPlotOptions(index)
        view.addPlotOptions(plot_options)
        
        return None
