from __future__ import (absolute_import, division, print_function)


import mantid.simpleapi as mantid

from Muon.GUI.Common import thread_model


class DummyCreatePresenter(object):

    """
    This class links the FFT model to the GUI
    """

    def __init__(self, view, alg,reporter):
        self.view = view
        self.alg = alg
        self.thread = None
        self.reporter = reporter
        # connect
        self.view.buttonSignal.connect(self.handleButton)

    @property
    def widget(self):
        return self.view

    def createThread(self):
        return thread_model.ThreadModel(self.alg, self.reporter)

    # constructs the inputs for the FFT algorithms
    # then executes them (see fft_model to see the order
    # of execution
    def handleButton(self):
        # put this on its own thread so not to freeze Mantid
        self.thread = self.createThread()
        self.thread.threadWrapperSetUp(self.deactivate,self.handleFinished)

        # make some inputs
        inputs = {}
        inputs["DataX"] = self.view.getDataX()
        inputs["DataY"] = self.view.getDataY()
        inputs["OutputWorkspace"] = self.view.getWS()
        self.thread.setInputs(inputs,None)
        self.thread.start()

    # kills the thread at end of execution
    def handleFinished(self):
        self.thread.threadWrapperTearDown(self.deactivate,self.handleFinished)
        self.thread.deleteLater()
        self.thread = None

    def deactivate(self):
        return
