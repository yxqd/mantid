from __future__ import (absolute_import, division, print_function)

from six import iteritems

import mantid.simpleapi as mantid


class DummyCreateWrapper(object):

    """
    A class to wrap the different parts
    of the FFT and its preprocessing.
    This keeps the main FFT class simple.
    """

    def __init__(self, FFT):
        self.name = "create"
        self.model = FFT

    def cancel(self):
        self.model.cancel()

    def setInputs(self, inputs,runName=None):
        """
        store the data in the wrapper for later
        """
        self.data = inputs

    def execute(self):
        """
        runs the relevant parts of the FFT and the preprocessing
        """
        self.model.make(self.data)

    def output(self):
        return

    def report(self):
        return self.model.report(self.data)


class DummyCreateModel(object):

    """
    A simple class which executes
    the relevant algorithms for
    the analysis.
    """

    def __init__(self):
        self.name = "create"
        self.alg = None

    def cancel(self):
        if self.alg is not None:
            self.alg.cancel()

    def make(self, inputs):
        """
        generates a phase table from CalMuonDetectorPhases
        """
        self.alg = mantid.AlgorithmManager.create("CreateWorkspace")
        self.alg.initialize()
        self.alg.setAlwaysStoreInADS(False)
        for name, value in iteritems(inputs):
            self.alg.setProperty(name, value)
        self.alg.execute()
        mantid.AnalysisDataService.addOrReplace(
            inputs["OutputWorkspace"],
            self.alg.getProperty("OutputWorkspace").value)
        self.alg = None

    def report(self,inputs):
        output = "CreateWorkspace("
        for name, value in iteritems(inputs):
            print (value)
            safeValue = self.convert(value)
            print (safeValue)
            output+=str(name )+ " = "+safeValue +", "
        output=output[:-2]
        output+=") \n"
        return output

    def convert(self,values_in):
        value = values_in.split(",")
        out = ""
        if len(value) >1 :
           out+='['
           for xx in value:
               out+=str(xx) +","
           out = out[:-1]
           out +="]"
        else: # if its not a list we know its output workspace
              # so need quotes
           out ='"'+str( value[0])+'"'
        return out
