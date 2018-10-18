.. _AlgorithmMPISupport:

===================
Saving Python GUI's 
===================

.. contents::
  :local:


Introduction
------------

The python interfaces can now be saved as part of project recovery. This is a brief guide on how to add this functionality into a python GUI. Most of the changes will be in the python 
code itself.

Project Saving 
--------------

For project saving to record a python interface it must know that the interface can be saved. This is done by adding the name of your interface to the ProjectSerialiser_.

.. _ProjectSerialiser: https://github.com/mantidproject/mantid/blob/master/MantidPlot/src/ProjectSerialiser.cpp#L133

The Interface Code
------------------

The code which creates the interface should have the following structure (pseudo code):

.. code:: python
   
    #imports
   
    Name = "name of GUI"
    class GUI(QtGui.QMainWindow):
	
        def __init__(self,parent =None):
	       #some stuff in here
		   
        def saveToProject(self):
	       # In this case; save the context = save the project
	       return self._context.save()
		   
        def loadFromProject(self,project):
            # code to load data
            self._context.loadFromProject(project)

    def main():
        # if GUI is open bring to front
        widget = getWidgetIfOpen()
            widget.raise_()
            return widget
        # creates the GUI window
        app = qapp()
        # some code
        myGUI = GUI()
        #need the follwoing line
        myGUI.setAccessibleName(Name)
        return myGUI
		
    def getWidgetIfOpen():
        allWidgets = QtGui.QApplication.allWidgets()
        for widget in allWidgets:
            if widget.accessibleName() == Name:
               return widget
        return None
	
    # the names of the following are fixed
    def saveToProject():
        widget = getWidgetIfOpen()
        if widget is None:
           return ""
        project = widget.saveToProject()
        return project
		
    def loadFromProject(project):
        myGUI = main()
        myGUI.loadFromProject(project)
        return myGUI

    if __name == '__main__':
       myGUI = main()
		
The important things to note are the load and save methods within the interface class, the inclusion of :code:`setAccessibleName` in the `main` method for the file, the :code:`saveToProject` 
and :code:`loadFromProject` methods. The :code:`setAccessibleName` is needed to retrieve the GUI when saving the data, see the :code:`getWidgetIfOpen` method. 

Code To Save Interace
---------------------

To save the project the TSV serialiser must be used. It is suggested that each bit of information has its own line with a name that can be used to identify the data when loading it
into the GUI. The simplest way to achieve this is to use the keys from a dictionary that stores the data. The `pythonTSV` file has some simple methods to help with passing data to 
a TSV serialiser. The :code:`pythonTSV.writeLine(TSV,name)` method takes a TSV serialiser (:code:`TSV`) and the name of the line you want to write, the advantage of using this method is that it 
will ensure that the name is compatible with the TSV serialiser (i.e. it will not contain spaces, underscores or dashes). The :code:`pythonTSV.saveToTSV(TSV,value)` method 
takes a TSV serialiser
and the value you want to store, it then tries to identify the correct data type for value and if it fails it will return a :code:`TypeError`. Below is some example code to save a context 
(dictionary):

.. code:: python

    #imports
    import pythonTSV
	
    # some code
    def saveToProject(self):
        TSVSec = MantidQT.API.TSVSerialiser()
        TSV = MantidQT.API.TSVSerialiser()
        # lets write the keys from the context
        TSV.writeLine("keys")
        keys = self.context.keys()
        # record how many items to read
        TSV.storeInt(len(keys))
        # record key values
        for key in keys:
            TSV.storeString(key)
        # lets record the context values
        for key in keys:
            pythonTSV.writeLine(TSV,key)
	        value = self.context[key]
	        try:
                pythonTSV.saveToTSV(TSV,value)
            except:
                # if it fails lets do the saving ourself
                self.saveCustom(TSV,key,value)
        lines = TSV.outputLines()
        TSVSec.writeSection("Demo",lines)
        return TSVSec.outputLines()
		
In the above code the data has been stored within a section called :code:`Demo`, the name of the section  has the same constraints as the line names. 
If this is a problem you can create
a safe name by using :code:`pythonTSV.makeLineNameSafe(unsafe_name)`. The function has to return a string, hence the :code:`outputLines` on the last line. 

Code To Load Interface
----------------------

The TSV serialiser is also required for loading data back into the GUI. The `pythonTSV` file also contains :code:`loadFromTSV(TSV,key,value)`, which takes a TSV serialiser, the key (line name),
and the current value (this is just needed to identify the data type). If the list of line names has been recorded then the load method can iterate over them, this would leave any 
key-value pairs unchanged that have been added since the GUI was last saved (maintains backwards compatibility). Below is an example piece of code:

.. code:: python

   #imports
   import pythonTSV

   def loadFromProject(self,project):
       #read section
       full_load = MantidQt.API.TSVSerialiser(project)
       secs = full_load("Demo")
	   
       # load in relevant section
       load = MantidQt.API.TSVSerialiser(secs[0])
       # read keys
       load.selectLine("keys")
       numKeys = load.readInt()
       keys = []
       for k in range(numKeys):
           tmp = load.readString()
           keys.append(tmp)
		   
        # read in values for keys    
        for key in keys:
            value = self.context[key]
            try:
	            self.context[key] = pythonTSV.loadFromTSV(load,key,value)
            except:
                self.customLoad(load,key,value)


If the name for the GUI contains unsafe characters the :code:`pythonTSV.makeLineNameSafe(unsafe_name)` method can be used to generate a safe name.