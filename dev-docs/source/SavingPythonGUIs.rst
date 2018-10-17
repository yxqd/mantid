.. _AlgorithmMPISupport:

===================
Saving Python GUI's 
===================

.. contents::
  :local:

Concept
#######

Introduction
------------

The python interfaces can now be saved as part of project recovery. This is a brief guide on how to add this functionality into a python GUI. Most of the changes will be in the python 
code itself.

Project Recovery 
----------------

For project recovery to record a python interface it must know that the interface can be saved. This is done by adding the name of your interface to the ProjectSerialiser_.

.. _ProjectSerialiser: https://github.com/mantidproject/mantid/blob/master/MantidPlot/src/ProjectSerialiser.cpp#L133

The Interface Code
------------------

The code which creates the interace should have the follwoing structure (pseudo code):

.. code:: python
   
    #imports
   
    Name = "name of GUI"
    class GUI(QtGui.QMainWindow):
	
        def __init__(self,parent =None):
	       #some stuff in here
        def saveToProject(self):
	       # In this case; save the context = save the project
	       return self._context.save()
	   
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
		
done. 