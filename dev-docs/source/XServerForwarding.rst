.. _XServerForwarding:

==================================
Run MantidPlot on a remote machine 
==================================

.. contents::
   :local:

Introduction
------------

Running MantidPlot with its graphical interface on a remote Linux machine requires an ssh client as well as an X server.

From a Linux computer
---------------------

Open a terminal and ...

::

  $ ssh -X RemoteMachineName
  $ MantidPlot

From a Windows computer
-----------------------

Choose Xming_ and PuTTY_ as X server and ssh client, respectively.

.. _PuTTY: https://www.putty.org/

.. _Xming: http://www.straightrunning.com/XmingNotes/

Fonts (and font sizes) between the Linux remote machine and your machine will likely differ.
Some of your options are:

- Xming-fonts can be downloaded from Xming_ (need to be installed in the same directory as Xming, e.g. `C:\\Program Files\\Xming`)
- Use a remote X font server and avoid installing extra fonts on your machine (start Xming with option `-fp`)
- Add fonts or use Windows System TrueType fonts (modify `C:\\Windows\\Fonts`)

Please consider as well the application Xming-portablePuTTY_.

.. _Xming-portablePuTTY: http://www.straightrunning.com/XmingNotes/portable.php

A summary from 2006 (fonts link does not work anymore) can be found here_.

.. _here: http://courses.cms.caltech.edu/cs11/misc/xwindows.html
