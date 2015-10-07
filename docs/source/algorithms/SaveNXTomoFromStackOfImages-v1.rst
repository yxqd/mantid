
.. algorithm::

.. summary::

.. alias::

.. properties::

Description
-----------

The aim of this algorithm is to produce an NXTomo output file as fast
as possible from a set of input files (stack of images, containing
sample and possibly flat field and dark field images). The output file
is in the NXTomo format, as specified in the `NeXus NXTomo application
definition
<http://download.nexusformat.org/sphinx/classes/applications/NXtomo.html>`__,

For file names, this algorithm in principle assumes the common
convention that files from a stack of images are named using a
trailing sequence number (with a varying number of digits before the
file format extension dot, if any), as for example in ImageJ. If you
have files sample000.fits, sample0001.fits and sample0002.fits, you
could use SampleFilename=sample000.fits, and unless any limit option
is given the algorithm will scan all the files named 'sampleZZZZ.fits'
where ZZZZ are the image sequence numbers. The maximum index that the
algorithm scans can be set by using the properties SampleMaxIndex,
FlatFieldMaxIndex, and DarkFieldMaxIndex.

An alternative to this algorithm, which operates only on workspaces is
:ref:`algm-SaveNXTomo`.

Usage
-----
..  Try not to use files in your examples,
    but if you cannot avoid it then the (small) files must be added to
    autotestdata\UsageData and the following tag unindented
    .. include:: ../usagedata-note.txt

**Example - SaveNXTomoFromStackOfImages**

.. testcode:: SaveNXTomoFromStackOfImagesExample

   # Create a host workspace
   ws = CreateWorkspace(DataX=range(0,3), DataY=(0,2))

   wsOut = SaveNXTomoFromStackOfImages()

   # Print the result
   print "The output workspace has %i spectra" % wsOut.getNumberHistograms()

Output:

.. testoutput:: SaveNXTomoFromStackOfImagesExample

  The output workspace has ?? spectra

.. categories::

