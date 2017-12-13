from __future__ import (absolute_import, print_function)

import unittest

from mantid.api import WorkspaceFactory
import numpy as np

class MatrixWorkspaceTest(unittest.TestCase):

    def test_setData_with_all_fields(self):
        nbins = 10
        nspec = 2
        xdata = np.arange((nbins+1)*nspec).reshape(nspec, nbins+1)
        ydata = 2*np.arange(nbins*nspec).reshape(nspec, nbins)
        edata = np.sqrt(ydata)
        dxdata = 0.5*np.sqrt(ydata)
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec, XLength=nbins+1, YLength=nbins)
        ws.setData(xdata, ydata, edata, dxdata)

        x_extracted, y_extracted = ws.extractX(), ws.extractY()
        e_extracted, dx_extracted = ws.extractE(), ws.extractDx()
        self.assertTrue(np.array_equal(xdata, x_extracted))
        self.assertTrue(np.array_equal(ydata, y_extracted))
        self.assertTrue(np.array_equal(edata, e_extracted))
        self.assertTrue(np.array_equal(dxdata, dx_extracted))

    def test_setData_accepts_single_field(self):
        nspec, nbins = 2, 10
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec, XLength=nbins+1, YLength=nbins)
        for field, extract_suffix in (('x', 'X'), ('y', 'Y'), ('e', 'E'), ('dx', 'Dx')):
            data = np.arange(nspec*(nbins+1)).reshape(nspec, nbins+1)
            ws.setData(y=data)
            extracted = getattr(ws, 'extract' + extract_suffix)()
            self.assertTrue(np.array_equal(data, extracted))

    def test_setData_raises_error_if_no_fields_supplied(self):
        nspec, nbins = 1, 1
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec, XLength=nbins+1, YLength=nbins)
        self.assertRaises(ValueError, ws.setData)

    def test_setData_raises_error_for_incorrect_shape(self):
        nbins = 10
        nspec = 2
        xdata = np.arange((nbins+1)*nspec)
        xdata = xdata.reshape(nspec, nbins+1)

        # bad shape
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec, XLength=nbins+1, YLength=nbins)
        self.assertRaises(ValueError, ws.setData, np.arange(nbins))
        # bad number of spectra
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec-1, XLength=nbins+1, YLength=nbins)
        self.assertRaises(ValueError, ws.setData, xdata)
        # bad number of bins
        ws = WorkspaceFactory.create("Workspace2D", NVectors=nspec, XLength=nbins-1, YLength=nbins-1)
        self.assertRaises(ValueError, ws.setData, xdata)


if __name__ == '__main__':
    unittest.main()
