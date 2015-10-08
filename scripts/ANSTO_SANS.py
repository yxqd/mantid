
#pylint: disable=invalid-name
"""
    Script used to start the SANS/ANSTO Batch GUI from MantidPlot
"""
from Interface.ui.sans.ansto.sans_ansto_gui import SANSBatchGui

ui = SANSBatchGui()
if ui.setup_layout():
    ui.show()
