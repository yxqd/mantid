include (UiToPy)

# List of UIs to Auto convert
set( UI_FILES
  IndirectQuickRun.ui
  EnergyWindowScan.ui
  SQWMomentScan.ui
  DiffractionScan.ui
  SampleChanger.ui
)

UiToPy( UI_FILES CompileUIIndirectQuickRun)
