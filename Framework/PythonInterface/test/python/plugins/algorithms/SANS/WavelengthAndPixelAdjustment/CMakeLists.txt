######################
# Tests for NormalizeToMonitor
######################

set ( TEST_PY_FILES
SANSNormalizeToMonitorTest.py
SANSCalculateTransmissionTest.py
SANSCreateAdjustmentWorkspacesTest.py
SANSCreateWavelengthAndPixelAdjustmentTest.py
)

check_tests_valid ( ${CMAKE_CURRENT_SOURCE_DIR} ${TEST_PY_FILES} )

# Prefix for test name=PythonAlgorithms_SANS
pyunittest_add_test ( ${CMAKE_CURRENT_SOURCE_DIR} PythonAlgorithms ${TEST_PY_FILES} )

