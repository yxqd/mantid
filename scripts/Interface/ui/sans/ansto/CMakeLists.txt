# UI files to convert
set( UI_FILES
ansto_sans.ui
 )

UiToPy( UI_FILES CompileUISANSANSTO)

set ( TEST_PY_FILES
      MainViewPresenterTest.py
      )


# python unit tests
if (PYUNITTEST_FOUND)
  pyunittest_add_test ( ${CMAKE_CURRENT_SOURCE_DIR} PythonUIScriptsTest ${TEST_PY_FILES} )
endif ()