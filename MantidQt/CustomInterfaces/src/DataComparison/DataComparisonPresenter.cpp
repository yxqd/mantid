#include "MantidQtCustomInterfaces/DataComparison/DataComparisonPresenter.h"
#include "MantidQtCustomInterfaces/DataComparison/IDataComparisonView.h"

using namespace MantidQt::CustomInterfaces;

/** Constructor
* @param view :: a pointer to the interface for the view we are managing
*/
DataComparisonPresenter::DataComparisonPresenter(IDataComparisonView *view)
    : m_view(view) {}

/** Destructor
*/
DataComparisonPresenter::~DataComparisonPresenter() {}
