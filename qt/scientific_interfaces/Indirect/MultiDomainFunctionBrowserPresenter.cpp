#include "MultiDomainFunctionBrowserPresenter.h"
#include "MultiDomainFunctionBrowser.h"
#include "MultiDomainFunctionModel.h"

#include "MDFEditLocalParameterPresenter.h"

namespace {
using MantidQt::MantidWidgets::MultiDomainFunctionModel;
using MantidQt::CustomInterfaces::MDF::EditLocalParameterModel;
using MantidQt::CustomInterfaces::MDF::EditLocalParameterPresenter;

bool executeEditParameterDialog(MultiDomainFunctionModel &model,
                                std::string const &parameter, QWidget *parent) {
  EditLocalParameterPresenter editPresenter(
      EditLocalParameterModel(model, parameter), parent);
  return editPresenter.executeDialog(model);
}
} // namespace

namespace MantidQt {
namespace MantidWidgets {

MultiDomainFunctionBrowserPresenter::MultiDomainFunctionBrowserPresenter(
    MultiDomainFunctionBrowser *browser, MultiDomainFunctionModel *model)
    : FunctionBrowserPresenter(browser, model) {
  browser->subscribeToMultiDomainBrowser(this);
}

void MultiDomainFunctionBrowserPresenter::globalChanged(
    std::string const &parameter, bool global) {
  if (global)
    m_multiDomainModel->addEqualityGlobalTie(parameter);
}

void MultiDomainFunctionBrowserPresenter::editParameter(
    std::string const &name) {
  if (auto const parent = m_multiDomainBrowser->parentWidget())
    executeEditParameterDialog(*m_multiDomainModel, name, parent);
  else
    executeEditParameterDialog(*m_multiDomainModel, name, m_multiDomainBrowser);
}

} // namespace MantidWidgets
} // namespace MantidQt