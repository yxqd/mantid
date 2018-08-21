#include "MultiDomainFunctionBrowser.h"
#include "MultiDomainFunctionBrowserSubscriber.h"

#include "MantidQtWidgets/Common/QtPropertyBrowser/CompositeEditorFactory.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/DoubleDialogEditor.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/qttreepropertybrowser.h"

namespace {
constexpr char *globalOptionName = "Global";
}

namespace MantidQt {
namespace MantidWidgets {

MultiDomainFunctionBrowser::MultiDomainFunctionBrowser() : FunctionBrowser() {}

MultiDomainFunctionBrowser::MultiDomainFunctionBrowser(QWidget *parent)
    : FunctionBrowser(parent) {}

void MultiDomainFunctionBrowser::subscribeToMultiDomainBrowser(
    MultiDomainFunctionBrowserSubscriber *subscriber) {
  m_multiDomainSubscriber = subscriber;
}

std::unique_ptr<QtTreePropertyBrowser>
MultiDomainFunctionBrowser::createNewBrowser() {
  auto browser = Mantid::Kernel::make_unique<QtTreePropertyBrowser>(
      nullptr, QStringList{globalOptionName});
  connect(browser.get(),
          SIGNAL(optionChanged(QtProperty *, const QString &, bool)), this,
          SLOT(globalChanged(QtProperty *, const QString &, bool)));
  return std::move(browser);
}

std::unique_ptr<QtAbstractEditorFactory<ParameterPropertyManager>>
MultiDomainFunctionBrowser::getParameterEditorFactory() {
  auto buttonFactory = new DoubleDialogEditorFactory(this);
  auto compositeFactory = Mantid::Kernel::make_unique<
      CompositeEditorFactory<ParameterPropertyManager>>(this, buttonFactory);
  compositeFactory->setSecondaryFactory(globalOptionName,
                                        new ParameterEditorFactory(this));
  connect(buttonFactory, SIGNAL(buttonClicked(QtProperty *)), this,
          SLOT(parameterButtonClicked(QtProperty *)));
  connectEditorCloseToBrowser(buttonFactory);
  return std::move(compositeFactory);
}

MultiDomainFunctionBrowser::AProperty
MultiDomainFunctionBrowser::addParameterProperty(QtProperty *prop,
                                                 QtProperty *parent,
                                                 QString const &name,
                                                 QString const &description,
                                                 double value) {
  prop->setOption(globalOptionName, false);
  return FunctionBrowser::addParameterProperty(prop, parent, name, description,
                                               value);
}

void MultiDomainFunctionBrowser::globalChanged(QtProperty *,
                                               QString const &parameter,
                                               bool global) {
  m_multiDomainSubscriber->globalChanged(parameter.toStdString(), global);
}

void MultiDomainFunctionBrowser::parameterButtonClicked(QtProperty *prop) {
  m_multiDomainSubscriber->editParameter(getParameterName(prop).toStdString());
}

} // namespace MantidWidgets
} // namespace MantidQt