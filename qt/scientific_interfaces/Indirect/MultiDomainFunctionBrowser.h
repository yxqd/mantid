#ifndef MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSER_H_
#define MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSER_H_

#include "FunctionBrowser.h"

#include <memory>

namespace MantidQt {
namespace MantidWidgets {

class MultiDomainFunctionBrowserSubscriber;

class MultiDomainFunctionBrowser : public FunctionBrowser {
public:
  MultiDomainFunctionBrowser();
  MultiDomainFunctionBrowser(QWidget *parent);

  void subscribeToMultiDomainBrowser(
      MultiDomainFunctionBrowserSubscriber *subscriber);

protected slots:
  void globalChanged(QtProperty *, QString const &, bool);
  void parameterButtonClicked(QtProperty *);

private:
  std::unique_ptr<QtTreePropertyBrowser> createNewBrowser() override;
  std::unique_ptr<QtAbstractEditorFactory<ParameterPropertyManager>>
  getParameterEditorFactory() override;
  AProperty addParameterProperty(QtProperty *prop, QtProperty *parent,
                                 QString const &name,
                                 QString const &description,
                                 double value) override;

  MultiDomainFunctionBrowserSubscriber *m_multiDomainSubscriber;
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif
