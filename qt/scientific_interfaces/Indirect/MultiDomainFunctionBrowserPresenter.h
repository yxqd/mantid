#ifndef MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSERPRESENTER_H_
#define MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSERPRESENTER_H_

#include "FunctionBrowserPresenter.h"
#include "MultiDomainFunctionBrowserSubscriber.h"

namespace MantidQt {
namespace MantidWidgets {

class MultiDomainFunctionModel;
class MultiDomainFunctionBrowser;

class MultiDomainFunctionBrowserPresenter
    : public FunctionBrowserPresenter,
      public MultiDomainFunctionBrowserSubscriber {
public:
  MultiDomainFunctionBrowserPresenter(MultiDomainFunctionBrowser *browser,
                                      MultiDomainFunctionModel *model);

  void globalChanged(std::string const &parameter, bool global) override;
  void editParameter(std::string const &name) override;

private:
  MultiDomainFunctionBrowser *m_multiDomainBrowser;
  MultiDomainFunctionModel *m_multiDomainModel;
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif