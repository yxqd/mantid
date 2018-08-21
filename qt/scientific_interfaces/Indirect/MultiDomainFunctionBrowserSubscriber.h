#ifndef MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSERSUBSCRIBER_H_
#define MANTIDWIDGETS_MULTIDOMAINFUNCTIONBROWSERSUBSCRIBER_H_

#include <string>

namespace MantidQt {
namespace MantidWidgets {

class MultiDomainFunctionBrowserSubscriber {
public:
  virtual void globalChanged(std::string const &parameter, bool global) = 0;
  virtual void editParameter(std::string const &name) = 0;
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif
