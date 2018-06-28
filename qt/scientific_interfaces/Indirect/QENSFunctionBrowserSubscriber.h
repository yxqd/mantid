#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSERSUBSCRIBER_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSERSUBSCRIBER_H_

#include <string>
#include <vector>

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

class QENSFunctionBrowserSubscriber {
public:
  virtual void addFunction(const std::string &name,
                           const std::vector<std::size_t> &position) = 0;
  virtual void removeFunction(const std::vector<std::size_t> &position) = 0;
  virtual void fixParameter(const std::string &name) = 0;
  virtual void removeTie(const std::string &name) = 0;
  virtual void addTie(const std::string &name,
                      const std::string &expression) = 0;
  virtual void addConstraints(const std::string &name, double upperBound,
                              double lowerBound) = 0;
  virtual void addConstraints10(const std::string &name) = 0;
  virtual void addConstraints50(const std::string &name) = 0;
  virtual void removeConstraint(const std::string &name) = 0;
  virtual void removeConstraints(const std::string &name) = 0;

  virtual void setStringAttribute(const std::string &name,
                                  const std::string &value,
                                  const std::vector<std::size_t> &position) = 0;
  virtual void setDoubleAttribute(const std::string &name, double value,
                                  const std::vector<std::size_t> &position) = 0;
  virtual void setIntAttribute(const std::string &name, int value,
                               const std::vector<std::size_t> &position) = 0;
  virtual void setBoolAttribute(const std::string &name, bool value,
                                const std::vector<std::size_t> &position) = 0;
  virtual void
  setVectorDoubleAttribute(const std::string &name,
                           const std::vector<double> &value,
                           const std::vector<std::size_t> &position) = 0;
};

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt

#endif