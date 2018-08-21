#ifndef MANTIDWIDGETS_IFUNCTIONMODEL_H_
#define MANTIDWIDGETS_IFUNCTIONMODEL_H_

#include "MantidAPI/IFunction.h"

#include <boost/optional.hpp>
#include <string>

namespace MantidQt {
namespace MantidWidgets {

class IFunctionModel {
public:
  virtual std::size_t numberOfParameters() const = 0;
  virtual std::string parameterName(std::size_t index) const = 0;
  virtual double parameterValue(std::string const &name) const = 0;
  virtual boost::optional<double>
  parameterError(std::string const &name) const = 0;
  virtual boost::optional<std::string>
  parameterTie(std::string const &name) const = 0;
  virtual boost::optional<double>
  parameterLowerBound(std::string const &name) const = 0;
  virtual boost::optional<double>
  parameterUpperBound(std::string const &name) const = 0;
  virtual std::vector<std::string> getAttributeNames() const = 0;
  virtual Mantid::API::IFunction::Attribute
  getAttribute(std::string const &name) = 0;
  virtual bool isComposite(std::vector<std::size_t> const &position) const = 0;
  virtual std::size_t
  numberOfFunctionsAt(std::vector<std::size_t> const &position) const = 0;
  virtual bool isParameterTied(std::string const &name) const = 0;
  virtual bool isParameterFixed(std::string const &name) const = 0;
  virtual bool isParameterConstrained(std::string const &name) const = 0;
  virtual std::string getLocalFunctionString() const = 0;
  virtual void setFunction(std::string const &functionString) = 0;
  virtual std::size_t addFunction(std::string const &name,
                                  std::vector<std::size_t> const &position) = 0;
  virtual void removeFunction(std::vector<std::size_t> const &position) = 0;
  virtual void setLocalParameterValue(std::string const &parameterName,
                                      double value) = 0;
  virtual void fixLocalParameter(std::string const &parameterName) = 0;
  virtual void unfixLocalParameter(std::string const &parameterName) = 0;
  virtual void setLocalTie(std::string const &parameterName,
                           std::string const &expression) = 0;
  virtual void removeLocalTie(std::string const &parameterName) = 0;
  virtual void removeLocalTies() = 0;
  virtual void addLocalUpperBound(std::string const &parameterName,
                                  double bound) = 0;
  virtual void addLocalLowerBound(std::string const &parameterName,
                                  double bound) = 0;
  virtual void addLocalBounds(std::string const &parameterName,
                              double lowerBound, double upperBound) = 0;
  virtual void addLocalBoundsWithinPercentile(std::string const &parameterName,
                                              double percentile) = 0;
  virtual void removeLocalConstraint(std::string const &parameterName,
                                     std::string const &type) = 0;
  virtual void removeLocalConstraints(std::string const &parameterName) = 0;
  virtual void setStringAttribute(std::string const &name,
                                  std::string const &value) = 0;
  virtual void setDoubleAttribute(std::string const &name, double value) = 0;
  virtual void setIntAttribute(std::string const &name, int value) = 0;
  virtual void setBoolAttribute(std::string const &name, bool value) = 0;
  virtual void setVectorAttribute(std::string const &name,
                                  std::vector<double> const &value) = 0;
  virtual void setVectorAttributeSize(std::string const &name,
                                      std::size_t size) = 0;
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif
