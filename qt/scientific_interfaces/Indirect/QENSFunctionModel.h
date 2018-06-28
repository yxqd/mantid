#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONMODEL_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONMODEl_H_

#include "MantidAPI/MultiDomainFunction.h"

#include <boost/optional.hpp>

#include <ostream>
#include <vector>

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

struct ParameterValue {
  explicit ParameterValue(double val) : value(val), error() {}
  ParameterValue(double value, double error) : value(value), error(error) {}
  double value;
  boost::optional<double> error;
};

struct ParameterData {
  explicit ParameterData(const std::string &aName, ParameterValue &&val,
                         const Mantid::API::ParameterTie *aTie)
      : name(aName), value(val), tie(aTie) {}
  const std::string name;
  const ParameterValue value;
  const Mantid::API::ParameterTie *tie;
};

class QENSFunctionModel {
public:
  QENSFunctionModel();

  Mantid::API::IFunction_sptr getFunction() const;
  Mantid::API::IFunction_sptr
  getConvolutionFunction(const std::vector<std::string> &resolution) const;
  std::size_t numberOfLocalParameters() const;

  void addFunction(const std::string &name,
                   const std::vector<std::size_t> &position);
  void removeFunction(const std::vector<std::size_t> &position);

  void setParameterValue(const std::string &name, const ParameterValue &value);

  template <typename T>
  void setAttribute(const std::string &name, const T &value,
                    const std::vector<std::size_t> &position);

  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace);
  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace,
                  const std::pair<std::size_t, std::size_t> &indexRange);
  void removeDataset(std::size_t index);
  void clearDatasets();

  template <typename BeginIterator, typename EndIterator>
  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace,
                  BeginIterator startIndexIt, EndIterator endIndexIt);

  void addEqualityGlobalTie(const std::string &parameterName);
  void addGlobalTie(const std::string &parameterName,
                    const std::string &expression);
  void removeGlobalTies(const std::string &parameterName);
  void clearTies();

  template <typename Functor> void forEachTie(const Functor &functor);

  void addGlobalUpperBound(const std::string &parameterName, double bound);
  void addGlobalLowerBound(const std::string &parameterName, double bound);
  void addGlobalBounds(const std::string &parameterName, double lowerBound,
                       double upperBound);

  void removeGlobalConstraints(const std::string &parameterName);
  void clearBoundConstraints();

protected:
  Mantid::API::IFunction_sptr getLocalFunction(std::size_t index) const;

  void fixParameterInDomain(const std::string &parameterName,
                            std::size_t domain);
  void unfixParameterInDomain(const std::string &parameterName,
                              std::size_t domain);

  void addLocalTieToDomain(const std::string &parameterName,
                           const std::string &expression, std::size_t domain);
  void removeLocalTieFromDomain(const std::string &parameterName,
                                std::size_t domain);
  void removeLocalTiesFromDomain(std::size_t domain);

  template <typename Predicate>
  void constrainLocalParametersInDomain(const std::string &expression,
                                        const Predicate &useInConstraint,
                                        std::size_t domain);

  template <typename Functor>
  void forEachParameterInDomain(const Functor &functor, std::size_t domain);

  void addUpperBoundToDomain(const std::string &parameterName, double bound,
                             std::size_t domain);
  void addLowerBoundToDomain(const std::string &parameterName, double bound,
                             std::size_t domain);
  void addBoundsToDomain(const std::string &parameterName, double lowerBound,
                         double upperBound, std::size_t domain);
  void addBoundsToDomainWithinPercentile(const std::string &parameterName,
                                         double percentile, std::size_t domain);

  void removeLocalConstraintFromDomain(const std::string &parameterName,
                                       std::size_t domain);
  void removeLocalConstraintsFromDomain(std::size_t domain);

private:
  Mantid::API::IFunction_sptr getNewLocalFunction() const;
  void addDatasetDomain(Mantid::API::MatrixWorkspace_sptr workspace,
                        std::size_t workspaceIndex);
  void addGlobalTie(const std::string &parameterName,
                    const std::string &expression, std::size_t fromDomain,
                    std::size_t toDomain);

  boost::shared_ptr<Mantid::API::MultiDomainFunction> m_function;
};

template <typename Predicate, typename OutputIterator>
OutputIterator findAllParametersWhich(const Mantid::API::IFunction &function,
                                      const Predicate &predicate,
                                      OutputIterator iterator) {
  for (auto i = 0; i < function.nParams(); ++i) {
    const auto parameter = function.getParameterName(i);
    if (predicate(parameter))
      *iterator++ = parameter;
  }
  return iterator;
}

ParameterValue getParameterValue(const Mantid::API::IFunction &function,
                                 std::size_t parameterIndex) {
  return ParameterValue(function.getParameter(parameterIndex),
                        function.getError(parameterIndex));
}

Mantid::API::CompositeFunction_sptr
getCompositeAt(Mantid::API::CompositeFunction_sptr composite,
               std::size_t index) {
  return boost::dynamic_pointer_cast<Mantid::API::CompositeFunction>(
      composite->getFunction(index));
}

template <typename BeginIterator, typename EndIterator>
Mantid::API::IFunction_sptr
getFunctionAt(const Mantid::API::CompositeFunction &composite,
              BeginIterator startIt, EndIterator endIt) {
  if (startIt >= endIt)
    return composite;

  for (auto i = startIt; i < endIt - 1; ++i) {
    if (!(composite = getCompositeAt(composite, *i)))
      return nullptr;
  }
  return composite->getFunction(*(endIt - 1));
}

Mantid::API::IFunction_sptr
getFunctionAt(const Mantid::API::CompositeFunction &function,
              const std::vector<std::size_t> &position) {
  getFunctionAt(function, position.begin(), position.end());
}

template <typename T>
void QENSFunctionModel::setAttribute(std::string &name, const T &value,
                                     const std::vector<std::size_t> &position) {
  if (auto function = getFunctionAt(*m_function, position))
    function->setAttribute(name, value);
}

template <typename BeginIterator, typename EndIterator>
void QENSFunctionModel::addDataset(Mantid::API::MatrixWorkspace_sptr workspace,
                                   BeginIterator startIndexIt,
                                   EndIterator endIndexIt) {
  for (auto i = startIndexIt; i < endIndexIt; ++i)
    addDatasetDomain(workspace, *i);
}

template <typename Predicate>
void QENSFunctionModel::constrainLocalParametersInDomain(
    const std::string &expression, const Predicate &useInConstraint,
    std::size_t domain) {
  const auto parameters = findAllParametersWhich(*m_function->getFunction(0));

  auto constraintExpression = expression;
  for (auto i = parameters.begin() + 1; i < parameters.end(); ++i)
    constraintExpression += "-" + *i;

  if (!parameters.empty())
    addLocalTie(parameters.front(), constraintExpression, domain);
}

template <typename Functor>
void QENSFunctionModel::forEachTie(const Functor &functor) {
  for (auto i = 0u; i < m_function->nParams(); ++i)
    functor(m_function->getTie(i));
}

template <typename Functor>
void QENSFunctionModel::forEachParameterInDomain(const Functor &functor,
                                                 std::size_t domain) {
  const auto function = getLocalFunction(domain);
  const auto names = function->getParameterNames();

  for (auto i = 0u; i < m_function->nParams(); ++i)
    functor(ParameterData(names[i], getParameterValue(*function, i),
                          function->getTie(i)));
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt

#endif
