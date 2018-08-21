#include "MultiDomainFunctionModel.h"

#include "MantidAPI/ConstraintFactory.h"
#include "MantidAPI/Expression.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/IConstraint.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/MultiDomainFunction.h"

#include "MantidKernel/make_unique.h"

using namespace Mantid::API;

namespace {
using namespace MantidQt::MantidWidgets;

Expression createExpression(std::string const &str) {
  Expression expression;
  expression.parse(str);
  return expression;
}

bool expressionContainsParameter(Expression const &expression,
                                 std::string const &parameter) {
  for (auto const &term : expression) {
    if (term.str() == parameter)
      return true;
  }
  return false;
}

bool expressionContainsParameter(std::string const &expression,
                                 std::string const &parameter) {
  return expressionContainsParameter(createExpression(expression), parameter);
}

void removeTiesWhichContainParameter(FunctionProperties &properties,
                                     std::string const &parameter) {
  properties.removeTieIf(
      [&](std::string const &tieParameter, std::string const &expression) {
        return tieParameter == parameter ||
               expressionContainsParameter(expression, parameter);
      });
}

void removeConstraintsWhichContainParameter(FunctionProperties &properties,
                                            std::string const &parameter) {
  properties.removeConstraintIf(
      [&](const std::pair<std::string, std::string> &constraint) {
        return constraint.first == parameter ||
               expressionContainsParameter(constraint.second, parameter);
      });
}

std::string functionIndex(std::size_t index) {
  return "f" + std::to_string(index) + ".";
}

template <typename T> IFunction::Attribute createAttribute(T const &value) {
  return IFunction::Attribute(value);
}

IFunction_sptr createFunction(std::string const &name) {
  return FunctionFactory::Instance().createFunction(name);
}

CompositeFunction_sptr getCompositeAt(CompositeFunction_sptr composite,
                                      std::size_t index) {
  return boost::dynamic_pointer_cast<CompositeFunction>(
      composite->getFunction(index));
}

template <typename BeginIterator, typename EndIterator>
IFunction_sptr getFunctionAt(CompositeFunction const &composite,
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
getFunctionAt(Mantid::API::CompositeFunction const &function,
              std::vector<std::size_t> const &position) {
  getFunctionAt(function, position.begin(), position.end());
}

std::size_t addFunctionAt(CompositeFunction_sptr composite,
                          std::string const &name,
                          std::vector<std::size_t> const &position) {
  auto const function = getFunctionAt(*composite, position);
  if (composite = boost::dynamic_pointer_cast<CompositeFunction>(function)) {
    composite->addFunction(createFunction(name));
    return composite->nFunctions() - 1;
  }
  throw std::runtime_error("Unable to add function to non-composite function.");
}

IFunction_sptr removeFunctionAt(CompositeFunction_sptr function,
                                std::vector<std::size_t> const &position) {
  function = getFunctionAt(*function, position.begin(), position.end() - 1);
  auto composite = boost::dynamic_pointer_cast<CompositeFunction>(function);
  if (composite && !position.empty()) {
    auto const removed = composite->getFunction(position.back());
    composite->removeFunction(position.back());
    return removed;
  }
  return nullptr;
}

std::unique_ptr<IConstraint> createConstraint(IFunction &function,
                                              std::string const &expression) {
  auto constraint =
      ConstraintFactory::Instance().createInitialized(&function, expression);
  return Mantid::Kernel::make_unique<IConstraint>(constraint);
}

void addTiesToFunction(IFunction &function,
                       FunctionProperties const &properties) {
  properties.forEachTie(
      [&function](std::string const &parameter, std::string const &expression) {
        function.tie(parameter, expression);
      });
}

void addConstraintsToFunction(IFunction &function,
                              FunctionProperties const &properties) {
  properties.forEachConstraint([&function](std::string const &constraint) {
    function.addConstraint(createConstraint(function, constraint));
  });
}

void setParameterValueInFunction(IFunction &function,
                                 std::string const &parameter,
                                 ParameterValue const &value) {
  auto const index = function.parameterIndex(parameter);
  function.setParameter(index, value.value);
  if (value.error)
    function.setError(index, *value.error);
}

void setParameterValuesInFunction(IFunction &function,
                                  LocalFunctionProperties const &properties) {
  properties.forEachParameter(
      [&function](std::string const &parameter, ParameterValue const &value) {
        setParameterValueInFunction(function, parameter, value);
      });
}

void setAttributeValuesInFunction(IFunction &function,
                                  LocalFunctionProperties const &properties) {
  properties.forEachAttribute([&function](std::string const &attribute,
                                          const IFunction::Attribute &value) {
    function.setAttribute(attribute, value);
  });
}

void addPropertiesToFunction(IFunction &function,
                             FunctionProperties const &properties) {
  addTiesToFunction(function, properties);
  addConstraintsToFunction(function, properties);
}

void addWorkspacePropertiesToFunction(
    IFunction &function, LocalFunctionProperties const &properties) {
  auto const workspace = properties.getWorkspace();
  function.setMatrixWorkspace(workspace, properties.getWorkspaceIndex(),
                              workspace->x(0).front(), workspace->x(0).back());
}

void addPropertiesToLocalFunction(IFunction &function,
                                  LocalFunctionProperties const &properties) {
  addPropertiesToFunction(function, properties);
  setParameterValuesInFunction(function, properties);
  setAttributeValuesInFunction(function, properties);
  addWorkspacePropertiesToFunction(function, properties);
}

IFunction_sptr createLocalFunction(IFunction const &function,
                                   LocalFunctionProperties const &properties) {
  auto newFunction = function.clone();
  addPropertiesToLocalFunction(*newFunction, properties);
  return newFunction;
}

boost::shared_ptr<MultiDomainFunction> createMultiDomainFunction(
    IFunction const &localFunction,
    std::vector<LocalFunctionProperties> const &localFunctionProperties,
    FunctionProperties const &globalFunctionProperties) {
  boost::shared_ptr<MultiDomainFunction> function(new MultiDomainFunction);
  for (auto i = 0u; i < localFunctionProperties.size(); ++i) {
    function->addFunction(
        createLocalFunction(localFunction, localFunctionProperties[i]));
    function->setDomainIndex(i, i);
  }
  addPropertiesToFunction(*function, globalFunctionProperties);
  return function;
}
} // namespace

namespace MantidQt {
namespace MantidWidgets {

MultiDomainFunctionModel::MultiDomainFunctionModel()
    : m_function(new CompositeFunction) {}

IFunction_sptr MultiDomainFunctionModel::getFitFunction() const {
  if (m_localFunctionProperties.size() == 1)
    return getLocalFunction(0);
  return createMultiDomainFunction(*m_function, m_localFunctionProperties,
                                   m_globalFunctionProperties);
}

std::size_t MultiDomainFunctionModel::numberOfParameters() const {
  return m_function->nParams();
}

std::size_t MultiDomainFunctionModel::numberOfDomains() const {
  return m_localFunctionProperties.size();
}

std::string MultiDomainFunctionModel::parameterName(std::size_t index) const {
  return m_function->parameterName(index);
}

double
MultiDomainFunctionModel::parameterValue(std::string const &parameter) const {
  return m_localFunctionProperties[m_activeDomain].getParameterValue(parameter);
}

boost::optional<double>
MultiDomainFunctionModel::parameterError(std::string const &parameter) const {
  return m_localFunctionProperties[m_activeDomain].getParameterError(parameter);
}

boost::optional<std::string>
MultiDomainFunctionModel::parameterTie(std::string const &parameter) const {
  try {
    return m_localFunctionProperties[m_activeDomain].getTie(parameter);
  } catch (std::runtime_error const &) {
    return boost::none;
  }
}

boost::optional<double>
MultiDomainFunctionModel::parameterLowerBound(std::string const &name) const {
  return m_localFunctionProperties[m_activeDomain].getParameterLowerBound(name);
}

boost::optional<double>
MultiDomainFunctionModel::parameterUpperBound(std::string const &name) const {
  return m_localFunctionProperties[m_activeDomain].getParameterUpperBound(name);
}

std::vector<std::string> MultiDomainFunctionModel::getAttributeNames() const {
  return m_function->getAttributeNames();
}

Mantid::API::IFunction::Attribute
MultiDomainFunctionModel::getAttribute(std::string const &name) {
  m_localFunctionProperties[m_activeDomain].getAttribute(name);
}

Mantid::API::MatrixWorkspace_sptr
MultiDomainFunctionModel::getWorkspace() const {
  m_localFunctionProperties[m_activeDomain].getWorkspace();
}

std::string MultiDomainFunctionModel::getWorkspaceName() const {
  return getWorkspace()->getName();
}

std::size_t MultiDomainFunctionModel::getWorkspaceIndex() const {
  m_localFunctionProperties[m_activeDomain].getWorkspaceIndex();
}

bool MultiDomainFunctionModel::isComposite(
    std::vector<std::size_t> const &position) const {
  return nullptr != boost::dynamic_pointer_cast<CompositeFunction>(
                        getFunctionAt(*m_function, position));
}

std::size_t MultiDomainFunctionModel::numberOfFunctionsAt(
    std::vector<std::size_t> const &position) const {
  auto const function = getFunctionAt(*m_function, position);
  if (auto const composite =
          boost::dynamic_pointer_cast<CompositeFunction>(function))
    return composite->nFunctions();
  return 0;
}

bool MultiDomainFunctionModel::isParameterTied(std::string const &name) const {
  return m_localFunctionProperties[m_activeDomain].isTied(name);
}

bool MultiDomainFunctionModel::isParameterFixed(std::string const &name) const {
  return m_localFunctionProperties[m_activeDomain].isFixed(name);
}

bool MultiDomainFunctionModel::isParameterConstrained(
    std::string const &name) const {
  return m_localFunctionProperties[m_activeDomain].isConstrained(name);
}

std::string MultiDomainFunctionModel::getLocalFunctionString() const {
  return getLocalFunction(m_activeDomain)->asString();
}

IFunction_sptr
MultiDomainFunctionModel::getLocalFunction(std::size_t domain) const {
  return createLocalFunction(*m_function, m_localFunctionProperties[domain]);
}

std::size_t MultiDomainFunctionModel::numberOfLocalParameters() const {
  return m_function->nParams();
}

std::size_t MultiDomainFunctionModel::getActiveDomain() const {
  return m_activeDomain;
}

void MultiDomainFunctionModel::setActiveDomain(std::size_t domain) {
  m_activeDomain = domain;
}

void MultiDomainFunctionModel::setFunction(std::string const &functionString) {
  m_function = Mantid::API::FunctionFactory::Instance().createInitialized(
      functionString);
}

std::size_t MultiDomainFunctionModel::addFunction(
    std::string const &name, std::vector<std::size_t> const &position) {
  return addFunctionAt(m_function, name, position);
}

void MultiDomainFunctionModel::removeFunction(
    std::vector<std::size_t> const &position) {
  auto const function = removeFunctionAt(m_function, position);
  for (auto i = 0u; i < function->nParams(); ++i)
    removeParameterProperties(function->parameterName(i));
}

void MultiDomainFunctionModel::removeTiesContainingParameter(
    std::string const &parameter) {
  for (auto &localProperties : m_localFunctionProperties)
    removeTiesWhichContainParameter(localProperties, parameter);
  removeTiesWhichContainParameter(m_globalFunctionProperties,
                                  functionIndex(m_activeDomain) + parameter);
}

void MultiDomainFunctionModel::removeConstraintsContainingParameter(
    std::string const &parameter) {
  for (auto &localProperties : m_localFunctionProperties)
    removeConstraintsWhichContainParameter(localProperties, parameter);
  removeConstraintsWhichContainParameter(
      m_globalFunctionProperties, functionIndex(m_activeDomain) + parameter);
}

void MultiDomainFunctionModel::removeParameterProperties(
    std::string const &parameterName) {
  removeTiesContainingParameter(parameterName);
  removeConstraintsContainingParameter(parameterName);
}

void MultiDomainFunctionModel::setParameterValue(std::size_t domain,
                                                 std::string const &name,
                                                 const ParameterValue &value) {
  m_localFunctionProperties[domain].setParameter(name, value);
}

void MultiDomainFunctionModel::setStringAttribute(std::string const &name,
                                                  std::string const &value) {
  m_localFunctionProperties[m_activeDomain].setAttribute(
      name, createAttribute(value));
}

void MultiDomainFunctionModel::setDoubleAttribute(std::string const &name,
                                                  double value) {
  m_localFunctionProperties[m_activeDomain].setAttribute(
      name, createAttribute(value));
}

void MultiDomainFunctionModel::setIntAttribute(std::string const &name,
                                               int value) {
  m_localFunctionProperties[m_activeDomain].setAttribute(
      name, createAttribute(value));
}

void MultiDomainFunctionModel::setBoolAttribute(std::string const &name,
                                                bool value) {
  m_localFunctionProperties[m_activeDomain].setAttribute(
      name, createAttribute(value));
}

void MultiDomainFunctionModel::setVectorAttribute(
    std::string const &name, const std::vector<double> &value) {
  m_localFunctionProperties[m_activeDomain].setAttribute(
      name, createAttribute(value));
}

void MultiDomainFunctionModel::setVectorAttributeSize(std::string const &name,
                                                      std::size_t size) {
  m_localFunctionProperties[m_activeDomain].resizeVectorAttribute(name, size);
}

void MultiDomainFunctionModel::addDataset(MatrixWorkspace_sptr workspace) {
  for (auto i = 0u; i < workspace->getNumberHistograms(); ++i)
    addDatasetDomain(workspace, i);
}

void MultiDomainFunctionModel::addDataset(
    MatrixWorkspace_sptr workspace,
    const std::pair<std::size_t, std::size_t> &indexRange) {
  for (auto i = indexRange.first; i <= indexRange.second; ++i)
    addDatasetDomain(workspace, i);
}

void MultiDomainFunctionModel::addDatasetDomain(MatrixWorkspace_sptr workspace,
                                                std::size_t workspaceIndex) {
  m_localFunctionProperties.emplace_back(
      LocalFunctionProperties(workspace, workspaceIndex));
}

void MultiDomainFunctionModel::removeDatasetDomain(std::size_t domain) {
  m_localFunctionProperties.erase(m_localFunctionProperties.begin() + domain);
}

void MultiDomainFunctionModel::clearDatasets() {
  m_function = CompositeFunction_sptr(new CompositeFunction);
  m_localFunctionProperties.clear();
  m_globalFunctionProperties = FunctionProperties();
}

void MultiDomainFunctionModel::setLocalParameterValue(
    std::string const &parameterName, double value) {
  m_localFunctionProperties[m_activeDomain].setParameterValue(parameterName,
                                                              value);
}

void MultiDomainFunctionModel::fixLocalParameter(
    std::string const &parameterName) {
  fixParameterInDomain(parameterName, m_activeDomain);
}

void MultiDomainFunctionModel::unfixLocalParameter(
    std::string const &parameterName) {
  unfixParameterInDomain(parameterName, m_activeDomain);
}

void MultiDomainFunctionModel::setLocalTie(std::string const &parameterName,
                                           std::string const &expression) {
  addLocalTieToDomain(parameterName, expression, m_activeDomain);
}

void MultiDomainFunctionModel::removeLocalTie(
    std::string const &parameterName) {
  removeLocalTieFromDomain(parameterName, m_activeDomain);
}

void MultiDomainFunctionModel::removeLocalTies() {
  removeLocalTiesFromDomain(m_activeDomain);
}

void MultiDomainFunctionModel::addLocalUpperBound(
    std::string const &parameterName, double bound) {
  addUpperBoundToDomain(parameterName, bound, m_activeDomain);
}

void MultiDomainFunctionModel::addLocalLowerBound(
    std::string const &parameterName, double bound) {
  addLowerBoundToDomain(parameterName, bound, m_activeDomain);
}

void MultiDomainFunctionModel::addLocalBounds(std::string const &parameterName,
                                              double lowerBound,
                                              double upperBound) {
  addBoundsToDomain(parameterName, lowerBound, upperBound, m_activeDomain);
}

void MultiDomainFunctionModel::addLocalBoundsWithinPercentile(
    std::string const &parameterName, double percentile) {
  addBoundsToDomainWithinPercentile(parameterName, percentile, m_activeDomain);
}

void MultiDomainFunctionModel::removeLocalConstraint(
    std::string const &parameterName, std::string const &type) {
  if (type == "LowerBound")
    m_localFunctionProperties[m_activeDomain].removeLowerBound(parameterName);
  else
    m_localFunctionProperties[m_activeDomain].removeUpperBound(parameterName);
}

void MultiDomainFunctionModel::removeLocalConstraints(
    std::string const &parameterName) {
  removeLocalConstraintsFromDomain(parameterName, m_activeDomain);
}

void MultiDomainFunctionModel::fixParameterInDomain(
    std::string const &parameterName, std::size_t domain) {
  m_localFunctionProperties[domain].fixParameter(parameterName);
}

void MultiDomainFunctionModel::unfixParameterInDomain(
    std::string const &parameterName, std::size_t domain) {
  m_localFunctionProperties[domain].removeTie(parameterName);
}

void MultiDomainFunctionModel::addLocalTieToDomain(
    std::string const &parameterName, std::string const &expression,
    std::size_t domain) {
  m_localFunctionProperties[domain].tie(parameterName, expression);
}

void MultiDomainFunctionModel::removeLocalTieFromDomain(
    std::string const &parameterName, std::size_t domain) {
  m_localFunctionProperties[domain].removeTie(parameterName);
}

void MultiDomainFunctionModel::removeLocalTiesFromDomain(std::size_t domain) {
  m_localFunctionProperties[domain].clearTies();
}

void MultiDomainFunctionModel::addEqualityGlobalTie(
    std::string const &parameterName) {
  auto const firstParameter = functionIndex(0) + parameterName;
  addGlobalTie(parameterName, firstParameter, 1,
               m_function->getNumberDomains());
}

void MultiDomainFunctionModel::addGlobalTie(std::string const &parameterName,
                                            std::string const &expression) {
  addGlobalTie(parameterName, expression, 0, m_function->getNumberDomains());
}

void MultiDomainFunctionModel::addGlobalTie(std::string const &parameterName,
                                            std::string const &expression,
                                            std::size_t fromDomain,
                                            std::size_t toDomain) {
  for (auto i = fromDomain; i < toDomain; ++i)
    m_globalFunctionProperties.tie(functionIndex(i) + parameterName,
                                   expression);
}

void MultiDomainFunctionModel::removeGlobalTies(
    std::string const &parameterName) {
  m_globalFunctionProperties.removeTie(parameterName);
}

void MultiDomainFunctionModel::removeLocalTies(
    std::string const &parameterName) {
  for (auto &properties : m_localFunctionProperties)
    properties.removeTie(parameterName);
}

void MultiDomainFunctionModel::clearTies() {
  m_globalFunctionProperties.clearTies();
  for (auto &localProperties : m_localFunctionProperties)
    localProperties.clearTies();
}

void MultiDomainFunctionModel::addUpperBoundToDomain(
    std::string const &parameterName, double bound, std::size_t domain) {
  m_localFunctionProperties[domain].setUpperBound(parameterName, bound);
}

void MultiDomainFunctionModel::addLowerBoundToDomain(
    std::string const &parameterName, double bound, std::size_t domain) {
  m_localFunctionProperties[domain].setLowerBound(parameterName, bound);
}

void MultiDomainFunctionModel::addBoundsToDomain(
    std::string const &parameterName, double lowerBound, double upperBound,
    std::size_t domain) {
  m_localFunctionProperties[domain].setConstraint(parameterName, lowerBound,
                                                  upperBound);
}

void MultiDomainFunctionModel::addBoundsToDomainWithinPercentile(
    std::string const &parameterName, double percentile, std::size_t domain) {
  auto const value = m_function->getParameter(parameterName);
  auto const lower = value * (1.0 - percentile);
  auto const upper = value * (1.0 + percentile);
  addBoundsToDomain(parameterName, lower, upper, domain);
}

void MultiDomainFunctionModel::addGlobalUpperBound(
    std::string const &parameterName, double bound) {
  m_globalFunctionProperties.setUpperBound(parameterName, bound);
}

void MultiDomainFunctionModel::addGlobalLowerBound(
    std::string const &parameterName, double bound) {
  m_globalFunctionProperties.setLowerBound(parameterName, bound);
}

void MultiDomainFunctionModel::addGlobalBounds(std::string const &parameterName,
                                               double lowerBound,
                                               double upperBound) {
  m_globalFunctionProperties.setConstraint(parameterName, lowerBound,
                                           upperBound);
}

void MultiDomainFunctionModel::removeLocalConstraintsFromDomain(
    std::string const &parameterName, std::size_t domain) {
  m_localFunctionProperties[domain].removeConstraints(parameterName);
}

void MultiDomainFunctionModel::clearLocalConstraintsFromDomain(
    std::size_t domain) {
  m_localFunctionProperties[domain].clearConstraints();
}

void MultiDomainFunctionModel::removeGlobalConstraints(
    std::string const &parameterName) {
  m_globalFunctionProperties.removeConstraints(parameterName);
}

void MultiDomainFunctionModel::removeLocalConstraints(
    std::string const &parameterName) {
  for (auto &properties : m_localFunctionProperties)
    properties.removeConstraints(parameterName);
}

} // namespace MantidWidgets
} // namespace MantidQt
