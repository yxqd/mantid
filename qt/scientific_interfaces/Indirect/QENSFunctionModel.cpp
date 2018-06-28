#include "QENSFunctionModel.h"

#include "MantidAPI/ConstraintFactory.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/IConstraint.h"
#include "MantidAPI/MatrixWorkspace.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>

using namespace Mantid::API;

namespace {
using namespace MantidQt::CustomInterfaces::IDA;

template <typename F>
void forEachFunctionEnumerated(const CompositeFunction &composite,
                               const F &functor) {
  for (auto i = 0u; i < composite.nFunctions(); ++i)
    functor(i, composite.getFunction(i));
}

template <typename F>
void forEachFunction(const CompositeFunction &composite, const F &functor) {
  for (auto i = 0u; i < composite.nFunctions(); ++i)
    functor(composite.getFunction(i));
}

std::unique_ptr<IConstraint> createConstraint(IFunction &function,
                                              const std::string &expression) {
  auto constraint =
      ConstraintFactory::Instance().createInitialized(&function, expression);
  return Mantid::Kernel::make_unique<IConstraint>(constraint);
}

IFunction_sptr createFunction(const std::string &name) {
  return FunctionFactory::Instance().createFunction(name);
}

CompositeFunction_sptr createComposite(const std::string &name) {
  return boost::dynamic_pointer_cast<CompositeFunction>(createFunction(name));
}

std::string functionIndex(std::size_t index) {
  return "f" + std::to_string(index) + ".";
}

std::string functionIndex(const std::vector<std::size_t> &indices) {
  std::string indexString = "";
  for (auto &&index : indices)
    indexString += functionIndex(index);
  return indexString;
}

std::string shortParameterName(const std::string &parameterName) {
  return parameterName.substr(parameterName.rfind("."));
}

bool isBackgroundFunction(IFunction_sptr function) {
  return function->category() == "Background";
}

bool isNotBackgroundFunction(IFunction_sptr function) {
  return !isBackgroundFunction(function);
}

bool isBackgroundTie(ParameterTie *tie) {
  return tie->getLocalFunction()->category() == "Background";
}

template <typename ExpressionCreator>
void addGlobalConstraints(MultiDomainFunction &function,
                          const ExpressionCreator &createExpression) {
  for (auto i = 0u; i < function.getNumberDomains(); ++i)
    function.addConstraint(createConstraint(function, createExpression(i)));
}

void addFunctionAt(CompositeFunction_sptr function, const std::string &name,
                   const std::vector<std::size_t> &position) {
  function = getFunctionAt(*function, position);
  if (auto composite = boost::dynamic_pointer_cast<CompositeFunction>(function))
    composite->addFunction(createFunction(name));
}

void removeFunctionAt(CompositeFunction_sptr function,
                      const std::vector<std::size_t> &position) {
  function = getFunctionAt(*function, position.begin(), position.end() - 1);
  auto composite = boost::dynamic_pointer_cast<CompositeFunction>(function);
  if (composite && !position.empty())
    composite->removeFunction(position.back());
}

CompositeFunction_sptr shallowCopy(CompositeFunction_sptr composite) {
  CompositeFunction_sptr copied(new CompositeFunction);
  forEachFunction(*composite, [&](IFunction_sptr function) {
    copied->addFunction(function);
  });
  copied->addConstraints(composite->writeConstraints());
  return copied;
}

IFunction_sptr getMinimalFormOfComposite(CompositeFunction_sptr composite) {
  const auto nFunctions = composite->nFunctions();
  if (nFunctions == 1)
    return composite->getFunction(0);
  else if (nFunctions > 1)
    return composite;
  return nullptr;
}

IFunction_sptr createResolutionFunction(const std::string &resolution) {
  auto function = FunctionFactory::Instance().createFunction("Resolution");
  IFunction::Attribute attribute(resolution);
  function->setAttribute("Workspace", attribute);
  return function;
}

CompositeFunction_sptr
createConvolutionFunction(IFunction_sptr function,
                          const std::string &resolution) {
  auto convolution = createComposite("Convolution");
  convolution->addFunction(createResolutionFunction(resolution));
  convolution->addFunction(function);
  return convolution;
}

IFunction_sptr
removeFunctions(CompositeFunction &composite,
                const std::vector<std::size_t> &functionIndices) {
  CompositeFunction_sptr removed(new CompositeFunction);
  for (auto &&index : functionIndices)
    removed->addFunction(composite.getFunction(index));

  for (auto i = functionIndices.rbegin(); i != functionIndices.rend(); ++i)
    composite.removeFunction(*i);
  return getMinimalFormOfComposite(removed);
}

template <typename Predicate>
std::vector<std::size_t> findFunctions(const CompositeFunction &composite,
                                       const Predicate &predicate) {
  std::vector<std::size_t> indices;
  forEachFunctionEnumerated(composite,
                            [&](std::size_t index, IFunction_sptr function) {
                              if (predicate(function))
                                functions.emplace_back(index);
                            });
  return indices;
}

template <typename Transform, typename OutputIterator>
void getTransformedParameterNames(const IFunction &function,
                                  const Transform &transform,
                                  OutputIterator outputIt) {
  for (auto i = 0u; i < function.nParams(); ++i)
    *outputIt++ = transform(function.parameterName(i));
}

template <typename OutputIterator>
void getParameterNamesWithPrefix(const IFunction &function,
                                 const std::string &prefix,
                                 OutputIterator outputIt) {
  getTransformedParameterNames(
      function, [&prefix](const std::string &name) { return prefix + name; },
      outputIt);
}

template <typename OutputIterator>
void getFunctionParameters(const CompositeFunction &composite,
                           std::size_t index, OutputIterator outputIt) {
  getParameterNamesWithPrefix(*composite.getFunction(index),
                              "f" + std::to_string(index) + ".", outputIt);
}

std::vector<std::vector<std::string>>
getFunctionParameters(const CompositeFunction &composite,
                      const std::vector<std::size_t> &functionIndices) {
  std::vector<std::vector<std::string>> parameters(functionIndices.size());
  for (auto i = 0u; i < functionIndices.size(); ++i)
    getFunctionParameters(composite, i, std::back_inserter(parameters[i]));
  return parameters;
}

std::vector<std::size_t>
getBackgroundIndices(const CompositeFunction &composite) {
  return findFunctions(composite, isBackgroundFunction);
}

std::vector<std::size_t> getModelIndices(const CompositeFunction &composite) {
  return findFunctions(composite, !isBackgroundFunction);
}

template <typename OutputIterator>
void mapParametersUsingPrefix(const IFunction &function,
                              const std::string &oldPrefix,
                              const std::string &newPrefix,
                              OutputIterator outputIt) {
  for (auto i = 0u; i < function.nParams(); ++i) {
    const auto name = function.parameterName(i);
    *outputIt++ = std::make_pair(oldPrefix + name, newPrefix + name);
  }
}

std::vector<std::pair<std::string, std::string>>
getParameterMapForPrefixChange(const IFunction &function,
                               const std::string &oldPrefix,
                               const std::string &newPrefix) {
  std::vector<std::pair<std::string, std::string>> parameterMap;
  parameterMap.reserve(function.nParams());
  mapParametersUsingPrefix(function, oldPrefix, newPrefix,
                           std::back_inserter(parameterMap));
  return parameterMap;
}

std::vector<std::pair<std::string, std::string>>
getParameterMapForMovedFunctions(const CompositeFunction &composite,
                                 const std::string &movedFrom,
                                 const std::vector<std::size_t> &functions,
                                 const std::string &movedTo) {
  if (functions.size() == 1)
    return getParameterMapForPrefixChange(
        *composite.getFunction(functions.front()), movedFrom, movedTo);

  std::vector<std::pair<std::string, std::string>> parameterMap;
  for (auto i = 0u; i < functions.size(); ++i) {
    const auto function = composite.getFunction(functions[i]);
    const auto oldPrefix = functionIndex(functions[i]) + movedFrom;
    const auto newPrefix = functionIndex(i) + movedTo;
    mapParametersUsingPrefix(*function, oldPrefix, newPrefix,
                             std::back_inserter(parameterMap));
  }
  return parameterMap;
}

std::vector<std::pair<std::string, std::string>>
getParameterMapForMovedFunctions(const CompositeFunction &composite,
                                 const std::vector<std::size_t> &movedFrom,
                                 const std::vector<std::size_t> &functions,
                                 const std::vector<std::size_t> &movedTo) {
  const auto function = getFunctionAt(composite, movedFrom);
  const auto innerComposite =
      boost::dynamic_pointer_cast<CompositeFunction>(function);

  if (innerComposite) {
    const auto initialPrefix = functionIndex(movedFrom);
    const auto movedPrefix = functionIndex(movedTo);
    return getParameterMapForMovedFunctions(*innerComposite, initialPrefix,
                                            functions, movedPrefix);
  }
  return std::vector<std::pair<std::string, std::string>>();
}

std::string replaceAll(
    const std::string &str,
    const std::vector<std::pair<std::string, std::string>> &replacementsMap) {
  auto result = str;
  for (auto &&replacementMap : replacementsMap)
    boost::replace_all(result, replacementMap.first, replacementMap.second);
  return result;
}

template <typename GetExpression>
std::vector<std::string> getModifiedParameterExpressions(
    const IFunction &function,
    const std::vector<std::pair<std::string, std::string>> &parameterMap,
    const GetExpression &getExpression) {
  std::vector<std::string> expressions;
  for (auto i = 0u; i < function.nParams(); ++i) {
    const auto originalExpression = getExpression(function, i);
    const auto expression = replaceAll(originalExpression, parameterMap);
    if (expression != originalExpression)
      expressions.emplace_back(expression);
  }
  return expressions;
}

std::string getTie(const IFunction &function, std::size_t index) {
  return function.getTie(index)->asString();
}

std::string getConstraint(const IFunction &function, std::size_t index) {
  return function.getConstraint(index)->asString();
}

std::string getModifiedTies(
    const IFunction &function,
    const std::vector<std::pair<std::string, std::string>> &parameterMap) {
  const auto modified =
      getModifiedParameterExpressions(function, parameterMap, getTie);
  return boost::algorithm::join(modified, ",");
}

std::string getModifiedConstraints(
    const IFunction &function,
    const std::vector<std::pair<std::string, std::string>> &parameterMap) {
  const auto modified =
      getModifiedParameterExpressions(function, parameterMap, getConstraint);
  return boost::algorithm::join(modified, ",");
}

void addModifiedTiesAndConstraints(
    IFunction &newFunction, const CompositeFunction &oldFunction,
    const std::vector<std::size_t> &movedFrom,
    const std::vector<std::size_t> &functionIndices,
    const std::vector<std::size_t> &movedTo) {
  const auto parameterMap = getParameterMapForMovedFunctions(
      oldFunction, movedFrom, functionIndices, movedTo);
  newFunction.addTies(getModifiedTies(oldFunction, parameterMap));
  newFunction.addConstraints(getModifiedTies(oldFunction, parameterMap));
}

void addModifiedBackgroundTiesAndConstraints(
    IFunction &newFunction, const CompositeFunction &oldFunction,
    const std::vector<std::size_t> &backgroundIndices) {
  addModifiedTiesAndConstraints(newFunction, oldFunction, {}, backgroundIndices,
                                {0});
}

CompositeFunction_sptr convolveWith(IFunction_sptr model,
                                    IFunction_sptr background,
                                    const std::string &resolution) {
  CompositeFunction_sptr composite(new CompositeFunction);
  composite->addFunction(background);
  composite->addFunction(createConvolutionFunction(model, resolution));
  return composite;
}

CompositeFunction_sptr
convolveWith(CompositeFunction_sptr composite, const std::string &resolution,
             const std::vector<std::size_t> &backgroundIndices) {
  if (backgroundIndices.empty())
    return createConvolutionFunction(composite, resolution);

  auto model = shallowCopy(composite);
  auto background = removeFunctions(*composite, backgroundIndices);
  auto function = convolveWith(composite, background, resolution);
  addModifiedBackgroundTiesAndConstraints(*function, *composite,
                                          backgroundIndices);
  return function;
}

CompositeFunction_sptr
convolveWith(IFunction_sptr function, const std::string &resolution,
             const std::vector<std::size_t> &backgroundIndices) {
  if (auto composite = boost::dynamic_pointer_cast<CompositeFunction>(function))
    return convolveWith(composite, resolution, backgroundIndices);
  return createConvolutionFunction(function, resolution);
}

CompositeFunction_sptr
convolveWith(boost::shared_ptr<MultiDomainFunction> function,
             const std::vector<std::string> &resolution,
             const std::vector<std::size_t> &backgroundIndices) {
  const auto composite = getCompositeAt(function, 0);
  const auto modelIndices = getModelIndices(*composite);

  boost::shared_ptr<MultiDomainFunction> result(new MultiDomainFunction);
  for (auto i = 0; i < function->nFunctions(); ++i) {
    result->addFunction(convolveWith(function->getFunction(i), resolution[i],
                                     backgroundIndices));

    if (backgroundIndices.empty())
      addModifiedTiesAndConstraints(*result, *function, {i}, modelIndices,
                                    {i, 1});
    else {
      addModifiedTiesAndConstraints(*result, *function, {i}, modelIndices,
                                    {i, 1, 1});
      addModifiedTiesAndConstraints(*result, *function, {i}, backgroundIndices,
                                    {i, 0});
    }
  }
  return result;
}

CompositeFunction_sptr
convolveWith(boost::shared_ptr<MultiDomainFunction> function,
             const std::vector<std::string> &resolution) {
  if (function->nFunctions() == 0)
    return function;
  else if (const auto composite = getCompositeAt(function, 0))
    return convolveWith(function, resolution, getBackgroundIndices(*composite));
}
} // namespace

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

QENSFunctionModel::QENSFunctionModel() : m_function(new MultiDomainFunction) {}

IFunction_sptr QENSFunctionModel::getFunction() const {
  return getMinimalFormOfComposite(m_function);
}

IFunction_sptr QENSFunctionModel::getConvolutionFunction(
    const std::vector<std::string> &resolution) const {
  const auto nDomains = m_function->getNumberDomains();
  if (nDomains != resolution.size())
    throw std::runtime_error("Number of provided resolution files must equal "
                             "the number of domains.");
  else if (nDomains == 1)
    return convolveWith(getLocalFunction(0), resolution.front());
}

IFunction_sptr QENSFunctionModel::getLocalFunction(std::size_t domain) const {
  return m_function->getFunction(domain);
}

Mantid::API::IFunction_sptr QENSFunctionModel::getNewLocalFunction() const {
  if (m_function->nFunctions() == 0)
    return boost::shared_ptr<IFunction>(new CompositeFunction);
  return getLocalFunction(0)->clone();
}

std::size_t QENSFunctionModel::numberOfLocalParameters() const {
  if (m_function->nFunctions() == 0)
    return 0;
  return getLocalFunction(0)->nParams();
}

void QENSFunctionModel::addFunction(const std::string &name,
                                    const std::vector<std::size_t> &position) {
  forEachFunction(*m_function, [&](IFunction_sptr function) {
    addFunctionAt(boost::dynamic_pointer_cast<CompositeFunction>(function),
                  name, position);
  });
}

void QENSFunctionModel::removeFunction(
    const std::vector<std::size_t> &position) {
  forEachFunction(*m_function, [&](IFunction_sptr function) {
    removeFunctionAt(boost::dynamic_pointer_cast<CompositeFunction>(function),
                     position);
  });
}

void QENSFunctionModel::setParameterValue(const std::string &name,
                                          const ParameterValue &value) {
  const auto index = m_function->parameterIndex(name);
  m_function->setParameter(index, value.value);
  if (value.error)
    m_function->setError(index, *value.error);
}

void QENSFunctionModel::addDataset(MatrixWorkspace_sptr workspace) {
  for (auto i = 0u; i < workspace->getNumberHistograms(); ++i)
    addDatasetDomain(workspace, i);
}

void QENSFunctionModel::addDataset(
    MatrixWorkspace_sptr workspace,
    const std::pair<std::size_t, std::size_t> &indexRange) {
  for (auto i = indexRange.first; i <= indexRange.second; ++i)
    addDatasetDomain(workspace, i);
}

void QENSFunctionModel::addDatasetDomain(MatrixWorkspace_sptr workspace,
                                         std::size_t workspaceIndex) {
  const auto index = m_function->nFunctions();
  const auto function = getNewLocalFunction();
  const auto &x = workspace->x(0);

  function->setMatrixWorkspace(workspace, workspaceIndex, x.front(), x.back());
  m_function->addFunction(function);
  m_function->setDomainIndex(index, index);
}

void QENSFunctionModel::removeDataset(std::size_t index) {
  m_function->removeFunction(index);
  for (auto i = index; i < m_function->nFunctions(); ++i)
    m_function->setDomainIndex(index, index);
}

void QENSFunctionModel::clearDatasets() { m_function->clear(); }

void QENSFunctionModel::fixParameterInDomain(const std::string &parameterName,
                                             std::size_t domain) {
  getLocalFunction(domain)->fixParameter(parameterName);
}

void QENSFunctionModel::unfixParameterInDomain(const std::string &parameterName,
                                               std::size_t domain) {
  getLocalFunction(domain)->unfixParameter(parameterName);
}

void QENSFunctionModel::addLocalTieToDomain(const std::string &parameterName,
                                            const std::string &expression,
                                            std::size_t domain) {
  getLocalFunction(domain)->tie(parameterName, expression);
}

void QENSFunctionModel::removeLocalTieFromDomain(
    const std::string &parameterName, std::size_t domain) {
  getLocalFunction(domain)->removeTie(parameterName);
}

void QENSFunctionModel::removeLocalTiesFromDomain(std::size_t domain) {
  getLocalFunction(domain)->clearTies();
}

void QENSFunctionModel::addEqualityGlobalTie(const std::string &parameterName) {
  const auto firstParameter = "f0." + parameterName;
  addGlobalTie(parameterName, firstParameter, 1,
               m_function->getNumberDomains());
}

void QENSFunctionModel::addGlobalTie(const std::string &parameterName,
                                     const std::string &expression) {
  addGlobalTie(parameterName, expression, 0, m_function->getNumberDomains());
}

void QENSFunctionModel::addGlobalTie(const std::string &parameterName,
                                     const std::string &expression,
                                     std::size_t fromDomain,
                                     std::size_t toDomain) {
  for (auto i = fromDomain; i < toDomain; ++i)
    m_function->tie("f" + std::to_string(i) + "." + parameterName, expression);
}

void QENSFunctionModel::removeGlobalTies(const std::string &parameterName) {
  const auto nParams = numberOfLocalParameters();
  if (nParams == 0)
    return;
  const auto start = getLocalFunction(0)->parameterIndex(parameterName);
  for (auto i = start; i < m_function->nParams(); i += nParams)
    m_function->removeTie(i);
}

void QENSFunctionModel::clearTies() { m_function->clearTies(); }

void QENSFunctionModel::addUpperBoundToDomain(const std::string &parameterName,
                                              double bound,
                                              std::size_t domain) {
  const auto function = getLocalFunction(domain);
  const auto expression = parameterName + "<" + std::to_string(bound);
  function->removeConstraint(parameterName);
  function->addConstraint(createConstraint(*function, expression));
}

void QENSFunctionModel::addLowerBoundToDomain(const std::string &parameterName,
                                              double bound,
                                              std::size_t domain) {
  const auto function = getLocalFunction(domain);
  const auto expression = std::to_string(bound) + "<" + parameterName;
  function->removeConstraint(parameterName);
  function->addConstraint(createConstraint(*function, expression));
}

void QENSFunctionModel::addBoundsToDomain(const std::string &parameterName,
                                          double lowerBound, double upperBound,
                                          std::size_t domain) {
  const auto function = getLocalFunction(domain);
  const auto expression = std::to_string(lowerBound) + "<" + parameterName +
                          "<" + std::to_string(upperBound);
  function->removeConstraint(parameterName);
  function->addConstraint(createConstraint(*function, expression));
}

void QENSFunctionModel::addBoundsToDomainWithinPercentile(
    const std::string &parameterName, double percentile, std::size_t domain) {
  const auto value = m_function->getParameter(parameterName);
  const auto lower = value * (1.0 - percentile);
  const auto upper = value * (1.0 + percentile);
  addBoundsToDomain(parameterName, lower, upper, domain);
}

void QENSFunctionModel::addGlobalUpperBound(const std::string &parameterName,
                                            double bound) {
  const auto expression = [&](std::size_t i) {
    return functionIndex(i) + parameterName + "<" + std::to_string(bound);
  };
  addGlobalConstraints(*m_function, expression);
}

void QENSFunctionModel::addGlobalLowerBound(const std::string &parameterName,
                                            double bound) {
  const auto expression = [&](std::size_t i) {
    return std::to_string(bound) + "<" + functionIndex(i) + parameterName;
  };
  addGlobalConstraints(*m_function, expression);
}

void QENSFunctionModel::addGlobalBounds(const std::string &parameterName,
                                        double lowerBound, double upperBound) {
  const auto expression = [&](std::size_t i) {
    return std::to_string(lowerBound) + "<" + functionIndex(i) + parameterName +
           "<" + std::to_string(upperBound);
  };
  addGlobalConstraints(*m_function, expression);
}

void QENSFunctionModel::removeLocalConstraintFromDomain(
    const std::string &parameterName, std::size_t domain) {
  getLocalFunction(domain)->removeConstraint(parameterName);
}

void QENSFunctionModel::removeLocalConstraintsFromDomain(std::size_t domain) {
  const auto function = getLocalFunction(domain);
  const auto parameters = function->getParameterNames();
  for (auto &&parameter : parameters)
    function->removeConstraint(parameter);
}

void QENSFunctionModel::removeGlobalConstraints(
    const std::string &parameterName) {
  for (auto i = 0u; i < m_function->nFunctions(); ++i)
    m_function->getFunction(i)->removeConstraint(parameterName);
}

void QENSFunctionModel::clearBoundConstraints() {
  m_function->clearConstraints();
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt
