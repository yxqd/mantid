#include "QENSIndexedFunctionModel.h"

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

QENSIndexedFunctionModel::QENSIndexedFunctionModel() : m_activeDomain(0) {}

void QENSIndexedFunctionModel::setActiveDomain(std::size_t activeDomain) {
  m_activeDomain = activeDomain;
}

void QENSIndexedFunctionModel::fixParameter(const std::string &parameterName) {
  fixParameterInDomain(parameterName, m_activeDomain);
}

void QENSIndexedFunctionModel::unfixParameter(
    const std::string &parameterName) {
  unfixParameterInDomain(parameterName, m_activeDomain);
}

void QENSIndexedFunctionModel::addLocalTie(const std::string &parameterName,
                                           const std::string &expression) {
  addLocalTieToDomain(parameterName, expression, m_activeDomain);
}

void QENSIndexedFunctionModel::removeLocalTie(
    const std::string &parameterName) {
  removeLocalTieFromDomain(parameterName, m_activeDomain);
}

void QENSIndexedFunctionModel::removeLocalTies() {
  removeLocalTiesFromDomain(m_activeDomain);
}

void QENSIndexedFunctionModel::addUpperBound(const std::string &parameterName,
                                             double bound) {
  addUpperBoundToDomain(parameterName, bound, m_activeDomain);
}

void QENSIndexedFunctionModel::addLowerBound(const std::string &parameterName,
                                             double bound) {
  addLowerBoundToDomain(parameterName, bound, m_activeDomain);
}

void QENSIndexedFunctionModel::addBounds(const std::string &parameterName,
                                         double lowerBound, double upperBound) {
  addBoundsToDomain(parameterName, lowerBound, upperBound, m_activeDomain);
}

void QENSIndexedFunctionModel::addBoundsWithinPercentile(
    const std::string &parameterName, double percentile) {
  addBoundsToDomainWithinPercentile(parameterName, percentile, m_activeDomain);
}

void QENSIndexedFunctionModel::removeLocalConstraint(
    const std::string &parameterName) {
  removeLocalConstraintFromDomain(parameterName, m_activeDomain);
}

void QENSIndexedFunctionModel::removeLocalConstraints() {
  removeLocalConstraintsFromDomain(m_activeDomain);
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt