#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSINDEXEDFUNCTIONMODEL_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSINDEXEDFUNCTIONMODEl_H_

#include "QENSFunctionModel.h"

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

class QENSIndexedFunctionModel : public QENSFunctionModel {
public:
  QENSIndexedFunctionModel();

  void setActiveDomain(std::size_t activeDomain);

  void fixParameter(const std::string &parameterName);
  void unfixParameter(const std::string &parameterName);

  void addLocalTie(const std::string &parameterName,
                   const std::string &expression);
  void removeLocalTie(const std::string &parameterName);
  void removeLocalTies();

  template <typename Predicate>
  void constrainLocalParametersTo(const std::string &expression,
                                  const Predicate &useInConstraint);

  template <typename Functor> void forEachParameter(const Functor &functor);

  void addUpperBound(const std::string &parameterName, double bound);
  void addLowerBound(const std::string &parameterName, double bound);
  void addBounds(const std::string &parameterName, double lowerBound,
                 double upperBound);
  void addBoundsWithinPercentile(const std::string &parameterName,
                                 double percentile);

  void removeLocalConstraint(const std::string &parameterName);
  void removeLocalConstraints();

private:
  std::size_t m_activeDomain;
};

template <typename Predicate>
void QENSIndexedFunctionModel::constrainLocalParametersTo(
    const std::string &expression, const Predicate &useInConstraint) {
  constrainLocalParametersInDomain(expression, useInConstraint, m_activeDomain);
}

template <typename Functor>
void QENSIndexedFunctionModel::forEachParameter(const Functor &functor) {
  forEachParameterInDomain(functor, m_activeDomain);
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt

#endif
