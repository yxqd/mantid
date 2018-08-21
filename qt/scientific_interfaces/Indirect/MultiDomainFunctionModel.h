#ifndef MANTIDWIDGETS_MULTIDOMAINFUNCTIONMODEL_H_
#define MANTIDWIDGETS_MULTIDOMAINFUNCTIONMODEL_H_

#include "FunctionProperties.h"
#include "IFunctionModel.h"

namespace Mantid {
namespace API {
class CompositeFunction;
} // namespace API
} // namespace Mantid

namespace MantidQt {
namespace MantidWidgets {

class MultiDomainFunctionModel : public IFunctionModel {
public:
  MultiDomainFunctionModel();

  virtual boost::shared_ptr<Mantid::API::IFunction> getFitFunction() const;

  std::size_t numberOfParameters() const override;
  std::size_t numberOfDomains() const;
  std::string parameterName(std::size_t index) const override;
  double parameterValue(std::string const &parameter) const override;
  boost::optional<double>
  parameterError(std::string const &parameter) const override;
  boost::optional<std::string>
  parameterTie(std::string const &parameter) const override;
  boost::optional<double>
  parameterLowerBound(std::string const &name) const override;
  boost::optional<double>
  parameterUpperBound(std::string const &name) const override;

  std::vector<std::string> getAttributeNames() const override;
  Mantid::API::IFunction::Attribute
  getAttribute(std::string const &name) override;

  Mantid::API::MatrixWorkspace_sptr getWorkspace() const;
  std::string getWorkspaceName() const;
  std::size_t getWorkspaceIndex() const;

  bool isComposite(std::vector<std::size_t> const &position) const override;
  std::size_t
  numberOfFunctionsAt(std::vector<std::size_t> const &position) const override;

  bool isParameterTied(std::string const &name) const override;
  bool isParameterFixed(std::string const &name) const override;
  bool isParameterConstrained(std::string const &name) const override;

  virtual std::string getLocalFunctionString() const override;
  std::size_t numberOfLocalParameters() const;
  std::size_t getActiveDomain() const;
  void setActiveDomain(std::size_t domain);

  void setFunction(std::string const &functionString) override;
  std::size_t addFunction(std::string const &name,
                          std::vector<std::size_t> const &position) override;
  void removeFunction(std::vector<std::size_t> const &position) override;

  void setStringAttribute(std::string const &name,
                          std::string const &value) override;
  void setDoubleAttribute(std::string const &name, double value) override;
  void setIntAttribute(std::string const &name, int value) override;
  void setBoolAttribute(std::string const &name, bool value) override;
  void setVectorAttribute(std::string const &name,
                          std::vector<double> const &value) override;
  void setVectorAttributeSize(std::string const &name,
                              std::size_t size) override;

  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace);
  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace,
                  std::pair<std::size_t, std::size_t> const &indexRange);
  virtual void removeDatasetDomain(std::size_t domain);
  void clearDatasets();

  template <typename BeginIterator, typename EndIterator>
  void addDataset(Mantid::API::MatrixWorkspace_sptr workspace,
                  BeginIterator startIndexIt, EndIterator endIndexIt);

  void addEqualityGlobalTie(std::string const &parameterName);
  void addGlobalTie(std::string const &parameterName,
                    std::string const &expression);
  void removeGlobalTies(std::string const &parameterName);
  void removeLocalTies(std::string const &parameterName);
  void clearTies();

  void addGlobalUpperBound(std::string const &parameterName, double bound);
  void addGlobalLowerBound(std::string const &parameterName, double bound);
  void addGlobalBounds(std::string const &parameterName, double lowerBound,
                       double upperBound);

  void removeGlobalConstraints(std::string const &parameterName);
  void removeLocalConstraints(std::string const &parameterName);

  boost::shared_ptr<Mantid::API::IFunction>
  getLocalFunction(std::size_t index) const;

  void setLocalParameterValue(std::string const &parameterName,
                              double value) override;
  void fixLocalParameter(std::string const &parameterName) override;
  void unfixLocalParameter(std::string const &parameterName) override;

  void setLocalTie(std::string const &parameterName,
                   std::string const &expression) override;
  void removeLocalTie(std::string const &parameterName) override;
  void removeLocalTies() override;

  void addLocalUpperBound(std::string const &parameterName,
                          double bound) override;
  void addLocalLowerBound(std::string const &parameterName,
                          double bound) override;
  void addLocalBounds(std::string const &parameterName, double lowerBound,
                      double upperBound) override;
  void addLocalBoundsWithinPercentile(std::string const &parameterName,
                                      double percentile) override;

  void removeLocalConstraint(std::string const &parameterName,
                             std::string const &type) override;
  void removeLocalConstraints(std::string const &parameterName) override;

protected:
  void setParameterValue(std::size_t domain, std::string const &name,
                         ParameterValue const &value);

  void fixParameterInDomain(std::string const &parameterName,
                            std::size_t domain);
  void unfixParameterInDomain(std::string const &parameterName,
                              std::size_t domain);

  void addLocalTieToDomain(std::string const &parameterName,
                           std::string const &expression, std::size_t domain);
  void removeLocalTieFromDomain(std::string const &parameterName,
                                std::size_t domain);
  void removeLocalTiesFromDomain(std::size_t domain);

  void addUpperBoundToDomain(std::string const &parameterName, double bound,
                             std::size_t domain);
  void addLowerBoundToDomain(std::string const &parameterName, double bound,
                             std::size_t domain);
  void addBoundsToDomain(std::string const &parameterName, double lowerBound,
                         double upperBound, std::size_t domain);
  void addBoundsToDomainWithinPercentile(std::string const &parameterName,
                                         double percentile, std::size_t domain);

  void removeLocalConstraintsFromDomain(std::string const &parameterName,
                                        std::size_t domain);
  void clearLocalConstraintsFromDomain(std::size_t domain);

protected:
  virtual void addDatasetDomain(Mantid::API::MatrixWorkspace_sptr workspace,
                                std::size_t workspaceIndex);

private:
  void removeTiesContainingParameter(std::string const &parameter);
  void removeConstraintsContainingParameter(std::string const &parameter);
  void removeParameterProperties(std::string const &parameterName);
  void addGlobalTie(std::string const &parameterName,
                    std::string const &expression, std::size_t fromDomain,
                    std::size_t toDomain);

  std::vector<LocalFunctionProperties> m_localFunctionProperties;
  FunctionProperties m_globalFunctionProperties;
  std::size_t m_activeDomain;
  boost::shared_ptr<Mantid::API::CompositeFunction> m_function;
};

template <typename BeginIterator, typename EndIterator>
void MultiDomainFunctionModel::addDataset(
    Mantid::API::MatrixWorkspace_sptr workspace, BeginIterator startIndexIt,
    EndIterator endIndexIt) {
  for (auto i = startIndexIt; i < endIndexIt; ++i)
    addDatasetDomain(workspace, *i);
}

} // namespace MantidWidgets
} // namespace MantidQt

#endif
