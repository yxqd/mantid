#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSERPRESENTER_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSERPRESENTER_H_

#include "QENSFunctionBrowserSubscriber.h"

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

class QENSFunctionBrowser;
class QENSIndexedFunctionModel;

class QENSFunctionBrowserPresenterSubscriber {
public:
  virtual void functionChanged() = 0;
  virtual void parameterValueChanged() = 0;
  virtual void attributeChanged() = 0;
};

class QENSFunctionBrowserPresenter : public QENSFunctionBrowserSubscriber {
public:
  QENSFunctionBrowserPresenter(QENSFunctionBrowser *m_browser,
                               QENSIndexedFunctionModel *m_model);

  void subscribe(QENSFunctionBrowserPresenterSubscriber *subscriber);

  void addFunction(const std::string &name,
                   const std::vector<std::size_t> &position) override;
  void removeFunction(const std::vector<std::size_t> &position) override;
  void fixParameter(const std::string &name) override;
  void removeTie(const std::string &name) override;
  void addTie(const std::string &name, const std::string &expression) override;
  void addConstraints(const std::string &name, double upperBound,
                      double lowerBound) override;
  void addConstraints10(const std::string &name) override;
  void addConstraints50(const std::string &name) override;
  void removeConstraint(const std::string &name) override;
  void removeConstraints(const std::string &name) override;

  void setStringAttribute(const std::string &name, const std::string &value,
                          const std::vector<std::size_t> &position) override;
  void setDoubleAttribute(const std::string &name, double value,
                          const std::vector<std::size_t> &position) override;
  void setIntAttribute(const std::string &name, int value,
                       const std::vector<std::size_t> &position) override;
  void setBoolAttribute(const std::string &name, bool value,
                        const std::vector<std::size_t> &position) override;
  void
  setVectorDoubleAttribute(const std::string &name,
                           const std::vector<double> &value,
                           const std::vector<std::size_t> &position) override;

private:
  QENSFunctionBrowserPresenterSubscriber *m_subscriber;
  QENSIndexedFunctionModel *m_model;
  QENSFunctionBrowser *m_browser;

  class EmptySubscriber : public QENSFunctionBrowserPresenterSubscriber {
  public:
    void functionChanged() override {}
    void parameterValueChanged() override {}
    void attributeChanged() override {}
  };
  static EmptySubscriber g_defaultSubscriber;
};

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt

#endif
