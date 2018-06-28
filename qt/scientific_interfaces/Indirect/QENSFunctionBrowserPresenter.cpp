#include "QENSFunctionBrowser.h"
#include "QENSFunctionBrowserPresenter.h"
#include "QENSIndexedFunctionModel.h"

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

QENSFunctionBrowserPresenter::QENSFunctionBrowserPresenter(
    QENSFunctionBrowser *browser, QENSIndexedFunctionModel *model)
    : m_browser(browser), m_model(model),
      m_subscriber(&QENSFunctionBrowserPresenter::g_defaultSubscriber) {
  browser->subscribe(this);
}

void QENSFunctionBrowserPresenter::subscribe(
    QENSFunctionBrowserPresenterSubscriber *subscriber) {
  m_subscriber = subscriber;
}

void QENSFunctionBrowserPresenter::addFunction(
    const std::string &name, const std::vector<std::size_t> &position) {
  m_model->addFunction(name, position);
  m_subscriber->functionChanged();
}

void QENSFunctionBrowserPresenter::removeFunction(
    const std::vector<std::size_t> &position) {
  m_model->removeFunction(position);
  m_subscriber->functionChanged();
}

void QENSFunctionBrowserPresenter::fixParameter(const std::string &name) {
  m_model->fixParameter(name);
}

void QENSFunctionBrowserPresenter::removeTie(const std::string &name) {
  m_model->removeLocalTie(name);
}

void QENSFunctionBrowserPresenter::addTie(const std::string &name,
                                          const std::string &expression) {
  m_model->addLocalTie(name, expression);
}

void QENSFunctionBrowserPresenter::addConstraints(const std::string &name,
                                                  double upperBound,
                                                  double lowerBound) {
  m_model->addBounds(name, upperBound, lowerBound);
}

void QENSFunctionBrowserPresenter::addConstraints10(const std::string &name) {
  m_model->addBoundsWithinPercentile(name, 0.1);
}

void QENSFunctionBrowserPresenter::addConstraints50(const std::string &name) {
  m_model->addBoundsWithinPercentile(name, 0.5);
}

void QENSFunctionBrowserPresenter::removeConstraint(const std::string &name) {}

void QENSFunctionBrowserPresenter::removeConstraints(const std::string &name) {
  m_model->removeLocalConstraint(name);
}

void QENSFunctionBrowserPresenter::setStringAttribute(
    const std::string &name, const std::string &value,
    const std::vector<std::size_t> &position) {
  m_model->setAttribute(name, value, position);
  m_subscriber->attributeChanged();
}

void QENSFunctionBrowserPresenter::setDoubleAttribute(
    const std::string &name, double value,
    const std::vector<std::size_t> &position) {
  m_model->setAttribute(name, value, position);
  m_subscriber->attributeChanged();
}

void QENSFunctionBrowserPresenter::setIntAttribute(
    const std::string &name, int value,
    const std::vector<std::size_t> &position) {
  m_model->setAttribute(name, value, position);
  m_subscriber->attributeChanged();
}

void QENSFunctionBrowserPresenter::setBoolAttribute(
    const std::string &name, bool value,
    const std::vector<std::size_t> &position) {
  m_model->setAttribute(name, value, position);
  m_subscriber->attributeChanged();
}

void QENSFunctionBrowserPresenter::setVectorDoubleAttribute(
    const std::string &name, const std::vector<double> &value,
    const std::vector<std::size_t> &position) {
  m_model->setAttribute(name, value, position);
  m_subscriber->attributeChanged();
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt
