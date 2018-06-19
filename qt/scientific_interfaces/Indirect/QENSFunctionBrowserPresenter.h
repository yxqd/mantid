#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSER_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSER_H_

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

class QENSFunctionBrowserPresenter {
public:
  QENSFunctionBrowserPresenter(FunctionBrowser *m_browser, QENSFunctionModel *m_model);

private:
  QENSFunctionModel *m_model;
  FunctionBrowser *m_browser;
};

#endif

}
}
}
