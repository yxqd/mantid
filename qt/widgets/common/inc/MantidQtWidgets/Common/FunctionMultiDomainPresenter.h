// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTIDWIDGETS_FUNCTIONMULTIDOMAINPRESENTER_H_
#define MANTIDWIDGETS_FUNCTIONMULTIDOMAINPRESENTER_H_

#include "DllOption.h"

#include "MantidAPI/IFunction_fwd.h"
#include "MantidQtWidgets/Common/FunctionModel.h"

#include <QObject>
#include <boost/optional.hpp>
#include <memory>
#include <vector>

namespace MantidQt {
namespace MantidWidgets {

class MultiDomainFunctionModel;
class IFunctionView;
class EditLocalParameterDialog;

using namespace Mantid::API;

class EXPORT_OPT_MANTIDQT_COMMON FunctionMultiDomainPresenter : public QObject {
  Q_OBJECT
public:
  FunctionMultiDomainPresenter(IFunctionView *view);
  void clear();
  void setFunction(IFunction_sptr fun);
  void setFunctionString(const QString &funStr);
  QString getFunctionString() const;
  IFunction_sptr getFunction() const;
  IFunction_sptr getFunctionByIndex(const QString &index);
  IFunction_sptr getFitFunction() const;
  QString getFitFunctionString() const;
  bool hasFunction() const;
  void setParameter(const QString &paramName, double value);
  void setParamError(const QString &paramName, double value);
  double getParameter(const QString &paramName);
  bool isParameterFixed(const QString &parName) const;
  QString getParameterTie(const QString &parName) const;
  void updateParameters(const IFunction &fun);
  void updateMultiDatasetParameters(const IFunction &fun);
  void clearErrors();
  boost::optional<QString> currentFunctionIndex() const;
  void setNumberOfDatasets(int);
  void setDatasetNames(const QStringList &names);
  QStringList getDatasetNames() const;
  int getNumberOfDatasets() const;
  int getCurrentDataset() const;
  void setCurrentDataset(int);
  void removeDatasets(QList<int> indices);
  double getLocalParameterValue(const QString &parName, int i) const;
  bool isLocalParameterFixed(const QString &parName, int i) const;
  QString getLocalParameterTie(const QString &parName, int i) const;
  void setLocalParameterValue(const QString &parName, int i, double value);
  void setLocalParameterValue(const QString &parName, int i, double value,
                              double error);
  void setLocalParameterFixed(const QString &parName, int i, bool fixed);
  void setLocalParameterTie(const QString &parName, int i, QString tie);
  QStringList getGlobalParameters() const;
  void setGlobalParameters(const QStringList &globals);
  QStringList getLocalParameters() const;

  void setColumnSizes(int s0, int s1, int s2);
  void setErrorsEnabled(bool enabled);
signals:
  void functionStructureChanged();
  void parameterChanged(const QString &funcIndex, const QString &paramName);
private slots:
  void viewChangedParameter(const QString &parName);
  void viewPastedFunction(const QString &funStr);
  void viewAddedFunction(const QString &funStr);
  void viewRemovedFunction(const QString &functionIndex);
  void viewChangedTie(const QString &parName, const QString &tie);
  void viewAddedConstraint(const QString &functionIndex,
                           const QString &constraint);
  void viewRemovedConstraint(const QString &parName);
  void viewRequestedCopyToClipboard();
  void viewChangedGlobals(const QStringList &globalParameters);
  void editLocalParameter(const QString &parName);
  void editLocalParameterFinish(int result);

private:
  IFunctionView *m_view;
  std::unique_ptr<MultiDomainFunctionModel> m_model;
  EditLocalParameterDialog *m_editLocalParameterDialog;

public:
  IFunctionView *view() const { return m_view; }
};

} // namespace MantidWidgets
} // namespace MantidQt

#endif // MANTIDWIDGETS_FUNCTIONMULTIDOMAINPRESENTER_H_
