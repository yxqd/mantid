#ifndef MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSER_H_
#define MANTIDQTCUSTOMINTERFACESIDA_QENSFUNCTIONBROWSER_H_

#include "QENSFunctionBrowserSubscriber.h"

#include <QWidget>

#include <boost/optional.hpp>

/* Forward declarations */

class QtTreePropertyBrowser;
class QtGroupPropertyManager;
class QtDoublePropertyManager;
class QtIntPropertyManager;
class QtBoolPropertyManager;
class QtStringPropertyManager;
class QtEnumPropertyManager;
class QtProperty;
class QtBrowserItem;
class ParameterPropertyManager;

class QPushButton;
class QLabel;
class QLineEdit;
class QComboBox;
class QSignalMapper;
class QMenu;
class QAction;
class QTreeWidget;

namespace Mantid {
namespace API {
class CompositeFunction;
class Workspace;
class ParameterTie;
} // namespace API
} // namespace Mantid

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

class DLLExport QENSFunctionBrowser : public QWidget {
  Q_OBJECT
public:
  /// To keep QtProperty and its QtBrowserItem in one place
  struct AProperty {
    QtProperty *prop;
    QtBrowserItem *item;
    QtProperty *parent;
  };
  /// Tie structure
  struct ATie {
    QtProperty *paramProp; ///< Parameter property
    QString paramName;     ///< Parameter name
    QtProperty *tieProp;   ///< Tie property
  };
  /// Constraint structure
  struct AConstraint {
    QtProperty *paramProp; ///< Parameter property
    QtProperty *lower;     ///< Constraint property
    QtProperty *upper;     ///< Constraint property
  };

  QENSFunctionBrowser(QWidget *parent = nullptr);

  void subscribe(QENSFunctionBrowserSubscriber *subscriber);

  /// Add a function property
  AProperty addFunctionProperty(QtProperty *parent, const QString &funName);
  /// Add a parameter property
  AProperty addParameterProperty(QtProperty *parent, const QString &paramName,
                                 const QString &paramDesc, double paramValue);
  /// Add a attribute property
  AProperty addAttributeProperty(QtProperty *parent, QString attName,
                                 const QString &value);
  /// Add property showing function's index in the composite function
  AProperty addIndexProperty(QtProperty *prop);
  /// Add a constraint property
  QList<AProperty> addConstraintProperties(QtProperty *prop,
                                           QString constraint);

signals:
  void localParameterButtonClicked(const QString &);

protected:
  /// Check if property is a function group
  bool isFunction(QtProperty *prop) const;
  /// Check if a property is a tie
  bool isTie(QtProperty *prop) const;
  /// Check if a property is a constraint
  bool isConstraint(QtProperty *prop) const;
  /// Check if property is a function attribute
  bool isAttribute(QtProperty *prop) const;
  /// Check if property is a string attribute
  bool isStringAttribute(QtProperty *prop) const;
  /// Check if property is a double attribute
  bool isDoubleAttribute(QtProperty *prop) const;
  /// Check if property is a int attribute
  bool isIntAttribute(QtProperty *prop) const;
  /// Check if property is a bool attribute
  bool isBoolAttribute(QtProperty *prop) const;
  /// Check if property is a vector attribute
  bool isVectorAttribute(QtProperty *prop) const;
  /// Check if property is a function paramater
  bool isParameter(QtProperty *prop) const;
  /// Check if a property is an index
  bool isIndex(QtProperty *prop) const;

  std::vector<std::size_t> getFunctionPosition(QtProperty *prop) const;
  /// Get the function index for a property
  QString getIndex(QtProperty *prop) const;

protected slots:
  /// Show the context menu
  void popupMenu(const QPoint &);
  /// Add a function
  void addFunction();
  /// Remove a function
  void removeFunction();
  /// Fix a parameter
  void fixParameter();
  /// Unfix a parameter
  void removeTie();
  /// Add a tie to a parameter
  void addTie();
  /// Copy function from the clipboard
  void copyFromClipboard();
  /// Copy the function to the clipboard
  void copyToClipboard();
  /// Add both constraints to current parameter
  void addConstraints();
  /// Remove both constraints from current parameter
  void removeConstraints();
  /// Add both constraints to current parameter
  void addConstraints10();
  /// Add both constraints to current parameter
  void addConstraints50();
  /// Remove one of the constraints
  void removeConstraint();
  /// Update current function index depending on currently selected item
  void updateCurrentFunctionIndex();

  void stringAttributeChanged(QtProperty *, const QString &);
  void intAttributeChanged(QtProperty *, int);
  void doubleAttributeChanged(QtProperty *, double);
  void boolAttributeChanged(QtProperty *, bool);
  void vectorDoubleAttributeChanged(QtProperty *);
  void vectorSizeAttributeChanged(QtProperty *);

  /// Called when button in local parameter editor was clicked
  void parameterButtonClicked(QtProperty *);

private:
  void createBrowser();
  void createBrowser(QStringList &&options);
  void createManagers();
  void createEditorFactories();
  void connectManagerSignals();
  void createActions();

  std::string getNameOfProperty(QtProperty *prop) const;
  QtProperty *getFirstProperty() const;
  QtProperty *getSelectedProperty() const;
  QtProperty *getContainingFunctionProperty(QtProperty *prop) const;

  /// Add a sub-property
  AProperty addProperty(QtProperty *parent, QtProperty *subproperty);
  /// Remove and delete property
  void removeProperty(QtProperty *prop);
  /// Update function index properties
  void updateFunctionIndices(QtProperty *prop = nullptr, QString index = "");
  void updateFunctionIndex(QtProperty *prop, QString prefix, std::size_t &i);

  boost::optional<std::string> getTieFromDialog();
  std::string getFunctionFromUserDialog() const;

  /// Qt property browser which displays properties
  QtTreePropertyBrowser *m_browser;

  QENSFunctionBrowserSubscriber *m_subscriber;

  /// Manager for function group properties
  QtGroupPropertyManager *m_functionManager;
  /// Manager for function parameter properties
  ParameterPropertyManager *m_parameterManager;
  /// Manager for function string attribute properties
  QtStringPropertyManager *m_attributeStringManager;
  /// Manager for function double attribute properties
  QtDoublePropertyManager *m_attributeDoubleManager;
  /// Manager for function int attribute properties
  QtIntPropertyManager *m_attributeIntManager;
  /// Manager for function bool attribute properties
  QtBoolPropertyManager *m_attributeBoolManager;
  /// Manager for function index properties
  QtStringPropertyManager *m_indexManager;
  /// Manager for function tie properties
  QtStringPropertyManager *m_tieManager;
  /// Manager for parameter constraint properties
  QtStringPropertyManager *m_constraintManager;
  /// Manager for file name attributes
  QtStringPropertyManager *m_filenameManager;
  /// Manager for Formula attributes
  QtStringPropertyManager *m_formulaManager;
  /// Manager for Workspace attributes
  QtStringPropertyManager *m_workspaceManager;
  /// Manager for vector attribute properties
  QtGroupPropertyManager *m_attributeVectorManager;
  /// Manager for vector attribute member properties
  QtDoublePropertyManager *m_attributeVectorDoubleManager;
  /// Manager for vector attribute size properties
  QtIntPropertyManager *m_attributeSizeManager;

  /// Add a function
  QAction *m_actionAddFunction;
  /// Remove a function
  QAction *m_actionRemoveFunction;
  /// Fix a parameter
  QAction *m_actionFixParameter;
  /// Unfix a parameter
  QAction *m_actionRemoveTie;
  /// Add a custom tie to a parameter
  QAction *m_actionAddTie;
  /// Copy a function from the clipboard
  QAction *m_actionFromClipboard;
  /// Copy a function to the clipboard
  QAction *m_actionToClipboard;
  /// Add both constraints to current parameter with 10% spread
  QAction *m_actionConstraints10;
  /// Add both constraints to current parameter with 50% spread
  QAction *m_actionConstraints50;
  /// Add both constraints to current parameter
  QAction *m_actionConstraints;
  /// Remove both constraints from current parameter
  QAction *m_actionRemoveConstraints;
  /// Remove one constraints from current parameter
  QAction *m_actionRemoveConstraint;

  /// Store all properties in a map for easy access
  QMap<QtProperty *, AProperty> m_properties;
  /// Store parameter ties. Keys are function properties.
  QMultiMap<QtProperty *, ATie> m_ties;
  /// Store parameter constraints. Keys are function properties.
  QMultiMap<QtProperty *, AConstraint> m_constraints;
  /// Store index of function properties in parent.
  QMap<QtProperty *, std::size_t> m_propertyToIndex;
};

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt

#endif
