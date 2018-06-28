#include "QENSFunctionBrowser.h"

#include "MantidAPI/Expression.h"

#include "MantidQtWidgets/Common/QtPropertyBrowser/DoubleDialogEditor.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/FilenameDialogEditor.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/FormulaDialogEditor.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/WorkspaceEditorFactory.h"
#include "MantidQtWidgets/Common/SelectFunctionDialog.h"
#include "MantidQtWidgets/Common/UserFunctionDialog.h"

#include "MantidQtWidgets/Common/QtPropertyBrowser/CompositeEditorFactory.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/DoubleEditorFactory.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/qteditorfactory.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/qtpropertymanager.h"
#include "MantidQtWidgets/Common/QtPropertyBrowser/qttreepropertybrowser.h"

#include <boost/lexical_cast.hpp>

#include <QInputDialog>

namespace {
const char *globalOptionName = "Global";

template <typename T, typename PropertyManager>
std::vector<T> getVectorFromProperty(QtProperty *prop,
                                     PropertyManager *manager) const {
  const auto members = prop->subProperties();
  if (members.empty())
    return std::vector<T>();

  std::vector<T> values;
  values.reserve(static_cast<std::size_t>(members.size()));
  for (auto &&member : members)
    values.emplace_back(manager->value(member));
  return values;
}
} // namespace

using namespace MantidQt::MantidWidgets;

namespace MantidQt {
namespace CustomInterfaces {
namespace IDA {

QENSFunctionBrowser::QENSFunctionBrowser(QWidget *parent) : QWidget(parent) {}

void QENSFunctionBrowser::subscribe(QENSFunctionBrowserSubscriber *subscriber) {
  m_subscriber = subscriber;
}

void QENSFunctionBrowser::createBrowser() {
  createBrowser({globalOptionName});
  createManagers();
  createEditorFactories();

  m_browser->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(m_browser, SIGNAL(customContextMenuRequested(const QPoint &)), this,
          SLOT(popupMenu(const QPoint &)));
  connect(m_browser, SIGNAL(optionChanged(QtProperty *, const QString &, bool)),
          this, SLOT(globalChanged(QtProperty *, const QString &, bool)));

  connectManagerSignals();

  connect(m_browser, SIGNAL(currentItemChanged(QtBrowserItem *)),
          SLOT(updateCurrentFunctionIndex()));

  m_browser->setFocusPolicy(Qt::StrongFocus);
}

void QENSFunctionBrowser::createBrowser(QStringList &&options) {
  m_browser = new QtTreePropertyBrowser(nullptr, options);
}

/**
 * Create property managers: they create, own properties, get and set values
 */
void QENSFunctionBrowser::createManagers() {
  m_functionManager = new QtGroupPropertyManager(this);
  m_parameterManager = new ParameterPropertyManager(this);
  m_attributeStringManager = new QtStringPropertyManager(this);
  m_attributeDoubleManager = new QtDoublePropertyManager(this);
  m_attributeIntManager = new QtIntPropertyManager(this);
  m_attributeBoolManager = new QtBoolPropertyManager(this);
  m_indexManager = new QtStringPropertyManager(this);
  m_tieManager = new QtStringPropertyManager(this);
  m_constraintManager = new QtStringPropertyManager(this);
  m_filenameManager = new QtStringPropertyManager(this);
  m_formulaManager = new QtStringPropertyManager(this);
  m_workspaceManager = new QtStringPropertyManager(this);
  m_attributeVectorManager = new QtGroupPropertyManager(this);
  m_attributeSizeManager = new QtIntPropertyManager(this);
  m_attributeVectorDoubleManager = new QtDoublePropertyManager(this);
}

/**
 * Creates the editor factories.
 */
void QENSFunctionBrowser::createEditorFactories() {
  QtSpinBoxFactory *spinBoxFactory = new QtSpinBoxFactory(this);
  DoubleEditorFactory *doubleEditorFactory = new DoubleEditorFactory(this);
  ParameterEditorFactory *paramEditorFactory = new ParameterEditorFactory(this);

  QtAbstractEditorFactory<ParameterPropertyManager> *parameterEditorFactory(
      nullptr);

  auto buttonFactory = new DoubleDialogEditorFactory(this);
  auto compositeFactory =
      new CompositeEditorFactory<ParameterPropertyManager>(this, buttonFactory);
  compositeFactory->setSecondaryFactory(globalOptionName, paramEditorFactory);
  parameterEditorFactory = compositeFactory;
  connect(buttonFactory, SIGNAL(buttonClicked(QtProperty *)), this,
          SLOT(parameterButtonClicked(QtProperty *)));
  connect(buttonFactory, SIGNAL(closeEditor()), m_browser, SLOT(closeEditor()));

  QtLineEditFactory *lineEditFactory = new QtLineEditFactory(this);
  QtCheckBoxFactory *checkBoxFactory = new QtCheckBoxFactory(this);
  FilenameDialogEditorFactory *filenameDialogEditorFactory =
      new FilenameDialogEditorFactory(this);
  FormulaDialogEditorFactory *formulaDialogEditFactory =
      new FormulaDialogEditorFactory(this);
  WorkspaceEditorFactory *workspaceEditorFactory =
      new WorkspaceEditorFactory(this);

  // assign factories to property managers
  m_browser->setFactoryForManager(m_parameterManager, parameterEditorFactory);
  m_browser->setFactoryForManager(m_attributeStringManager, lineEditFactory);
  m_browser->setFactoryForManager(m_attributeDoubleManager,
                                  doubleEditorFactory);
  m_browser->setFactoryForManager(m_attributeIntManager, spinBoxFactory);
  m_browser->setFactoryForManager(m_attributeBoolManager, checkBoxFactory);
  m_browser->setFactoryForManager(m_indexManager, lineEditFactory);
  m_browser->setFactoryForManager(m_tieManager, lineEditFactory);
  m_browser->setFactoryForManager(m_constraintManager, lineEditFactory);
  m_browser->setFactoryForManager(m_filenameManager,
                                  filenameDialogEditorFactory);
  m_browser->setFactoryForManager(m_formulaManager, formulaDialogEditFactory);
  m_browser->setFactoryForManager(m_workspaceManager, workspaceEditorFactory);
  m_browser->setFactoryForManager(m_attributeSizeManager, spinBoxFactory);
  m_browser->setFactoryForManager(m_attributeVectorDoubleManager,
                                  doubleEditorFactory);
}

void QENSFunctionBrowser::connectManagerSignals() {
  connect(m_attributeStringManager,
          SIGNAL(valueChanged(QtProperty *, const QString &)), this,
          SLOT(stringAttributeChanged(QtProperty *, const QString &)));
  connect(m_attributeDoubleManager, SIGNAL(valueChanged(QtProperty *, double)),
          this, SLOT(doubleAttributeChanged(QtProperty *, double)));
  connect(m_attributeIntManager, SIGNAL(valueChanged(QtProperty *, int)), this,
          SLOT(intAttributeChanged(QtProperty *, int)));
  connect(m_attributeBoolManager, SIGNAL(valueChanged(QtProperty *, bool)),
          this, SLOT(boolAttributeChanged(QtProperty *, bool)));
  connect(m_formulaManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(attributeChanged(QtProperty *)));
  connect(m_filenameManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(attributeChanged(QtProperty *)));
  connect(m_workspaceManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(attributeChanged(QtProperty *)));
  connect(m_attributeVectorDoubleManager, SIGNAL(propertyChanged(QtProperty *)),
          this, SLOT(vectorAttributeChanged(QtProperty *)));
  connect(m_attributeSizeManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(attributeVectorSizeChanged(QtProperty *)));
  connect(m_tieManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(tieChanged(QtProperty *)));
  connect(m_constraintManager, SIGNAL(propertyChanged(QtProperty *)), this,
          SLOT(constraintChanged(QtProperty *)));
  connect(m_parameterManager, SIGNAL(valueChanged(QtProperty *, double)),
          SLOT(parameterChanged(QtProperty *)));
}

void QENSFunctionBrowser::createActions() {
  m_actionAddFunction = new QAction("Add function", this);
  connect(m_actionAddFunction, SIGNAL(triggered()), this, SLOT(addFunction()));

  m_actionRemoveFunction = new QAction("Remove function", this);
  connect(m_actionRemoveFunction, SIGNAL(triggered()), this,
          SLOT(removeFunction()));

  m_actionFixParameter = new QAction("Fix", this);
  connect(m_actionFixParameter, SIGNAL(triggered()), this,
          SLOT(fixParameter()));

  m_actionRemoveTie = new QAction("Remove tie", this);
  connect(m_actionRemoveTie, SIGNAL(triggered()), this, SLOT(removeTie()));

  m_actionAddTie = new QAction("Add tie", this);
  connect(m_actionAddTie, SIGNAL(triggered()), this, SLOT(addTie()));

  m_actionFromClipboard = new QAction("Copy from clipboard", this);
  connect(m_actionFromClipboard, SIGNAL(triggered()), this,
          SLOT(copyFromClipboard()));

  m_actionToClipboard = new QAction("Copy to clipboard", this);
  connect(m_actionToClipboard, SIGNAL(triggered()), this,
          SLOT(copyToClipboard()));

  m_actionConstraints = new QAction("Custom", this);
  connect(m_actionConstraints, SIGNAL(triggered()), this,
          SLOT(addConstraints()));

  m_actionConstraints10 = new QAction("10%", this);
  connect(m_actionConstraints10, SIGNAL(triggered()), this,
          SLOT(addConstraints10()));

  m_actionConstraints50 = new QAction("50%", this);
  connect(m_actionConstraints50, SIGNAL(triggered()), this,
          SLOT(addConstraints50()));

  m_actionRemoveConstraints = new QAction("Remove constraints", this);
  connect(m_actionRemoveConstraints, SIGNAL(triggered()), this,
          SLOT(removeConstraints()));

  m_actionRemoveConstraint = new QAction("Remove", this);
  connect(m_actionRemoveConstraint, SIGNAL(triggered()), this,
          SLOT(removeConstraint()));
}

void QENSFunctionBrowser::popupMenu(const QPoint &) {}

void QENSFunctionBrowser::addFunction() {
  auto prop = getSelectedProperty();
  if (isFunction(prop)) {
    const auto position = getFunctionPosition(prop);
    m_subscriber->addFunction(getFunctionFromUserDialog(), position);
  }
}

void QENSFunctionBrowser::removeFunction() {
  auto prop = getSelectedProperty();
  if (isFunction(prop)) {
    removeProperty(prop);
    updateFunctionIndices();
    m_subscriber->removeFunction(getFunctionPosition(prop));
  }
}

void QENSFunctionBrowser::fixParameter() {
  m_subscriber->fixParameter(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::removeTie() {
  const auto prop = getSelectedProperty();
  if (prop && isParameter(prop))
    m_subscriber->removeTie(getNameOfProperty(prop));
}

void QENSFunctionBrowser::addTie() {
  const auto prop = getSelectedProperty();
  if (prop && isParameter(prop)) {
    const auto name = getNameOfProperty(prop);
    if (auto tie = getTieFromDialog())
      m_subscriber->addTie(name, *tie);
  }
}

void QENSFunctionBrowser::addConstraints() {
  m_subscriber->addConstraints(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::removeConstraints() {
  m_subscriber->removeConstraints(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::addConstraints10() {
  m_subscriber->addConstraints10(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::addConstraints50() {
  m_subscriber->addConstraints50(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::removeConstraint() {
  m_subscriber->removeConstraint(getNameOfSelectedProperty());
}

void QENSFunctionBrowser::copyFromClipboard() {}

void QENSFunctionBrowser::copyToClipboard() {}

void QENSFunctionBrowser::updateCurrentFunctionIndex() {}

/**
 * Add a sub-property to a parent property
 * @param parent :: The parent property
 * @param subproperty :: New sub-property
 */
QENSFunctionBrowser::AProperty
QENSFunctionBrowser::addProperty(QtProperty *parent, QtProperty *subproperty) {
  AProperty ap;
  ap.prop = subproperty;
  if (parent == nullptr) {
    ap.item = m_browser->addProperty(subproperty);
  } else {
    parent->addSubProperty(subproperty);
    auto items = m_browser->items(subproperty);
    if (items.isEmpty()) {
      throw std::runtime_error("Unexpected error in FunctionBrowser [1]");
    }
    ap.item = items[0];
  }
  ap.parent = parent;
  m_properties[subproperty] = ap;
  return ap;
}

/**
 * Remove and delete property
 * @param prop :: Property to remove.
 */
void QENSFunctionBrowser::removeProperty(QtProperty *prop) {
  auto p = m_properties.find(prop);
  if (p == m_properties.end())
    return;
  AProperty ap = *p;

  // remove references to the children
  auto children = prop->subProperties();
  foreach (QtProperty *child, children) { removeProperty(child); }
  m_properties.erase(p);

  if (isFunction(prop)) {
    m_ties.remove(prop);
  }

  if (isTie(prop)) { //
    for (auto it = m_ties.begin(); it != m_ties.end(); ++it) {
      if (it.value().tieProp == prop) {
        m_ties.erase(it);
        break;
      }
    }
  }

  if (isConstraint(prop)) {
    for (auto it = m_constraints.begin(); it != m_constraints.end(); ++it) {
      auto &cp = it.value();
      if (cp.lower == prop) {
        if (!cp.upper) {
          m_constraints.erase(it);
        } else {
          cp.lower = nullptr;
        }
        break;
      } else if (cp.upper == prop) {
        if (!cp.lower) {
          m_constraints.erase(it);
        } else {
          cp.upper = nullptr;
        }
        break;
      }
    }
  }

  // remove property from Qt browser
  if (ap.parent) {
    ap.parent->removeSubProperty(prop);
  } else {
    m_browser->removeProperty(prop);
  }
  delete prop;
}

/**
 * Add a function property
 * @param parent :: Parent function property or NULL
 * @param funName :: Function name
 * @return :: A set AProperty struct
 */
QENSFunctionBrowser::AProperty
QENSFunctionBrowser::addFunctionProperty(QtProperty *parent,
                                         const QString &funName) {
  // check that parent is a function property
  if (parent && dynamic_cast<QtAbstractPropertyManager *>(m_functionManager) !=
                    parent->propertyManager()) {
    throw std::runtime_error("Unexpected error in FunctionBrowser [2]");
  }
  QtProperty *prop = m_functionManager->addProperty(funName);
  if (parent)
    m_propertyToIndex[prop] = parent->subProperties().size() - 1;
  return addProperty(parent, prop);
}

/**
 * Add a parameter property
 * @param parent :: Parent function property
 * @param paramName :: Parameter name
 * @param paramDesc :: Parameter description
 * @param paramValue :: Parameter value
 */
QENSFunctionBrowser::AProperty QENSFunctionBrowser::addParameterProperty(
    QtProperty *parent, const QString &paramName, const QString &paramDesc,
    double paramValue) {
  // check that parent is a function property
  if (!parent || dynamic_cast<QtAbstractPropertyManager *>(m_functionManager) !=
                     parent->propertyManager()) {
    throw std::runtime_error("Unexpected error in FunctionBrowser [3]");
  }
  QtProperty *prop = m_parameterManager->addProperty(paramName);
  m_parameterManager->blockSignals(true);
  m_parameterManager->setDecimals(prop, 6);
  m_parameterManager->setValue(prop, paramValue);
  m_parameterManager->setDescription(prop, paramDesc.toStdString());
  m_parameterManager->blockSignals(false);
  prop->setOption(globalOptionName, false);
  return addProperty(parent, prop);
}

/**
 * Add property showing function's index in the composite function
 * @param prop :: A function property
 * @return :: AProperty struct for added property. If all fields are NULL -
 * property wasn't added
 *  because it is the top function
 */
QENSFunctionBrowser::AProperty
QENSFunctionBrowser::addIndexProperty(QtProperty *prop) {
  AProperty ap;
  ap.item = nullptr;
  ap.parent = nullptr;
  ap.prop = nullptr;
  if (!prop)
    return ap;
  if (!isFunction(prop))
    return ap;
  if (!m_properties[prop].parent)
    return ap;

  QString index = "fff";
  QtProperty *ip = m_indexManager->addProperty("Index");
  ip->setEnabled(false);
  m_indexManager->setValue(ip, index);
  auto retval = addProperty(prop, ip);
  updateFunctionIndices();
  return retval;
}

/**
 * Add a constraint property
 * @param prop :: Parent parameter property
 * @param constraint :: A constraint string
 */
QList<QENSFunctionBrowser::AProperty>
QENSFunctionBrowser::addConstraintProperties(QtProperty *prop,
                                             QString constraint) {
  if (!isParameter(prop))
    return QList<AProperty>();
  QString lowerBoundStr = "";
  QString upperBoundStr = "";
  Mantid::API::Expression expr;
  expr.parse(constraint.toStdString());
  if (expr.name() != "==")
    return QList<AProperty>();
  if (expr.size() == 3) { // lower < param < upper
    try {
      // check that the first and third terms are numbers
      double d1 = boost::lexical_cast<double>(expr[0].name());
      (void)d1;
      double d2 = boost::lexical_cast<double>(expr[2].name());
      (void)d2;
      if (expr[1].operator_name() == "<" && expr[2].operator_name() == "<") {
        lowerBoundStr = QString::fromStdString(expr[0].name());
        upperBoundStr = QString::fromStdString(expr[2].name());
      } else // assume that the operators are ">"
      {
        lowerBoundStr = QString::fromStdString(expr[2].name());
        upperBoundStr = QString::fromStdString(expr[0].name());
      }
    } catch (...) { // error in constraint
      return QList<AProperty>();
    }
  } else if (expr.size() == 2) { // lower < param or param > lower etc
    size_t paramPos = 0;
    try // find position of the parameter name in expression
    {
      double d = boost::lexical_cast<double>(expr[1].name());
      (void)d;
    } catch (...) {
      paramPos = 1;
    }
    std::string op = expr[1].operator_name();
    if (paramPos == 0) { // parameter goes first
      if (op == "<") {   // param < number
        upperBoundStr = QString::fromStdString(expr[1].name());
      } else { // param > number
        lowerBoundStr = QString::fromStdString(expr[1].name());
      }
    } else {           // parameter is second
      if (op == "<") { // number < param
        lowerBoundStr = QString::fromStdString(expr[0].name());
      } else { // number > param
        upperBoundStr = QString::fromStdString(expr[0].name());
      }
    }
  }

  // add properties
  QList<AProperty> plist;
  AConstraint ac;
  ac.paramProp = prop;
  ac.lower = ac.upper = nullptr;
  if (!lowerBoundStr.isEmpty()) {
    auto ap = addProperty(prop, m_constraintManager->addProperty("LowerBound"));
    plist << ap;
    ac.lower = ap.prop;
    m_constraintManager->setValue(ac.lower, lowerBoundStr);
  }
  if (!upperBoundStr.isEmpty()) {
    auto ap = addProperty(prop, m_constraintManager->addProperty("UpperBound"));
    plist << ap;
    ac.upper = ap.prop;
    m_constraintManager->setValue(ac.upper, upperBoundStr);
  }
  return plist;
}

/**
 * Update function index properties
 * @param prop :: A function property
 * @param index :: The parent function's index
 */
void QENSFunctionBrowser::updateFunctionIndices(QtProperty *prop,
                                                QString index) {
  prop = prop ? getFirstProperty() : prop;
  if (!prop)
    return;

  auto children = prop->subProperties();
  std::size_t i = 0;
  foreach (QtProperty *child, children) { updateFunctionIndex(prop, index, i); }
}

void QENSFunctionBrowser::updateFunctionIndex(QtProperty *prop, QString prefix,
                                              std::size_t &i) {
  if (isFunction(prop)) {
    m_propertyToIndex[prop] = i;
    updateFunctionIndices(prop, prefix + "f" + QString::number(i) + ".");
    ++i;
  } else if (isIndex(prop)) {
    m_indexManager->setValue(prop, prefix);
  }
}

boost::optional<std::string> QENSFunctionBrowser::getTieFromDialog() {
  bool ok;
  QString tie = QInputDialog::getText(this, "Add a tie",
                                      "Tie:", QLineEdit::Normal, "", &ok);
  return ok ? boost::optional<std::string>(tie.toStdString()) : boost::none;
}

std::string QENSFunctionBrowser::getFunctionFromUserDialog() const {
  SelectFunctionDialog dlg(m_browser);
  if (dlg.exec() == QDialog::Accepted)
    return dlg.getFunction().toStdString();
  return "";
}

/**
 * Check if property is a function group
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isFunction(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(m_functionManager) ==
                     prop->propertyManager();
}

/**
 * Check if a property is a tie
 * @param prop :: A property
 */
bool QENSFunctionBrowser::isTie(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(m_tieManager) ==
                     prop->propertyManager();
}

/**
 * Check if a property is a constraint
 * @param prop :: Property to check.
 */
bool QENSFunctionBrowser::isConstraint(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_constraintManager) == prop->propertyManager();
}

/**
 * Check if property is any of the string attributes
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isStringAttribute(QtProperty *prop) const {
  return prop &&
         (dynamic_cast<QtAbstractPropertyManager *>(m_attributeStringManager) ==
              prop->propertyManager() ||
          dynamic_cast<QtAbstractPropertyManager *>(m_formulaManager) ==
              prop->propertyManager() ||
          dynamic_cast<QtAbstractPropertyManager *>(m_filenameManager) ==
              prop->propertyManager() ||
          dynamic_cast<QtAbstractPropertyManager *>(m_workspaceManager) ==
              prop->propertyManager());
}

/**
 * Check if property is a function attribute
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isDoubleAttribute(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_attributeDoubleManager) == prop->propertyManager();
}

/**
 * Check if property is a function attribute
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isIntAttribute(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_attributeIntManager) == prop->propertyManager();
}

/**
 * Check if property is a function bool attribute
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isBoolAttribute(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_attributeBoolManager) == prop->propertyManager();
}

/**
 * Check if property is a function vector attribute
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isVectorAttribute(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_attributeVectorManager) == prop->propertyManager();
}

/**
 * Check if property is a function attribute
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isAttribute(QtProperty *prop) const {
  return isStringAttribute(prop) || isDoubleAttribute(prop) ||
         isIntAttribute(prop) || isBoolAttribute(prop) ||
         isVectorAttribute(prop);
}

/**
 * Check if property is a function parameter
 * @param prop :: Property to check
 */
bool QENSFunctionBrowser::isParameter(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(
                     m_parameterManager) == prop->propertyManager();
}

/**
 * Check if a property is an index
 * @param prop :: A property
 */
bool QENSFunctionBrowser::isIndex(QtProperty *prop) const {
  return prop && dynamic_cast<QtAbstractPropertyManager *>(m_indexManager) ==
                     prop->propertyManager();
}

std::string QENSFunctionBrowser::getNameOfProperty(QtProperty *prop) const {
  const auto name = getIndex(prop) + prop->propertyName();
  return name.toStdString();
}

QtProperty *QENSFunctionBrowser::getFirstProperty() const {
  auto top = m_browser->properties();
  if (!top.isEmpty())
    return top.front();
  return nullptr;
}

QtProperty *QENSFunctionBrowser::getSelectedProperty() const {
  auto item = m_browser->currentItem();
  if (item)
    return item->property();
  return getFirstProperty();
}

QtProperty *
QENSFunctionBrowser::getContainingFunctionProperty(QtProperty *prop) const {
  auto aPropIt = m_properties.find(prop);
  while (aPropIt != m_properties.end()) {
    auto parent = aPropIt->parent;
    if (isFunction(parent))
      return parent;
    aPropIt = m_properties.find(parent);
  }
  return nullptr;
}

std::vector<std::size_t>
QENSFunctionBrowser::getFunctionPosition(QtProperty *prop) const {
  if (!prop)
    return std::vector<std::size_t>();

  std::vector<std::size_t> position;
  QMap<QtProperty *, std::size_t>::const_iterator indexIt;
  while (prop &&
         (indexIt = m_propertyToIndex.find(prop)) != m_propertyToIndex.end()) {
    position.emplace_back(*indexIt);
    prop = getContainingFunctionProperty(prop);
  }
  return position;
}

/**
 * Get the function index for a property
 * @param prop :: A property
 */
QString QENSFunctionBrowser::getIndex(QtProperty *prop) const {
  if (!prop)
    return "";
  if (isFunction(prop)) {
    auto props = prop->subProperties();
    if (props.isEmpty())
      return "";
    for (auto it = props.begin(); it != props.end(); ++it) {
      if (isIndex(*it)) {
        return m_indexManager->value(*it);
      }
    }
    return "";
  }
  return getIndex(m_properties[prop].parent);
}

void QENSFunctionBrowser::stringAttributeChanged(QtProperty *prop,
                                                 const QString &value) {
  const auto name = prop->propertyName().toStdString();
  const auto position = getFunctionPosition(prop);
  m_subscriber->setStringAttribute(name, value.toStdString(), position);
}

void QENSFunctionBrowser::intAttributeChanged(QtProperty *prop, int value) {
  const auto name = prop->propertyName().toStdString();
  const auto position = getFunctionPosition(prop);
  m_subscriber->setIntAttribute(name, value, position);
}

void QENSFunctionBrowser::doubleAttributeChanged(QtProperty *prop,
                                                 double value) {
  const auto name = prop->propertyName().toStdString();
  const auto position = getFunctionPosition(prop);
  m_subscriber->setDoubleAttribute(name, value, position);
}

void QENSFunctionBrowser::boolAttributeChanged(QtProperty *prop, bool value) {
  const auto name = prop->propertyName().toStdString();
  const auto position = getFunctionPosition(prop);
  m_subscriber->setBoolAttribute(name, value, position);
}

void QENSFunctionBrowser::vectorDoubleAttributeChanged(QtProperty *prop) {
  const auto name = prop->propertyName().toStdString();
  const auto position = getFunctionPosition(prop);
  const auto value =
      getVectorFromProperty<double>(prop, m_attributeVectorDoubleManager);
  m_subscriber->setVectorDoubleAttribute(name, value, position);
}

void QENSFunctionBrowser::parameterButtonClicked(QtProperty *prop) {
  emit localParameterButtonClicked(getIndex(prop) + prop->propertyName());
}

} // namespace IDA
} // namespace CustomInterfaces
} // namespace MantidQt
