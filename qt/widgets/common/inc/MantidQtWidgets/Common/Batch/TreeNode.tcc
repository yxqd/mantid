#include <cassert>
#include <algorithm>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {
template <typename T>
TreeNode<T>::TreeNode(T const &value)
    : m_value(value), m_parent(nullptr), m_children() {}

template <typename T>
auto TreeNode<T>::erase(const_iterator position) -> iterator {
  return m_children.erase(position);
}

template <typename T>
auto TreeNode<T>::erase(const_iterator begin, const_iterator end) -> iterator {
  return m_children.erase(begin, end);
}

template <typename T>
void TreeNode<T>::updateParentPointers() {
  for (auto &child : m_children) {
    child.setParent(this);
    child.updateParentPointers();
  }
}

template <typename T>
TreeNode<T>::TreeNode(T const &value, std::vector<TreeNode<T>> children)
    : m_value(value), m_parent(nullptr), m_children(children) {
  for (auto &child : m_children)
    child.setParent(this);
}

template <typename T> bool TreeNode<T>::hasParent() const {
  return m_parent != nullptr;
}

template <typename T> TreeNode<T> &TreeNode<T>::parent() const {
  assert(hasParent());
  return *m_parent;
}

template <typename T> TreeNode<T> *TreeNode<T>::parentOrNullptr() const {
  return m_parent;
}

template <typename T> void TreeNode<T>::setParent(TreeNode<T> *parent) {
  m_parent = parent;
}

template <typename T> T const &TreeNode<T>::value() const { return m_value; }
template <typename T> T &TreeNode<T>::value() { return m_value; }

template <typename T> void TreeNode<T>::setValue(T const &value) {
  m_value = value;
}

template <typename T>
TreeNode<T> &TreeNode<T>::push_back(TreeNode<T> newChild) {
  newChild.setParent(this);
  m_children.emplace_back(std::move(newChild));
  return m_children.back();
}

template <typename T>
TreeNode<T> &TreeNode<T>::insert(const_iterator position,
                                 TreeNode<T> newChild) {
  newChild.setParent(this);
  auto it = m_children.insert(position, std::move(newChild));
  return *it;
}

template <typename T> std::size_t TreeNode<T>::size() const {
  return m_children.size();
}

template <typename T> auto TreeNode<T>::begin() -> iterator {
  return m_children.begin();
}

template <typename T> auto TreeNode<T>::end() -> iterator {
  return m_children.end();
}

template <typename T> auto TreeNode<T>::begin() const -> const_iterator {
  return cbegin();
}

template <typename T> auto TreeNode<T>::end() const -> const_iterator {
  return cend();
}

template <typename T> auto TreeNode<T>::cbegin() const -> const_iterator {
  return m_children.cbegin();
}

template <typename T> auto TreeNode<T>::cend() const -> const_iterator {
  return m_children.cend();
}

template <typename T> bool isLeaf(TreeNode<T> const &node) {
  return node.size() == 0;
}

template <typename T> bool isRoot(TreeNode<T> const &node) {
  return !node.hasParent();
}

template <typename T>
int indexOfNodeWithinParent(const TreeNode<T> &nodeWithParent) {
  auto &parent = nodeWithParent.parent();
  auto iteratorAtNode =
      std::find(parent.cbegin(), parent.cend(), nodeWithParent);
  auto index = std::distance(parent.cbegin(), iteratorAtNode);
  return static_cast<int>(index);
}

template <typename T> TreeNode<T> &TreeNode<T>::operator[](std::size_t index) {
  return m_children[index];
}

template <typename T>
TreeNode<T> const &TreeNode<T>::operator[](std::size_t index) const {
  return m_children[index];
}

template <typename T>
bool operator==(TreeNode<T> const &lhs, TreeNode<T> const &rhs) {
  return &lhs == &rhs;
}

template <typename T>
bool operator!=(TreeNode<T> const &lhs, TreeNode<T> const &rhs) {
  return !(lhs == rhs);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, TreeNode<T> const& node) {
  os << "n\n";
  for (const auto& c : node) {
    os << " c\n";
  }
  return os;
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
