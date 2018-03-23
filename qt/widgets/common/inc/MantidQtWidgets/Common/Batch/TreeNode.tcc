namespace MantidQt {
namespace MantidWidgets {
namespace Batch {
template <typename T>
TreeNode<T>::TreeNode(T const &value)
    : m_value(value), m_parent(nullptr), m_children() {}

template <typename T>
TreeNode<T>::TreeNode(T const &value, std::vector<TreeNode<T>> children)
    : m_value(value), m_parent(nullptr), m_children(children) {
  for (auto &child : m_children)
    child->setParent(this);
}

template <typename T> bool TreeNode<T>::hasParent() const {
  return m_parent != nullptr;
}

template <typename T> TreeNode<T> &TreeNode<T>::parent() const {
  return *m_parent;
}

template <typename T> void TreeNode<T>::setParent(TreeNode<T> *parent) {
  m_parent = parent;
}

template <typename T> void TreeNode<T>::setValue(T const &value) {
  m_value = value;
}

template <typename T> TreeNode<T> &TreeNode<T>::push_back(T const &childValue) {
  auto newNode = TreeNode<T>(childValue);
  newNode.setParent(this);
  m_children.emplace_back(newNode);
}

template <typename T>
TreeNode<T> &TreeNode<T>::insert(const_iterator pos, T const &value) {
  auto newNode = TreeNode<T>(childValue);
  newNode.setParent(this);
}

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
