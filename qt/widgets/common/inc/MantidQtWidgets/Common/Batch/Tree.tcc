namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

template <typename T>
Tree<T>::Tree(Tree const &other)
  : m_root(other.m_root) {
  m_root.updateParentPointers();
}

template <typename T>
TreeNode<T>* Tree<T>::at(NodeLocation location) {
  TreeNode<T>* currentNode = &m_root;
  for (auto const& pathComponent : location.path()) {
    currentNode = &(*currentNode[pathComponent]);
  }
  return currentNode;
}

TreeNode<T> const* Tree<T>::at(NodeLocation location) const {
  TreeNode<T>* currentNode = &m_root;
  for (auto const& pathComponent : location.path()) {
    currentNode = &(*currentNode[pathComponent]);
  }
  return currentNode;
}

TreeNode<T>& addBefore(NodeLocation const& location, T value) {
  auto* node = at(location);
  node->push_back()
}
  void addAllAt(NodeLocation const& location, std::vector<T> value);

  void removeAt(NodeLocation location);
  void removeAll(std::vector<NodeLocation> const& locations);

  Tree& operator=(Tree const& other);
private:
  TreeNode<T> m_root;
};

}
}
}
