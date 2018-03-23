#ifndef MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
#define MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
#include <vector>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

template<typename T>
class TreeNode {
  using ChildrenStorage = std::vector<TreeNode<T>>;
  using iterator = typename ChildrenStorage::iterator;
  using const_iterator = typename ChildrenStorage::const_iterator;

  TreeNode();
  explicit TreeNode(std::vector<TreeNode<T>> children);

  bool hasParent() const;
  TreeNode<T>& parent();
  TreeNode<T> const& parent() const;

  void setParent(TreeNode<T>& parent);
  void setValue(T const& value);

  TreeNode<T>& push_back(T const&);
  TreeNode<T>& insert(const_iterator pos, T const& value);
  iterator erase(const_iterator pos);
  iterator erase(const_iterator begin, const_iterator end);

  iterator begin();
  const_iterator cbegin() const;
  const_iterator begin() const;

  iterator end();
  const_iterator cend() const;
  const_iterator end() const;

  std::size_t size() const;

  T value;
  TreeNode<T>* m_parent;
  ChildrenStorage m_children;
};

template <typename T>
bool isLeaf(TreeNode<T>& node);

}
}
}
#include "TreeNode.tcc"
#endif // MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
