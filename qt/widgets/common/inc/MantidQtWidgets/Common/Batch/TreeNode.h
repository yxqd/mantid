#ifndef MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
#define MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
#include <ostream>
#include <vector>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

template <typename T> class TreeNode {
public:
  using ChildrenStorage = std::vector<TreeNode<T>>;
  using iterator = typename ChildrenStorage::iterator;
  using const_iterator = typename ChildrenStorage::const_iterator;

  explicit TreeNode(T const &value);
  explicit TreeNode(T const &value, std::vector<TreeNode<T>> children);

  bool hasParent() const;
  TreeNode<T> &parent() const;
  TreeNode<T> *parentOrNullptr() const;
  void updateParentPointers();

  void setParent(TreeNode<T> *parent);
  void setValue(T const &value);
  T const &value() const;
  T &value();

  TreeNode<T> &push_back(TreeNode<T> newChild);
  TreeNode<T> &insert(const_iterator pos, TreeNode<T> newChild);
  iterator erase(const_iterator pos);
  iterator erase(const_iterator begin, const_iterator end);

  iterator begin();
  const_iterator cbegin() const;
  const_iterator begin() const;

  iterator end();
  const_iterator cend() const;
  const_iterator end() const;

  std::size_t size() const;

  TreeNode<T> &operator[](std::size_t index);
  TreeNode<T> const &operator[](std::size_t index) const;

private:
  T m_value;
  TreeNode<T> *m_parent;
  ChildrenStorage m_children;
};

template <typename T>
bool operator==(TreeNode<T> const &lhs, TreeNode<T> const &rhs);

template <typename T>
bool operator!=(TreeNode<T> const &lhs, TreeNode<T> const &rhs);

template <typename T>
std::ostream &operator<<(std::ostream &os, TreeNode<T> const &node);

template <typename T> bool isLeaf(TreeNode<T> const &node);
template <typename T> bool isRoot(TreeNode<T> const &node);

template <typename T>
int indexOfNodeWithinParent(const TreeNode<T> &nodeWithParent);

} // namespace Batch
} // namespace MantidWidgets
} // namespace MantidQt
#include "TreeNode.tcc"
#endif // MANTIDQTMANTIDWIDGETS_BATCHWIDGETROW_H_
