#ifndef MANTIDQTMANTIDWIDGETS_BATCHWIDGETTREE_H_
#define MANTIDQTMANTIDWIDGETS_BATCHWIDGETTREE_H_
#include "TreeNode.h"
#include <vector>
namespace MantidQt {
namespace MantidWidgets {
namespace Batch {

class NodeLocation {
public:
  NodeLocation(std::vector<std::size_t> path);
  std::vector<std::size_t> const& path() const;
private:
  std::vector<std::size_t> m_path;
};

template <typename T> class Tree {
public:
  Tree(Tree const &other);

  TreeNode<T>* at(NodeLocation location);
  TreeNode<T> const* at(NodeLocation location) const;

  NodeLocation appendChild(NodeLocation const& parent, T value);
  void appendManyChildren(NodeLocation const& parent, std::vector<T> value);

  TreeNode<T>& addBeforeSibling(NodeLocation const& sibling, T value);
  void addManyBeforeSibling(NodeLocation const& sibling, std::vector<T> value);

  void removeAt(NodeLocation location);
  void removeMany(std::vector<NodeLocation> const& locations);

  Tree& operator=(Tree const& other);
private:
  TreeNode<T> m_root;
};

}
}
}
#include "Tree.tcc"
#endif // MANTIDQTMANTIDWIDGETS_BATCHWIDGETTREE_H_
