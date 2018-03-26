#include "MantidQtWidgets/Common/Batch/TreeNode.h"
#include <cxxtest/TestSuite.h>
using namespace MantidQt::MantidWidgets;
using namespace MantidQt::MantidWidgets::Batch;
class TreeNodeTest : public CxxTest::TestSuite {
public:
  void testHasCorrectValueWhenDefaultConstructed() {
    auto constexpr VALUE = 0;
    auto node = TreeNode<int>(VALUE);
    TS_ASSERT_EQUALS(VALUE, node.value());
  }

  void testHasNoParentWhenDefaultConstructed() {
    auto node = TreeNode<int>(0);
    TS_ASSERT(!node.hasParent());
  }

  void testIsNotLeafWhenConstructedWithChildren() {
    auto firstChild = TreeNode<int>(1);
    auto secondChild = TreeNode<int>(2);
    auto rootNode = TreeNode<int>(0, {firstChild, secondChild});
    TS_ASSERT(!isLeaf(rootNode));
  }

  void testUpdatesParentsOnCopy() {
    auto parent = TreeNode<int>(0);
    auto child = TreeNode<int>(1);
    parent.push_back(std::move(child));

    TS_ASSERT_EQUALS(0, parent[0].parent().value());

    auto parentCopy = parent;
    parentCopy.setValue(2);
    parentCopy.updateParentPointers();
    TS_ASSERT_EQUALS(2, parentCopy[0].parent().value());
  }

  void testSizeIncreasesWhenNodeAppended() {
    auto rootNode = TreeNode<int>(0);
    auto firstChild = TreeNode<int>(2);
    rootNode.push_back(std::move(firstChild));
    TS_ASSERT_EQUALS(1, rootNode.size());
  }

  void testParentNodeSetWhenNodeAppended() {
    auto rootNode = TreeNode<int>(0);
    auto firstChild = TreeNode<int>(2);
    rootNode.push_back(std::move(firstChild));
    TS_ASSERT_EQUALS(0, (*rootNode.begin()).parent().value());
  }

  void testCanIterateOverChildrenWhenConst() {
    auto firstChild = TreeNode<int>(1);
    auto secondChild = TreeNode<int>(2);
    auto rootNode = TreeNode<int>(0, {firstChild, secondChild});

    TS_ASSERT_EQUALS(2, rootNode.size());

    auto it = rootNode.cbegin();
    TS_ASSERT_EQUALS(1, (*it).value());
    ++it;
    TS_ASSERT_EQUALS(2, (*it).value());
    ++it;
    TS_ASSERT_EQUALS(rootNode.cend(), it);
  }

  void testCanAppendChildNode() {
    auto rootNode = TreeNode<int>(0);
    auto firstChild = TreeNode<int>(1);
    rootNode.push_back(std::move(firstChild));
    auto it = rootNode.cbegin();
    TS_ASSERT_EQUALS(1, (*it).value());
  }

  void testCanInsertNodeInBetweenOtherChildren() {
    auto rootNode = TreeNode<int>(0, {TreeNode<int>(1), TreeNode<int>(3)});
    rootNode.insert(rootNode.cend() - 1, TreeNode<int>(2));

    TS_ASSERT_EQUALS(3, rootNode.size());
    TS_ASSERT_EQUALS(2, rootNode[1].value());
  }

  void testCanEraseNode() {
    auto rootNode = TreeNode<int>(0, {TreeNode<int>(1), TreeNode<int>(3)});
    rootNode.erase(rootNode.cbegin());

    TS_ASSERT_EQUALS(1, rootNode.size());
    TS_ASSERT_EQUALS(3, rootNode[0].value());
  }

  void testCanEraseRangeOfNodes() {
    auto rootNode = TreeNode<int>(0, {TreeNode<int>(1), TreeNode<int>(2),
                                      TreeNode<int>(3), TreeNode<int>(4)});
    rootNode.erase(rootNode.cbegin() + 1, rootNode.cend() - 1);

    TS_ASSERT_EQUALS(2, rootNode.size());
    TS_ASSERT_EQUALS(1, rootNode[0].value());
    TS_ASSERT_EQUALS(4, rootNode[1].value());
  }

  void testParentNodeSetWhenFirstNodeInserted() {
    auto rootNode = TreeNode<int>(0);
    rootNode.insert(rootNode.begin(), TreeNode<int>(2));

    TS_ASSERT_EQUALS(0, rootNode[0].parent().value());
  }

  void testParentNodeSetWhenNodeInserted() {
    auto rootNode = TreeNode<int>(0, {TreeNode<int>(1), TreeNode<int>(3)});
    rootNode.insert(rootNode.cend() - 1, TreeNode<int>(2));

    TS_ASSERT_EQUALS(0, rootNode[1].parent().value());
  }

  void testSubscriptOperatorReturnsNthChild() {
    auto rootNode = TreeNode<int>(
        0, {TreeNode<int>(1), TreeNode<int>(2), TreeNode<int>(3)});

    TS_ASSERT_EQUALS(1, rootNode[0].value());
    TS_ASSERT_EQUALS(2, rootNode[1].value());
    TS_ASSERT_EQUALS(3, rootNode[2].value());
  }

  void testDoublyNestedTreeSubscript() {
    auto rootNode = TreeNode<int>(0, {TreeNode<int>(1),
                                      TreeNode<int>(2, {TreeNode<int>(4)}),
                                      TreeNode<int>(3)});

    TS_ASSERT_EQUALS(1, rootNode[0].value());
    TS_ASSERT_EQUALS(2, rootNode[1].value());
    TS_ASSERT_EQUALS(4, rootNode[1][0].value());
    TS_ASSERT_EQUALS(3, rootNode[2].value());
  }

  void testParentPointersOnThreeLevelDeepTree() {
    // clang-format off
    auto tree = TreeNode<int>(0,
      {
        TreeNode<int>(1,
          {
            TreeNode<int>(2)
          }),
        TreeNode<int>(3)
      });
    // clang-format on

    tree.updateParentPointers();
    TS_ASSERT_EQUALS(2, tree[0][0].value());
    TS_ASSERT_EQUALS(1, tree[0][0].parent().value());
    TS_ASSERT_EQUALS(0, tree[0][0].parent().parent().value());
  }
};
