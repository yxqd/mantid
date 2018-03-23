#ifndef MANTID_MANTIDWIDGETS_BATCHTREENODETEST_H
#define MANTID_MANTIDWIDGETS_BATCHTREENODETEST_H
#include <cxxtest/TestSuite.h>
#include "MantidQtWidgets/Common/Batch/TreeNode.h"
using namespace MantidQt::MantidWidgets;
using namespace MantidQt::MantidWidgets::DataProcessor;
class TreeNodeTest : public CxxTest::TestSuite {
public:
  void testHasNoParentWhenDefaultConstructed() {
    auto node = TreeNode<int>(0);
    TS_ASSERT(!node.hasParent());
  }

  void testIsNotLeafWhenConstructedWithChildren() {
    auto firstChild = TreeNode<int>(1);
    auto secondChild = TreeNode<int>(2);
    auto rootNode = TreeNode<int>({ firstChild, secondChild });
    TS_ASSERT(!isLeaf(rootNode));
  }

};
#endif // MANTID_MANTIDWIDGETS_BATCHTREENODETEST_H
