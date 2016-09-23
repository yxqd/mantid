#ifndef MANTID_CUSTOMINTERFACES_DATACOMPARISONPRESENTERTEST_H_
#define MANTID_CUSTOMINTERFACES_DATACOMPARISONPRESENTERTEST_H_

#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>

#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"

#include "MantidQtCustomInterfaces/DataComparison/IDataComparisonView.h"
#include "MantidQtCustomInterfaces/DataComparison/DataComparisonPresenter.h"
#include "qwt_data.h"
#include <vector>

using namespace Mantid::API;
using namespace MantidQt::CustomInterfaces;
using namespace testing;

/** Mock class, we'll use it to test the presenter
*
*/
class MockDataComparisonView : public IDataComparisonView {

public:
  MOCK_CONST_METHOD0(getSelectedWorkspaceName, std::string());
  MOCK_CONST_METHOD0(getWorkspaceNames, std::vector<std::string>());
  MOCK_CONST_METHOD0(getSelectedWorkspaceNames, std::vector<std::string>());
  MOCK_CONST_METHOD0(getWorkspaceColors, std::vector<std::string>());
  MOCK_CONST_METHOD0(getAvailableColors, std::vector<std::string>());

  MOCK_METHOD2(addWorkspace, void(const std::string &, int));
  MOCK_METHOD1(removeWorkspace, void(const std::string &));
  MOCK_CONST_METHOD1(containsWorkspace, bool(const std::string &));

  MOCK_METHOD1(blockTableSignals, void(bool));

  MOCK_METHOD1(detachCurve, void(const std::string &));
  MOCK_METHOD3(plotCurve, void(const std::string &, const QwtArrayData &,
                               const std::string &));

  MOCK_CONST_METHOD0(getSelectedWorkspaceIndex, int());

  MOCK_METHOD1(printError, void(const std::string &));
  MOCK_METHOD1(printInformation, void(const std::string &));
  MOCK_METHOD1(printDebug, void(const std::string &));

  // We don't necessarily have to mock all the virtual methods
  // If there's a method we are not interested in we can just
  // implement it here

  // We can also define new methods if needed
};

MATCHER_P3(QwtDataX, i, value, delta, "") {
  return fabs(arg.x(i) - value) < delta;
}
MATCHER_P3(QwtDataY, i, value, delta, "") {
  return fabs(arg.y(i) - value) < delta;
}

class DataComparisonPresenterTest : public CxxTest::TestSuite {

public:
  // The mock view
  NiceMock<MockDataComparisonView> m_view;
  // NiceMock<>: Avoids GMOCK Warnings (Uninteresting mock function call...)
  // when function returning default (see setUp() below)

  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static DataComparisonPresenterTest *createSuite() {
    return new DataComparisonPresenterTest();
  }
  static void destroySuite(DataComparisonPresenterTest *suite) { delete suite; }

  DataComparisonPresenterTest() { FrameworkManager::Instance(); }

  void setUp() override {

    /// We can set here some default return values for the view mock object
    /// getters

    std::vector<std::string> colors = {"Red", "Green", "Blue", "Black", "Pink"};
    ON_CALL(m_view, getAvailableColors()).WillByDefault(Return(colors));
  }

  void tearDown() override {}

  void test_addWorkspace_workspace_does_not_exist() {
    DataComparisonPresenter presenter(&m_view);

    // Test what happens when we try to add a workspace that does not exist in
    // the ADS

    EXPECT_CALL(m_view, getSelectedWorkspaceName())
        .Times(1)
        .WillOnce(Return("TestWorkspace"));
    EXPECT_CALL(m_view, printError(_));

    presenter.notify(IDataComparisonPresenter::AddWorkspace);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));

    // We can even be more specific about the error message:
    // EXPECT_CALL(m_view, printError("Cannot add workspace TestWorkspace"));
  }

  void test_addWorkspace_wrong_workspace_type() {
    DataComparisonPresenter presenter(&m_view);

    // We shouldn't be able to load a table ws into the interface

    ITableWorkspace_sptr ws = WorkspaceFactory::Instance().createTable();

    EXPECT_CALL(m_view, printError(_));

    presenter.addWorkspace(ws);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
  }

  void test_addWorkspace_existing_workspace() {
    DataComparisonPresenter presenter(&m_view);

    // We shouldn't be able to load a ws that already exists in the UI

    MatrixWorkspace_sptr ws =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 1, 1);

    EXPECT_CALL(m_view, containsWorkspace(_)).Times(1).WillOnce(Return(true));
    EXPECT_CALL(m_view,
                printInformation("Workspace already shown in comparison"));

    presenter.addWorkspace(ws);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
  }

  void test_addWorkspace_index_out_of_range() {
    DataComparisonPresenter presenter(&m_view);

    // We shouldn't be able to load a ws that already exists in the UI

    MatrixWorkspace_sptr ws =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 1, 1);

    EXPECT_CALL(m_view, containsWorkspace(_)).Times(1).WillOnce(Return(false));
    EXPECT_CALL(m_view, getSelectedWorkspaceIndex())
        .Times(1)
        .WillOnce(Return(100));
    EXPECT_CALL(m_view, printError(_));
    EXPECT_CALL(m_view, getAvailableColors()).Times(0);
    EXPECT_CALL(m_view, getWorkspaceColors()).Times(0);

    presenter.addWorkspace(ws);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
  }

  void test_addWorkspace_matrix_workspace_OK() {
    DataComparisonPresenter presenter(&m_view);

    // A valid ws should be added to the table and to the plot widget

    MatrixWorkspace_sptr ws =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 2, 2);
    ws->mutableX(0)[0] = 0.0;
    ws->mutableX(0)[1] = 1.0;
    ws->mutableY(0)[0] = 2.0;
    ws->mutableY(0)[1] = 3.0;

    EXPECT_CALL(m_view, containsWorkspace(_)).Times(1).WillOnce(Return(false));
    EXPECT_CALL(m_view, getWorkspaceColors())
        .Times(1)
        .WillOnce(Return(std::vector<std::string>()));
    EXPECT_CALL(m_view, addWorkspace(_, 0)).Times(1);
    EXPECT_CALL(m_view,
                plotCurve(_,
                          AllOf(Property(&QwtData::size, 2),
                                QwtDataX(0, 0, 1E-12), QwtDataX(1, 1, 1E-12),
                                QwtDataY(0, 2, 1E-12), QwtDataY(1, 3, 1E-12)),
                          "Black")).Times(1);

    presenter.addWorkspace(ws);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
  }

  void test_addWorkspace_workspace_group_OK() {

    // We could test the same but with workspace groups:
    // We could test that all ws in the group are added
    // We could test values and colors passed to the view for plotting
  }

  void test_plotWorkspaces_workspace_not_in_ADS() {

    // We could test that nothing is plotted if workpsace is not in the ADS
  }

  void test_plotWorkspaces_index_out_of_range() {

    // We could test that nothing is plotted if workpsace index is out of range
  }

  void test_plotDiffWorkspace_more_than_two_ws_selected() {

    // Test that if there are more than two ws selected in the table the view
    // prints an error message
  }

  void test_plotDiffWorkspace_workspaces_not_in_ADS() {

    // Test that if selected workspaces are not in the ADS the view prints an
    // error message
  }

  void test_plotDiffWorkspace_OK() {
    DataComparisonPresenter presenter(&m_view);

    MatrixWorkspace_sptr ws1 =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 3, 2);
    ws1->mutableX(0)[0] = 0.0;
    ws1->mutableX(0)[1] = 1.0;
    ws1->mutableX(0)[2] = 2.0;
    ws1->mutableY(0)[0] = 2.0;
    ws1->mutableY(0)[1] = 3.0;
    AnalysisDataService::Instance().add("ws1", ws1);
    MatrixWorkspace_sptr ws2 =
        WorkspaceFactory::Instance().create("Workspace2D", 1, 3, 2);
    ws2->mutableX(0)[0] = 0.0;
    ws2->mutableX(0)[1] = 1.0;
    ws2->mutableX(0)[2] = 2.0;
    ws2->mutableY(0)[0] = 1.0;
    ws2->mutableY(0)[1] = 2.0;
    AnalysisDataService::Instance().add("ws2", ws2);

    EXPECT_CALL(m_view, getSelectedWorkspaceNames())
        .Times(1)
        .WillOnce(Return(std::vector<std::string>{"ws1", "ws2"}));
    EXPECT_CALL(m_view, getSelectedWorkspaceIndex())
        .Times(1)
        .WillOnce(Return(0));
    EXPECT_CALL(m_view, detachCurve("ws1ws2"));
    EXPECT_CALL(m_view,
                plotCurve("ws1ws2",
                          AllOf(QwtDataX(0, 0, 1E-12), QwtDataX(1, 1, 1E-12),
                                QwtDataY(0, 1, 1E-12), QwtDataY(1, 1, 1E-12)),
                          "")).Times(1);

    presenter.notify(IDataComparisonPresenter::PlotDiffWorkspace);

    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));

    AnalysisDataService::Instance().clear();
  }

  // More tests...

};

#endif /* MANTID_CUSTOMINTERFACES_DATACOMPARISONPRESENTERTEST_H_ */
