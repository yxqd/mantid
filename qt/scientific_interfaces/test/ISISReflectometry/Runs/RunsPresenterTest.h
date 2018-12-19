// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H
#define MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H

#include "../../../ISISReflectometry/GUI/Runs/RunsPresenter.h"
#include "../../ReflMockObjects.h"
#include "../RunsTable/MockRunsTablePresenterFactory.h"
#include "../RunsTable/MockRunsTableView.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/make_unique.h"
#include "MantidQtWidgets/Common/Batch/MockJobTreeView.h"
#include "MantidQtWidgets/Common/MockProgressableView.h"
#include "MockRunsView.h"
#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace MantidQt::CustomInterfaces;
using namespace testing;

//=====================================================================================
// Functional tests
//=====================================================================================
class RunsPresenterTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static RunsPresenterTest *createSuite() { return new RunsPresenterTest(); }
  static void destroySuite(RunsPresenterTest *suite) { delete suite; }

  RunsPresenterTest()
      : m_thetaTolerance(0.01), m_instruments{"INTER", "SURF", "CRISP",
                                              "POLREF", "OFFSPEC"},
        m_view(), m_runsTableView(), m_progressView(), m_messageHandler(),
        m_autoreduction(new MockReflAutoreduction),
        m_searcher(new MockReflSearcher) {
    ON_CALL(m_view, table()).WillByDefault(Return(&m_runsTableView));
    ON_CALL(m_runsTableView, jobs()).WillByDefault(ReturnRef(m_jobs));
  }

  void testMakePresenter() { auto presenter = makePresenter(); }

  void testSettingsChanged() {
    // TODO
  }

  void testSearchWithEmptyString() {
    auto presenter = makePresenter();
    auto searchString = std::string("");
    EXPECT_CALL(m_view, getSearchString())
        .Times(1)
        .WillOnce(Return(searchString));
    expectSearchFailed();
    presenter.notify(IRunsPresenter::SearchFlag);
    verifyAndClear();
  }

  void testSearchCatalogLoginFails() {
    auto presenter = makePresenter();
    auto searchString = std::string("test string");
    EXPECT_CALL(m_view, getSearchString())
        .Times(1)
        .WillOnce(Return(searchString));
    // TODO: add expected call to python runner when implemented
    EXPECT_CALL(m_view, noActiveICatSessions()).Times(1);
    expectSearchFailed();
    presenter.notify(IRunsPresenter::SearchFlag);
    verifyAndClear();
  }

  void testSearchSucceeds() {
    // TODO: add this test when python runner implemented
  }

private:
  RunsPresenter makePresenter() {
    auto defaultInstrumentIndex = 0;
    EXPECT_CALL(m_view, subscribe(_)).Times(1);
    EXPECT_CALL(m_runsTableView, subscribe(_)).Times(1);
    EXPECT_CALL(m_view, table()).Times(1).WillOnce(Return(&m_runsTableView));
    EXPECT_CALL(m_runsTableView, jobs()).Times(1).WillOnce(ReturnRef(m_jobs));
    EXPECT_CALL(m_view,
                setInstrumentList(m_instruments, defaultInstrumentIndex))
        .Times(1);
    expectUpdateViewWhenMonitorStopped();

    auto presenter = RunsPresenter(
        &m_view, &m_progressView,
        MockRunsTablePresenterFactory(m_instruments, m_thetaTolerance),
        m_thetaTolerance, m_instruments, defaultInstrumentIndex,
        &m_messageHandler, m_autoreduction, m_searcher);

    verifyAndClear();
    return presenter;
  }

  void verifyAndClear() {
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_view));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_runsTableView));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_progressView));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_messageHandler));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_autoreduction));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_searcher));
    TS_ASSERT(Mock::VerifyAndClearExpectations(&m_jobs));
  }

  void expectUpdateViewWhenMonitorStarting() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(false));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(false));
  }

  void expectUpdateViewWhenMonitorStarted() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(false));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(true));
  }

  void expectUpdateViewWhenMonitorStopped() {
    EXPECT_CALL(m_view, setStartMonitorButtonEnabled(true));
    EXPECT_CALL(m_view, setStopMonitorButtonEnabled(false));
  }

  void expectStopAutoreduction() {
    EXPECT_CALL(m_view, stopTimer()).Times(1);
    EXPECT_CALL(*m_autoreduction, stop()).Times(1);
  }

  void expectSearchFailed() {
    EXPECT_CALL(m_view, getAlgorithmRunner()).Times(0);
    expectStopAutoreduction();
  }

  double m_thetaTolerance;
  std::vector<std::string> m_instruments;
  MockRunsView m_view;
  MockRunsTableView m_runsTableView;
  MockProgressableView m_progressView;
  MockMessageHandler m_messageHandler;
  boost::shared_ptr<MockReflAutoreduction> m_autoreduction;
  boost::shared_ptr<MockReflSearcher> m_searcher;
  NiceMock<MantidQt::MantidWidgets::Batch::MockJobTreeView> m_jobs;
};

#endif /* MANTID_CUSTOMINTERFACES_RUNSPRESENTERTEST_H */
