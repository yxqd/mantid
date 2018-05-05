#include "MantidAlgorithms/Warp.h"

#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/MuParserUtils.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidDataObjects/WorkspaceCreation.h"

#include <boost/geometry.hpp>



namespace {
using namespace Mantid;
using Mantid::Algorithms::Warp;

Warp::Ring makeQuadrilateral(const HistogramData::HistogramX &Xs, const std::vector<double> &twoThetaEdges, const size_t binIndex, const size_t wsIndex) {
  const auto y1 = twoThetaEdges[wsIndex];
  const auto y2 = twoThetaEdges[wsIndex + 1];
  const auto top = std::max(y1, y2);
  const auto bottom = std::min(y1, y2);
  const auto x1 = Xs[binIndex];
  const auto x2 = Xs[binIndex + 1];
  const auto left = std::min(x1, x2);
  const auto right = std::max(x1, x2);
  Warp::Ring q;
  q.reserve(4);
  q.emplace_back(left, top);
  q.emplace_back(right, top);
  q.emplace_back(right, bottom);
  q.emplace_back(left, bottom);
  return q;
}

class BinCache {
public:
  BinCache(const std::vector<double> &twoThetaEdges, const API::MatrixWorkspace &ws);

  double area(const size_t wsIndex, const size_t binIndex) const;
  const Warp::Ring & ring(const size_t wsIndex, const size_t binIndex) const;
private:
  std::vector<std::vector<double>> m_binAreas{};
  std::vector<std::vector<Warp::Ring>> m_binRings{};
};

BinCache::BinCache(const std::vector<double> &twoThetaEdges, const API::MatrixWorkspace &ws) {
  const auto nHist = ws.getNumberHistograms();
  m_binAreas.reserve(nHist);
  m_binRings.reserve(nHist);
  for (size_t i = 0; i < nHist; ++i) {
    const auto &Xs = ws.x(i);
    const auto nBins = Xs.size() - 1;
    m_binAreas.emplace_back();
    auto &areas = m_binAreas.back();
    areas.reserve(nBins);
    m_binRings.emplace_back();
    auto &rings = m_binRings.back();
    rings.reserve(nBins);
    for (size_t j = 0; j < nBins; ++j) {
      rings.emplace_back(makeQuadrilateral(Xs, twoThetaEdges, j, i));
      const auto &ring = rings.back();
      areas.emplace_back((ring[1].x() - ring[0].x()) * (ring[1].y() - ring[2].y()));
    }
  }
}

double BinCache::area(const size_t wsIndex, const size_t binIndex) const {
  return m_binAreas[wsIndex][binIndex];
}

const Warp::Ring & BinCache::ring(const size_t wsIndex, const size_t binIndex) const {
  return m_binRings[wsIndex][binIndex];
}

void fillNaNs(API::MatrixWorkspace &ws) {
	for (size_t i = 0; i < ws.getNumberHistograms(); ++i) {
		ws.mutableY(i) = std::nan("");
	}
}

std::vector<double> makeTwoThetaEdges(const API::SpectrumInfo &spectrumInfo) {
	const auto n = spectrumInfo.size();
	std::vector<double> edges(n + 1);
  auto t1 = spectrumInfo.signedTwoTheta(0);
  auto t2 = spectrumInfo.signedTwoTheta(1);
	edges.front() = t1 - (t2 - t1) / 2.;
	for (size_t i = 1; i < n - 1; ++i) {
		edges[i] = (t2 + t1) / 2.;
		t1 = t2;
    t2 = spectrumInfo.signedTwoTheta(i + 1);
	}
	edges[n - 1] = (t2 + t1) / 2.;
	edges.back() = t2 + (t2 - t1) / 2.;
	return edges;
}

void redistributeCounts(const double counts, const Warp::Ring &warpedBin, const std::vector<double> &twoThetaEdges, API::MatrixWorkspace &ws, const BinCache &cache) {
  const auto warpedArea = boost::geometry::area(warpedBin);
  Warp::Box boundingBox;
  boost::geometry::envelope(warpedBin, boundingBox);
  const auto &minCorner = boundingBox.min_corner();
  auto twoThetaIter = std::lower_bound(twoThetaEdges.cbegin(), twoThetaEdges.cend(), minCorner.y());
  if (twoThetaIter == twoThetaEdges.cend()) {
    return;
  }
  if (twoThetaIter != twoThetaEdges.cbegin()) {
    --twoThetaIter;
  }
  auto i = static_cast<size_t>(std::distance(twoThetaEdges.cbegin(), twoThetaIter));
  ++twoThetaIter;
  const auto &maxCorner = boundingBox.max_corner();
  const auto twoThetaEnd = std::lower_bound(twoThetaIter, twoThetaEdges.cend(), maxCorner.y());
  const auto nHist = ws.getNumberHistograms();
  const auto endHist = std::min(nHist, static_cast<decltype(nHist)>(std::distance(twoThetaEdges.cbegin(), twoThetaEnd)));
  for (; i < endHist; ++i) {
    const auto &Xs = ws.x(i);
		auto &Ys = ws.mutableY(i);
    auto xIter = std::lower_bound(Xs.cbegin(), Xs.cend(), minCorner.x());
    if (xIter == Xs.cend()) {
      continue;
    }
    if (xIter != Xs.cbegin()) {
      --xIter;
    }
    auto j = static_cast<size_t>(std::distance(Xs.cbegin(), xIter));
    ++xIter;
    const auto xEnd = std::lower_bound(xIter, Xs.cend(), maxCorner.x());
    const auto nBins = Ys.size();
    const auto endBin = std::min(nBins, static_cast<decltype(nBins)>(std::distance(Xs.cbegin(), xEnd)));
    for (; j < endBin; ++j) {
      const auto &bin = cache.ring(i, j);
      Warp::PolygonVector intersections;
      if (boost::geometry::intersection(warpedBin, bin, intersections)) {
        double intersectionArea{0.};
        for (const auto &intersection : intersections) {
          intersectionArea += boost::geometry::area(intersection);
        }
        const auto f = intersectionArea / warpedArea * counts;
        if (std::isnan(Ys[j])) {
          Ys[j] = f;
        } else {
          Ys[j] += f;
        }
      }
    }
	}
}

}

namespace Mantid {
namespace Algorithms {
// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(Warp)

Warp::Warp() {
	m_tExpressionParser.DefineVar("t", &m_t);
  m_tExpressionParser.DefineVar("x", &m_x);
  m_xExpressionParser.DefineVar("t", &m_t);
	m_xExpressionParser.DefineVar("x", &m_x);
	API::MuParserUtils::addDefaultConstants(m_tExpressionParser);
	API::MuParserUtils::addDefaultConstants(m_xExpressionParser);
	API::MuParserUtils::extraOneVarFunctions(m_tExpressionParser);
	API::MuParserUtils::extraOneVarFunctions(m_xExpressionParser);
}

/// Algorithms name for identification. @see Algorithm::name
const std::string Warp::name() const {
  return "Warp";
}

/// Algorithm's version for identification. @see Algorithm::version
int Warp::version() const { return 1; }

/// Algorithm's category for identification. @see Algorithm::category
const std::string Warp::category() const {
  return "Transforms";
}

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string Warp::summary() const {
  return "Transforms workspace bins according to a given expression.";
}

/** Initialize the algorithm's properties.
 */
void Warp::init() {
	declareProperty(
		Kernel::make_unique<API::WorkspaceProperty<API::MatrixWorkspace>>(
			"InputWorkspace", "", Kernel::Direction::Input),
		"A workspace to warp.");
	declareProperty(
		Kernel::make_unique<API::WorkspaceProperty<API::MatrixWorkspace>>(
			"OutputWorkspace", "", Kernel::Direction::Output),
		"The warped workspace.");
	declareProperty("TExpression", "", "A muParser expression t' = T(x,t) where x is the horizontal and t is the vertical axis coordinate.");
	declareProperty("XExpression", "", "A muParser expression x' = X(x,t) where x is the horizontal and t is the vertical axis coordinate.");
}

/** Execute the algorithm.
 */
void Warp::exec() {
	m_tExpressionParser.SetExpr(getProperty("TExpression"));
	m_xExpressionParser.SetExpr(getProperty("XExpression"));
	API::MatrixWorkspace_sptr inWS = getProperty("InputWorkspace");
	API::MatrixWorkspace_sptr outWS = DataObjects::create<DataObjects::Workspace2D>(*inWS);
	fillNaNs(*outWS);
	const auto &spectrumInfo = inWS->spectrumInfo();
  const auto nHist = spectrumInfo.size();
	const auto twoThetaEdges = makeTwoThetaEdges(spectrumInfo);
  BinCache cache(twoThetaEdges, *inWS);
  API::Progress progress(this, 0., 1., nHist);
  for (size_t i = 0; i < nHist; ++i) {
    progress.report("Processing histogram " + std::to_string(i));
    const auto &Xs = inWS->x(i);
		const auto &Ys = inWS->y(i);
		for (size_t j = 0; j < Xs.size() - 1; ++j) {
      interruption_point();
      const auto &bin = cache.ring(i, j);
			const auto warped = warp(bin);
      redistributeCounts(Ys[j], warped, twoThetaEdges, *outWS, cache);
    }
  }
	setProperty("OutputWorkspace", outWS);
}

double Warp::evaluate(mu::Parser &parser, const Point &p) {
	m_x = boost::geometry::get<0>(p);
	m_t = boost::geometry::get<1>(p);
	try {
		return parser.Eval();
	}
	catch (mu::Parser::exception_type &) {
		throw std::runtime_error("Error in evaluating expression '" + parser.GetExpr() + "'");
	}
}

Warp::Ring Warp::warp(const Ring &q) {
  Ring warped;
	for (size_t i = 0; i < q.size(); ++i) {
		const auto warpedT = evaluate(m_tExpressionParser, q[i]);
		const auto warpedX = evaluate(m_xExpressionParser, q[i]);
    warped.emplace_back(warpedX, warpedT);
	}
	return warped;
}

} // namespace Algorithms
} // namespace Mantid
