#include "Algorithms.h"

#include <map>
#include <iostream>
#include <algorithm>

namespace Geometry
{
	using namespace Eigen;

	bool isInside(const std::vector<Point>& polygon, Point point)
	{
		Point pivot = (polygon[0] + polygon[1] + polygon[2]) / 3;

		std::map<float, Point> vertices;
		for (auto vertex : polygon) {
			vertices[angleOf(vertex - pivot)] = vertex;
		}

		Point first, second;
		auto upperBound = vertices.upper_bound(angleOf(point - pivot));
		if (upperBound == vertices.begin() || upperBound == vertices.end()) {
			first = vertices.rbegin()->second;
			second = vertices.begin()->second;
		}
		else {
			first = prev(upperBound)->second;
			second = upperBound->second;
		}

		Matrix3f matrix;
		matrix <<
			first.transpose(), 1,
			second.transpose(), 1,
			point.transpose(), 1;

		auto determinant =
			(matrix <<
				first.transpose(), 1,
				second.transpose(), 1,
				point.transpose(), 1)
			.finished().determinant();

		return !signbit(determinant);
	}

	std::unique_ptr<GraphHelper> regular(const Graph& _graph)
	{
		auto graph = std::make_unique<GraphHelper>(_graph);
		auto& coords = _graph.coords;

		auto intersection = [&] (float y, Edge e)
		{
			auto& begin = coords[e.begin];
			auto& end = coords[e.end];
			auto yDelta = end.y() - begin.y();
			if (abs(yDelta) < 1e-8f) {
				constexpr auto inf = std::numeric_limits<float>::infinity();
				return begin.x() < end.x() ? inf : -inf;
			}
			return (end.x() - begin.x()) * (y - begin.y()) / yDelta + begin.x();
		};

		float sweepLine;
		auto comparator = [&] (Edge left, Edge right)
		{
			bool commonBegin = left.begin == right.begin,
				commonEnd = left.end == right.end;
			if (commonBegin || commonEnd)
			{
				auto leftVector = coords[left.end] - coords[left.begin],
					rightVector = coords[right.end] - coords[right.begin];
				auto crossProduct = cross(leftVector, rightVector);

				return commonBegin ?
					(crossProduct < 0 ? true : false) :
					(crossProduct > 0 ? true : false);
			}

			return intersection(sweepLine, left) < intersection(sweepLine, right);
		};

		auto vertexCount = coords.size();
		std::set<Edge, decltype(comparator)> edges{ comparator };
		std::vector<std::vector<decltype(edges)::iterator>> associatedEdges(vertexCount);
		std::vector<int> associatedLeft(vertexCount, 0);
		std::vector<int> associatedRight(vertexCount, 0);

		for (int vertex = 0; vertex < coords.size() - 1; vertex++)
		{
			sweepLine = coords[vertex].y();
			for (auto edge : associatedEdges[vertex]) {
				edges.erase(edge);
			}

			auto edgeIter = edges.lower_bound(Edge{ 0, vertex }),
				leftBound = edgeIter != edges.begin() ? std::prev(edgeIter) : edges.end(),
				rightBound = edgeIter;
			for (auto edge : graph->outgoingEdges[vertex])
			{
				edgeIter = edges.insert(edgeIter, graph->edges[edge]);
				associatedEdges[graph->edges[edge].end].push_back(edgeIter);
				rightBound = std::next(edgeIter);
			}

			int left = leftBound != edges.end() ? leftBound->begin : 0,
				right = rightBound != edges.end() ? rightBound->begin : 0;
			if (graph->ingoingEdges[vertex].empty() && vertex)
			{
				int candidate = std::max({ left, right, associatedRight[left], associatedLeft[right] });
				graph->insertFakeEdge(candidate, vertex);
			}
			else if (graph->outgoingEdges[vertex].empty())
			{
				associatedRight[left] = std::max(associatedRight[left], vertex);
				associatedLeft[right] = std::max(associatedLeft[right], vertex);
			}
		}

		for (auto& vertexEdges : associatedEdges) {
			vertexEdges.clear();
		}
		edges.clear();
		std::fill(associatedLeft.begin(), associatedLeft.end(), vertexCount - 1);
		std::fill(associatedRight.begin(), associatedRight.end(), vertexCount - 1);

		int lastVertex = coords.size() - 1;
		for (int vertex = lastVertex; vertex > 0; vertex--)
		{
			sweepLine = coords[vertex].y();
			for (auto edge : associatedEdges[vertex]) {
				edges.erase(edge);
			}

			auto edgeIter = edges.lower_bound(Edge{ 0, vertex }),
				leftBound = edgeIter != edges.begin() ? std::prev(edgeIter) : edges.end(),
				rightBound = edgeIter;
			for (auto edge : graph->ingoingEdges[vertex])
			{
				edgeIter = edges.insert(edgeIter, graph->edges[edge]);
				associatedEdges[graph->edges[edge].begin].push_back(edgeIter);
				rightBound = std::next(edgeIter);
			}

			int left = leftBound != edges.end() ? leftBound->end : vertexCount - 1,
				right = rightBound != edges.end() ? rightBound->end : vertexCount - 1;
			if (graph->outgoingEdges[vertex].empty() && vertex != lastVertex)
			{
				int candidate = std::min({ left, right, associatedRight[left], associatedLeft[right] });
				graph->insertFakeEdge(vertex, candidate);
			}
			else if (graph->ingoingEdges.empty())
			{
				associatedRight[left] = std::min(associatedRight[left], vertex);
				associatedLeft[right] = std::min(associatedLeft[right], vertex);
			}
		}

		return graph;
	}

	OrderedChains monotonousChainsSet(const GraphHelper& graph)
	{
		int vertexCount = graph.coords.size();
		int edgeCount = graph.edges.size();

		auto& outgoingEdges = graph.outgoingEdges;
		auto& ingoingEdges = graph.ingoingEdges;

		std::vector<int>
			weights(edgeCount, 1),
			ingoingWeights(vertexCount),
			outgoingWeights(vertexCount);

		for (int vertex = 1; vertex < vertexCount - 1; vertex++)
		{
			for (int edge : ingoingEdges[vertex]) {
				ingoingWeights[vertex] += weights[edge];
			}

			int outgoingEdgesCount = outgoingEdges[vertex].size();
			int leftmostOutgoingEdge = *(outgoingEdges[vertex].begin());
			if (ingoingWeights[vertex] > outgoingEdgesCount) {
				weights[leftmostOutgoingEdge] = ingoingWeights[vertex] - outgoingEdgesCount + 1;
			}
		}

		for (int vertex = vertexCount - 2; vertex > 0; vertex--)
		{
			for (int edge : outgoingEdges[vertex]) {
				outgoingWeights[vertex] += weights[edge];
			}

			int leftmostIngoingEdge = *(ingoingEdges[vertex].begin());
			if (outgoingWeights[vertex] > ingoingWeights[vertex]) {
				weights[leftmostIngoingEdge] = outgoingWeights[vertex] - ingoingWeights[vertex] + weights[leftmostIngoingEdge];
			}
		}

		int chainsCount = 0;
		for (int edge : outgoingEdges[0]) {
			chainsCount += weights[edge];
		}

		OrderedChains chains(chainsCount);
		std::vector<int> leftmostChains(vertexCount, chainsCount - 1); leftmostChains[0] = 0;
		for (int vertex = 0; vertex < vertexCount; vertex++)
		{
			int leftmostChain = leftmostChains[vertex];
			for (int edge : outgoingEdges[vertex])
			{
				int endpoint = graph.edges[edge].end;
				if (leftmostChain < leftmostChains[endpoint]) {
					leftmostChains[endpoint] = leftmostChain;
				}

				for (int weight = 0; weight < weights[edge]; weight++, leftmostChain++) {
					chains[leftmostChain].push_back(edge);
				}
			}
		}

		return chains;
	}

	Facet locatePoint(const Graph& graph, Point point)
	{
		auto& coords = graph.coords;
		if (coords.empty() || point.y() <= coords.front().y() || point.y() >= coords.back().y()) 
			return std::vector<int> { };

		auto regularGraph = regular(graph);
		auto chains = monotonousChainsSet(*regularGraph);
		auto rightBound = std::upper_bound(chains.begin(), chains.end(), point, [&] (auto point, const auto& chain)
		{
			int edgeIndex = *std::upper_bound(chain.begin(), chain.end(), point, [&] (auto point, int edge) {
				return point.y() < coords[regularGraph->edges[edge].end].y();
			});

			auto edge = regularGraph->edges[edgeIndex];
			auto left = coords[edge.end] - coords[edge.begin],
				right = point - coords[edge.begin];

			return cross(left, right) > 0 ? true : false;
		});

		if (rightBound == chains.end() || rightBound == chains.begin()) 
			return std::vector<int> { };

		int edge = *std::upper_bound(rightBound->begin(), rightBound->end(), point, [&](auto point, int edge) {
			return point.y() < coords[regularGraph->edges[edge].end].y();
		});

		return regularGraph->facetLeftToEdge(edge);
	}
}