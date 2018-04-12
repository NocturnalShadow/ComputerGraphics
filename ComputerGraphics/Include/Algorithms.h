#include <set>
#include <vector>
#include <memory>
#include <functional>

#include <Eigen/Dense>

namespace Geometry
{
	typedef Eigen::Vector2f Point;
	typedef std::vector<int> Chain;				// sequence of edge indices
	typedef std::vector<Chain> OrderedChains;	
	typedef std::vector<int> Facet;				// sequence of vertices

	struct Graph
	{
		std::vector<Point> coords;
		std::vector<std::vector<int>> endpoints;
	};

	struct Edge { int begin, end; };			// pair of vertex indices

	inline float cross(Eigen::Vector2f first, Eigen::Vector2f second) {
		return first.x() * second.y() - first.y() * second.x();
	}

	struct GraphHelper
	{
		std::vector<Point> coords;
		std::vector<Edge> edges;
		std::vector<std::set<int, std::function<bool(int, int)>>> ingoingEdges, outgoingEdges;

		int realEdgesCount;

		float cross(int left, int right)
		{
			auto leftVector = coords[edges[left].end] - coords[edges[left].begin],
				rightVector = coords[edges[right].end] - coords[edges[right].begin];
			return Geometry::cross(leftVector, rightVector);
		}

		GraphHelper(const Graph& graph)
			: coords{ graph.coords }
		{
			auto ingoingComparator = [this] (int leftEdge, int rightEdge) {
				return cross(leftEdge, rightEdge) > 0 ? true : false;
			};

			auto outgoingComparator = [this] (int leftEdge, int rightEdge) {
				return cross(leftEdge, rightEdge) < 0 ? true : false;
			};

			for (int vertex = 0; vertex < graph.coords.size(); vertex++)
			{
				for (int endpoint : graph.endpoints[vertex]) {
					edges.push_back({ vertex, endpoint });
				}

				ingoingEdges.emplace_back(ingoingComparator);
				outgoingEdges.emplace_back(outgoingComparator);
			}

			for (int e = 0; e < edges.size(); e++) 
			{
				outgoingEdges[edges[e].begin].insert(e);
				ingoingEdges[edges[e].end].insert(e);
			}

			realEdgesCount = edges.size();
		}

		void insertFakeEdge(int begin, int end)
		{
			int edge = edges.size();
			edges.push_back({ begin, end });
			outgoingEdges[begin].insert(edge);
			ingoingEdges[end].insert(edge);
		}

		bool isFake(int edge) {
			return edge >= realEdgesCount;
		}

		Facet facetLeftToEdge(int edge)
		{
			Facet facet;
			
			int initEdge = edge;
			int vertex;
			while (true)
			{
				vertex = edges[edge].end;
				facet.push_back(vertex);
				if (*ingoingEdges[vertex].begin() != edge) break;
				edge = *outgoingEdges[vertex].begin();
			}

			edge = *std::prev(ingoingEdges[vertex].find(edge));

			while (true)
			{
				vertex = edges[edge].begin;
				facet.push_back(vertex);
				if (*outgoingEdges[vertex].rbegin() != edge) break;
				edge = *ingoingEdges[vertex].rbegin();
			}

			edge = *std::next(outgoingEdges[vertex].find(edge));

			while (edge != initEdge)
			{
				vertex = edges[edge].end;
				facet.push_back(vertex);
				edge = *outgoingEdges[vertex].begin();
			}

			return facet;
		}
	};

	inline float angleOf(Eigen::Vector2f vector) {
		return atan2(vector.y(), vector.x());
	}

	bool isInside(const std::vector<Point>& polygon, Point point);

	OrderedChains monotonousChainsSet(const GraphHelper& graph);
	std::unique_ptr<GraphHelper> regular(const Graph& graph);
	Facet locatePoint(const Graph& graph, Point point);
}