#include <Catch>

#include "Algorithms.h"

TEST_CASE(
	"Function isInside(polygon, point) tells if the point is inside of the given convex polygon", 
	"[2d][location][point-in-polygon][convex]") 
{
	using namespace Eigen;
	using namespace Geometry;

	std::vector<Point> polygon { 
		{ 0.f, 0.f }, 
		{ 1.0f, 0.0f },
		{ 1.5f, 1.5f },
		{ 0.5f, 1.5f },
		{ -0.5f, 0.5f }
	};

	SECTION("Should return true if the point is inside")
	{
		REQUIRE(isInside(polygon, Point { 0.5f, 0.5f }) == true);
	}

	SECTION("Should return false if the point is outside")
	{
		REQUIRE(isInside(polygon, Point { 2.f, 2.f }) == false);
	}

	SECTION("Should return true if the point is on the boundary")
	{
		REQUIRE(isInside(polygon, Point { 0.5f, 0.0f }) == true);
		REQUIRE(isInside(polygon, Point { 1.0f, 1.5f }) == true);
	}
}

TEST_CASE(
	"Function regular(graph) produces a regularized graph.",
	"[2d][graph][regular]")
{
	SECTION("Simple graph")
	{
		Geometry::Graph graph {
			{
				{ 3, 1 },{ 9, 2 },{ 2, 3 },
				{ 4, 4 },{ 8, 5 },{ 6, 6 },
				{ 7, 7 },{ 1, 8 },{ 5, 9 }
			},
			{
				{ 2, 4, 1 },
				{ 4 },
				{ 7, 5 },
				{ 5, 4 },
				{ 5 },
				{ 8 },
				{ 8 },
				{ },
				{ },
			}
		};

		CHECK(graph.endpoints[0].size() == 3);
		CHECK(graph.endpoints[2].size() == 2);
		CHECK(graph.endpoints[4].size() == 1);
		CHECK(graph.endpoints[5].size() == 1);
		CHECK(graph.endpoints[7].size() == 0);

		auto regularGraph = Geometry::regular(graph);
		REQUIRE(regularGraph->outgoingEdges[0].size() == 3);
		REQUIRE(regularGraph->outgoingEdges[2].size() == 3);
		REQUIRE(regularGraph->outgoingEdges[4].size() == 1);
		REQUIRE(regularGraph->outgoingEdges[5].size() == 2);
		REQUIRE(regularGraph->outgoingEdges[7].size() == 1);
	}

	SECTION("Complex graph")
	{
		Geometry::Graph graph {
			{
				{ 165, 10 },{ 329, 17 },{ 235, 143 },
				{ 160, 175 },{ 322, 185 },{ 252, 240 },
				{ 9, 327 },{ 462, 451 }
			},
			{
				{ 6, 2 },
				{ 2, 7 },
				{ },
				{ 6, 5 },
				{ 5, 7 },
				{ },
				{ },
				{ },
			}
		};

		CHECK(graph.endpoints[0].size() == 2);
		CHECK(graph.endpoints[2].size() == 0);
		CHECK(graph.endpoints[3].size() == 2);
		CHECK(graph.endpoints[5].size() == 0);
		
		auto regularGraph = Geometry::regular(graph);
		REQUIRE(regularGraph->outgoingEdges[0].size() == 3);
		REQUIRE(regularGraph->outgoingEdges[2].size() == 1);
		REQUIRE(regularGraph->outgoingEdges[3].size() == 3);
		REQUIRE(regularGraph->outgoingEdges[5].size() == 1);
		REQUIRE(regularGraph->ingoingEdges[3].size() == 1);
		REQUIRE(regularGraph->ingoingEdges[4].size() == 1);
	}
}

TEST_CASE(
	"Function monotonousChainsSet(graph) returns a full set of monotonous chains of the given graph",
	"[2d][location]")
{
	Geometry::Graph graph {
		{
			{ 3, 1 },{ 9, 2 },{ 2, 3 },
			{ 4, 4 },{ 8, 5 },{ 6, 6 },
			{ 7, 7 },{ 1, 8 },{ 5, 9 }
		},
		{
			{ 2, 3, 4, 1 },
			{ 4 },
			{ 7, 5, 3 },
			{ 5, 4 },
			{ 5, 6 },
			{ 8, 6 },
			{ 8 },
			{ 8 },
			{ },
		}
	};

	auto chains = Geometry::monotonousChainsSet(graph);
	REQUIRE(chains.size() == 6);
	REQUIRE(chains[0].size() == 3);
	REQUIRE(chains[1].size() == 3);
	REQUIRE(chains[2].size() == 4);
	REQUIRE(chains[3].size() == 4);
	REQUIRE(chains[4].size() == 4);
	REQUIRE(chains[5].size() == 4);
}

TEST_CASE(
	"Function localizePoint(graph, point) returns a facet of regularized graph where the point is located.",
	"[2d][location]")
{
	Geometry::Graph graph {
		{
			{ 165, 10 },	{ 329, 17 },	{ 235, 143 },
			{ 160, 175 },	{ 322, 185 },	{ 252, 240 },
			{ 9, 327 },		{ 462, 451 }
		},
		{
			{ 6, 2 },
			{ 2, 7 },
			{ },
			{ 6, 5 },
			{ 5, 7 },
			{ },
			{ },
			{ },
		}
	};

	SECTION("Point is outside of the regularized graph.")
	{
		auto actual = Geometry::locatePoint(graph, { 200, 0 });
		REQUIRE(actual.empty());
	}
	
	SECTION("Point is outside of the initial graph, but inside of a regularized area.")
	{
		auto actual = Geometry::locatePoint(graph, { 200, 40 });
		auto expected = std::vector<int> { 2, 0, 1 };
		REQUIRE(expected == actual);
	}

	SECTION("Point is inside of the initial graph.")
	{
		auto actual = Geometry::locatePoint(graph, { 155, 115 });
		auto expected = std::vector<int> { 2, 3, 6, 0 };
		REQUIRE(expected == actual);
	}
}