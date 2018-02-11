#include <Catch>

#include "Algorithms.h"

TEST_CASE("isInside() function tells if the point is inside of the given convex polygon", "[2d][location][point-in-polygon][convex]") 
{
	using namespace Eigen;
	using namespace Geometry;

	std::vector<Vector2f> polygon { 
		{ 0.f, 0.f }, 
		{ 1.0f, 0.0f },
		{ 1.5f, 1.5f },
		{ 0.5f, 1.5f },
		{ -0.5f, 0.5f }
	};

	SECTION("Should return true if the point is inside")
	{
		Vector2f point { 0.5f, 0.5f };
		REQUIRE(isInside(polygon, point) == true);
	}

	SECTION("Should return false if the point is outside")
	{
		Vector2f point { 2.f, 2.f };
		REQUIRE(isInside(polygon, point) == false);
	}

	SECTION("Should return true if the point is on the boundary")
	{

		Vector2f point1{ 0.5f, 0.0f };
		REQUIRE(isInside(polygon, point1) == true);

		Vector2f point3{ 1.0f, 1.5f };
		REQUIRE(isInside(polygon, point3) == true);
	}
}