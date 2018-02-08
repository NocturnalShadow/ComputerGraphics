#include "Algorithms.h"

#include <map>
#include <algorithm>

namespace Geometry
{
	using namespace std;
	using namespace Eigen;

	bool isInside(const vector<Vector2f>& polygon, Vector2f point)
	{
		Vector2f pivot = (polygon[0] + polygon[1] + polygon[2]) / 3;

		// TODO: replace with vector
		map<float, Vector2f> vertices;
		for (auto vertex : polygon) {
			vertices[angleOf(vertex - pivot)] = vertex;
		}

		Vector2f first, second;
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
}