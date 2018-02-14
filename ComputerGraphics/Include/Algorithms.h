#include <vector>

#include <Eigen/Dense>

namespace Geometry
{
	inline float angleOf(Eigen::Vector2f vector) {
		return atan2(vector.y(), vector.x());
	}

	bool isInside(const std::vector<Eigen::Vector2f>& polygon, Eigen::Vector2f point);
}