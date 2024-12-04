#include "interpolation.hpp"

glm::vec3
interpolation::evalLERP(glm::vec3 const& p0, glm::vec3 const& p1, float const x)
{
	//! \todo Implement this function

	return (1 - x) * p0 + x * p1;
}

glm::vec3
interpolation::evalCatmullRom(glm::vec3 const& p0, glm::vec3 const& p1,
                              glm::vec3 const& p2, glm::vec3 const& p3,
                              float const t, float const x)
{
	//! \todo Implement this function

	glm::vec4 xs = glm::vec4(1, x, x * x, x * x * x);

	glm::mat4 M = glm::mat4(
		glm::vec4(0, -t, 2 * t, -t),
		glm::vec4(1, 0, t - 3, 2 - t),
		glm::vec4(0, t, 3 - 2 * t, t - 2),
		glm::vec4(0, 0, -t, t)
	);

	glm::mat3x4 points = glm::transpose(glm::mat4x3(p0, p1, p2, p3));

	return xs * M * points;
}
