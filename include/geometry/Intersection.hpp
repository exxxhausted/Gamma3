#ifndef G3_INTERSECTION_HPP
#define G3_INTERSECTION_HPP

#include <limits>

#include <glm/glm.hpp>

namespace gamma3::geometry {

using Point = glm::dvec3;
using Vector = glm::dvec3;

struct Intersection {

    double t = std::numeric_limits<double>::max();

    Point point;

    Vector normal;

    uint32_t face_id = 0;

};

}

#endif // G3_INTERSECTION_HPP
