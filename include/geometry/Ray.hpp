#ifndef G3_RAY_HPP
#define G3_RAY_HPP

#include <glm/glm.hpp>
#include <nanort.h>

namespace gamma3::geometry {

using Point = glm::dvec3;
using Vector = glm::dvec3;

class Ray {
public:

    Ray(const Point& source, const Vector& direction) {
        nanort_ray_.org[0] = source.x;
        nanort_ray_.org[1] = source.y;
        nanort_ray_.org[2] = source.z;

        Vector normalized_dir = glm::normalize(direction);
        nanort_ray_.dir[0] = normalized_dir.x;
        nanort_ray_.dir[1] = normalized_dir.y;
        nanort_ray_.dir[2] = normalized_dir.z;

        nanort_ray_.min_t = 0.0;
        nanort_ray_.max_t = std::numeric_limits<double>::max();
    }

    Ray(const Point& source, const Vector& direction, double min_t, double max_t) {
        nanort_ray_.org[0] = source.x;
        nanort_ray_.org[1] = source.y;
        nanort_ray_.org[2] = source.z;

        Vector normalized_dir = glm::normalize(direction);
        nanort_ray_.dir[0] = normalized_dir.x;
        nanort_ray_.dir[1] = normalized_dir.y;
        nanort_ray_.dir[2] = normalized_dir.z;

        nanort_ray_.min_t = min_t;
        nanort_ray_.max_t = max_t;
    }

    explicit Ray(const nanort::Ray<double>& nanort_ray) : nanort_ray_(nanort_ray) {}

    Point pointAt(double t) const {
        return Point(
            nanort_ray_.org[0] + t * nanort_ray_.dir[0],
            nanort_ray_.org[1] + t * nanort_ray_.dir[1],
            nanort_ray_.org[2] + t * nanort_ray_.dir[2]
            );
    }

    Point source() const { return Point(nanort_ray_.org[0], nanort_ray_.org[1], nanort_ray_.org[2]); }

    Vector direction() const { return Vector(nanort_ray_.dir[0], nanort_ray_.dir[1], nanort_ray_.dir[2]); }

    void moveSource(double t) {
        nanort_ray_.org[0] += t * nanort_ray_.dir[0];
        nanort_ray_.org[1] += t * nanort_ray_.dir[1];
        nanort_ray_.org[2] += t * nanort_ray_.dir[2];
    }

    void setSource(const Point& source) {
        nanort_ray_.org[0] = source.x;
        nanort_ray_.org[1] = source.y;
        nanort_ray_.org[2] = source.z;
    }

    void setDirection(const Vector& direction) {
        Vector normalized_dir = glm::normalize(direction);
        nanort_ray_.dir[0] = normalized_dir.x;
        nanort_ray_.dir[1] = normalized_dir.y;
        nanort_ray_.dir[2] = normalized_dir.z;
    }

    double minT() const { return nanort_ray_.min_t; }

    double maxT() const { return nanort_ray_.max_t; }

    void setMinT(double min_t) { nanort_ray_.min_t = min_t; }

    void setMaxT(double max_t) { nanort_ray_.max_t = max_t; }

    bool isInRange(double t) const { return t >= minT() && t <= maxT(); }

    const nanort::Ray<double>& getNanortRay() const { return nanort_ray_; }

    nanort::Ray<double>& getNanortRay() { return nanort_ray_; }

    operator const nanort::Ray<double>() const { return nanort_ray_; }

    operator nanort::Ray<double>() { return nanort_ray_; }

    Ray reversed() const { return Ray(source(), -direction(), minT(), maxT()); }

private:

    nanort::Ray<double> nanort_ray_;

};

}

#endif //G3_RAY_HPP
