#ifndef G3_INTERSECTOR_H
#define G3_INTERSECTOR_H

#include "Ray.hpp"
#include "Intersection.hpp"

#include <nanort.h>

namespace gamma3::geometry {

class Intersector {
public:

    Intersector() = default;

    explicit Intersector(const std::vector<Point>& vertices,
                                   const std::vector<glm::u32vec3>& faces) {
        buildBVH(vertices, faces);
    }

    void buildBVH(const std::vector<Point>& vertices,
                  const std::vector<glm::u32vec3>& faces) {
        setupGeometry(vertices, faces);
        buildAccelerationStructure();
    }

    std::optional<Intersection> intersect(const Ray& ray) const {
        if (!bvh_) return std::nullopt;

        nanort::TriangleIntersector<double> intersector(vertices_.data(), faces_.data(), sizeof(double) * 3);
        nanort::TriangleIntersection<double> isect;

        bool hit = bvh_->Traverse(ray.getNanortRay(), intersector, &isect);

        if (hit && ray.isInRange(isect.t)) {
            Intersection result;
            result.t = isect.t;
            result.face_id = isect.prim_id;
            result.point = ray.pointAt(result.t);
            result.normal = computeNormal(isect.prim_id);
            return result;
        }

        return std::nullopt;
    }

private:

    void setupGeometry(const std::vector<Point>& vertices,
                       const std::vector<glm::u32vec3>& faces) {
        vertices_.resize(vertices.size() * 3);
        for (size_t i = 0; i < vertices.size(); ++i) {
            vertices_[i * 3] = vertices[i].x;
            vertices_[i * 3 + 1] = vertices[i].y;
            vertices_[i * 3 + 2] = vertices[i].z;
        }

        faces_.resize(faces.size() * 3);
        for (size_t i = 0; i < faces.size(); ++i) {
            faces_[i * 3] = faces[i].x;
            faces_[i * 3 + 1] = faces[i].y;
            faces_[i * 3 + 2] = faces[i].z;
        }
    }

    void buildAccelerationStructure() {
        nanort::TriangleMesh<double> triangle_mesh(vertices_.data(), faces_.data(), sizeof(double) * 3);
        nanort::TriangleSAHPred<double> triangle_pred(vertices_.data(), faces_.data(), sizeof(double) * 3);

        bvh_ = std::make_shared<nanort::BVHAccel<double>>();

        nanort::BVHBuildOptions<double> build_options;
        build_options.cache_bbox = false;

        bool success = bvh_->Build(static_cast<unsigned int>(faces_.size() / 3),
                                   triangle_mesh, triangle_pred, build_options);

        if (!success) {
            throw std::runtime_error("Failed to build BVH");
        }
    }

    Vector computeNormal(uint32_t face_id) const {
        unsigned int idx = face_id * 3;
        unsigned int v0_idx = faces_[idx] * 3;
        unsigned int v1_idx = faces_[idx + 1] * 3;
        unsigned int v2_idx = faces_[idx + 2] * 3;

        Point v0(vertices_[v0_idx], vertices_[v0_idx + 1], vertices_[v0_idx + 2]);
        Point v1(vertices_[v1_idx], vertices_[v1_idx + 1], vertices_[v1_idx + 2]);
        Point v2(vertices_[v2_idx], vertices_[v2_idx + 1], vertices_[v2_idx + 2]);

        Vector edge1 = v1 - v0;
        Vector edge2 = v2 - v0;
        return glm::normalize(glm::cross(edge1, edge2));
    }

    std::vector<double> vertices_;
    std::vector<unsigned int> faces_;
    std::shared_ptr<nanort::BVHAccel<double>> bvh_;

};

}

#endif // G3_INTERSECTOR_H
