#ifndef G3_SURFACE_HPP
#define G3_SURFACE_HPP

#include <vector>
#include <utility>
#include <string>
#include <fstream>

#include "Ray.hpp"
#include "Intersector.h"

#include <glm/gtc/matrix_transform.hpp>

namespace gamma3::geometry {

class Surface {
public:

    Surface(const std::vector<Point>& vertices, const std::vector<glm::u32vec3>& faces)
        : vertices_(vertices), faces_(faces) {
        validateMesh();
        ensureIntersector();
    }

    // ==================== ТРАНСФОРМАЦИИ ====================
    void translate(const Vector& translation) {
        glm::dmat4 transform = glm::translate(glm::dmat4(1.0), translation);
        applyTransform(transform);
        invalidateIntersector();
    }

    void rotate(double angle_rad, const Vector& axis) {
        glm::dmat4 transform = glm::rotate(glm::dmat4(1.0), angle_rad, glm::normalize(axis));
        applyTransform(transform);
        invalidateIntersector();
    }

    void scale(double factor) {
        glm::dmat4 transform = glm::scale(glm::dmat4(1.0), glm::dvec3(factor));
        applyTransform(transform);
        invalidateIntersector();
    }

    void scale(const glm::dvec3& factors) {
        glm::dmat4 transform = glm::scale(glm::dmat4(1.0), factors);
        applyTransform(transform);
        invalidateIntersector();
    }

    // ==================== ЭКСПОРТ ====================
    void saveOBJ(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        writeOBJ(file);
    }

    std::size_t writeOBJ(std::ofstream& os, std::size_t offset = 0) const {
        for (const auto& vertex : vertices_) {
            os << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
        }

        for (const auto& face : faces_) {
            os << "f " << (face.x + offset) << " " << (face.y + offset) << " " << (face.z + offset) << "\n";
        }

        return offset + vertices_.size();
    }

    // ==================== ГЕОМЕТРИЧЕСКИЕ СВОЙСТВА ====================
    double area() const {
        double total_area = 0.0;
        for (const auto& face : faces_) {
            const Point& v0 = vertices_[face.x];
            const Point& v1 = vertices_[face.y];
            const Point& v2 = vertices_[face.z];

            Vector edge1 = v1 - v0;
            Vector edge2 = v2 - v0;
            total_area += 0.5 * glm::length(glm::cross(edge1, edge2));
        }
        return total_area;
    }

    double volume() const {
        if (!isClosed()) {
            throw std::runtime_error("Cannot compute volume for open surface");
        }

        double vol = 0.0;
        for (const auto& face : faces_) {
            const Point& v0 = vertices_[face.x];
            const Point& v1 = vertices_[face.y];
            const Point& v2 = vertices_[face.z];

            double triple_product = glm::dot(v0, glm::cross(v1, v2));
            vol += triple_product;
        }

        return std::abs(vol) / 6.0;
    }

    bool contains(const Point& point) const {
        if (!isClosed()) {
            throw std::runtime_error("Point containment check requires closed surface");
        }

        Ray ray(point, Vector(1.0, 0.0, 0.0));

        int intersection_count = 0;
        bool continue_search = true;

        while (continue_search) {
            auto intersection = intersect(ray);
            if (intersection) {
                intersection_count++;
                Ray new_ray(intersection->point + intersection->normal * 1e-6, ray.direction());
                ray = new_ray;
            } else {
                continue_search = false;
            }
        }

        return (intersection_count % 2) == 1;
    }

    bool isClosed() const {
        std::unordered_map<std::pair<uint32_t, uint32_t>, int, PairHash> edge_count;

        for (const auto& face : faces_) {
            addEdge(edge_count, face.x, face.y);
            addEdge(edge_count, face.y, face.z);
            addEdge(edge_count, face.z, face.x);
        }

        for (const auto& [edge, count] : edge_count) {
            if (count != 2) return false;
        }

        return true;
    }

    // ==================== ТРАССИРОВКА ЛУЧЕЙ ====================
    std::optional<Intersection> intersect(const Ray& ray) const {
        return intersector_->intersect(ray);
    }


    // ==================== ДОСТУП К ДАННЫМ ====================
    const std::vector<Point>& vertices() const { return vertices_; }
    const std::vector<glm::u32vec3>& faces() const { return faces_; }

    Point centroid() const {
        Point center(0.0);
        for (const auto& vertex : vertices_) {
            center += vertex;
        }
        return center / static_cast<double>(vertices_.size());
    }

private:
    std::vector<Point> vertices_;
    std::vector<glm::u32vec3> faces_;
    std::shared_ptr<Intersector> intersector_;

    struct PairHash {
        template <typename T, typename U>
        std::size_t operator()(const std::pair<T, U>& p) const {
            return std::hash<T>()(p.first) ^ std::hash<U>()(p.second);
        }
    };

    void validateMesh() {
        if (vertices_.empty()) {
            throw std::invalid_argument("Surface must have vertices");
        }
        if (faces_.empty()) {
            throw std::invalid_argument("Surface must have faces");
        }

        for (const auto& face : faces_) {
            if (face.x >= vertices_.size() || face.y >= vertices_.size() || face.z >= vertices_.size()) {
                throw std::invalid_argument("Face index out of bounds");
            }
        }
    }

    void applyTransform(const glm::dmat4& transform) {
        for (auto& vertex : vertices_) {
            glm::dvec4 transformed = transform * glm::dvec4(vertex, 1.0);
            vertex = glm::dvec3(transformed);
        }
    }

    void invalidateIntersector() {
        intersector_.reset();
    }

    void ensureIntersector() {
        if (!intersector_) {
            intersector_ = std::make_shared<Intersector>(vertices_, faces_);
        }
    }

    void addEdge(std::unordered_map<std::pair<uint32_t, uint32_t>, int, PairHash>& edge_count,
                 uint32_t a, uint32_t b) const {
        auto edge = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
        edge_count[edge]++;
    }
};

}

#endif // G3_SURFACE_HPP
