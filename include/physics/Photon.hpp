#ifndef G3_PHOTON_HPP
#define G3_PHOTON_HPP

#include "geometry/Ray.hpp"

namespace gamma3::physics {

class Photon {
public:

    Photon(geometry::Ray ray,
           double MeV) :
        ray_(std::move(ray)),
        photon_energy_MeV_(MeV) {}

    double energy() const { return photon_energy_MeV_; }

    geometry::Ray ray() const { return ray_; }

    void move(double l) { ray_.moveSource(l); }

    void setEnergy(double MeV) { photon_energy_MeV_ = MeV; }

    void setDirection(geometry::Vector dir) { ray_ = geometry::Ray(ray_.source(), dir); }

private:

    double photon_energy_MeV_;
    geometry::Ray ray_;

};

}

#endif // G3_PHOTON_HPP
