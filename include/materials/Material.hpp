#ifndef G3_MATERIAL_HPP
#define G3_MATERIAL_HPP

#include "Molecule.hpp"

namespace gamma3::materials {

class Material {
public:

    Material(std::string name,
             double density,
             std::vector<std::pair<Molecule, double>> components) :
        name_(std::move(name)),
        density_(density),
        components_(std::move(components)) {}

    double density() const { return density_; }

    const std::vector<std::pair<Molecule, double>>& components() const { return components_; }

    const std::string& name() const { return name_; }

private:

    std::string name_;
    double density_;
    std::vector<std::pair<Molecule, double>> components_;

};



}

#endif // G3_MATERIAL_HPP
