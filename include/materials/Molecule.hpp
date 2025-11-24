#ifndef G3_MOLECULE_HPP
#define G3_MOLECULE_HPP

#include "Atom.hpp"

#include <vector>
#include <utility>

namespace gamma3::materials {

class Molecule {
public:

    Molecule(std::string name,
             double mass,
             std::vector<std::pair<Atom, std::size_t>> components) :
        name_(std::move(name)),
        mass_(mass),
        components_(std::move(components)) {}

    const std::string& name() const { return name_; }

    double mass() const { return mass_; }

    const std::vector<std::pair<Atom, std::size_t>>& components() const { return components_; }

private:

    std::string name_;
    double mass_;
    std::vector<std::pair<Atom, std::size_t>> components_;

};



}

#endif // G3_MOLECULE_HPP
