#ifndef G3_IREPOSITORY_HPP
#define G3_IREPOSITORY_HPP

#include "materials/Material.hpp"
#include "etc/Function2.hpp"

namespace gamma3::io {

class IRepository {
public:

    virtual materials::Atom loadAtom(const std::string& symbol) = 0;
    virtual materials::Molecule loadMolecule(const std::string& name) = 0;
    virtual materials::Material loadMaterial(const std::string& name) = 0;
    virtual etc::Function2<double> loadFormula2(const std::string& name) = 0;

};

}

#endif // G3_IREPOSITORY_HPP
