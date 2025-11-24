#ifndef G3_JSONREPOSITORY_HPP
#define G3_JSONREPOSITORY_HPP

#include "IRepository.hpp"

#include <nlohmann/json.hpp>

namespace gamma3::io {

class JsonRepository : public IRepository {
public:

    JsonRepository(std::string path) : path_(std::move(path)) {
        std::ifstream f(path_);
        if (!f) throw std::runtime_error("JsonRepository: cannot open " + path);
        f >> j_;
    }

    virtual materials::Atom loadAtom(const std::string& symbol) override {
        int Z = j_.at("Atoms").at(symbol).at("Z");
        double A = j_.at("Atoms").at(symbol).at("A");

        return materials::Atom(symbol, Z, A);
    }

    virtual materials::Molecule loadMolecule(const std::string& name) override {
        double mass = 0.0;
        std::vector<std::pair<materials::Atom, std::size_t>> comps;
        for (auto& c : j_.at("Molecules").at(name).at("composition")) {
            std::string elem = c.at("element");
            int count = c.at("count");

            materials::Atom atom = loadAtom(elem);
            mass += atom.A();

            comps.push_back({std::move(atom), std::move(count)});
        }

        return materials::Molecule(name, mass, std::move(comps));
    }

    virtual materials::Material loadMaterial(const std::string& name) override {
        const auto j = j_.at("Materials").at(name);

        double density = j.at("density");

        double fraction_sum = 0.0;
        std::vector<std::pair<materials::Molecule, double>> comps;
        for (auto& c : j.at("composition")) {

            std::string molecule_name = c.at("molecule");
            double mass_fraction = c.at("mass_fraction");

            fraction_sum += mass_fraction;
            if(fraction_sum > 1.0) throw std::runtime_error("MaterialLoader: mass_fraction > 1.0");

            materials::Molecule molecule = loadMolecule(molecule_name);
            comps.push_back({std::move(molecule), std::move(mass_fraction)});
        }

        return materials::Material(name, density, std::move(comps));

    }

    virtual etc::Function2<double> loadFormula2(const std::string& name) override {
        auto j = j_.at("Formulas").at(name);
        std::string x1 = j.at("vars")[0];
        std::string x2 = j.at("vars")[1];
        std::string formula = j.at("formula");
        return etc::Function2<double>(formula, {x1, x2});
    }

private:

    const std::string path_;
    nlohmann::json j_;

};

}

#endif // G3_JSONREPOSITORY_HPP
