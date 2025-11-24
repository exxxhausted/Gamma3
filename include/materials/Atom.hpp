#ifndef G3_ATOM_HPP
#define G3_ATOM_HPP

#include <string>

namespace gamma3::materials {

class Atom {
public:

    explicit Atom(std::string symbol,
                  int z,
                  double a) :
        symbol_(std::move(symbol)),
        Z_(z),
        A_(a) {}

    const std::string& symbol() const { return symbol_; }

    int Z() const { return Z_; }

    double A() const { return A_; }

private:

    const std::string symbol_;
    const int Z_;
    const double A_;

};


}

#endif // G3_ATOM_HPP
