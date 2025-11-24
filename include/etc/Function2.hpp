#ifndef G3_FUNCTION2_HPP
#define G3_FUNCTION2_HPP

#include "exprtk.hpp"

namespace gamma3::etc {

template<typename T>
class Function2 {
public:

    Function2(const std::string& expression_str, const std::array<std::string, 2>& vars) {
        symbol_table.add_variable(vars[0], x);
        symbol_table.add_variable(vars[1], y);
        symbol_table.add_constants();

        expression.register_symbol_table(symbol_table);

        if (!parser.compile(expression_str, expression)) {
            throw std::runtime_error("Ошибка компиляции: " + parser.error());
        }
    }

    T operator()(T x_val, T y_val) {
        x = x_val;
        y = y_val;
        return expression.value();
    }

private:

    exprtk::symbol_table<T> symbol_table;
    exprtk::expression<T> expression;
    exprtk::parser<T> parser;

    T x{}, y{};

};

}

#endif // G3_FUNCTION2_HPP
