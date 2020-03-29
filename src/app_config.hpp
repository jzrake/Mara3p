/**
 ==============================================================================
 Copyright 2019 - 2020, Jonathan Zrake

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 ==============================================================================
*/




#pragma once
#include <iomanip>
#include <map>
#include <string>
#include <variant>




//=============================================================================
namespace mara
{
    using config_parameter_t     = std::variant<int, double, std::string>;
    using config_parameter_map_t = std::map<std::string, config_parameter_t>;
    using config_string_map_t    = std::map<std::string, std::string>;

    //=========================================================================
    class config_t;
    class config_template_t;

    //=========================================================================
    inline auto argv_to_string_map(int argc, const char* argv[]);
    inline auto config_template();
    inline void pretty_print(std::ostream& os, std::string header, const config_t& parameters);
};




//=============================================================================
class mara::config_t
{
public:

    //=========================================================================
    config_t() {}
    config_t(config_parameter_map_t parameters, config_string_map_t usage_desc)
    : parameters(parameters)
    , template_items(parameters)
    , usage_desc(usage_desc) {}

    auto get_int(std::string key) const
    {
        return get<int>(key);
    }

    auto get_double(std::string key) const
    {
        return get<double>(key);
    }
    
    auto get_string(std::string key) const
    {
        return get<std::string>(key);
    }

    const auto& at(std::string key) const
    {
        return parameters.at(key);
    }

    const auto& help_at(std::string key) const
    {
        return usage_desc.at(key);
    }

    template<typename ValueType>
    ValueType get(std::string key) const
    {
        if (! template_items.count(key))
        {
            throw std::invalid_argument("config has no option " + key);
        }
        try {
            return std::get<ValueType>(parameters.at(key));
        }
        catch (const std::bad_variant_access&)
        {
            throw std::invalid_argument("mara::config_t (wrong type for key '" + key + "')");
        }
    }

    template<typename Mapping>
    config_t& update(const Mapping& parameters)
    {
        for (auto item : parameters)
        {
            set(item.first, item.second);
        }
        return *this;
    }

    config_t& set(std::string key, std::string value)
    {
        if (! template_items.count(key))
        {
            throw std::invalid_argument("config has no option " + key);
        }
        switch (template_items.at(key).index())
        {
            case 0: return set(key, config_parameter_t(std::stoi(value)));
            case 1: return set(key, config_parameter_t(std::stod(value)));
            case 2: return set(key, config_parameter_t(value));
        }
        return *this;
    }

    config_t& set(std::string key, config_parameter_t value)
    {
        if (! template_items.count(key))
        {
            throw std::invalid_argument("config has no option " + key);
        }
        if (value.index() != template_items.at(key).index())
        {
            throw std::invalid_argument("config got wrong data type for option " + key);
        }
        parameters[key] = value;

        return *this;
    }

    auto begin() const
    {
        return parameters.begin();
    }
    
    auto end() const
    {
        return parameters.end();
    }

private:
    //=========================================================================
    config_parameter_map_t parameters;
    config_parameter_map_t template_items;
    config_string_map_t usage_desc;
};




//=============================================================================
class mara::config_template_t
{
public:

    //=========================================================================
    config_template_t() {}
    config_template_t(config_parameter_map_t parameters) : parameters(parameters) {}

    config_template_t& item(std::string key, config_parameter_t default_value, std::string usage="")
    {
        if (parameters.count(key))
        {
            throw std::invalid_argument("mara::config_template_t::item (option already exists in config template)");
        }
        parameters[key] = default_value;
        usage_desc[key] = usage;
        return *this;
    }

    config_template_t& import(const config_template_t& other)
    {
        for (const auto& [k, v] : other.parameters)
        {
            item(k, v, other.usage_desc.at(k));
        }
        return *this;
    }

    config_t create() const
    {
        return config_t(parameters, usage_desc);
    }

    const auto& at(std::string key) const
    {
        return parameters.at(key);
    }

private:
    //=========================================================================
    config_parameter_map_t parameters;
    config_string_map_t usage_desc;
};




//=============================================================================
auto mara::argv_to_string_map(int argc, const char* argv[])
{
    config_string_map_t items;

    for (int n = 0; n < argc; ++n)
    {
        std::string arg = argv[n];
        std::string::size_type eq_index = arg.find('=');

        if (eq_index != std::string::npos)
        {
            std::string key = arg.substr(0, eq_index);
            std::string val = arg.substr(eq_index + 1);

            if (items.count(key))
            {
                throw std::invalid_argument("mara::argv_to_string_map (duplicate parameter " + key + ")");
            }
            items[key] = val;
        }
    }
    return items;
}

auto mara::config_template()
{
    return config_template_t();
}

void mara::pretty_print(std::ostream& os, std::string header, const config_t& parameters)
{
    using std::left;
    using std::setw;
    using std::setfill;

    os << std::string(52, '=') << "\n";
    os << header << ":\n\n";

    std::ios orig(nullptr);
    orig.copyfmt(os);

    for (const auto& [k, v] : parameters)
    {
        auto put_value = [&os] (auto v) { os << v; };

        os << '\t' << left << setw(24) << setfill('.') << k << ' ';
        os << setw(12) << setfill(' ');

        std::visit(put_value, v);

        os << parameters.help_at(k);
        os << '\n';
    }

    os << '\n';
    os.copyfmt(orig);
}
