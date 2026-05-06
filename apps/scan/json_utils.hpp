#ifndef JSON_UTILS_HPP
#define JSON_UTILS_HPP

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "json.hpp"

namespace json_utils
{
    using json = nlohmann::json;

    struct IntRange {
        int min;
        int max;
    };

    struct DoubleRange {
        double min;
        double max;
    };

    inline std::string path_to_string(const std::vector<std::string>& path)
    {
        std::string out;

        for (size_t i = 0; i < path.size(); ++i) {
            out += "/";
            out += path[i];
        }
        return out.empty() ? "/" : out;
    }

    inline std::string make_updated_filename(const std::string& input_filename)
    {
        std::string::size_type slash_pos = input_filename.find_last_of("/\\");
        std::string dir;
        std::string file;

        if (slash_pos == std::string::npos) {
            dir = "";
            file = input_filename;
        } else {
            dir = input_filename.substr(0, slash_pos + 1);
            file = input_filename.substr(slash_pos + 1);
        }

        std::string::size_type dot_pos = file.find_last_of('.');
        if (dot_pos == std::string::npos || dot_pos == 0) {
            return dir + file + "_updated";
        }

        std::string base = file.substr(0, dot_pos);
        std::string ext = file.substr(dot_pos);
        return dir + base + "_updated" + ext;
    }

    /* ============================================================
    ** LOAD / WRITE
    ** ============================================================ */

    inline json load_json(const std::string& filename)
    {
        std::ifstream file(filename.c_str());
        if (!file)
            throw std::runtime_error("cannot open input file: " + filename);

        json j;
        file >> j;
        return j;
    }

    inline void write_json(const json& j, const std::string& filename)
    {
        std::ofstream file(filename.c_str());
        if (!file)
            throw std::runtime_error("cannot open output file: " + filename);

        file << j.dump(4) << '\n';
    }

    inline bool maybe_write_json(const json& j,
                                 const std::string& input_filename,
                                 bool write_enabled,
                                 std::string* output_filename = NULL)
    {
        if (!write_enabled)
            return false;

        std::string out = make_updated_filename(input_filename);
        write_json(j, out);

        if (output_filename)
            *output_filename = out;

        return true;
    }

    /* ============================================================
    ** PATH ACCESS
    ** ============================================================ */

    inline const json* get_path(const json& j, const std::vector<std::string>& path)
    {
        const json* cur = &j;

        for (size_t i = 0; i < path.size(); ++i) {
            const std::string& key = path[i];

            if (cur->is_object()) {
                json::const_iterator it = cur->find(key);
                if (it == cur->end())
                    return NULL;
                cur = &(*it);
            }
            else if (cur->is_array()) {
                size_t idx = 0;
                try {
                    idx = static_cast<size_t>(std::stoul(key));
                }
                catch (...) {
                    return NULL;
                }

                if (idx >= cur->size())
                    return NULL;
                cur = &(*cur)[idx];
            }
            else {
                return NULL;
            }
        }
        return cur;
    }

    inline json* get_path_mut(json& j, const std::vector<std::string>& path)
    {
        json* cur = &j;

        for (size_t i = 0; i < path.size(); ++i) {
            const std::string& key = path[i];

            if (cur->is_object()) {
                json::iterator it = cur->find(key);
                if (it == cur->end())
                    return NULL;
                cur = &(*it);
            }
            else if (cur->is_array()) {
                size_t idx = 0;
                try {
                    idx = static_cast<size_t>(std::stoul(key));
                }
                catch (...) {
                    return NULL;
                }

                if (idx >= cur->size())
                    return NULL;
                cur = &(*cur)[idx];
            }
            else {
                return NULL;
            }
        }
        return cur;
    }

    /* ============================================================
    ** FIND
    ** ============================================================ */

    template <typename T>
    T find_value(const json& j, const std::vector<std::string>& path)
    {
        const json* node = get_path(j, path);

        if (!node)
            throw std::runtime_error("missing JSON path: " + path_to_string(path));

        try {
            return node->get<T>();
        }
        catch (const std::exception&) {
            throw std::runtime_error("type mismatch at JSON path: " + path_to_string(path));
        }
    }

    template <typename T>
    T find_value_or(const json& j,
                    const std::vector<std::string>& path,
                    const T& default_value)
    {
        const json* node = get_path(j, path);

        if (!node)
            return default_value;

        try {
            return node->get<T>();
        }
        catch (const std::exception&) {
            return default_value;
        }
    }

    /* ============================================================
    ** UPDATE
    ** ============================================================ */

    template <typename T>
    void update_value(json& j, const std::vector<std::string>& path, const T& value)
    {
        json* node = get_path_mut(j, path);

        if (!node)
            throw std::runtime_error("cannot update missing JSON path: " + path_to_string(path));

        *node = value;
    }

    /* ============================================================
    ** ADD
    ** ============================================================ */

    inline json& ensure_object_path(json& j, const std::vector<std::string>& path)
    {
        json* cur = &j;
        std::vector<std::string> partial_path;

        for (size_t i = 0; i < path.size(); ++i) {
            const std::string& key = path[i];
            partial_path.push_back(key);

            if (cur->is_null())
                *cur = json::object();

            if (!cur->is_object()) {
                throw std::runtime_error(
                    "cannot create key inside non-object at: " + path_to_string(partial_path)
                );
            }

            cur = &((*cur)[key]);
        }
        return *cur;
    }
    template <typename T>
    void add_value(json& j, const std::vector<std::string>& path, const T& value)
    {
        json& node = ensure_object_path(j, path);
        node = value;
    }

    inline void add_type(json& j,
                         const std::string& name,
                         int molecules,
                         int atoms_per_molecule)
    {
        if (!j.contains("types") || !j["types"].is_array())
            j["types"] = json::array();

        json entry;
        entry["name"] = name;
        entry["molecules"] = molecules;
        entry["atoms_per_molecule"] = atoms_per_molecule;

        j["types"].push_back(entry);
    }

    /* ============================================================
    ** DOMAIN HELPERS
    ** ============================================================ */

    inline DoubleRange find_channel_x_range(const json& j)
    {
        const json* min_arr = get_path(j, {"channel", "min"});
        const json* max_arr = get_path(j, {"channel", "max"});

        if (!min_arr || !max_arr)
            throw std::runtime_error("missing channel/min or channel/max");

        if (!min_arr->is_array() || !max_arr->is_array() ||
            min_arr->size() < 1 || max_arr->size() < 1)
        {
            throw std::runtime_error("channel.min/max are invalid");
        }

        return DoubleRange{
            (*min_arr).at(0).get<double>(),
            (*max_arr).at(0).get<double>()
        };
    }

    inline IntRange find_molecule_range(const json& j, const std::string& mol_name)
    {
        const json* types = get_path(j, {"types"});

        if (!types || !types->is_array())
            throw std::runtime_error("\"types\" is not an array");

        int start = 0;

        for (size_t i = 0; i < types->size(); ++i) {
            const json& t = (*types)[i];

            std::string name = find_value<std::string>(t, {"name"});
            int molecules = find_value<int>(t, {"molecules"});

            if (name == mol_name) {
                if (molecules <= 0)
                    throw std::runtime_error("molecule count is invalid for " + mol_name);
                return IntRange{start, start + molecules - 1};
            }

            start += molecules;
        }

        throw std::runtime_error("molecule name not found: " + mol_name);
    }

    inline IntRange find_atom_range(const json& j, const std::string& mol_name)
    {
        const json* types = get_path(j, {"types"});

        if (!types || !types->is_array())
            throw std::runtime_error("\"types\" is not an array");

        int start = 0;

        for (size_t i = 0; i < types->size(); ++i) {
            const json& t = (*types)[i];

            std::string name = find_value<std::string>(t, {"name"});
            int molecules = find_value<int>(t, {"molecules"});
            int atoms_per_molecule = find_value<int>(t, {"atoms_per_molecule"});

            int natoms = molecules * atoms_per_molecule;

            if (name == mol_name) {
                if (natoms <= 0)
                    throw std::runtime_error("atom count is invalid for " + mol_name);
                return IntRange{start, start + natoms - 1};
            }

            start += natoms;
        }

        throw std::runtime_error("molecule name not found: " + mol_name);
    }

} // namespace json_utils

#endif