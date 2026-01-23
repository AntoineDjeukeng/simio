#pragma once

#include <stdexcept>
#include <string>

namespace mdconv {

struct ParseError : public std::runtime_error {
    explicit ParseError(const std::string& msg) : std::runtime_error("ParseError: " + msg) {}
};

struct ValidationError : public std::runtime_error {
    explicit ValidationError(const std::string& msg) : std::runtime_error("ValidationError: " + msg) {}
};

struct WriteError : public std::runtime_error {
    explicit WriteError(const std::string& msg) : std::runtime_error("WriteError: " + msg) {}
};

} // namespace mdconv
