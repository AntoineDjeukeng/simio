#ifndef UTILS_PATH_H
#define UTILS_PATH_H
#include <stddef.h>
int path_make_subset_gro(const char *gro_full, const char *suffix,
                         int preserve_ext, char *out, size_t outsz);
#endif