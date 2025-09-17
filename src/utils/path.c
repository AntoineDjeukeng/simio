#include <string.h>
#include <stdio.h>
#include "../../include/utils/path.h"

static int ends_with(const char *s, const char *suf) {
    size_t n=strlen(s), m=strlen(suf);
    if (m>n) return 0;
    return strcmp(s+n-m, suf)==0;
}

int path_make_subset_gro(const char *gro_full, const char *suffix,
                         int preserve_ext, char *out, size_t outsz)
{
    if (!gro_full || !suffix || !out || outsz==0) return -1;
    const char *ext = NULL;
    size_t n = strlen(gro_full);
    if (preserve_ext) {
        if (ends_with(gro_full, ".gro.gz")) ext = ".gro.gz";
        else if (ends_with(gro_full, ".gro")) ext = ".gro";
    }
    char stem[1024];
    if (ext) {
        size_t stem_len = n - strlen(ext);
        if (stem_len >= sizeof(stem)) return -1;
        memcpy(stem, gro_full, stem_len);
        stem[stem_len] = '\0';
        int r = snprintf(out, outsz, "%s%s%s", stem, suffix, ext);
        return (r>0 && (size_t)r < outsz) ? 0 : -1;
    } else {
        int r = snprintf(out, outsz, "%s%s.gro", gro_full, suffix);
        return (r>0 && (size_t)r < outsz) ? 0 : -1;
    }
}