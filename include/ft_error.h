#ifndef FT_ERROR_H
#define FT_ERROR_H
enum {
    FT_OK = 0,
    FT_EOF = 1,
    FT_ERR = -1,
    FT_ENOMEM = -2,
    FT_EFORMAT = -3,
    FT_EINVAL = -4,
    FT_ENOTSUP = -5
};
#endif