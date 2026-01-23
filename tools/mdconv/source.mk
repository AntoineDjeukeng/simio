SRC := src/cli/main.cpp \
       src/gmx/parser.cpp \
       src/backend/writer.cpp \
       src/validate/validator.cpp \
       src/util/dump.cpp \
       src/gmx/gro.cpp \
       src/gmx/expander.cpp \
       src/normalize/resolver.cpp

LIB_SRC := src/gmx/parser.cpp \
           src/backend/writer.cpp \
           src/validate/validator.cpp \
           src/util/dump.cpp \
           src/gmx/gro.cpp \
           src/gmx/expander.cpp \
           src/normalize/resolver.cpp

TEST_SRC := tests/test_main.cpp $(LIB_SRC)
