NAME        := mdconv
TEST_NAME   := mdconv_test
BUILD_DIR   := build

CXX         := c++
CXXFLAGS    := -std=c++17 -Wall -Wextra -Werror -pedantic
CPPFLAGS    := -Iinclude

include source.mk

OBJ         := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRC))
TEST_OBJ    := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRC))

DEBUG_FLAGS := -O0 -g3 -DDEBUG -fsanitize=address,undefined -fno-omit-frame-pointer
RELEASE_FLAGS := -O2 -DNDEBUG

.PHONY: all debug test clean fclean re

all: CXXFLAGS += $(RELEASE_FLAGS)
all: $(NAME)

debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(NAME)

$(NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@

test: CXXFLAGS += $(DEBUG_FLAGS)
test: $(TEST_NAME)
	LSAN_OPTIONS=detect_leaks=0 ./$(TEST_NAME)

$(TEST_NAME): $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $(TEST_OBJ) -o $@

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(NAME) $(TEST_NAME)

fclean: clean

re: fclean all
