NAME := ../build/reader/libreader.a
CC := gcc
AR := ar
RM := rm -f
MKDIR_P := mkdir -p

INCS := -I../include -Iinclude
CFLAGS := -O2 -Wall -Wextra -Wpedantic $(INCS)

# ---------- Colors ----------
GREEN  := \033[1;32m
YELLOW := \033[1;33m
RED    := \033[1;31m
RESET  := \033[0m

CFILES := \
	src/trr.c \
	src/xtc.c

OBJDIR := ../build/reader
OBJS := $(CFILES:%.c=$(OBJDIR)/%.o)
DEPS := $(OBJS:.o=.d)

.PHONY: all clean fclean re

all: $(NAME)
	@printf "%bMake success: %s%b\n" "$(GREEN)" "$(NAME)" "$(RESET)"

$(NAME): $(OBJS)
	@printf "%b[AR]%b %s\n" "$(YELLOW)" "$(RESET)" "$@"
	$(AR) rcs $@ $^

$(OBJDIR)/%.o: %.c
	@$(MKDIR_P) $(dir $@)
	@printf "%b[CC]%b %s\n" "$(YELLOW)" "$(RESET)" "$<"
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	@printf "%b[CLEAN]%b reader objects\n" "$(YELLOW)" "$(RESET)"
	$(RM) $(OBJS) $(DEPS)

fclean: clean
	@printf "%b[FCLEAN]%b %s\n" "$(YELLOW)" "$(RESET)" "$(NAME)"
	$(RM) $(NAME)

re: fclean all

-include $(DEPS)
