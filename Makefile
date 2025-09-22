CC       := cc
CFLAGS   := -Wall -Wextra -Werror -O2
CPPFLAGS :=
LDFLAGS  :=
LDLIBS   :=

INCDIR   := include
CPPFLAGS += -I$(INCDIR)

OBJDIR   := build/obj
BINDIR   := build/bin
SRC      := src/gro_parse.c src/gro_channel.c src/gro_print.c src/main.c
OBJ      := $(SRC:src/%.c=$(OBJDIR)/%.o)
BIN      := $(BINDIR)/app

all: $(BIN)

$(BIN): $(OBJ)
	@mkdir -p $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(OBJDIR)/%.o: src/%.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(OBJDIR)

fclean: clean
	rm -rf $(BINDIR)

re: fclean all

-include $(OBJ:.o=.d)
.PHONY: all clean fclean re
