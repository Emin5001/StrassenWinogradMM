# CC = icc
# CFLAGS = -I.
# LDFLAGS = -lpapi
# OPT_FLAGS = -O3 -fp-model precise
# BUILD = build
# SOURCES = w_strassen.c matrix.c
# OBJ = $(SOURCES:.c=.o)

# $(BUILD)/main: $(patsubst %,obj/%,$(OBJ)) 
# 	$(CC) $(CFLAGS) $(LDFLAGS) $(OPT_FLAGS) -o $@ $^

# .PHONY: clean
# clean:
# 	rm -f $(BUILD)/main $(patsubst %,obj/%,$(OBJ))

# .PHONY: compile
# compile: $(BUILD)/main 

# obj/%.o: %.c
# 	$(CC) $(CFLAGS) $(OPT_FLAGS) -c -o $@ $<

# VALGRIND = valgrind
# VALGRIND_FLAGS = --leak-check=yes

# .PHONY: valgrind
# valgrind: $(BUILD)/main
# 	$(VALGRIND) $(VALGRIND_FLAGS) $(BUILD)/main

CC = gcc
CFLAGS = -I.
BUILD = build
SOURCES = w_strassen.c matrix.c
OBJ = $(SOURCES:.c=.o)

$(BUILD)/main: $(patsubst %,obj/%,$(OBJ)) 
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(BUILD)/main $(patsubst %,obj/%,$(OBJ))

.PHONY: compile
compile: $(BUILD)/main 

obj/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

VALGRIND = valgrind
VALGRIND_FLAGS = --leak-check=yes

.PHONY: valgrind
valgrind: $(BUILD)/main
	$(VALGRIND) $(VALGRIND_FLAGS) $(BUILD)/main