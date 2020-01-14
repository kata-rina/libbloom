SRC += main.c 2sidekmer.c mutate.c best_fit.c
SRC += bloom.c MurmurHash2.c
INC += $(wildcard ./*.h) ../bloom.h ../murmur2/murmurhash2.h
OBJ = $(addsuffix .o, $(basename $(SRC)))
TARGET = test
MUT_TARGET = mutate
MUT_OBJ = build/bloom.o build/MurmurHash2.o build/mutate.o
BUILDIR = build
OBJS = $(addprefix $(BUILDIR)/, $(OBJ))
SRC_DIR = .
BLOOM_DIR = ..
MUR_DIR = ../murmur2

CC = gcc

INC_DIRS = -I../murmur2 -I.. -I.
CFLAGS = -Wall -Wextra

all: $(TARGET) $(INC)

# $(MUT_TARGET)
#
# $(MUT_TARGET): $(MUT_OBJ)
# 	# mkdir -p $(BUILDIR)
	# $(CC) -g $^ -o $@ -lm

$(TARGET): $(OBJS)
	# mkdir -p $(BUILDIR)
	$(CC) -g $^ -o $@ -lm

$(BUILDIR)/%.o: %.c
	mkdir -p $(BUILDIR)
	$(CC) -g $(CFLAGS) $(INC_DIRS) -c $< -o $@

$(BUILDIR)/%.o: $(BLOOM_DIR)/%.c
	mkdir -p $(BUILDIR)
	$(CC) -g $(CFLAGS) $(INC_DIRS) -c $< -o $@

$(BUILDIR)/%.o: $(MUR_DIR)/%.c
	mkdir -p $(BUILDIR)
	$(CC) -g $(CFLAGS) $(INC_DIRS) -c $< -o $@

clean:
	rm -rf build
