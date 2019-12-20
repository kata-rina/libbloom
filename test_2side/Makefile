SRC += main.c 2sidekmer.c
SRC += bloom.c MurmurHash2.c
INC += $(wildcard ./*.h) ../bloom.h ../murmur2/murmurhash2.h
OBJ = $(addsuffix .o, $(basename $(SRC)))
TARGET = test
BUILDIR = build
OBJS = $(addprefix $(BUILDIR)/, $(OBJ))
SRC_DIR = .
BLOOM_DIR = ..
MUR_DIR = ../murmur2

CC = gcc

INC_DIRS = -I../murmur2 -I..
CFLAGS = -Wall -Wextra

all: $(TARGET) $(INC)

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
