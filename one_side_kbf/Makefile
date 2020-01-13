CC = gcc

PROJ_DIR = /mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/libbloom/one_side_kbf
INC_DIR = $(PROJ_DIR)/inc
SRC_DIR = $(PROJ_DIR)/src

TEST_FILE_DIR = /mnt/Jupiter/FAKS/Diplomski/3_semestar/Bioinformarika/genom

CFLAGS = -Wall -Wextra -c -g
INC = -I$(INC_DIR)

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(addsuffix .o, $(basename $(SRC)))
OBJ += test_prog.o

TARGET = test

all : $(TARGET)

test : $(OBJ)
				$(CC) -g $^ -o $@ -lm
				rm -f $(SRC_DIR)/*.o
				rm -f test_prog.o

test_prog.o : test_prog.c
				$(CC) $(CFLAGS) $(INC) $< -o $@ -lm

$(SRC_DIR)/%.o : $(SRC_DIR)/%.c
				$(CC) $(CFLAGS) $(INC) $< -o $@ -lm


clean:
	rm -f $(SRC_DIR)/*.o
	rm -f test_prog.o
	rm -f $(TEST_FILE_DIR)/mutate.txt
	rm -f $(TEST_FILE_DIR)/nonmutate.txt
	rm -f test
