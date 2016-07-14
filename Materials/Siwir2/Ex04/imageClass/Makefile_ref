CC = g++
CFLAGS = -ansi -Wall -pedantic -O3 $(INCLUDES)
LIBS =
INCLUDES =
TARGET = lbm

OBJS = Fileread.o \
	main.o \
	LBM.o 

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

Fileread.o: Fileread.cpp Fileread.h
	$(CC) -c $(CFLAGS) Fileread.cpp

LBM.o: LBM.cpp LBM.h Grid.h
	$(CC) -c $(CFLAGS) LBM.cpp

main.o: main.cpp LBM.h Timing.h
	$(CC) -c $(CFLAGS) main.cpp

.PHONY : clean
clean: rm -f *.o

