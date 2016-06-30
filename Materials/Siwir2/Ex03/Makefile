CC = g++
CFLAGS = -Wall -Wshadow -std=c++11 -pedantic
LIBS = 
INCLUDES = 
TARGET = mdsim

SRC = mdsim.cpp
OBJ = $(patsubst %.c,%.o,$(SRC))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(LIBS)

debug: $(SRC)
	$(CC) -g $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(LIBS)

.PHONY : clean depend 

clean:
	rm -f *.vtk
	@/bin/rm -f $(OBJS)
	@/bin/rm -f $(TARGET)
	@/bin/rm -f mdsim.output

depend:
	@makedepend -- $(CFLAGS) -- $(SRC)
