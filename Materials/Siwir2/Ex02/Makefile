CC = g++
CFLAGS =  -Wall  -Wshadow  -pedantic -std=c++11
LIBS = 
INCLUDES = 
TARGET = waveguide

SRC = waveguide.cpp
OBJ = $(patsubst %.c,%.o,$(SRC))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(LIBS)

debug: $(SRC)
	$(CC) -g $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SRC) $(LIBS)

.PHONY : clean depend 

clean:
	@/bin/rm -f $(OBJS)
	@/bin/rm -f $(TARGET)
	@/bin/rm -f mgsolve.output

depend:
	@makedepend -- $(CFLAGS) -- $(SRC)
