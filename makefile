LIBS  = -lRmath -lm
CC = g++

HBOOST = /home/duyv/workspace/source/boost_1_50_0
HRLIB = /usr/share/R/include
RLIB = /usr/local/lib
HMYDIR = /home/duyv/workspace/cpp/mixrerg/src

CFLAGS = -O0 -g3 -Wall

COMPILE = $(CC) -I$(HBOOST) -I$(HRLIB) -I$(HMYDIR) $(CFLAGS) -c

OBJFILES := $(patsubst src/util/%.cpp,src/util/%.o,$(wildcard src/util/*.cpp))	$(patsubst src/model/%.cpp,src/model/%.o,$(wildcard src/model/*.cpp))	$(patsubst src/model/binary/%.cpp,src/model/binary/%.o,$(wildcard src/model/binary/*.cpp))	$(patsubst src/model/signed/%.cpp,src/model/signed/%.o,$(wildcard src/model/signed/*.cpp))	$(patsubst src/var/%.cpp,src/var/%.o,$(wildcard src/var/*.cpp))	$(patsubst src/var/binary/%.cpp,src/var/binary/%.o,$(wildcard src/var/binary/*.cpp))	$(patsubst src/var/signed/%.cpp,src/var/signed/%.o,$(wildcard src/var/signed/*.cpp))	$(patsubst src/em/%.cpp,src/em/%.o,$(wildcard src/em/*.cpp))	$(patsubst src/%.cpp,src/%.o,$(wildcard src/*.cpp))

all: mixrerg

mixrerg: $(OBJFILES)

	$(CC) -L$(RLIB) -o mixrerg $(OBJFILES) $(LIBS)

%.o: %.cpp

	$(COMPILE) -o $@ $<  

clean:

	rm mixrerg $(OBJFILES)

