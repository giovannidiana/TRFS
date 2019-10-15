CC	=	g++
CFLAGS	=	-ansi -pedantic -Wall -g
LIBS    =       -lgsl -lblas -I/usr/include/boost -L/usr/lib/x86_64-linux-gnu -larmadillo
##############################################################
# platform specific redefines (to use, set ARCH appropriately)

## ARCH = X86_64

#EXECUTABLES = genNFL genXTF genCTM 
EXECUTABLES = gen 
OBJECT_FILES= 
#HPP=include/functions.hpp include/model_*.hpp 
HPP=include/functions.hpp include/solveCTM.hpp include/constants.hpp

all: $(EXECUTABLES)

$(EXECUTABLES):	% : %.cpp $(OBJECT_FILES) $(HPP)
	$(CC) $(CFLAGS) -o $@ $< $(OBJECT_FILES) $(LIBS)

$(OBJECT_FILES): %.o : %.cpp 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm gen
	rm *.o


