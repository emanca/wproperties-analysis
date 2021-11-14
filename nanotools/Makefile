GCC=g++
CXXFLAGS=`root-config --libs --cflags` -O3 -fPIC -Wall 
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin
HPPDIR=interface
AUXDIR=aux
COREDIR=../RDFprocessor/framework
PATHFLAGS=-I$(COREDIR)/ -I$(COREDIR)/interface/ -I./interface/ -I/scratchnvme/wmass/eigen-3.3.7/ -I. -I/scratchnvme/wmass/CMSJMECalculators/include -I/scratchnvme/wmass/CMSJMECalculators/interface 
CORELIB=$(COREDIR)/lib/libRDFProcessorCore.so
JMELIB=/scratchnvme/wmass/CMSJMECalculators/lib/libCMSJMECalculatorsDict.so /scratchnvme/wmass/CMSJMECalculators/lib/libCMSJMECalculators.so
SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJ=$(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o , $(SRC)  )
HPPLINKDEF=$(patsubst $(SRCDIR)/%.cpp, ../interface/%.hpp , $(SRC)  )

.PHONY: all
all:
	$(info, "--- Full compilation --- ")	
	$(info, "-> if you want just to recompile something use 'make fast' ")	
	$(info, "------------------------ ")	
	$(MAKE) clean
	$(MAKE) libNanoTools.so

.PHONY: fast
fast:
	$(MAKE) libNanoTools.so

libNanoTools.so: $(OBJ) Dict | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(PATHFLAGS) $(RPATH) $(SOFLAGS) -o $(BINDIR)/$@ $(OBJ) $(BINDIR)/dict.o $(CORELIB)  $(JMELIB)

$(OBJ) : $(BINDIR)/%.o : $(SRCDIR)/%.cpp interface/%.hpp | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(PATHFLAGS) $(RPATH) -c -o $(BINDIR)/$*.o $<

.PHONY: Dict
Dict: $(BINDIR)/dict.o

$(BINDIR)/dict.o: $(SRC) | $(BINDIR)
	genreflex $(SRCDIR)/classes.h -s $(SRCDIR)/classes_def.xml -o $(BINDIR)/dict.cc --fail_on_warnings --rootmap=$(BINDIR)/dict.rootmap --rootmap-lib=libNanoTools.so $(PATHFLAGS)
	$(GCC) -c -o $(BINDIR)/dict.o $(CXXFLAGS) $(PATHFLAGS) $(RPATH) $(BINDIR)/dict.cc

$(BINDIR):
	mkdir -p $(BINDIR)
	mkdir -p $(AUXDIR)

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm $(BINDIR)/dict*
	-rm $(BINDIR)/NanoTools.so
	-rmdir $(BINDIR)
