CXXFLAGS ?= -O3 --std=c++11 # -fdiagnostics-absolute-paths
#CXXFLAGS ?= -O3 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -lomp

all: distocc

DOL_HEADERS = DistanceFields.hpp OcclusionFields.hpp MaskUtilities.hpp DOFImage.hpp
DOL_SOURCES = DistanceFields.cpp OcclusionFields.cpp MaskUtilities.cpp DOFImage.cpp
DOL_DEPS    = $(DOL_HEADERS) Makefile
DOL_OBJS    = $(DOL_SOURCES:.cpp=.o)
DOL_DOBJS   = $(DOL_SOURCES:.cpp=D.o)

libdistocc.a: $(DOL_OBJS)
	$(AR) rcs $@ $^

libdistoccd.a: $(DOL_DOBJS)
	$(AR) rcs $@ $^

%.o: %.cpp $(DOL_DEPS)
	$(CXX) $(CXXFLAGS) $(OPTS) -c $< -o $@

%D.o: %.cpp $(DOL_DEPS)
	$(CXX) $(CXXFLAGS) $(DBG_OPTS) -c $< -o $@

distocc: $(DOL_DEPS) libdistocc.a MeshSupport.hpp MeshSupport.cpp DistanceOcclusionTool.cpp
	$(CXX) $(CXXFLAGS) -o $@ DistanceOcclusionTool.cpp -L. -ldistocc

clean:
	$(RM) distocc *.o *.a
