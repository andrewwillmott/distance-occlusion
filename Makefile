TARGET = distocc
COMPILE_FLAGS = -O3

all: $(TARGET)

$(TARGET): DistanceOcclusionLib.cpp DistanceOcclusionLib.h DistanceOcclusionTool.cpp
	$(CXX) $(COMPILE_FLAGS) -o $(TARGET) DistanceOcclusionLib.cpp DistanceOcclusionTool.cpp

clean:
	$(RM) $(TARGET)
