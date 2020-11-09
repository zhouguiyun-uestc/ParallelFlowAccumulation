CC=mpic++

#头文件1
h1=/ParallelFlowAccumulation/include/

#头文件2
h2=/ParallelFlowAccumulation/src/flowdir/

#头文件3
h3=/ParallelFlowAccumulation/
CFLAGS=-I$(h1) -I$(h2) -I$(h3) -std=c++11 -fpermissive -O3
LIBS=-lm -lgdal 
#.h文件
DEPS=$(shell find /ParallelFlowAccumulation/src/flowdir/ -name "*.h") $(shell find /ParallelFlowAccumulation/include/paradem/ -name "*.h") $(shell find /ParallelFlowAccumulation/ -name "*.hpp")

#不同文件目录下的.cpp文件
src=$(shell find /ParallelFlowAccumulation/src/common/ -name "*.cpp") $(shell find /ParallelFlowAccumulation/src/flowdir/ -name "*.cpp")

#所有的.o文件
OBJ = $(src:%.cpp=%.o)

%.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ParaFlowAccum:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(OBJ) ParaFlowAccum
