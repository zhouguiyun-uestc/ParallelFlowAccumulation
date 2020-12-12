# ParallelFlowAccumulation  
**Manuscript Title**: Parallel Multi-Flow Accumulation for Massive Digital ELevation Models   
**Author**: Guiyun Zhou, Lihui Song, Wenyan Dong    
**Corresponding Author**: Guiyun Zhou(zhouguiyun@uestc.edu.cn)
# Abstract
Flow accumulation matrices are important input parameters in many hydrological and geospatial applications. With the increasing volume of digital elevation models (DEMs), there is a growing need to speed up the computation of flow accumulations and to meet the challenges in processing massive DEMs. This study proposes a new sequential multi-flow accumulation algorithm and a parallel algorithm for multi-flow accumulation calculation based on the three-step parallel framework. When a cell has multiple immediate downstream cells are encountered and all of its upstream cells have been processed, the proposed sequential algorithm conducts a breadth-first search of its downstream cells. The proposed parallel algorithm uses the proposed sequential algorithm to process each tile, sends necessary data to the producer to construct a global multi-flow graph, and obtains the final flow accumulation values of each tile based on the global offsets calculated in the global graph. Compared with the existing sequential algorithm, the calculation speed of our sequential algorithm is improved by approximately 20%. For the proposed parallel algorithm, the speed-up ratios generally increase with more consumer processes and the scaling efficiencies decrease with more consumer processes. Our proposed parallel algorithm can process massive DEMs that cannot be successfully processed using existing sequential algorithms

# Prerequisite  
GDAL, MPI and cereal  
The packages can be installed using standard install procedures. For example, to install cereal, the following command can be used on Ubuntu:  
```
sudo apt install libcereal-dev
```

# Compilation  
To compile and programs run:  
```
mkdir build  
cd build  
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..  
make
```  
The result is a program called `ParallelFlowAccum`.  
# Program arguments  
The program can be run with four modes: `parallel`, `recursive`, `wang`, and `proposed_sequntial`. In `parallel` mode, the program derives flow accumulation in parallel. In `recursive` mode, the program derives flow accumulation using recursive algorithm. In `wang` mode, the program derives flow accumulation using wang's algorithm. In `proposed_sequential` mode, the program derives flow accumulation using our proposed sequential algorithm.  
In `parallel` mode, the program has the following arguments:
```
mpirun -np <PROCESSES_NUMBER>  ParaFlowAccum parallel <INPUT_DEM> <INPUT_DIR> <NAME> <OUTPUT_ACCUM>
```
The `<INPUT_DEM>` argument is a folder path which contains the DEM files.  
The `<INPUT_DIR>` argument is a folder path which contains the flow directions files.  
The `<NAME>` argument is a string that is the longest string in DEM and flow directions files with the same name.  
The `<OUTPUT_ACCUM>` argument specifes the output folder.  

An example command is:  
```
mpirun -np 3 ParaFlowAccum ./test_data/freeborn/dem ./test_data/freeborn/dir "freeborn9" ./test_data/freeborn/accum  
```
`-np 3` indicates that the program is run in parallel over 3 processes, which includes 1 producer process and 2 consumer processes.     
In `recursive` mode, the program has the following arguments: 
```
ParaFlowAccum recursive <INPUT_DEM_FILE> <INPUT_DIR_FILE> <OUTPUT_ACCUM_FILE>
```
In `wang` mode, the program has the following arguments: 
```
ParaFlowAccum wang <INPUT_DEM_FILE> <INPUT_DIR_FILE> <OUTPUT_ACCUM_FILE>
```
In `proposed_sequential` mode, the program has the following arguments: 
```
ParaFlowAccum proposed_sequential <INPUT_DEM_FILE> <INPUT_DIR_FILE> <OUTPUT_ACCUM_FILE>
```
The `<INPUT_DEM_FILE>` argument specifes the input DEM.  
The `<INPUT_DIR_FILE>` argument specifes the input flow_directions.  
The `<OUTPUT_ACCUM_FILE>` argument specifes the output flow_accum.  


# gridInfo.txt  
In fact, in `parallel` mode, the program reads the input data through the text file. In DEM folder, there is a `gridInfo.txt` which contains the information of tiles of DEM. In the girdInfo.txt   
the first line specifies the height of the tile.  
The second line specifies the width of the tile.   
The third line specifies the height of the grid.  
The fourth line specifies the width of the grid.  
The fifth line specifies the total height of the whole DEM.   
The sixth line specifies the total width of the whole DEM.   
The seventh line specifies the cell size.  
similarity, there is a gridInfo.txt which contains the information of tiles of the flow directions in flow directions folder.  

# Test data
The test_data folder contains one test datasets that can be used with the program.  



