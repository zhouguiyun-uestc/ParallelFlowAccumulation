# ParallelFlowAccumulation  
**Manuscript Title**: Parallel extraction algorithm on multiple flow accumulation in large-scale digital elevation model  
**Author**: Guiyun Zhou, Lihui Song, Wenyan Dong, Junjie Zhou  
**Corresponding Author**: Guiyun Zhou(zhouguiyun@uestc.edu.cn)
# Prerequisite  
GDAL, MPI and cereal  
# Compilation  
To compile and programs run:  
```
make
```  
The result is a program called `ParaFlowAccum`.  
# Program argumens
```
mpirun -np <PROCESSES_NUMBER>  ParaFlowAccum "evict"  <INPUT_DEM> <INPUT_DIR> <NAME> <OUTPUT_ACCUM>
```
The `<INPUT_DEM>` argument is a folder path which contains the DEM files.  
The `<INPUT_DIR>` argument is a folder path which contains the flow directions files.  
The `<NAME>` argument is a string that is the longest string in DEM and flow directions files with the same name.  
The `<OUTPUT_ACCUM>` argument specifes the output folder.  

An example command is: `mpirun -np 3 ParaFlowAccum "evict" ./test_data/freeborn/dem ./test_data/freeborn/dir "freeborn9" ./test_data/freeborn/accum  
`-np 3` indicates that the program is run in parallel over 3 processes, which includes 1 producer process and 2 consumer processes.   
# gridInfo.txt  
In fact, the program reads the input data through the text file. In DEM folder, there is a gridInfo.txt which contains the information of tiles of DEM. In the girdInfo.txt   
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



