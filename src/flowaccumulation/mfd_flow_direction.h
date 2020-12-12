#pragma once

#include <paradem/gdal.h>
#include <paradem/grid_info.h>
#include <paradem/raster.h>
#include <paradem/tool.h>

#include <deque>
#include <iostream>
#include <queue>
#include <string>

Raster<uint8_t> computeDir(Raster<float>& dem, Raster<uint8_t>& flowDir);
bool haveData(int row, int col, int gridHeight, int gridWidth);
void setborder(const int i, Raster<uint8_t>& bigDir, Raster<int>& iDir);
bool getBigData(Raster<int>& dir, Raster<uint8_t>& bigDir, const int tileRow, const int tileCol, const int gridHeight, const int gridWidth, const std::string pathDir, const std::string filenameDir);
void setborder(const int i, Raster<float>& bigDEM, Raster<float>& iDEM);
bool getBigData(Raster<float>& dem, Raster<float>& bigDEM, const int tileRow, const int tileCol, const int gridHeight, const int gridWidth, const std::string pathDEM, const std::string filenameDEM);
bool preprocess(const std::string& pathDir, const std::string& filenameDir, const std::string& pathDEM, const std::string& filenameDEM, const std::string& saveDirPath, const std::string& saveDEMPath,
                const std::string& savefilename);
bool checkDir(const std::string& pathDir, const std::string& filenameDir, const std::string& pathDEM);