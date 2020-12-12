#pragma once
#include "Consumer.h"
#include "Node.h"
#include "accumulation.h"

#include <paradem/gdal.h>
#include <paradem/raster.h>

#include <chrono>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

const double SQRT2 = 1.4142135623730950488016887242097;
const double MAX_SINGLE = -1;

void wang(Raster<float>& dem, Raster<uint8_t> flowdDir, Raster<double>& accum, const double& cellSize, const std::vector<std::vector<double>>& index);
extern int BitCount(uint8_t n, int& k);
void improve(Raster<float> dem, Raster<uint8_t> flowdDir, Raster<double>& accum, const double& cellSize, const std::vector<std::vector<double>>& index);
double checkNeighbor(int row, int col, Raster<float>& dem, Raster<uint8_t>& flowdirs, const double& cellSize, const std::vector<std::vector<double>>& index,
                     const std::vector<std::vector<double>>& fractionSum, Raster<double>& accum);
void recursive(Raster<float> dem, Raster<uint8_t> flowdDir, const double& cellSize, Raster<double>& accum, const std::vector<std::vector<double>>& index);
bool computeAccum(std::string inputDEMFile, std::string inputDirFile, std::string outputAccuFile, char* method);
