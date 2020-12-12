#ifndef PARADEM_ACCUMULATION_H
#define PARADEM_ACCUMULATION_H

#include "Node.h"
#include "atype.h"
#include "producer_2_consumer.h"

#include <paradem/raster.h>

#include <queue>
#include <stdint.h>

typedef double fraction_t;
typedef int link_t;

bool isOuterEdges(int row, int col, Raster<float>& dem);
bool isInnerEdges(int row, int col, Raster<float>& dem);
double getOffset1(int row, int col, int height, int width, Producer2Consumer* p2c);
bool is2RowCol(int row, int col, int i, int iRow, int iCol, uint8_t iDir);
void FlowAccumulationInnerTradition(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum);
void FlowAccumulationInnerImprove(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum);
double CheckNeighbor(int row, int col, Raster<float>& dem, Raster<uint8_t>& flowdirs, const double& cellSize, const std::vector<std::vector<double>>& index,
                     const std::vector<std::vector<double>>& fractionSum, Raster<double>& accum);
void FlowAccumulationInnerRecursion(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum);
double calculateIndex(int row, int col, int cellSize, Raster<float>& dem, Raster<uint8_t>& flowdirs);
std::vector<Node> calculateFra(int row, int col, double value, int cellSize, double index, Raster<float>& dem, Raster<uint8_t>& flowdirs);
bool handle(const int row, const int col, const std::vector<std::vector<bool>>& flag, Raster<uint8_t>& flowdirs, Raster<float>& dem);
void FollowPath(const int row0, const int col0, const int tileHeight, const int tileWidth, Raster<uint8_t>& flowdirs, Raster<float>& dem, std::vector<link_t>& links,
                std::vector<fraction_t>& fractions, int cellSize);
void FollowPath(Raster<uint8_t>& flowdirs, Raster<float>& dem, const int tileHeight, const int tileWidth, const int row, const int col, std::vector<link_t>& links);
extern int BitCount2(uint8_t n, int& k);

void FlowAccumulationTradition(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c);
void FlowAccumulationImprove(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c);
double CheckNeighbor(int row, int col, Raster<float>& dem, Raster<uint8_t>& flowDir, const double& cellSize, const std::vector<std::vector<double>>& index,
                     const std::vector<std::vector<double>>& fractionSum, Raster<double>& accum, Producer2Consumer* p2c);
void FlowAccumulationRecursion(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c);

#endif