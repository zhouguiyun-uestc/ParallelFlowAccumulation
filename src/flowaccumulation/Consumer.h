#ifndef PARADEM_CONSUMER_H
#define PARADEM_CONSUMER_H

#include "accumulation.h"
#include "consumer_2_producer.h"
#include "paradem/timer.h"
#include "producer_2_consumer.h"

#include <paradem/gdal.h>
#include <paradem/grid.h>
#include <paradem/grid_info.h>
#include <paradem/i_consumer.h>
#include <paradem/raster.h>
#include <paradem/tile_info.h>

#include <iostream>
#include <memory>

class Consumer : public IConsumer {
public:
    Timer timer_io;
    Timer timer_calc;

public:
    virtual bool processRound1(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                               IConsumer2Producer* pIC2P);
    virtual bool processRound2(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                               IProducer2Consumer* pIP2C);
    virtual ~Consumer() = default;

public:
    virtual void free();
    void getEdgeLink(const int tileHeight, const int tileWidth, Raster<uint8_t>& flowdirs, Raster<float>& dem, const int cellSize, Consumer2Producer* pC2P);
    bool isInnerEdge(int row, int col, int height, int width);
    void calculateStatistics(Raster<double>& accum, double* min, double* max, double* mean, double* stdDev);
};

#endif
