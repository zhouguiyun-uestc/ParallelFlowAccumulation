
#ifndef PARADEM_CONSUMER_I_H
#define PARADEM_CONSUMER_I_H

#include <paradem/grid_info.h>
#include <paradem/i_consumer_2_producer.h>
#include <paradem/i_producer_2_consumer.h>
#include <paradem/tile_info.h>

class IConsumer : public IObject {
public:
    virtual bool processRound1(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                               IConsumer2Producer* pIC2P) = 0;
    virtual bool processRound2(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                               IProducer2Consumer* pIP2C) = 0;
};

#endif
