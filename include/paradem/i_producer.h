#ifndef PARADEM_PRODUCER_I_H
#define PARADEM_PRODUCER_I_H

#include <paradem/grid.h>
#include <paradem/grid_info.h>
#include <paradem/i_consumer_2_producer.h>
#include <paradem/i_producer_2_consumer.h>
#include <paradem/tile_info.h>

#include <memory>

class IProducer : public IObject {
public:
    virtual void process(const GridInfo& gridInfo, const std::vector<TileInfo>& tileInfos, Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) = 0;
    virtual std::shared_ptr<IProducer2Consumer> toConsumer(const IConsumer2Producer* ic2p) = 0;
};

#endif
