#ifndef PARADEM_PRODUCER_H
#define PARADEM_PRODUCER_H

#include "atype.h"
#include "consumer_2_producer.h"
#include "producer_2_consumer.h"

#include <paradem/grid.h>
#include <paradem/grid_info.h>
#include <paradem/i_producer.h>
#include <paradem/tile_info.h>
#include <paradem/timer.h>

#include <memory>
class Producer : public IProducer {
public:
    Timer timer_io, timer_calc;

public:
    virtual void process( const GridInfo& gridInfo, const std::vector< TileInfo >& tileInfos, Grid< std::shared_ptr< IConsumer2Producer > >& gridIConsumer2Producer );
    virtual std::shared_ptr< IProducer2Consumer > toConsumer( const IConsumer2Producer* ic2p );
    virtual void free();

public:
    void globalAccumOffset_tradition( const GridInfo& gridInfo, const std::vector< TileInfo >& tileInfos, Grid< std::shared_ptr< IConsumer2Producer > >& gridIConsumer2Producer );
    void globalAccumOffset_improve( const GridInfo& gridInfo, const std::vector< TileInfo >& tileInfos, Grid< std::shared_ptr< IConsumer2Producer > >& gridIConsumer2Producer );
};

#endif
