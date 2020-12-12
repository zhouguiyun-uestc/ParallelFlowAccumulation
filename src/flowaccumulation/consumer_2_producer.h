#ifndef PARADEM_CONSUMER2PRODUCER_H
#define PARADEM_CONSUMER2PRODUCER_H

#include "Node.h"
#include "atype.h"
#include "globalPoint.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <paradem/i_consumer_2_producer.h>
#include <paradem/tile_info.h>
#include <paradem/timeInfo.h>
#include <paradem/timer.h>

#include <iterator>
#include <map>
#include <stdint.h>
#include <vector>

typedef double accum_t;
typedef float elev_t;
typedef uint8_t flowdir_t;
typedef int link_t;
typedef int dependency_t;
typedef double fraction_t;

class Consumer2Producer : public IConsumer2Producer {
private:
    friend class cereal::access;
    template <class Archive> void serialize(Archive& ar) {
        ar(elevation, flowdirs, links, fraction, inner_accum, gridRow, gridCol, time_info);
    }

public:
    std::vector<elev_t> elevation;
    std::vector<flowdir_t> flowdirs;
    std::vector<std::vector<link_t>> links;
    std::vector<std::vector<fraction_t>> fraction;

    std::vector<dependency_t> dependencies;

public:
    std::vector<accum_t> inner_accum;

    std::vector<accum_t> inner_accum_setoff1;

    std::vector<accum_t> inner_accum_setoff2;

public:
    std::vector<std::vector<atype>> outerGlobal_links;
    std::vector<std::vector<fraction_t>> outerGlobal_fraction;
    std::vector<accum_t> outer_accum;

public:
    int gridRow, gridCol;
    TimeInfo time_info;

public:
    Consumer2Producer() = default;
    virtual ~Consumer2Producer() = default;

public:
    virtual void free();
};

#endif
