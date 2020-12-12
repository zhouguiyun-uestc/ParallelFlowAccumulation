#ifndef PARADEM_PRODUCER2CONSUMER_H
#define PARADEM_PRODUCER2CONSUMER_H

#include <cereal/archives/binary.hpp>
#include <paradem/i_producer_2_consumer.h>

#include <vector>

typedef double accum_t;

class Producer2Consumer : public IProducer2Consumer {
private:
    friend class cereal::access;
    template <class Archive> void serialize(Archive& ar) {
        ar(inner_accum_setoff1, inner_accum_setoff2, outer_accum);
    }

public:
    Producer2Consumer() = default;
    virtual ~Producer2Consumer() = default;

public:
    virtual void free();

public:
    std::vector<accum_t> inner_accum_setoff1;
    std::vector<accum_t> inner_accum_setoff2;
    std::vector<accum_t> outer_accum;
};

#endif
