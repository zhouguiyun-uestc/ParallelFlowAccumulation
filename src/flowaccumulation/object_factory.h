#ifndef PARADEM_OBJECT_FACTORY_H
#define PARADEM_OBJECT_FACTORY_H

#include "Consumer.h"
#include "Producer.h"
#include "consumer_2_producer.h"
#include "producer_2_consumer.h"

#include <memory>

class ObjectFactory {
public:
    virtual std::shared_ptr<IConsumer> createConsumer();
    virtual std::shared_ptr<IProducer> createProducer();
    virtual std::shared_ptr<IConsumer2Producer> createConsumer2Producer();
};

#endif
