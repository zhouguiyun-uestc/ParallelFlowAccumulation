#include <paradem/grid.h>
#include <paradem/i_consumer_2_producer.h>
#include <paradem/i_producer_2_consumer.h>

#include <memory>
#include <stdint.h>

template < class T > Grid< T >::Grid() {}

template < class T > Grid< T >::~Grid() {}

template < class T > bool Grid< T >::init( int height, int width ) {
    this->height = height;
    this->width  = width;
    data.resize( height * width );
    return true;
}

template < class T > void Grid< T >::setAllValues( T val ) {
    for ( uint32_t i = 0; i < data.size(); i++ ) {
        data[ i ] = val;
    }
}

template < class T > int Grid< T >::getHeight() {
    return height;
}

template < class T > int Grid< T >::getWidth() {
    return width;
}

template < class T > T& Grid< T >::at( int row, int col ) {
    return data[ row * width + col ];
}

template < class T > T* Grid< T >::operator&() {
    return &data[ 0 ];
}

template class Grid< float >;
template class Grid< int32_t >;
template class Grid< std::shared_ptr< IConsumer2Producer > >;
template class Grid< double >;
template class Grid< uint8_t >;