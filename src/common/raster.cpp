#include <paradem/raster.h>

template < class T > Raster< T >::Raster() {

    NoDataValue = -9999;
}

template < class T > Raster< T >::~Raster() {}

template < class T > bool Raster< T >::isNoData( int row, int col ) {
    if ( fabs( at( row, col ) - NoDataValue ) < 0.00001 )
        return true;
    return false;
}

template < class T > int Raster< T >::getRow( int dir, int row ) {
    return row + dRow[ dir ];
}

template < class T > int Raster< T >::getCol( int dir, int col ) {
    return ( col + dCol[ dir ] );
}

template < class T > bool Raster< T >::isInGrid( int row, int col ) {
    if ( ( row >= 0 && row < height ) && ( col >= 0 && col < width ) )
        return true;

    return false;
}
template < class T > T Raster< T >::getNoDataValue() {
    return NoDataValue;
}

template < class T > double* Raster< T >::getGeoTransformsPtr() {
    if ( geoTransforms == nullptr )
        return nullptr;
    return &geoTransforms->at( 0 );
}

template < class T > std::vector< T > Raster< T >::getRowData( int row ) {
    return std::vector< T >( &data[ 0 ] + row * width, &data[ 0 ] + ( row + 1 ) * width );
}

template < class T > std::vector< T > Raster< T >::getColData( int col ) {
    std::vector< T > vec( height );
    for ( int row = 0; row < height; row++ ) {
        vec[ row ] = at( row, col );
    }
    return vec;
}

template class Raster< float >;
template class Raster< int32_t >;
template class Raster< double >;
template class Raster< uint8_t >;