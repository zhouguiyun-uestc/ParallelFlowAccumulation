#include "mfd_serial_accum.h"

void traditon( Raster< float >& dem, Raster< uint8_t > flowdDir, Raster< double >& accum, const double& cellSize, const std::vector< std::vector< double > >& index ) {
    int height = dem.getHeight();
    int width  = dem.getWidth();

    std::cout << "Using tradition algorithm to compute the flow accumulation matrix ... " << std::endl;
    auto start = std::chrono::system_clock::now();

    std::cout << "Computing dependencies.." << std::endl;
    std::vector< std::vector< int > > dependencies( height, std::vector< int >( width, 0 ) );
    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( dem.isNoData( row, col ) )
                continue;
            accum.at( row, col ) = 1;
            if ( flowdDir.isNoData( row, col ) )
                continue;
            uint8_t dirs = flowdDir.at( row, col );
            for ( int i = 0; i < 8; i++ ) {
                if ( dirs & ( 1 << i ) ) {
                    int iRow = dem.getRow( i, row );
                    int iCol = dem.getCol( i, col );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                        continue;
                    dependencies[ iRow ][ iCol ]++;
                }
            }
        }
    }

    std::queue< Node > sources;
    for ( int row = 0; row < height; row++ )
        for ( int col = 0; col < width; col++ )
            if ( ( !dem.isNoData( row, col ) ) && dependencies[ row ][ col ] == 0 ) {
                Node node( row, col, dem.at( row, col ) );
                sources.push( node );
            }

    while ( !sources.empty() ) {
        Node node = sources.front();
        sources.pop();
        if ( dem.isNoData( node.row, node.col ) )
            continue;
        std::vector< double > slopes( 8, 0 );
        double fracSum = 0;
        if ( flowdDir.isNoData( node.row, node.col ) )
            continue;
        uint8_t dir = flowdDir.at( node.row, node.col );
        for ( int k = 0; k < 8; k++ ) {
            if ( dir & ( 1 << k ) ) {
                int iRow, iCol;
                iRow = dem.getRow( k, node.row );
                iCol = dem.getCol( k, node.col );
                if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                    continue;
                double dSlope = ( dem.at( node.row, node.col ) - dem.at( iRow, iCol ) ) / cellSize;
                if ( k % 2 == 1 ) {
                    dSlope      = dSlope / SQRT2;
                    slopes[ k ] = pow( dSlope, index[ node.row ][ node.col ] ) * SQRT2 / 4;
                }
                else
                    slopes[ k ] = pow( dSlope, index[ node.row ][ node.col ] ) / 2;
                fracSum += slopes[ k ];
            }
        }

        if ( fracSum == 0 ) {
            int count = 0;
            for ( int k = 0; k < 8; k++ ) {
                if ( dir & ( 1 << k ) )
                    count++;
            }
            double temp = accum.at( node.row, node.col ) / count;
            for ( int k = 0; k < 8; k++ ) {
                if ( dir & ( 1 << k ) ) {
                    int iRow, iCol;
                    iRow = dem.getRow( k, node.row );
                    iCol = dem.getCol( k, node.col );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                        continue;
                    accum.at( iRow, iCol ) += temp;
                    dependencies[ iRow ][ iCol ]--;
                    if ( dependencies[ iRow ][ iCol ] == 0 ) {
                        Node tempNode( iRow, iCol, dem.at( iRow, iCol ) );
                        sources.push( tempNode );
                    }
                }
            }
        }
        else {
            for ( int k = 0; k < 8; k++ ) {
                if ( slopes[ k ] != 0 ) {
                    int iRow, iCol;
                    iRow        = dem.getRow( k, node.row );
                    iCol        = dem.getCol( k, node.col );
                    double temp = slopes[ k ] / fracSum;
                    accum.at( iRow, iCol ) += accum.at( node.row, node.col ) * temp;
                    dependencies[ iRow ][ iCol ]--;
                    if ( dependencies[ iRow ][ iCol ] == 0 ) {
                        Node tempNode( iRow, iCol, dem.at( iRow, iCol ) );
                        sources.push( tempNode );
                    }
                }
            }
        }
    }

    auto end      = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast< std::chrono::microseconds >( end - start );
    std::cout << "it takes£º" << double( duration.count() ) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "seconds" << std::endl;
}

extern int BitCount( uint8_t n, int& k ) {
    int cnt = 0;
    for ( ; n; n >>= 1 ) {
        cnt += n & 1;
        k++;
    }
    return cnt;
}
void improve( Raster< float > dem, Raster< uint8_t > flowdDir, Raster< double >& accum, const double& cellSize, const std::vector< std::vector< double > >& index ) {
    int height = dem.getHeight();
    int width  = dem.getWidth();
    int iRow, iCol;
    double dSlope;
    std::cout << "Using improved algorithm to compute the flow accumulation matrix ... " << std::endl;
    auto start = std::chrono::system_clock::now();

    std::cout << "Computing dependencies.." << std::endl;
    std::vector< std::vector< int > > dependencies( height, std::vector< int >( width, 0 ) );
    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( dem.isNoData( row, col ) )
                continue;
            accum.at( row, col ) = 1;
            if ( flowdDir.isNoData( row, col ) )
                continue;
            uint8_t dirs = flowdDir.at( row, col );
            for ( int i = 0; i < 8; i++ ) {
                if ( dirs & ( 1 << i ) ) {
                    int iRow = dem.getRow( i, row );
                    int iCol = dem.getCol( i, col );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                        continue;
                    dependencies[ iRow ][ iCol ]++;
                }
            }
        }
    }

    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( !dem.isNoData( row, col ) && dependencies[ row ][ col ] == 0 ) {
                std::queue< Node > sources;
                uint8_t dir0 = flowdDir.at( row, col );
                dependencies[ row ][ col ]--;
                int tempRow0   = row;
                int tempCol0   = col;
                int k0         = -1;
                bool visited0  = false;
                bool isPushed0 = true;
                while ( BitCount( dir0, k0 ) == 1 ) {
                    visited0 = true;
                    int iRow = dem.getRow( k0, tempRow0 );
                    int iCol = dem.getCol( k0, tempCol0 );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) ) {
                        isPushed0 = false;
                        break;
                    }
                    accum.at( iRow, iCol ) += accum.at( tempRow0, tempCol0 );
                    dependencies[ iRow ][ iCol ]--;
                    if ( dependencies[ iRow ][ iCol ] == 0 ) {
                        dependencies[ iRow ][ iCol ]--;
                        dir0     = flowdDir.at( iRow, iCol );
                        tempRow0 = iRow;
                        tempCol0 = iCol;
                        k0       = -1;
                    }
                    else {
                        isPushed0 = false;
                        break;
                    }
                }
                if ( visited0 ) {
                    if ( isPushed0 ) {
                        Node tempNode( tempRow0, tempCol0, dem.at( tempRow0, tempCol0 ) );
                        sources.push( tempNode );
                    }
                }
                else {
                    Node node( tempRow0, tempCol0, dem.at( tempRow0, tempCol0 ) );
                    sources.push( node );
                }
                while ( !sources.empty() ) {
                    Node tempNode = sources.front();
                    sources.pop();
                    int tempRow = tempNode.row;
                    int tempCol = tempNode.col;
                    dependencies[ tempRow ][ tempCol ]--;
                    if ( dem.isNoData( tempRow, tempCol ) )
                        continue;
                    std::vector< double > slopes( 8, 0 );
                    double fracSum = 0;
                    uint8_t dir    = flowdDir.at( tempRow, tempCol );
                    int k          = -1;
                    bool visited   = false;
                    bool isPushed  = true;
                    while ( BitCount( dir, k ) == 1 ) {
                        visited  = true;
                        int iRow = dem.getRow( k, tempRow );
                        int iCol = dem.getCol( k, tempCol );
                        if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) ) {
                            isPushed = false;
                            break;
                        }
                        accum.at( iRow, iCol ) += accum.at( tempRow, tempCol );
                        dependencies[ iRow ][ iCol ]--;
                        if ( dependencies[ iRow ][ iCol ] == 0 ) {
                            dependencies[ iRow ][ iCol ]--;
                            dir     = flowdDir.at( iRow, iCol );
                            tempRow = iRow;
                            tempCol = iCol;
                            k       = -1;
                        }
                        else {
                            isPushed = false;
                            break;
                        }
                    }
                    if ( visited ) {
                        if ( isPushed ) {
                            Node tempNode( tempRow, tempCol, dem.at( tempRow, tempCol ) );
                            sources.push( tempNode );
                        }
                        continue;
                    }
                    for ( int k = 0; k < 8; k++ ) {

                        if ( dir & ( 1 << k ) ) {
                            int iRow, iCol;
                            iRow = dem.getRow( k, tempRow );
                            iCol = dem.getCol( k, tempCol );
                            if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                                continue;
                            dSlope = ( dem.at( tempRow, tempCol ) - dem.at( iRow, iCol ) ) / cellSize;
                            if ( k % 2 == 1 ) {
                                dSlope      = dSlope / SQRT2;
                                slopes[ k ] = pow( dSlope, index[ tempRow ][ tempCol ] ) * SQRT2 / 4;
                            }
                            else {
                                slopes[ k ] = pow( dSlope, index[ tempRow ][ tempCol ] ) / 2;
                            }
                            fracSum += slopes[ k ];
                        }
                    }
                    if ( fracSum == 0 ) {
                        int count = 0;
                        for ( int k = 0; k < 8; k++ ) {
                            if ( dir & ( 1 << k ) ) {
                                count++;
                            }
                        }
                        if ( accum.at( tempRow, tempCol ) < 0 ) {
                            int less0 = 0;
                        }
                        double temp = accum.at( tempRow, tempCol ) / count;
                        for ( int k = 0; k < 8; k++ ) {
                            if ( dir & ( 1 << k ) ) {
                                int iRow, iCol;
                                iRow = dem.getRow( k, tempRow );
                                iCol = dem.getCol( k, tempCol );
                                if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                                    continue;
                                accum.at( iRow, iCol ) += temp;
                                dependencies[ iRow ][ iCol ]--;
                                if ( dependencies[ iRow ][ iCol ] == 0 ) {
                                    // dependencies[iRow][iCol]--;
                                    Node tempNode( iRow, iCol, dem.at( iRow, iCol ) );
                                    sources.push( tempNode );
                                }
                            }
                        }
                    }
                    else {
                        for ( int k = 0; k < 8; k++ ) {
                            if ( slopes[ k ] == 0 )
                                continue;
                            int iRow, iCol;
                            iRow        = dem.getRow( k, tempRow );
                            iCol        = dem.getCol( k, tempCol );
                            double temp = slopes[ k ] / fracSum;
                            accum.at( iRow, iCol ) += accum.at( tempRow, tempCol ) * temp;
                            dependencies[ iRow ][ iCol ]--;
                            if ( dependencies[ iRow ][ iCol ] == 0 ) {
                                Node tempNode( iRow, iCol, dem.at( iRow, iCol ) );
                                sources.push( tempNode );
                            }
                        }
                    }
                }
            }
        }
    }

    auto end      = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast< std::chrono::microseconds >( end - start );
    std::cout << "it takes£º" << double( duration.count() ) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "seconds" << std::endl;
}

double checkNeighbor( int row, int col, Raster< float >& dem, Raster< uint8_t >& flowdirs, const double& cellSize, const std::vector< std::vector< double > >& index,
                      const std::vector< std::vector< double > >& fractionSum, Raster< double >& accum ) {
    if ( dem.isNoData( row, col ) )
        return 0;
    if ( accum.isNoData( row, col ) ) {
        accum.at( row, col ) = 1;
        for ( int k = 0; k < 8; k++ ) {
            int iRow, iCol;
            iRow = dem.getRow( k, row );
            iCol = dem.getCol( k, col );
            if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                continue;
            double temp = 0, dSlope = 0;
            if ( is2RowCol( row, col, k, iRow, iCol, flowdirs.at( iRow, iCol ) ) ) {
                if ( fractionSum[ iRow ][ iCol ] == MAX_SINGLE )  //Æ½µØ
                {
                    uint8_t dir = flowdirs.at( iRow, iCol );
                    int count   = 0;
                    for ( int k = 0; k < 8; k++ ) {
                        if ( dir & ( 1 << k ) ) {
                            count++;
                        }
                    }
                    temp = 1.0 / count;
                }
                else {
                    dSlope = ( dem.at( iRow, iCol ) - dem.at( row, col ) ) / cellSize;
                    if ( k % 2 == 1 ) {
                        dSlope = dSlope / SQRT2;
                        temp   = pow( dSlope, index[ iRow ][ iCol ] ) * SQRT2 / ( 4 * fractionSum[ iRow ][ iCol ] );
                    }
                    else {
                        temp = pow( dSlope, index[ iRow ][ iCol ] ) / ( 2 * fractionSum[ iRow ][ iCol ] );
                    }
                }

                if ( temp > 0 ) {
                    double value         = accum.at( row, col ) + temp * checkNeighbor( iRow, iCol, dem, flowdirs, cellSize, index, fractionSum, accum );
                    accum.at( row, col ) = value;
                }
            }
        }
    }
    return accum.at( row, col );
}

void recursive( Raster< float > dem, Raster< uint8_t > flowdDir, const double& cellSize, Raster< double >& accum, const std::vector< std::vector< double > >& index ) {

    int height = dem.getHeight();
    int width  = dem.getWidth();

    double dSlope;
    double a        = 8.9;
    double b        = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;

    std::cout << "Computing accuGrid by Recursive method..." << std::endl;
    auto start = std::chrono::system_clock::now();

    std::vector< std::vector< double > > fractionSum( height, std::vector< double >( width, 0 ) );

    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( dem.isNoData( row, col ) || flowdDir.isNoData( row, col ) )
                continue;
            double dSum = 0;
            uint8_t dir = flowdDir.at( row, col );
            for ( int k = 0; k < 8; k++ ) {
                if ( dir & ( 1 << k ) ) {
                    int iRow, iCol;
                    iRow = dem.getRow( k, row );
                    iCol = dem.getCol( k, col );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                        continue;
                    dSlope = ( dem.at( row, col ) - dem.at( iRow, iCol ) ) / cellSize;
                    if ( k % 2 == 1 ) {
                        dSlope = dSlope / SQRT2;
                        dSum   = dSum + pow( dSlope, index[ row ][ col ] ) * SQRT2 / 4;
                    }
                    else {
                        dSum = dSum + pow( dSlope, index[ row ][ col ] ) / 2;
                    }
                }
            }
            fractionSum[ row ][ col ] = dSum == 0 ? MAX_SINGLE : dSum;  //¸Ä£¿£¿£¿£¿£¿£¿
        }
    }

    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( !dem.isNoData( row, col ) )
                checkNeighbor( row, col, dem, flowdDir, cellSize, index, fractionSum, accum );
        }
    }

    auto end      = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast< std::chrono::microseconds >( end - start );
    std::cout << "it takes£º" << double( duration.count() ) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "seconds" << std::endl;
}

bool computeAccum( std::string inputDEMFile, std::string inputDirFile, std::string outputAccuFile ) {
    Raster< float > dem;
    if ( !readTif( inputDEMFile.data(), dem ) ) {
        std::cout << "Error occured when reading DEM!" << std::endl;
        return false;
    }
    int width  = dem.getWidth();
    int height = dem.getHeight();

    Raster< uint8_t > dir;
    if ( !readTif( inputDirFile.data(), dir ) ) {
        std::cout << "Error occured when reading Dir!" << std::endl;
        return false;
    }

    double a        = 8.9;
    double b        = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    double cellSize = abs( dem.geoTransforms->at( 1 ) );
    std::cout << "Computing index..." << std::endl;
    std::vector< std::vector< double > > index( height, std::vector< double >( width, 0 ) );
    int iRow, iCol;
    double dSlope;
    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( dem.isNoData( row, col ) )
                continue;
            double dMax  = 0;
            uint8_t dirs = dir.at( row, col );
            for ( int i = 0; i < 8; i++ ) {
                if ( dirs & ( 1 << i ) ) {
                    iRow = dem.getRow( i, row );
                    iCol = dem.getCol( i, col );
                    if ( !dem.isInGrid( iRow, iCol ) || dem.isNoData( iRow, iCol ) )
                        continue;
                    dSlope = ( dem.at( row, col ) - dem.at( iRow, iCol ) ) / cellSize;
                    if ( i % 2 == 1 )
                        dSlope = dSlope / SQRT2;
                    if ( dSlope > dMax )
                        dMax = dSlope;
                }
            }
            if ( dMax >= slopeMax )
                index[ row ][ col ] = a + b;
            else if ( dMax <= slopeMin )
                index[ row ][ col ] = b;
            else
                index[ row ][ col ] = a * dMax + b;
        }
    }

    Raster< double > accum;
    accum.init( dem.getHeight(), dem.getWidth() );
    accum.setAllValues( accum.NoDataValue );
    // traditon(dem, dir, accum, cellSize, index);
    improve( dem, dir, accum, cellSize, index );
    // recursive(dem, dir, cellSize, accum, index);
    double min, max, mean, stdDev;
    Consumer con;
    con.calculateStatistics( accum, &min, &max, &mean, &stdDev );
    WriteGeoTIFF( outputAccuFile.data(), accum.getHeight(), accum.getWidth(), &accum, GDALDataType::GDT_Float64, dem.getGeoTransformsPtr(), &min, &max, &mean, &stdDev, accum.NoDataValue );
    return true;
}
