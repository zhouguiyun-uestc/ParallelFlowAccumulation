
#include "communication.h"
#include "mfd_flow_direction.h"
#include "mfd_serial_accum.h"
#include "object_factory.h"
#include "perimeters.h"

#include <paradem/Timer.h>
#include <paradem/gdal.h>
#include <paradem/grid_info.h>
#include <paradem/tile_info.h>
#include <paradem/tool.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

bool processTileGrid( GridInfo& gridDEMInfo, std::vector< TileInfo >& tileDEMInfos, GridInfo& gridDirInfo, std::vector< TileInfo >& tileDirInfos, ObjectFactory* pIObjFactory,
                      const std::string filename ) {
    Grid< std::shared_ptr< IConsumer2Producer > > gridIConsumer2Producer;
    gridIConsumer2Producer.init( gridDirInfo.gridHeight, gridDirInfo.gridWidth );
    std::cout << "tile number= " << tileDirInfos.size() << std::endl;
    for ( int i = 0; i < tileDirInfos.size(); i++ ) {
        std::cout << "Consumer process job1: " << i << std::endl;
        TileInfo& tileDirInfo                                     = tileDirInfos[ i ];
        TileInfo& tileDEMInfo                                     = tileDEMInfos[ i ];
        std::shared_ptr< IConsumer2Producer > pIConsumer2Producer = pIObjFactory->createConsumer2Producer();
        std::shared_ptr< IConsumer > pIConsumer                   = pIObjFactory->createConsumer();
        pIConsumer->processRound1( gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, pIConsumer2Producer.get() );
        gridIConsumer2Producer.at( tileDirInfo.gridRow, tileDirInfo.gridCol ) = pIConsumer2Producer;
    }
    std::cout << "Producer process global solution. " << std::endl;
    std::shared_ptr< IProducer > pIProducer = pIObjFactory->createProducer();
    pIProducer->process( gridDirInfo, tileDirInfos, gridIConsumer2Producer );
    for ( int i = 0; i < tileDirInfos.size(); i++ ) {
        std::cout << "Consumer process job2: " << i << std::endl;
        TileInfo& tileDirInfo                     = tileDirInfos[ i ];
        TileInfo& tileDEMInfo                     = tileDEMInfos[ i ];
        std::shared_ptr< IProducer2Consumer > p2c = pIProducer->toConsumer( gridIConsumer2Producer.at( tileDirInfo.gridRow, tileDirInfo.gridCol ).get() );
        std::shared_ptr< IConsumer > pIConsumer   = pIObjFactory->createConsumer();
        pIConsumer->processRound2( gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, p2c.get() );  //三种方法
    }

    std::string txtPath = gridDirInfo.outputFolder + "\\" + "gridInfo.txt";
    std::ofstream fout;
    fout.open( txtPath, std::ofstream::app );
    if ( fout.fail() ) {
        std::cout << "Open " << txtPath << " error!" << std::endl;
        return false;
    }
    fout << gridDirInfo.tileHeight << std::endl;
    fout << gridDirInfo.tileWidth << std::endl;
    fout << gridDirInfo.gridHeight << std::endl;
    fout << gridDirInfo.gridWidth << std::endl;
    fout << gridDirInfo.grandHeight << std::endl;
    fout << gridDirInfo.grandWidth << std::endl;
    fout << gridDirInfo.cellSize << std::endl;
    fout.close();

    return true;
}

void readTXTInfo( const std::string path, const std::string savePath, GridInfo& gridInfo, std::vector< TileInfo >& tileInfos ) {
    std::string txt = path + "\\gridInfo.txt";
    if ( !readGridInfo( txt.data(), gridInfo ) )
        return;
    gridInfo.inputFolder  = path;
    gridInfo.outputFolder = savePath;
    createTileInfoArray( gridInfo, tileInfos );
}

bool check( const std::string inputDEMFile, const std::string inputDirFile ) {
    Raster< double > accumSerial;
    if ( !readTif( inputDEMFile.data(), accumSerial ) ) {
        std::cout << "Error occured when reading accumSerial!" << std::endl;
        return false;
    }
    int width  = accumSerial.getWidth();
    int height = accumSerial.getHeight();

    Raster< double > accumParallel;
    if ( !readTif( inputDirFile.data(), accumParallel ) ) {
        std::cout << "Error occured when reading accumParallel!" << std::endl;
        return false;
    }

    for ( int row = 0; row < height; row++ ) {
        for ( int col = 0; col < width; col++ ) {
            if ( accumParallel.isNoData( row, col ) )
                continue;
            double accP = accumParallel.at( row, col );
            double accS = accumSerial.at( row, col );
            if ( std::abs( accP - accS ) / accP > 0.0001 )
                return false;
        }
    }
    return true;
}

#pragma region Pararllel
int main( int argc, char** argv ) {
    int rank, size;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    if ( rank == 0 ) {
        std::string retention = argv[ 1 ];
        std::string pathDEM   = argv[ 2 ];
        std::string pathDir   = argv[ 3 ];
        std::string filename  = argv[ 4 ];
        std::string savePath  = argv[ 5 ];
        std::cerr << "c Processes = " << size << std::endl;
        std::cerr << "c Input DEM file = " << pathDEM << std::endl;
        std::cerr << "c Input Dir file = " << pathDir << std::endl;
        std::cerr << "c Input file name = " << filename << std::endl;
        std::cerr << "c Retention strategy = " << retention << std::endl;
        std::cerr << "c Save path = " << savePath << std::endl;

        Timer timer_master;
        timer_master.start();

        Timer timer_overall;
        timer_overall.start();

        GridInfo gridDEMInfo;
        std::vector< TileInfo > tileDEMInfos;
        readTXTInfo( pathDEM, savePath, gridDEMInfo, tileDEMInfos );

        GridInfo gridDirInfo;
        std::vector< TileInfo > tileDirInfos;
        readTXTInfo( pathDir, savePath, gridDirInfo, tileDirInfos );

        std::cout << "c Total rows and cols in the grid = " << gridDEMInfo.grandHeight << "*" << gridDEMInfo.grandWidth << std::endl;
        std::cout << "c Grid has been divided into " << gridDEMInfo.gridHeight << "*" << gridDEMInfo.gridWidth << " tiles." << std::endl;
        std::cout << "c Tile dimensions = " << gridDEMInfo.tileHeight << "*" << gridDEMInfo.tileWidth << std::endl;

        timer_overall.stop();
        std::cerr << "t Preparer time = " << timer_overall.elapsed() << "s" << std::endl;

        ObjectFactory objectFactory;
        Preparer( gridDEMInfo, tileDEMInfos, gridDirInfo, tileDirInfos, &objectFactory, filename );

        timer_master.stop();
        std::cerr << "t Total wall-time=" << timer_master.elapsed() << "s" << std::endl;
    }
    else {
        int good_to_go;
        CommBroadcast( &good_to_go, 0 );
        if ( good_to_go ) {

            ObjectFactory pIObjFactory;
            GridInfo gridDEMInfo, gridDirInfo;
            while ( true ) {
                int the_job = CommGetTag( 0 );

                if ( the_job == TagFirst ) {
                    Timer timer_overall;
                    timer_overall.start();
                    TileInfo tileDEMInfo, tileDirInfo;
                    std::string filename;
                    CommRecv( &gridDEMInfo, &gridDirInfo, &tileDEMInfo, &tileDirInfo, &filename, 0 );
                    std::shared_ptr< IConsumer2Producer > pIConsumer2Producer = pIObjFactory.createConsumer2Producer();
                    std::shared_ptr< IConsumer > pIConsumer                   = pIObjFactory.createConsumer();
                    pIConsumer->processRound1( gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, pIConsumer2Producer.get() );
                    timer_overall.stop();
                    long vmpeak, vmhwm;
                    ProcessMemUsage( vmpeak, vmhwm );
                    Consumer2Producer* p = ( Consumer2Producer* )pIConsumer2Producer.get();
                    p->time_info = TimeInfo( ( ( Consumer* )pIConsumer.get() )->timer_calc.elapsed(), timer_overall.elapsed(), ( ( Consumer* )pIConsumer.get() )->timer_io.elapsed(), vmpeak, vmhwm );
                    CommSend( p, nullptr, 0, ObjectFirst );
                }
                else if ( the_job == TagSecond ) {
                    Timer timer_overall;
                    timer_overall.start();
                    TileInfo tileDEMInfo, tileDirInfo;
                    Producer2Consumer p2c;
                    std::string filename;
                    CommRecv( &gridDEMInfo, &gridDirInfo, &tileDEMInfo, &tileDirInfo, &filename, &p2c, 0 );
                    std::shared_ptr< IConsumer > pIConsumer = pIObjFactory.createConsumer();
                    pIConsumer->processRound2( gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, &p2c );
                    timer_overall.stop();
                    long vmpeak, vmhwm;
                    ProcessMemUsage( vmpeak, vmhwm );
                    TimeInfo temp( ( ( Consumer* )pIConsumer.get() )->timer_calc.elapsed(), timer_overall.elapsed(), ( ( Consumer* )pIConsumer.get() )->timer_io.elapsed(), vmpeak, vmhwm );
                    CommSend( &temp, nullptr, 0, ObjectSecond );
                }
                else if ( the_job == TagNull )
                    break;
            }
        }
    }
    MPI_Finalize();
    return 0;
}
#pragma endregion

std::string itostring( int x ) {
    if ( x == 1 )
        return "1";
    if ( x == 2 )
        return "2";
    if ( x == 3 )
        return "3";
    if ( x == 4 )
        return "4";
    if ( x == 5 )
        return "5";
    if ( x == 6 )
        return "2";
    if ( x == 7 )
        return "3";
    if ( x == 8 )
        return "4";
    if ( x == 9 )
        return "5";
    if ( x == 10 )
        return "10";
    if ( x == 11 )
        return "11";
    if ( x == 12 )
        return "12";
    if ( x == 13 )
        return "13";
    if ( x == 14 )
        return "14";
    if ( x == 15 )
        return "15";
    if ( x == 16 )
        return "16";
    if ( x == 17 )
        return "17";
    if ( x == 18 )
        return "18";
    if ( x == 19 )
        return "19";
    if ( x == 20 )
        return "20";
}

//#pragma region serial
//
// int main(int argc, char **argv) {
//
//	//std::string doJob="creatDEM";
//	//std::string doJob = "preprocess";
//	//std::string doJob = "big data preprocess";
//	std::string doJob = "toAccum";
//	//std::string doJob = "merge";
//	//std::string doJob = "merge dem";
//	//std::string doJob = "serial compute accum";
//	//std::string doJob = "check";
//	//std::string doJob = "test";
//	if (doJob == "creatDEM") {
//		std::string path = "E:\\dem\\grid0\\dem";
//		createNewTif(path);
//	}
//	else if (doJob == "preprocess")
//	{
//		for (int i = 11; i < 31; i++)
//		{
//			std::string path = "F:\\dem\\grid01";
//			std::string pathDEM = path + "\\dem\\dem.tif";
//			std::string pathDir = path + "\\dir\\flowDir.tif";
//
//			int tileHieght = 100 * i, tileWidth = 100 * i;
//			std::string ouputPath = "F:\\dem\\grid0";
//			std::string outputDEMtiles = ouputPath + "\\" + itostring(i) + "\\dem";
//			std::string filename = "freeborn";
//			std::string type = "dem";
//			generateTiles(pathDEM.data(), tileHieght, tileWidth, outputDEMtiles.data(), filename, type);
//
//			std::string outputDirtiles = ouputPath + "\\" + itostring(i) + "\\dir";
//			type = "dir";
//			generateTiles(pathDir.data(), tileHieght, tileWidth, outputDirtiles.data(), filename, type);
//
//		}
//		//std::string path = "F:\\dem\\grid01";
//		//std::string pathDEM = path + "\\dem\\dem.tif";
//		//std::string pathDir = path + "\\dir\\flowDir.tif";
//		//
//		//int tileHieght = 400, tileWidth = 400;
//		//std::string ouputPath = "F:\\dem\\grid0";
//		//std::string outputDEMtiles = ouputPath + "\\04\\dem";
//		//std::string filename = "freeborn";
//		//std::string type = "dem";
//		//generateTiles(pathDEM.data(), tileHieght, tileWidth, outputDEMtiles.data(), filename,type);
//
//		//std::string outputDirtiles = ouputPath + "\\04\\dir";
//		//type = "dir";
//		//generateTiles(pathDir.data(), tileHieght, tileWidth, outputDirtiles.data(),filename,type);
//
//	}
//	else if (doJob == "big data preprocess") {
//		std::string pathDir = "F:\\dem\\anoka\\anokaflow";
//		std::string filenameDir = "anokafilling";
//		std::string pathDEM = "F:\\dem\\anoka\\anoka";
//		std::string filenameDEM = "anokafilling";
//		std::string saveDirPath = "F:\\dem\\anoka\\dir";
//		std::string saveDEMPath = "F:\\dem\\anoka\\dem";
//		std::string savefilename = "anoka";
//		preprocess(pathDir, filenameDir, pathDEM, filenameDEM, saveDirPath, saveDEMPath, savefilename);
//		if (!checkDir(saveDirPath, savefilename, saveDEMPath))
//			std::cout << "dir is error!" << std::endl;
//		else
//			std::cout << "dir is ok!" << std::endl;
//
//	}
//	else if (doJob == "toAccum")
//	{
//		GridInfo gridDEMInfo;
//		std::string pathDEM = "F:\\dem\\dodge\\weak_scaling\\weak_file\\dem\\gridInfo.txt";
//		if (!readGridInfo(pathDEM.data(), gridDEMInfo)) return 0;
//		gridDEMInfo.inputFolder = "F:\\dem\\dodge\\weak_scaling\\weak_file\\dem";
//		std::vector<TileInfo> tileDEMInfos;
//		createTileInfoArray(gridDEMInfo, tileDEMInfos);
//		GridInfo gridDirInfo;
//		std::string pathDir = "F:\\dem\\dodge\\weak_scaling\\weak_file\\dir\\gridInfo.txt";
//		if (!readGridInfo(pathDir.data(), gridDirInfo)) return 0;
//		gridDirInfo.inputFolder = "F:\\dem\\dodge\\weak_scaling\\weak_file\\dir";
//		std::vector<TileInfo> tileDirInfos;
//		createTileInfoArray(gridDirInfo, tileDirInfos);
//
//		gridDEMInfo.outputFolder = "F:\\dem\\grid0\\accum";
//		gridDirInfo.outputFolder = "F:\\dem\\grid0\\accum";
//
//		std::string filename = "Ddodge-1m-filling";
//		ObjectFactory objectFactory;
//		processTileGrid(gridDEMInfo, tileDEMInfos, gridDirInfo, tileDirInfos, &objectFactory, filename);
//		/*std::string pathAccum = "F:\\dem\\grid0\\accum\\gridInfo.txt";
//		GridInfo gridAccumInfo;
//		readGridInfo(pathAccum.data(), gridAccumInfo);
//		gridAccumInfo.inputFolder = "F:\\dem\\grid0\\accum";
//		gridAccumInfo.outputFolder = "F:\\dem\\grid0\\accum";
//		mergeTiles(gridAccumInfo);*/
//	}
//	else if (doJob == "merge dem") {
//		std::string pathDEM = "F:\\dem\\aitkin\\dem\\gridInfo.txt";
//		GridInfo gridDEMInfo;
//		readGridInfo(pathDEM.data(), gridDEMInfo);
//		gridDEMInfo.inputFolder = "F:\\dem\\aitkin\\dem";
//		gridDEMInfo.outputFolder = "F:\\dem\\aitkin";
//		mergeDEM(gridDEMInfo);
//	}
//	else if (doJob == "serial compute accum") {
//		Timer timer_master;
//		timer_master.start();
//		/*std::string inputDEMFile = argv[1];
//		std::string inputDirFile = argv[2];
//		std::string outputAccuFile = argv[3];*/
//
//		std::string inputDEMFile = "F:\\dem\\aitkin\\aitkin_mergeDEM_song.tif";
//		std::string inputDirFile = "F:\\dem\\aitkin\\aitkin-3mmergeDIR_song.tif";
//		std::string outputAccuFile = "F:\\dem\\aitkin\\aitkiaccum_Improved_0918.tif";
//		if (!computeAccum(inputDEMFile, inputDirFile, outputAccuFile))
//			return 0;
//		timer_master.stop();
//		std::cout << timer_master.elapsed() << std::endl;
//
//	}
//	else if (doJob == "check") {
//		std::string inputSeFile = "F:\\dem\\grid0\\accum\\accum_recur.tif";
//		std::string inputParaFile = "F:\\dem\\grid0\\accum\\merge.tif";
//		check(inputSeFile, inputParaFile);
//	}
//	else if (doJob == "test")
//	{
//		std::string path = "F:\\dem\\grid00";
//		for (int i = 1; i <= 20; i++)
//		{
//			for (int j = 1; j < 4; j++)
//			{
//				std::string inputDEMFile = path + "\\" + itostring(i) + "\\dem\\"+itostring(j)+".tif";
//				std::string inputDirFile =path + "\\" + itostring(i) + "\\dir\\" + itostring(j) + ".tif";
//				std::string outputAccuFile = path + "\\" + itostring(i) + "\\accum\\" + itostring(j) + ".tif";
//				if (!computeAccum(inputDEMFile, inputDirFile, outputAccuFile))
//					return 0;
//			}
//		}
//	}
//	return 0;
//}
//
//#pragma endregion