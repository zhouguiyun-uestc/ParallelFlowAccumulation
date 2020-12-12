
#include "communication.h"
#include "mfd_flow_direction.h"
#include "mfd_serial_accum.h"
#include "object_factory.h"
#include "perimeters.h"

#include <CLI/CLI.hpp>
#include <paradem/gdal.h>
#include <paradem/grid_info.h>
#include <paradem/tile_info.h>
#include <paradem/timer.h>
#include <paradem/tool.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

bool processTileGrid(GridInfo& gridDEMInfo, std::vector<TileInfo>& tileDEMInfos, GridInfo& gridDirInfo, std::vector<TileInfo>& tileDirInfos, ObjectFactory* pIObjFactory, const std::string filename) {
    Grid<std::shared_ptr<IConsumer2Producer>> gridIConsumer2Producer;
    gridIConsumer2Producer.init(gridDirInfo.gridHeight, gridDirInfo.gridWidth);
    std::cout << "tile number= " << tileDirInfos.size() << std::endl;
    for (size_t i = 0; i < tileDirInfos.size(); i++) {
        std::cout << "Consumer process job1: " << i << std::endl;
        TileInfo& tileDirInfo = tileDirInfos[i];
        TileInfo& tileDEMInfo = tileDEMInfos[i];
        std::shared_ptr<IConsumer2Producer> pIConsumer2Producer = pIObjFactory->createConsumer2Producer();
        std::shared_ptr<IConsumer> pIConsumer = pIObjFactory->createConsumer();
        pIConsumer->processRound1(gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, pIConsumer2Producer.get());
        gridIConsumer2Producer.at(tileDirInfo.gridRow, tileDirInfo.gridCol) = pIConsumer2Producer;
    }
    std::cout << "Producer process global solution. " << std::endl;
    std::shared_ptr<IProducer> pIProducer = pIObjFactory->createProducer();
    pIProducer->process(gridDirInfo, tileDirInfos, gridIConsumer2Producer);
    for (size_t i = 0; i < tileDirInfos.size(); i++) {
        std::cout << "Consumer process job2: " << i << std::endl;
        TileInfo& tileDirInfo = tileDirInfos[i];
        TileInfo& tileDEMInfo = tileDEMInfos[i];
        std::shared_ptr<IProducer2Consumer> p2c = pIProducer->toConsumer(gridIConsumer2Producer.at(tileDirInfo.gridRow, tileDirInfo.gridCol).get());
        std::shared_ptr<IConsumer> pIConsumer = pIObjFactory->createConsumer();
        pIConsumer->processRound2(gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, p2c.get());  //ÈýÖÖ·½·¨
    }

    std::string txtPath = gridDirInfo.outputFolder + "/" + "gridInfo.txt";
    std::ofstream fout;
    fout.open(txtPath, std::ofstream::app);
    if (fout.fail()) {
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

void readTXTInfo(const std::string path, const std::string savePath, GridInfo& gridInfo, std::vector<TileInfo>& tileInfos) {
    std::string txt = path + "/gridInfo.txt";
    if (!readGridInfo(txt.data(), gridInfo))
        return;
    gridInfo.inputFolder = path;
    gridInfo.outputFolder = savePath;
    createTileInfoArray(gridInfo, tileInfos);
}

bool check(const std::string inputDEMFile, const std::string inputDirFile) {
    Raster<double> accumSerial;
    if (!readTif(inputDEMFile.data(), accumSerial)) {
        std::cout << "Error occured when reading accumSerial!" << std::endl;
        return false;
    }
    int width = accumSerial.getWidth();
    int height = accumSerial.getHeight();

    Raster<double> accumParallel;
    if (!readTif(inputDirFile.data(), accumParallel)) {
        std::cout << "Error occured when reading accumParallel!" << std::endl;
        return false;
    }

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            if (accumParallel.isNoData(row, col))
                continue;
            double accP = accumParallel.at(row, col);
            double accS = accumSerial.at(row, col);
            if (std::abs(accP - accS) / accP > 0.0001)
                return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cout << "Parameter error! There are five arguments (method,the folder path of DEMs, the folder path of flow direction file, fileName, outputPath.)" << std::endl;
    }
    char* method = argv[1];
    if (strcmp(method, "parallel") == 0) {
        int rank, size;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if (size < 2) {
            std::cout << "Specify at least 2 processes in the parallel algorithm." << std::endl;
            return -1;
        }
        if (rank == 0) {
            CLI::App app("Parallel-Multi-Flow-Accumulation Program");
            std::string pathDEM, pathDir, prefixName, outputPath;
            app.add_option("method", method, "Parallel model")->required();
            app.add_option("pathDEM", pathDEM, "Path of DEMs folder")->required();
            app.add_option("pathDir", pathDir, "Path of flow-directions folder")->required();
            app.add_option("prefixName", prefixName, "Files with the same prefix")->required();
            app.add_option("outputPath", outputPath, "Path of flow-accumulation output folder")->required();
            CLI11_PARSE(app, argc, argv);

            std::cerr << "Processes = " << size << std::endl;
            std::cerr << "Input DEM file = " << pathDEM << std::endl;
            std::cerr << "Input Dir file = " << pathDir << std::endl;
            std::cerr << "Input prefix name of file = " << prefixName << std::endl;
            std::cerr << "Save path = " << outputPath << std::endl;

            Timer timer_master;
            timer_master.start();

            Timer timer_overall;
            timer_overall.start();

            GridInfo gridDEMInfo;
            std::vector<TileInfo> tileDEMInfos;
            readTXTInfo(pathDEM, outputPath, gridDEMInfo, tileDEMInfos);

            GridInfo gridDirInfo;
            std::vector<TileInfo> tileDirInfos;
            readTXTInfo(pathDir, outputPath, gridDirInfo, tileDirInfos);

            std::cout << "Total rows and cols in the grid = " << gridDEMInfo.grandHeight << "*" << gridDEMInfo.grandWidth << std::endl;
            std::cout << "Grid has been divided into " << gridDEMInfo.gridHeight << "*" << gridDEMInfo.gridWidth << " tiles." << std::endl;
            std::cout << "Tile dimensions = " << gridDEMInfo.tileHeight << "*" << gridDEMInfo.tileWidth << std::endl;

            timer_overall.stop();
            std::cerr << "Preparer time = " << timer_overall.elapsed() << "s" << std::endl;

            ObjectFactory objectFactory;
            Preparer(gridDEMInfo, tileDEMInfos, gridDirInfo, tileDirInfos, &objectFactory, prefixName);

            timer_master.stop();
            std::cerr << "Total wall-time= " << timer_master.elapsed() << "s" << std::endl;
        }
        else {
            int good_to_go;
            CommBroadcast(&good_to_go, 0);
            if (good_to_go) {

                ObjectFactory pIObjFactory;
                GridInfo gridDEMInfo, gridDirInfo;
                while (true) {
                    int the_job = CommGetTag(0);

                    if (the_job == TagFirst) {
                        Timer timer_overall;
                        timer_overall.start();
                        TileInfo tileDEMInfo, tileDirInfo;
                        std::string filename;
                        CommRecv(&gridDEMInfo, &gridDirInfo, &tileDEMInfo, &tileDirInfo, &filename, 0);
                        std::shared_ptr<IConsumer2Producer> pIConsumer2Producer = pIObjFactory.createConsumer2Producer();
                        std::shared_ptr<IConsumer> pIConsumer = pIObjFactory.createConsumer();
                        pIConsumer->processRound1(gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, pIConsumer2Producer.get());
                        timer_overall.stop();
                        long vmpeak, vmhwm;
                        ProcessMemUsage(vmpeak, vmhwm);
                        Consumer2Producer* p = (Consumer2Producer*)pIConsumer2Producer.get();
                        p->time_info = TimeInfo(((Consumer*)pIConsumer.get())->timer_calc.elapsed(), timer_overall.elapsed(), ((Consumer*)pIConsumer.get())->timer_io.elapsed(), vmpeak, vmhwm);
                        CommSend(p, nullptr, 0, ObjectFirst);
                    }
                    else if (the_job == TagSecond) {
                        Timer timer_overall;
                        timer_overall.start();
                        TileInfo tileDEMInfo, tileDirInfo;
                        Producer2Consumer p2c;
                        std::string filename;
                        CommRecv(&gridDEMInfo, &gridDirInfo, &tileDEMInfo, &tileDirInfo, &filename, &p2c, 0);
                        std::shared_ptr<IConsumer> pIConsumer = pIObjFactory.createConsumer();
                        pIConsumer->processRound2(gridDEMInfo, gridDirInfo, tileDEMInfo, tileDirInfo, filename, &p2c);
                        timer_overall.stop();
                        long vmpeak, vmhwm;
                        ProcessMemUsage(vmpeak, vmhwm);
                        TimeInfo temp(((Consumer*)pIConsumer.get())->timer_calc.elapsed(), timer_overall.elapsed(), ((Consumer*)pIConsumer.get())->timer_io.elapsed(), vmpeak, vmhwm);
                        CommSend(&temp, nullptr, 0, ObjectSecond);
                    }
                    else if (the_job == TagNull)
                        break;
                }
            }
        }
        MPI_Finalize();
    }
    else if (strcmp(method, "recursive") == 0 || strcmp(method, "wang") == 0 || strcmp(method, "proposed_sequential") == 0) {
        CLI::App app("Sequential-Multi-Flow-Accumulation Program");
        std::string inputDEMFile, inputDirFile, outputAccumFile;
        app.add_option("method", method, "Sequential model")->required();
        app.add_option("inputDEMFile", inputDEMFile, "Path of DEM input file")->required();
        app.add_option("inputDirFile", inputDirFile, "Path of flow-directions input file")->required();
        app.add_option("outputAccumFile", outputAccumFile, "Path of flow-accumulation output file")->required();
        CLI11_PARSE(app, argc, argv);
        std::cout << "Sequential model = " << method << std::endl;
        std::cout << "Input DEM tif = " << inputDEMFile << std::endl;
        std::cout << "Input Dir tif = " << inputDirFile << std::endl;
        std::cout << "Output Accum tif = " << outputAccumFile << std::endl;

        computeAccum(inputDEMFile, inputDirFile, outputAccumFile, method);
    }
    else {
        std::cout << "Parameters error! You should specify 'parallel', 'recursive','wang',or 'proposed_sequential' mode." << std::endl;
        return -1;
    }

    return 0;
}
