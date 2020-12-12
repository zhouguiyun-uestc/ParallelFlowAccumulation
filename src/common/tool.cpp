#include <paradem/gdal.h>
#include <paradem/grid.h>
#include <paradem/i_consumer.h>
#include <paradem/i_object_factory.h>
#include <paradem/i_producer.h>
#include <paradem/object_deleter.h>
#include <paradem/tool.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

bool generateTiles(const char* filePath, int tileHeight, int tileWidth, const char* outputFolder, const std::string filename, const std::string dojob) {

    std::string inputFilePath = filePath;
    std::string output = outputFolder;
    GDALAllRegister();
    GDALDataset* fin = (GDALDataset*)GDALOpen(inputFilePath.c_str(), GA_ReadOnly);
    if (fin == NULL)
        throw std::runtime_error("Could not open file '" + inputFilePath + "' to get dimensions.");
    GDALRasterBand* band = fin->GetRasterBand(1);

    GDALDataType type = band->GetRasterDataType();
    int grandHeight = band->GetYSize();
    int grandWidth = band->GetXSize();
    std::vector<double> geotransform(6);
    fin->GetGeoTransform(&geotransform[0]);
    int height, width;
    int gridHeight = std::ceil((double)grandHeight / tileHeight);
    int gridWidth = std::ceil((double)grandWidth / tileWidth);
    double cellSize = geotransform[1];  //+cell size info
    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            std::string outputFileName = filename + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            std::string path = output + "/" + outputFileName;
            height = (grandHeight - tileHeight * tileRow >= tileHeight) ? tileHeight : (grandHeight - tileHeight * tileRow);
            width = (grandWidth - tileWidth * tileCol >= tileWidth) ? tileWidth : (grandWidth - tileWidth * tileCol);

            if (dojob == "dem") {
                Raster<float> tile;
                if (!tile.init(height, width)) {
                    GDALClose((GDALDatasetH)fin);
                    return false;
                }
                band->RasterIO(GF_Read, tileWidth * tileCol, tileHeight * tileRow, tile.getWidth(), tile.getHeight(), (void*)&tile, tile.getWidth(), tile.getHeight(), type, 0, 0);
                std::vector<double> tileGeotransform(geotransform);
                tileGeotransform[0] = geotransform[0] + tileWidth * tileCol * geotransform[1] + tileHeight * tileRow * geotransform[2];
                tileGeotransform[3] = geotransform[3] + tileHeight * tileRow * geotransform[5] + tileWidth * tileCol * geotransform[4];
                WriteGeoTIFF(path.data(), tile.getHeight(), tile.getWidth(), &tile, type, &tileGeotransform[0], nullptr, nullptr, nullptr, nullptr, tile.NoDataValue);
            }
            else if (dojob == "dir") {
                Raster<uint8_t> tile;
                if (!tile.init(height, width)) {
                    GDALClose((GDALDatasetH)fin);
                    return false;
                }
                tile.NoDataValue = (uint8_t)band->GetNoDataValue();
                band->RasterIO(GF_Read, tileWidth * tileCol, tileHeight * tileRow, tile.getWidth(), tile.getHeight(), (void*)&tile, tile.getWidth(), tile.getHeight(), type, 0, 0);
                std::vector<double> tileGeotransform(geotransform);
                tileGeotransform[0] = geotransform[0] + tileWidth * tileCol * geotransform[1] + tileHeight * tileRow * geotransform[2];
                tileGeotransform[3] = geotransform[3] + tileHeight * tileRow * geotransform[5] + tileWidth * tileCol * geotransform[4];
                WriteGeoTIFF(path.data(), tile.getHeight(), tile.getWidth(), &tile, type, &tileGeotransform[0], nullptr, nullptr, nullptr, nullptr, tile.NoDataValue);
            }
        }
    }
    GDALClose((GDALDatasetH)fin);
    std::string txtPath = output + "/" + "gridInfo.txt";
    std::ofstream fout;
    fout.open(txtPath, std::ofstream::app);
    if (fout.fail()) {
        std::cout << "Open " << txtPath << " error!" << std::endl;
        return false;
    }
    fout << tileHeight << std::endl;
    fout << tileWidth << std::endl;
    fout << gridHeight << std::endl;
    fout << gridWidth << std::endl;
    fout << grandHeight << std::endl;
    fout << grandWidth << std::endl;
    fout << cellSize << std::endl;  //+cell size info
    fout.close();
    return true;
}

bool readGridInfo(const char* tileFolder, GridInfo& gridInfo) {
    std::ifstream infile;
    infile.open(tileFolder);
    if (!infile.is_open()) {
        std::cout << "Open gridInfo error!" << std::endl;
        return false;
    }
    std::string s;
    std::vector<std::string> input(10);
    int k = 0;
    while (getline(infile, s))
        input[k++] = s;

    size_t pos;
    gridInfo.tileHeight = stod(input[0], &pos);
    gridInfo.tileWidth = stod(input[1], &pos);
    gridInfo.gridHeight = stod(input[2], &pos);
    gridInfo.gridWidth = stod(input[3], &pos);
    gridInfo.grandHeight = stod(input[4], &pos);
    gridInfo.grandWidth = stod(input[5], &pos);
    gridInfo.cellSize = stod(input[6], &pos);  //+cell size
    infile.close();
    return true;
}

void createTileInfoArray(GridInfo& gridInfo, std::vector<TileInfo>& tileInfos) {
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    tileInfos.resize(gridHeight * gridWidth);

    int grandHeight = gridInfo.grandHeight;
    int grandWidth = gridInfo.grandWidth;
    int tileHeight = gridInfo.tileHeight;
    int tileWidth = gridInfo.tileWidth;

    int height, width;
    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            height = (grandHeight - tileHeight * tileRow >= tileHeight) ? tileHeight : (grandHeight - tileHeight * tileRow);
            width = (grandWidth - tileWidth * tileCol >= tileWidth) ? tileWidth : (grandWidth - tileWidth * tileCol);
            TileInfo tileInfo;
            tileInfo.gridRow = tileRow;
            tileInfo.gridCol = tileCol;
            tileInfo.height = height;
            tileInfo.width = width;
            tileInfos[tileRow * gridWidth + tileCol] = tileInfo;
        }
    }
}

bool mergeTiles(GridInfo& gridInfo) {
    int grandHeight = gridInfo.grandHeight;
    int grandWidth = gridInfo.grandWidth;
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    int tileHeight = gridInfo.tileHeight;
    int tileWidth = gridInfo.tileWidth;
    Raster<double> tiles;
    if (!tiles.init(grandHeight, grandWidth))
        return false;
    std::string inputFolder = gridInfo.inputFolder;
    std::vector<double> geotransform(6);
    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            std::string fileName = inputFolder + "/anoka-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            Raster<double> tile;
            if (!readGeoTIFF2(fileName.data(), GDALDataType::GDT_Float64, tile))
                return false;
            if (tileRow == 0 && tileCol == 0) {
                for (int i = 0; i < 6; i++) {
                    geotransform[i] = tile.geoTransforms->at(i);
                }
            }
            //写数据
            int height = tile.getHeight();
            int width = tile.getWidth();
            int startRow = tileHeight * tileRow;
            int startCol = tileWidth * tileCol;
            for (int row = 0; row < height; row++) {
                for (int col = 0; col < width; col++) {
                    tiles.at(startRow + row, startCol + col) = tile.at(row, col);
                }
            }
        }
    }
    //保存tiles
    std::string output = gridInfo.outputFolder;
    std::string path = output + "/mergedem.tif";
    WriteGeoTIFF(path.data(), tiles.getHeight(), tiles.getWidth(), &tiles, GDALDataType::GDT_Float64, &geotransform[0], nullptr, nullptr, nullptr, nullptr, tiles.NoDataValue);
    // delete tiles?
    return true;
}

void createNewTif(const std::string path) {
    Raster<float> tif;
    int height = 21;
    int width = 21;
    tif.init(height, width);
    tif.at(0, 0) = 9;
    tif.at(0, 1) = 9;
    tif.at(0, 2) = 7;
    tif.at(0, 3) = 6;
    tif.at(0, 4) = 7;
    tif.at(0, 5) = 6;
    tif.at(0, 6) = 4;
    tif.at(0, 7) = 2;
    tif.at(0, 8) = 1;
    tif.at(0, 9) = 2;
    tif.at(0, 10) = 4;
    tif.at(0, 11) = 1;
    tif.at(0, 12) = 0;
    tif.at(0, 13) = 1;
    tif.at(0, 14) = 3;
    tif.at(0, 15) = 4;
    tif.at(0, 16) = 4;
    tif.at(0, 17) = 5;
    tif.at(0, 18) = 5;
    tif.at(0, 19) = 6;
    tif.at(0, 20) = 7;
    tif.at(1, 0) = 6;
    tif.at(1, 1) = 7;
    tif.at(1, 2) = 6;
    tif.at(1, 3) = 5;
    tif.at(1, 4) = 5;
    tif.at(1, 5) = 4;
    tif.at(1, 6) = 4;
    tif.at(1, 7) = 4;
    tif.at(1, 8) = 4;
    tif.at(1, 9) = 4;
    tif.at(1, 10) = 5;
    tif.at(1, 11) = 2;
    tif.at(1, 12) = 2;
    tif.at(1, 13) = 5;
    tif.at(1, 14) = 6;
    tif.at(1, 15) = 6;
    tif.at(1, 16) = 5;
    tif.at(1, 17) = 4;
    tif.at(1, 18) = 4;
    tif.at(1, 19) = 6;
    tif.at(1, 20) = 8;
    tif.at(2, 0) = 2;
    tif.at(2, 1) = 5;
    tif.at(2, 2) = 5;
    tif.at(2, 3) = 4;
    tif.at(2, 4) = 4;
    tif.at(2, 5) = 4;
    tif.at(2, 6) = 4;
    tif.at(2, 7) = 4;
    tif.at(2, 8) = 4;
    tif.at(2, 9) = 4;
    tif.at(2, 10) = 5;
    tif.at(2, 11) = 4;
    tif.at(2, 12) = 5;
    tif.at(2, 13) = 6;
    tif.at(2, 14) = 6;
    tif.at(2, 15) = 6;
    tif.at(2, 16) = 5;
    tif.at(2, 17) = 4;
    tif.at(2, 18) = 4;
    tif.at(2, 19) = 5;
    tif.at(2, 20) = 6;
    tif.at(3, 0) = 1;
    tif.at(3, 1) = 2;
    tif.at(3, 2) = 4;
    tif.at(3, 3) = 4;
    tif.at(3, 4) = 4;
    tif.at(3, 5) = 4;
    tif.at(3, 6) = 4;
    tif.at(3, 7) = 4;
    tif.at(3, 8) = 4;
    tif.at(3, 9) = 4;
    tif.at(3, 10) = 5;
    tif.at(3, 11) = 4;
    tif.at(3, 12) = 5;
    tif.at(3, 13) = 6;
    tif.at(3, 14) = 6;
    tif.at(3, 15) = 6;
    tif.at(3, 16) = 5;
    tif.at(3, 17) = 4;
    tif.at(3, 18) = 4;
    tif.at(3, 19) = 5;
    tif.at(3, 20) = 6;
    tif.at(4, 0) = 5;
    tif.at(4, 1) = 4;
    tif.at(4, 2) = 4;
    tif.at(4, 3) = 4;
    tif.at(4, 4) = 4;
    tif.at(4, 5) = 4;
    tif.at(4, 6) = 4;
    tif.at(4, 7) = 5;
    tif.at(4, 8) = 6;
    tif.at(4, 9) = 6;
    tif.at(4, 10) = 7;
    tif.at(4, 11) = 5;
    tif.at(4, 12) = 5;
    tif.at(4, 13) = 6;
    tif.at(4, 14) = 6;
    tif.at(4, 15) = 5;
    tif.at(4, 16) = 4;
    tif.at(4, 17) = 4;
    tif.at(4, 18) = 4;
    tif.at(4, 19) = 4;
    tif.at(4, 20) = 4;
    tif.at(5, 0) = 6;
    tif.at(5, 1) = 4;
    tif.at(5, 2) = 4;
    tif.at(5, 3) = 4;
    tif.at(5, 4) = 4;
    tif.at(5, 5) = 5;
    tif.at(5, 6) = 6;
    tif.at(5, 7) = 7;
    tif.at(5, 8) = 8;
    tif.at(5, 9) = 8;
    tif.at(5, 10) = 6;
    tif.at(5, 11) = 5;
    tif.at(5, 12) = 5;
    tif.at(5, 13) = 6;
    tif.at(5, 14) = 7;
    tif.at(5, 15) = 6;
    tif.at(5, 16) = 4;
    tif.at(5, 17) = 4;
    tif.at(5, 18) = 4;
    tif.at(5, 19) = 4;
    tif.at(5, 20) = 4;
    tif.at(6, 0) = 7;
    tif.at(6, 1) = 4;
    tif.at(6, 2) = 4;
    tif.at(6, 3) = 4;
    tif.at(6, 4) = 4;
    tif.at(6, 5) = 5;
    tif.at(6, 6) = 7;
    tif.at(6, 7) = 9;
    tif.at(6, 8) = 9;
    tif.at(6, 9) = 8;
    tif.at(6, 10) = 5;
    tif.at(6, 11) = 5;
    tif.at(6, 12) = 5;
    tif.at(6, 13) = 6;
    tif.at(6, 14) = 8;
    tif.at(6, 15) = 7;
    tif.at(6, 16) = 4;
    tif.at(6, 17) = 4;
    tif.at(6, 18) = 4;
    tif.at(6, 19) = 4;
    tif.at(6, 20) = 4;

    tif.at(7, 0) = 8;
    tif.at(7, 1) = 6;
    tif.at(7, 2) = 5;
    tif.at(7, 3) = 5;
    tif.at(7, 4) = 7;
    tif.at(7, 5) = 6;
    tif.at(7, 6) = 6;
    tif.at(7, 7) = 8;
    tif.at(7, 8) = 7;
    tif.at(7, 9) = 5;
    tif.at(7, 10) = 5;
    tif.at(7, 11) = 5;
    tif.at(7, 12) = 5;
    tif.at(7, 13) = 6;
    tif.at(7, 14) = 8;
    tif.at(7, 15) = 8;
    tif.at(7, 16) = 6;
    tif.at(7, 17) = 4;
    tif.at(7, 18) = 4;
    tif.at(7, 19) = 4;
    tif.at(7, 20) = 6;
    tif.at(8, 0) = 6;
    tif.at(8, 1) = 7;
    tif.at(8, 2) = 8;
    tif.at(8, 3) = 7;
    tif.at(8, 4) = 8;
    tif.at(8, 5) = 7;
    tif.at(8, 6) = 6;
    tif.at(8, 7) = 5;
    tif.at(8, 8) = 5;
    tif.at(8, 9) = 5;
    tif.at(8, 10) = 5;
    tif.at(8, 11) = 5;
    tif.at(8, 12) = 5;
    tif.at(8, 13) = 6;
    tif.at(8, 14) = 7;
    tif.at(8, 15) = 7;
    tif.at(8, 16) = 6;
    tif.at(8, 17) = 4;
    tif.at(8, 18) = 4;
    tif.at(8, 19) = 5;
    tif.at(8, 20) = 6;
    tif.at(9, 0) = 5;
    tif.at(9, 1) = 7;
    tif.at(9, 2) = 8;
    tif.at(9, 3) = 8;
    tif.at(9, 4) = 7;
    tif.at(9, 5) = 7;
    tif.at(9, 6) = 6;
    tif.at(9, 7) = 5;
    tif.at(9, 8) = 5;
    tif.at(9, 9) = 5;
    tif.at(9, 10) = 5;
    tif.at(9, 11) = 5;
    tif.at(9, 12) = 5;
    tif.at(9, 13) = 6;
    tif.at(9, 14) = 7;
    tif.at(9, 15) = 7;
    tif.at(9, 16) = 6;
    tif.at(9, 17) = 4;
    tif.at(9, 18) = 4;
    tif.at(9, 19) = 5;
    tif.at(9, 20) = 6;
    tif.at(10, 0) = 6;
    tif.at(10, 1) = 6;
    tif.at(10, 2) = 6;
    tif.at(10, 3) = 6;
    tif.at(10, 4) = 6;
    tif.at(10, 5) = 5;
    tif.at(10, 6) = 5;
    tif.at(10, 7) = 5;
    tif.at(10, 8) = 5;
    tif.at(10, 9) = 5;
    tif.at(10, 10) = 5;
    tif.at(10, 11) = 5;
    tif.at(10, 12) = 5;
    tif.at(10, 13) = 7;
    tif.at(10, 14) = 6;
    tif.at(10, 15) = 6;
    tif.at(10, 16) = 6;
    tif.at(10, 17) = 6;
    tif.at(10, 18) = 5;
    tif.at(10, 19) = 6;
    tif.at(10, 20) = 8;
    tif.at(11, 0) = 4;
    tif.at(11, 1) = 4;
    tif.at(11, 2) = 4;
    tif.at(11, 3) = 4;
    tif.at(11, 4) = 6;
    tif.at(11, 5) = 7;
    tif.at(11, 6) = 8;
    tif.at(11, 7) = 7;
    tif.at(11, 8) = 6;
    tif.at(11, 9) = 5;
    tif.at(11, 10) = 5;
    tif.at(11, 11) = 5;
    tif.at(11, 12) = 5;
    tif.at(11, 13) = 7;
    tif.at(11, 14) = 6;
    tif.at(11, 15) = 7;
    tif.at(11, 16) = 8;
    tif.at(11, 17) = 7;
    tif.at(11, 18) = 7;
    tif.at(11, 19) = 6;
    tif.at(11, 20) = 6;
    tif.at(12, 0) = 4;
    tif.at(12, 1) = 4;
    tif.at(12, 2) = 4;
    tif.at(12, 3) = 5;
    tif.at(12, 4) = 6;
    tif.at(12, 5) = 7;
    tif.at(12, 6) = 7;
    tif.at(12, 7) = 8;
    tif.at(12, 8) = 7;
    tif.at(12, 9) = 5;
    tif.at(12, 10) = 5;
    tif.at(12, 11) = 5;
    tif.at(12, 12) = 5;
    tif.at(12, 13) = 7;
    tif.at(12, 14) = 6;
    tif.at(12, 15) = 7;
    tif.at(12, 16) = 7;
    tif.at(12, 17) = 7;
    tif.at(12, 18) = 6;
    tif.at(12, 19) = 5;
    tif.at(12, 20) = 4;
    tif.at(13, 0) = 7;
    tif.at(13, 1) = 5;
    tif.at(13, 2) = 5;
    tif.at(13, 3) = 7;
    tif.at(13, 4) = 7;
    tif.at(13, 5) = 6;
    tif.at(13, 6) = 6;
    tif.at(13, 7) = 7;
    tif.at(13, 8) = 7;
    tif.at(13, 9) = 6;
    tif.at(13, 10) = 6;
    tif.at(13, 11) = 5;
    tif.at(13, 12) = 5;
    tif.at(13, 13) = 6;
    tif.at(13, 14) = 6;
    tif.at(13, 15) = 5;
    tif.at(13, 16) = 5;
    tif.at(13, 17) = 5;
    tif.at(13, 18) = 4;
    tif.at(13, 19) = 4;
    tif.at(13, 20) = 3;

    tif.at(14, 0) = 8;
    tif.at(14, 1) = 8;
    tif.at(14, 2) = 8;
    tif.at(14, 3) = 7;
    tif.at(14, 4) = 6;
    tif.at(14, 5) = 6;
    tif.at(14, 6) = 6;
    tif.at(14, 7) = 7;
    tif.at(14, 8) = 8;
    tif.at(14, 9) = 8;
    tif.at(14, 10) = 8;
    tif.at(14, 11) = 7;
    tif.at(14, 12) = 5;
    tif.at(14, 13) = 4;
    tif.at(14, 14) = 4;
    tif.at(14, 15) = 4;
    tif.at(14, 16) = 3;
    tif.at(14, 17) = 3;
    tif.at(14, 18) = 3;
    tif.at(14, 19) = 3;
    tif.at(14, 20) = 4;
    tif.at(15, 0) = 8;
    tif.at(15, 1) = 8;
    tif.at(15, 2) = 7;
    tif.at(15, 3) = 7;
    tif.at(15, 4) = 6;
    tif.at(15, 5) = 6;
    tif.at(15, 6) = 6;
    tif.at(15, 7) = 7;
    tif.at(15, 8) = 8;
    tif.at(15, 9) = 8;
    tif.at(15, 10) = 8;
    tif.at(15, 11) = 7;
    tif.at(15, 12) = 5;
    tif.at(15, 13) = 4;
    tif.at(15, 14) = 5;
    tif.at(15, 15) = 4;
    tif.at(15, 16) = 3;
    tif.at(15, 17) = 3;
    tif.at(15, 18) = 3;
    tif.at(15, 19) = 3;
    tif.at(15, 20) = 4;
    tif.at(16, 0) = 8;
    tif.at(16, 1) = 7;
    tif.at(16, 2) = 6;
    tif.at(16, 3) = 6;
    tif.at(16, 4) = 6;
    tif.at(16, 5) = 6;
    tif.at(16, 6) = 6;
    tif.at(16, 7) = 7;
    tif.at(16, 8) = 7;
    tif.at(16, 9) = 7;
    tif.at(16, 10) = 7;
    tif.at(16, 11) = 6;
    tif.at(16, 12) = 6;
    tif.at(16, 13) = 6;
    tif.at(16, 14) = 5;
    tif.at(16, 15) = 4;
    tif.at(16, 16) = 3;
    tif.at(16, 17) = 3;
    tif.at(16, 18) = 4;
    tif.at(16, 19) = 5;
    tif.at(16, 20) = 5;
    tif.at(17, 0) = 9;
    tif.at(17, 1) = 7;
    tif.at(17, 2) = 5;
    tif.at(17, 3) = 4;
    tif.at(17, 4) = 6;
    tif.at(17, 5) = 7;
    tif.at(17, 6) = 7;
    tif.at(17, 7) = 6;
    tif.at(17, 8) = 5;
    tif.at(17, 9) = 4;
    tif.at(17, 10) = 4;
    tif.at(17, 11) = 5;
    tif.at(17, 12) = 6;
    tif.at(17, 13) = 6;
    tif.at(17, 14) = 6;
    tif.at(17, 15) = 5;
    tif.at(17, 16) = 5;
    tif.at(17, 17) = 3;
    tif.at(17, 18) = 5;
    tif.at(17, 19) = 7;
    tif.at(17, 20) = 8;
    tif.at(18, 0) = 8;
    tif.at(18, 1) = 6;
    tif.at(18, 2) = 5;
    tif.at(18, 3) = 4;
    tif.at(18, 4) = 6;
    tif.at(18, 5) = 7;
    tif.at(18, 6) = 7;
    tif.at(18, 7) = 6;
    tif.at(18, 8) = 5;
    tif.at(18, 9) = 5;
    tif.at(18, 10) = 6;
    tif.at(18, 11) = 4;
    tif.at(18, 12) = 5;
    tif.at(18, 13) = 6;
    tif.at(18, 14) = 6;
    tif.at(18, 15) = 5;
    tif.at(18, 16) = 5;
    tif.at(18, 17) = 6;
    tif.at(18, 18) = 6;
    tif.at(18, 19) = 6;
    tif.at(18, 20) = 6;
    tif.at(19, 0) = 4;
    tif.at(19, 1) = 4;
    tif.at(19, 2) = 4;
    tif.at(19, 3) = 5;
    tif.at(19, 4) = 5;
    tif.at(19, 5) = 6;
    tif.at(19, 6) = 6;
    tif.at(19, 7) = 6;
    tif.at(19, 8) = 6;
    tif.at(19, 9) = 5;
    tif.at(19, 10) = 6;
    tif.at(19, 11) = 4;
    tif.at(19, 12) = 4;
    tif.at(19, 13) = 6;
    tif.at(19, 14) = 5;
    tif.at(19, 15) = 5;
    tif.at(19, 16) = 5;
    tif.at(19, 17) = 8;
    tif.at(19, 18) = 7;
    tif.at(19, 19) = 5;
    tif.at(19, 20) = 3;
    tif.at(20, 0) = 0;
    tif.at(20, 1) = 1;
    tif.at(20, 2) = 2;
    tif.at(20, 3) = 5;
    tif.at(20, 4) = 5;
    tif.at(20, 5) = 6;
    tif.at(20, 6) = 6;
    tif.at(20, 7) = 4;
    tif.at(20, 8) = 6;
    tif.at(20, 9) = 6;
    tif.at(20, 10) = 3;
    tif.at(20, 11) = 3;
    tif.at(20, 12) = 4;
    tif.at(20, 13) = 5;
    tif.at(20, 14) = 3;
    tif.at(20, 15) = 3;
    tif.at(20, 16) = 5;
    tif.at(20, 17) = 9;
    tif.at(20, 18) = 7;
    tif.at(20, 19) = 4;
    tif.at(20, 20) = 1;

    std::vector<double> geotransform(6, 0);  //需要添加信息
    geotransform[1] = 1;
    geotransform[5] = -1;
    geotransform[0] = 1000;
    geotransform[3] = 1000;
    std::string tifPath = path + "/dem000.tif";
    WriteGeoTIFF(tifPath.data(), height, width, &tif, GDALDataType::GDT_Float32, &geotransform[0], nullptr, nullptr, nullptr, nullptr, tif.NoDataValue);
}

bool mergeDEM(GridInfo& gridInfo) {
    int grandHeight = gridInfo.grandHeight;
    int grandWidth = gridInfo.grandWidth;
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    int tileHeight = gridInfo.tileHeight;
    int tileWidth = gridInfo.tileWidth;
    // Raster<uint8_t> tiles;
    Raster<float> tiles;
    if (!tiles.init(grandHeight, grandWidth))
        return false;
    tiles.setAllValues(tiles.NoDataValue);
    std::string inputFolder = gridInfo.inputFolder;
    std::vector<double> geotransform(6);

    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            std::cout << tileRow << "," << tileCol << std::endl;
            std::string fileName = inputFolder + "/aitkin-3m-filling-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            Raster<float> tile;
            if (!readGeoTIFF(fileName.data(), GDALDataType::GDT_Float32, tile))
                return false;

            if (tileRow == 0 && tileCol == 0) {
                for (int i = 0; i < 6; i++) {
                    geotransform[i] = tile.geoTransforms->at(i);
                }
            }
            //写数据
            int height = tile.getHeight();
            int width = tile.getWidth();
            int startRow = tileHeight * tileRow;
            int startCol = tileWidth * tileCol;
            for (int row = 0; row < height; row++) {
                for (int col = 0; col < width; col++) {
                    if (tile.isNoData(row, col))
                        continue;
                    tiles.at(startRow + row, startCol + col) = tile.at(row, col);
                }
            }
        }
    }

    //保存tiles
    std::string output = gridInfo.outputFolder;
    std::string path;
    WriteGeoTIFF(path.data(), tiles.getHeight(), tiles.getWidth(), &tiles, GDALDataType::GDT_Float32, &geotransform[0], nullptr, nullptr, nullptr, nullptr, tiles.NoDataValue);
    // WriteGeoTIFF(path.data(), tiles.getHeight(), tiles.getWidth(), &tiles, GDALDataType::GDT_Byte, &geotransform[0], nullptr, nullptr, nullptr, nullptr, tiles.NoDataValue);
    // delete tiles?
    return true;
}
