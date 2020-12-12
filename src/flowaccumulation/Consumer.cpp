#include "consumer.h"
#include "perimeters.h"

const double SQRT2 = 1.4142135623730950488016887242097;
const double MAX_SINGLE = -1;
const double NODATA = -9999;

void Consumer::getEdgeLink(const int tileHeight, const int tileWidth, Raster<uint8_t>& flowdirs, Raster<float>& dem, const int cellSize, Consumer2Producer* pC2P) {
    int length = (2 * tileHeight + 2 * tileWidth) * 2 - 16;
    pC2P->links.resize(length);
    pC2P->fraction.resize(length);
    int mid = 2 * tileHeight + 2 * tileWidth - 4;

    for (int i = 0; i < mid; i++) {
        int row, col;
        serialToXYAll(i, row, col, tileWidth, tileHeight);
        std::vector<link_t> links;
        FollowPath(flowdirs, dem, tileHeight, tileWidth, row, col, links);
        pC2P->links[i] = links;
    }

    for (int i = mid; i < length; i++) {
        int row, col;
        serialToXYAll(i, row, col, tileWidth, tileHeight);
        std::vector<link_t> links;
        std::vector<fraction_t> fractions;
        FollowPath(row, col, tileHeight, tileWidth, flowdirs, dem, links, fractions, cellSize);
        pC2P->links[i] = links;
        pC2P->fraction[i] = fractions;
    }
}

bool Consumer::processRound1(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                             IConsumer2Producer* pIC2P) {
    int gridRow = tileDirInfo.gridRow;
    int gridCol = tileDirInfo.gridCol;
    Raster<float> dem;
    std::string pathDEM = gridDEMInfo.inputFolder + "/" + filename + "-" + std::to_string(gridCol) + "_" + std::to_string(gridRow) + ".tif";
    if (!readGeoTIFF(pathDEM.data(), GDALDataType::GDT_Float32, dem)) {
        std::cout << "read " << pathDEM << " error." << std::endl;
        return false;
    }
    Raster<uint8_t> flowDir;
    std::string pathDir = gridDirInfo.inputFolder + "/" + filename + "-" + std::to_string(gridCol) + "_" + std::to_string(gridRow) + ".tif";
    if (!readGeoTIFF3(pathDir.data(), GDALDataType::GDT_Byte, flowDir)) {
        std::cout << "read " << pathDir << " error." << std::endl;
        return false;
    }
    Raster<double> accum;
    FlowAccumulationInnerImprove(dem, flowDir, accum);

    Consumer2Producer* pC2P = (Consumer2Producer*)pIC2P;
    pC2P->gridRow = tileDEMInfo.gridRow;
    pC2P->gridCol = tileDEMInfo.gridCol;
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    GridInnerPerimToArray(accum, pC2P->inner_accum);
    GridAllPerimToArray(dem, pC2P->elevation);
    GridAllPerimToArray(flowDir, pC2P->flowdirs);
    getEdgeLink(tileDirInfo.height, tileDirInfo.width, flowDir, dem, cellSize, pC2P);
    return true;
}

bool Consumer::isInnerEdge(int row, int col, int height, int width) {
    if (row == 1 || col == 1 || row == height - 2 || col == width - 2)
        return true;
    return false;
}

void Consumer::calculateStatistics(Raster<double>& accum, double* min, double* max, double* mean, double* stdDev) {
    int width = accum.getWidth();
    int height = accum.getHeight();

    int validElements = 0;
    double minValue = 0, maxValue = 0;
    double sum = 0.0;
    double sumSqurVal = 0.0;
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            if (!accum.isNoData(row, col)) {
                double value = accum.at(row, col);
                // Initialize minValue and maxValue using the first valid value
                if (validElements == 0) {
                    minValue = maxValue = value;
                }
                validElements++;
                if (minValue > value) {
                    minValue = value;
                }
                if (maxValue < value) {
                    maxValue = value;
                }

                sum += value;
                sumSqurVal += (value * value);
            }
        }
    }

    double meanValue = sum / validElements;
    double stdDevValue = sqrt((sumSqurVal / validElements) - (meanValue * meanValue));
    *min = minValue;
    *max = maxValue;
    *mean = meanValue;
    *stdDev = stdDevValue;
}

bool Consumer::processRound2(const GridInfo& gridDEMInfo, const GridInfo& gridDirInfo, const TileInfo& tileDEMInfo, const TileInfo& tileDirInfo, const std::string& filename,
                             IProducer2Consumer* pIP2C) {
    int gridRow = tileDirInfo.gridRow;
    int gridCol = tileDirInfo.gridCol;
    Raster<float> dem;
    std::string pathDEM = gridDEMInfo.inputFolder + "/" + filename + "-" + std::to_string(gridCol) + "_" + std::to_string(gridRow) + ".tif";
    if (!readGeoTIFF(pathDEM.data(), GDALDataType::GDT_Float32, dem)) {
        std::cout << "read " << pathDEM << " error." << std::endl;
        return false;
    }
    Raster<uint8_t> flowDir;
    std::string pathDir = gridDirInfo.inputFolder + "/" + filename + "-" + std::to_string(gridCol) + "_" + std::to_string(gridRow) + ".tif";
    if (!readGeoTIFF3(pathDir.data(), GDALDataType::GDT_Byte, flowDir)) {
        std::cout << "read " << pathDir << " error." << std::endl;
        return false;
    }
    Raster<double> accum;
    Producer2Consumer* p2c = (Producer2Consumer*)pIP2C;
    FlowAccumulationImprove(dem, flowDir, accum, p2c);
    double min, max, mean, stdDev;
    double noDataValue = -9999.0;
    calculateStatistics(accum, &min, &max, &mean, &stdDev);

    std::string outpath = gridDirInfo.outputFolder + "/" + std::to_string(gridRow) + "_" + std::to_string(gridCol) + ".tif";

    WriteGeoTIFF(outpath.data(), accum.getHeight(), accum.getWidth(), &accum, GDALDataType::GDT_Float64, dem.getGeoTransformsPtr(), &min, &max, &mean, &stdDev, noDataValue);

    return false;
}

void Consumer::free() {
    delete this;
}
