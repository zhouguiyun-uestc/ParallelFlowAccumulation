#include "mfd_flow_direction.h"

int dirs[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

Raster<uint8_t> computeDir(Raster<float>& dem, Raster<uint8_t>& flowDir) {
    int height = dem.getHeight();
    int width = dem.getWidth();
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (flowDir.isNoData(row, col)) {
                continue;
            }
            uint8_t dir = 0;
            for (int i = 0; i < 8; i++) {
                int iRow, iCol;
                iRow = dem.getRow(i, row);
                iCol = dem.getCol(i, col);
                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                    continue;
                if (dem.at(row, col) > dem.at(iRow, iCol))
                    dir += dirs[i];
            }
            if (dir != 0 && dir != flowDir.at(row, col))
                flowDir.at(row, col) = dir;
        }
    }

    Raster<uint8_t> newDir;
    newDir.init(flowDir.getHeight() - 2, flowDir.getWidth() - 2);
    newDir.NoDataValue = 0;
    newDir.geoTransforms = flowDir.geoTransforms;
    for (int row = 0; row < newDir.getHeight(); row++) {
        for (int col = 0; col < newDir.getWidth(); col++) {
            if (flowDir.isNoData(row + 1, col + 1))
                continue;
            newDir.at(row, col) = flowDir.at(row + 1, col + 1);
        }
    }
    return newDir;
}
bool haveData(int row, int col, int gridHeight, int gridWidth) {
    if ((row >= 0 && row < gridHeight) && (col >= 0 && col < gridWidth))
        return true;

    return false;
}

void setborder(const int i, Raster<uint8_t>& bigDir, Raster<int>& iDir) {
    if (i == 0) {
        int col = bigDir.getWidth() - 1;
        for (int row = 1; row < bigDir.getHeight() - 1; row++) {
            if (iDir.isNoData(row - 1, 0))
                continue;
            bigDir.at(row, col) = 2 ^ iDir.at(row - 1, 0);
        }
    }
    else if (i == 1) {
        if (iDir.isNoData(0, 0))
            return;
        bigDir.at(bigDir.getWidth() - 1, bigDir.getHeight() - 1) = pow(2, iDir.at(0, 0));
    }
    else if (i == 2) {
        int row = bigDir.getHeight() - 1;
        for (int col = 1; col < bigDir.getWidth() - 1; col++) {
            if (iDir.isNoData(0, col - 1))
                continue;
            bigDir.at(row, col) = pow(2, iDir.at(0, col - 1));
        }
    }
    else if (i == 3) {
        if (iDir.isNoData(0, iDir.getWidth() - 1))
            return;
        bigDir.at(bigDir.getHeight() - 1, 0) = pow(2, iDir.at(0, iDir.getWidth() - 1));
    }
    else if (i == 4) {
        int col = 0;
        for (int row = 1; row < bigDir.getHeight() - 1; row++) {
            if (iDir.isNoData(row - 1, iDir.getWidth() - 1))
                continue;

            bigDir.at(row, col) = pow(2, iDir.at(row - 1, iDir.getWidth() - 1));
        }
    }
    else if (i == 5) {
        if (iDir.isNoData(iDir.getHeight() - 1, iDir.getWidth() - 1))
            return;
        bigDir.at(0, 0) = pow(2, iDir.at(iDir.getHeight() - 1, iDir.getWidth() - 1));
    }
    else if (i == 6) {
        int row = 0;
        for (int col = 1; col < bigDir.getWidth() - 1; col++) {
            if (iDir.isNoData(iDir.getHeight() - 1, col - 1))
                continue;
            bigDir.at(row, col) = pow(2, iDir.at(iDir.getHeight() - 1, col - 1));
        }
    }
    else if (i == 7) {
        if (iDir.isNoData(iDir.getHeight() - 1, 0))
            return;
        bigDir.at(0, bigDir.getWidth() - 1) = pow(2, iDir.at(iDir.getHeight() - 1, 0));
    }
}

bool getBigData(Raster<int>& dir, Raster<uint8_t>& bigDir, const int tileRow, const int tileCol, const int gridHeight, const int gridWidth, const std::string pathDir, const std::string filenameDir) {

    bigDir.init(dir.getHeight() + 2, dir.getWidth() + 2);
    bigDir.NoDataValue = 0;
    bigDir.geoTransforms = dir.geoTransforms;
    bigDir.setAllValues(0);
    for (int row = 1; row < bigDir.getHeight() - 1; row++) {
        for (int col = 1; col < bigDir.getWidth() - 1; col++) {
            if (dir.isNoData(row - 1, col - 1))
                continue;
            bigDir.at(row, col) = pow(2, dir.at(row - 1, col - 1));
        }
    }
    for (int i = 0; i < 8; i++) {
        int iTileRow, iTileCol;
        iTileRow = bigDir.getRow(i, tileRow);
        iTileCol = bigDir.getCol(i, tileCol);
        if (!haveData(iTileRow, iTileCol, gridHeight, gridWidth))
            continue;
        std::string iPath = pathDir + "/" + filenameDir + "-" + std::to_string(iTileCol) + "_" + std::to_string(iTileRow) + "flowdir.tif";
        Raster<int> iDir;
        if (!readTif(iPath.data(), iDir))
            return false;
        setborder(i, bigDir, iDir);
    }

    return true;
}

void setborder(const int i, Raster<float>& bigDEM, Raster<float>& iDEM) {
    if (i == 0) {
        int col = bigDEM.getWidth() - 1;
        for (int row = 1; row < bigDEM.getHeight() - 1; row++) {
            if (iDEM.isNoData(row - 1, 0))
                continue;
            bigDEM.at(row, col) = iDEM.at(row - 1, 0);
        }
    }
    else if (i == 1) {
        if (iDEM.isNoData(0, 0))
            return;
        bigDEM.at(bigDEM.getWidth() - 1, bigDEM.getHeight() - 1) = iDEM.at(0, 0);
    }
    else if (i == 2) {
        int row = bigDEM.getHeight() - 1;
        for (int col = 1; col < bigDEM.getWidth() - 1; col++) {
            if (iDEM.isNoData(0, col - 1))
                continue;
            bigDEM.at(row, col) = iDEM.at(0, col - 1);
        }
    }
    else if (i == 3) {
        if (iDEM.isNoData(0, iDEM.getWidth() - 1))
            return;
        bigDEM.at(bigDEM.getHeight() - 1, 0) = iDEM.at(0, iDEM.getWidth() - 1);
    }
    else if (i == 4) {
        int col = 0;
        for (int row = 1; row < bigDEM.getHeight() - 1; row++) {
            if (iDEM.isNoData(row - 1, iDEM.getWidth() - 1))
                continue;

            bigDEM.at(row, col) = iDEM.at(row - 1, iDEM.getWidth() - 1);
        }
    }
    else if (i == 5) {
        if (iDEM.isNoData(iDEM.getHeight() - 1, iDEM.getWidth() - 1))
            return;
        bigDEM.at(0, 0) = iDEM.at(iDEM.getHeight() - 1, iDEM.getWidth() - 1);
    }
    else if (i == 6) {
        int row = 0;
        for (int col = 1; col < bigDEM.getWidth() - 1; col++) {
            if (iDEM.isNoData(iDEM.getHeight() - 1, col - 1))
                continue;
            bigDEM.at(row, col) = iDEM.at(iDEM.getHeight() - 1, col - 1);
        }
    }
    else if (i == 7) {
        if (iDEM.isNoData(iDEM.getHeight() - 1, 0))
            return;
        bigDEM.at(0, bigDEM.getWidth() - 1) = iDEM.at(iDEM.getHeight() - 1, 0);
    }
}

bool getBigData(Raster<float>& dem, Raster<float>& bigDEM, const int tileRow, const int tileCol, const int gridHeight, const int gridWidth, const std::string pathDEM, const std::string filenameDEM) {

    bigDEM.init(dem.getHeight() + 2, dem.getWidth() + 2);
    bigDEM.geoTransforms = dem.geoTransforms;
    bigDEM.setAllValues(bigDEM.NoDataValue);
    for (int row = 1; row < bigDEM.getHeight() - 1; row++) {
        for (int col = 1; col < bigDEM.getWidth() - 1; col++) {
            if (dem.isNoData(row - 1, col - 1))
                continue;
            bigDEM.at(row, col) = dem.at(row - 1, col - 1);
        }
    }
    for (int i = 0; i < 8; i++) {
        int iTileRow, iTileCol;
        iTileRow = bigDEM.getRow(i, tileRow);
        iTileCol = bigDEM.getCol(i, tileCol);
        if (!haveData(iTileRow, iTileCol, gridHeight, gridWidth))
            continue;
        std::string iPath = pathDEM + "/" + filenameDEM + "-" + std::to_string(iTileCol) + "_" + std::to_string(iTileRow) + ".tif";
        Raster<float> iDEM;
        if (!readTif(iPath.data(), iDEM))
            return false;
        setborder(i, bigDEM, iDEM);
    }

    return true;
}

bool preprocess(const std::string& pathDir, const std::string& filenameDir, const std::string& pathDEM, const std::string& filenameDEM, const std::string& saveDirPath, const std::string& saveDEMPath,
                const std::string& savefilename) {
    std::string gridInfoPath = pathDir + "/gridInfo.txt";
    GridInfo gridInfo;
    if (!readGridInfo(gridInfoPath.data(), gridInfo))
        return false;
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            std::string fileNameDir = pathDir + "/" + filenameDir + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + "flowdir.tif";
            Raster<int> dir;
            if (!readTif(fileNameDir.data(), dir))
                return false;
            Raster<uint8_t> bigDir;
            if (!getBigData(dir, bigDir, tileRow, tileCol, gridHeight, gridWidth, pathDir, filenameDir))
                return false;
            std::string fileNameDEM = pathDEM + "/" + filenameDEM + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            Raster<float> dem;
            if (!readTif(fileNameDEM.data(), dem))
                return false;
            Raster<float> bigDEM;
            if (!getBigData(dem, bigDEM, tileRow, tileCol, gridHeight, gridWidth, pathDEM, filenameDEM))
                return false;
            Raster<uint8_t> newDir = computeDir(bigDEM, bigDir);
            std::string path = saveDirPath + "/" + savefilename + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            if (!WriteGeoTIFF(path.data(), newDir.getHeight(), newDir.getWidth(), &newDir, GDALDataType::GDT_Byte, newDir.getGeoTransformsPtr(), nullptr, nullptr, nullptr, nullptr,
                              newDir.getNoDataValue()))
                return false;
            float nodata = -9999;
            for (int row = 0; row < dem.getHeight(); row++) {
                for (int col = 0; col < dem.getWidth(); col++) {
                    if (dem.isNoData(row, col))
                        dem.at(row, col) = nodata;
                }
            }
            dem.NoDataValue = nodata;
            std::string pathDEM = saveDEMPath + "/" + savefilename + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            if (!WriteGeoTIFF(pathDEM.data(), dem.getHeight(), dem.getWidth(), &dem, GDALDataType::GDT_Float32, dem.getGeoTransformsPtr(), nullptr, nullptr, nullptr, nullptr, dem.NoDataValue))
                return false;
        }
    }

    return true;
}

bool checkDir(const std::string& pathDir, const std::string& filenameDir, const std::string& pathDEM) {
    std::string gridInfoPath = pathDir + "/gridInfo.txt";
    GridInfo gridInfo;
    if (!readGridInfo(gridInfoPath.data(), gridInfo))
        return false;
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    for (int tileRow = 0; tileRow < gridHeight; tileRow++) {
        for (int tileCol = 0; tileCol < gridWidth; tileCol++) {
            std::string fileNameDir = pathDir + "/" + filenameDir + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            Raster<uint8_t> flowdir;
            if (!readTif(fileNameDir.data(), flowdir))
                return false;

            std::string fileNameDEM = pathDEM + "/" + filenameDir + "-" + std::to_string(tileCol) + "_" + std::to_string(tileRow) + ".tif";
            Raster<float> dem;
            if (!readTif(fileNameDEM.data(), dem))
                return false;

            for (int row = 0; row < flowdir.getHeight(); row++) {
                for (int col = 0; col < flowdir.getWidth(); col++) {
                    if (dem.isNoData(row, col))
                        continue;
                    uint8_t dirs = flowdir.at(row, col);
                    for (int i = 0; i < 8; i++) {
                        if (dirs & (1 << i)) {
                            int iRow = flowdir.getRow(i, row);
                            int iCol = flowdir.getCol(i, col);
                            if (!dem.isInGrid(iRow, iCol))
                                continue;
                            if (dem.at(iRow, iCol) > dem.at(row, col))
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}
