#include "accumulation.h"
#include "Node.h"
#include "perimeters.h"
#include "producer_2_consumer.h"
#include <paradem/gdal.h>
#include <paradem/grid.h>
#include <paradem/raster.h>

#include <iostream>
#include <queue>
#include <stdint.h>

const double SQRT2 = 1.4142135623730950488016887242097;
const double MAX_SINGLE = -1;

bool isOuterEdges(int row, int col, Raster<float>& dem) {
    int height = dem.getHeight();
    int width = dem.getWidth();
    if (row == 0 || col == 0 || row == height - 1 || col == width - 1)
        return true;
    return false;
}
bool isInnerEdges(int row, int col, Raster<float>& dem) {
    int height = dem.getHeight();
    int width = dem.getWidth();
    if (row == 1 || col == 1 || row == height - 2 || col == width - 2)
        return true;
    return false;
}

double getOffset1(int row, int col, int height, int width, Producer2Consumer* p2c) {
    if ((row == 1) || (col == 1) || (row == height - 2) || (col == width - 2)) {
        int s = xyToSerialInner(row, col, width, height);
        return p2c->inner_accum_setoff1[s];
    }

    else
        return 0;
}

double calculateIndex(int row, int col, int cellSize, Raster<float>& dem, Raster<uint8_t>& flowdirs) {
    double dSlope, index;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    double dMax = 0;
    uint8_t dirs = flowdirs.at(row, col);
    for (int i = 0; i < 8; i++) {
        if (dirs & (1 << i)) {
            int iRow, iCol;
            iRow = dem.getRow(i, row);
            iCol = dem.getCol(i, col);
            dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
            if (i % 2 == 1)
                dSlope = dSlope / SQRT2;
            if (dSlope > dMax)
                dMax = dSlope;
        }
    }
    if (dMax > slopeMax)
        index = a + b;
    else if (dMax <= slopeMin)
        index = b;
    else
        index = a * dMax + b;
    return index;
}

std::vector<Node> calculateFra(int row, int col, double value, int cellSize, double index, Raster<float>& dem, Raster<uint8_t>& flowdirs) {
    std::vector<Node> link;
    double dSum = 0;
    double dSlope;
    uint8_t dir = flowdirs.at(row, col);
    if (dir == 0)
        return link;
    for (int k = 0; k < 8; k++) {
        if (dir & (1 << k)) {
            int iRow, iCol;
            iRow = dem.getRow(k, row);
            iCol = dem.getCol(k, col);
            dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
            if (k % 2 == 1) {
                dSlope = dSlope / SQRT2;
                dSlope = pow(dSlope, index) * SQRT2 / 4;
            }
            else {
                dSlope = pow(dSlope, index) / 2;
            }
            Node node(iRow, iCol, dSlope);
            link.push_back(node);
            dSum = dSum + dSlope;
        }
    }

    if (dSum == 0) {
        int count = link.size();
        double temp = 1.0 / count;
        for (int i = 0; i < count; i++) {
            link[i].elevation = temp * value;
        }
    }
    else {
        for (size_t i = 0; i < link.size(); i++) {
            link[i].elevation = link[i].elevation / dSum * value;
        }
    }
    return link;
}

class Point {
public:
    int row, col;
    Point(int row, int col) {
        this->row = row;
        this->col = col;
    }
};

bool handle(const int row, const int col, const std::vector<std::vector<bool>>& flag, Raster<uint8_t>& flowdirs, Raster<float>& dem) {
    unsigned char redir = ~(unsigned char)flowdirs.at(row, col);
    std::vector<int> redirections = { 4, 5, 6, 7, 0, 1, 2, 3 };
    for (int i = 0; i < 8; i++) {
        if (redir & (1 << i)) {
            int iRow = dem.getRow(i, row);
            int iCol = dem.getCol(i, col);
            if (!flowdirs.isInGrid(iRow, iCol))
                continue;
            uint8_t iDir = flowdirs.at(iRow, iCol);
            if (iDir & (1 << redirections[i])) {
                if (flag[iRow][iCol])
                    return false;
            }
        }
    }
    return true;
}

void FollowPath(const int row0, const int col0, const int tileHeight, const int tileWidth, Raster<uint8_t>& flowdirs, Raster<float>& dem, std::vector<link_t>& links,
                std::vector<fraction_t>& fractions, int cellSize) {
    if (dem.isNoData(row0, col0))
        return;
    std::vector<std::vector<double>> fraction(tileHeight, std::vector<double>(tileWidth, 0));
    std::vector<std::vector<bool>> flag(tileHeight, std::vector<bool>(tileWidth, false));
    int row = row0;
    int col = col0;
    Point p(row, col);
    std::queue<Point> flagQ;
    flagQ.push(p);
    while (!flagQ.empty()) {
        Point temp = flagQ.front();
        flagQ.pop();
        int tempRow = temp.row;
        int tempCol = temp.col;
        if (flag[tempRow][tempCol])
            continue;

        flag[tempRow][tempCol] = true;
        int8_t dir = flowdirs.at(tempRow, tempCol);
        for (int i = 0; i < 8; i++) {
            if (dir & (1 << i)) {
                int iRow, iCol;
                iRow = dem.getRow(i, tempRow);
                iCol = dem.getCol(i, tempCol);
                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                    continue;
                if (!(iRow > 1 && iRow < tileHeight - 2 && iCol > 1 && iCol < tileWidth - 2))
                    continue;
                if (flag[iRow][iCol])
                    continue;
                Point tempP(iRow, iCol);
                flagQ.push(tempP);
            }
        }
    }
    std::queue<Point> fracQ;
    fraction[row0][col0] = 1;
    fracQ.push(p);
    while (!fracQ.empty()) {
        Point temp = fracQ.front();
        fracQ.pop();
        int tempRow = temp.row;
        int tempCol = temp.col;

        flag[tempRow][tempCol] = false;

        double index = calculateIndex(tempRow, tempCol, cellSize, dem, flowdirs);
        std::vector<Node> link = calculateFra(tempRow, tempCol, fraction[tempRow][tempCol], cellSize, index, dem, flowdirs);
        for (size_t i = 0; i < link.size(); i++) {
            int iRow = link[i].row;
            int iCol = link[i].col;
            fraction[iRow][iCol] += link[i].elevation;
            if (!flag[iRow][iCol])
                continue;
            if (handle(iRow, iCol, flag, flowdirs, dem)) {
                Point inputCell(iRow, iCol);
                fracQ.push(inputCell);
            }
        }
    }

    std::vector<fraction_t> perimeters;
    GridAllPerimToArray(fraction, perimeters);
    int i0 = xyToSerialAll(row0, col0, tileWidth, tileHeight);
    for (int i = 0; i < (int)perimeters.size(); i++) {
        if (perimeters[i] != 0 && i != i0) {
            links.push_back(i);
            fractions.push_back(perimeters[i]);
        }
    }
}

void FollowPath(Raster<uint8_t>& flowdirs, Raster<float>& dem, const int tileHeight, const int tileWidth, const int row, const int col, std::vector<link_t>& links) {
    if (!dem.isNoData(row, col)) {
        uint8_t dir = flowdirs.at(row, col);
        for (int k = 0; k < 8; k++) {
            if (dir & (1 << k)) {
                int iRow, iCol;
                iRow = dem.getRow(k, row);
                iCol = dem.getCol(k, col);
                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                    continue;
                links.push_back(xyToSerialAll(iRow, iCol, tileWidth, tileHeight));
            }
        }
    }
}

extern int BitCount2(uint8_t n, int& k) {
    int cnt = 0;
    for (; n; n >>= 1) {
        cnt += n & 1;
        k++;
    }
    return cnt;
}

void FlowAccumulationInnerTradition(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum) {

    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }

    std::vector<std::vector<int>> dependencies(height, std::vector<int>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            accum.at(row, col) = 1;
            if (flowdirs.isNoData(row, col))
                continue;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dependencies[iRow][iCol]++;
                }
            }
        }
    }

    std::queue<Node> sources;
    for (int row = 1; row < height - 1; row++)
        for (int col = 1; col < width - 1; col++)
            if (dependencies[row][col] == 0 && !dem.isNoData(row, col)) {
                Node node(row, col, dem.at(row, col));
                sources.push(node);
            }

    while (!sources.empty()) {
        Node node = sources.front();
        sources.pop();
        if (dem.isNoData(node.row, node.col))
            continue;
        std::vector<double> slopes(8, 0);
        double fracSum = 0;
        if (flowdirs.isNoData(node.row, node.col))
            continue;  //+-?
        uint8_t dir = flowdirs.at(node.row, node.col);
        for (int k = 0; k < 8; k++) {
            if (dir & (1 << k)) {
                int iRow, iCol;
                iRow = dem.getRow(k, node.row);
                iCol = dem.getCol(k, node.col);
                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                    continue;
                dSlope = (dem.at(node.row, node.col) - dem.at(iRow, iCol)) / cellSize;
                if (k % 2 == 1) {
                    dSlope = dSlope / SQRT2;
                    slopes[k] = pow(dSlope, index[node.row][node.col]) * SQRT2 / 4;
                }
                else
                    slopes[k] = pow(dSlope, index[node.row][node.col]) / 2;
                fracSum += slopes[k];
            }
        }

        if (fracSum == 0) {
            //平地
            int count = 0;
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    count++;
                }
            }
            double temp = accum.at(node.row, node.col) / count;
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, node.row);
                    iCol = dem.getCol(k, node.col);
                    if (isOuterEdges(iRow, iCol, dem))
                        continue;
                    accum.at(iRow, iCol) += temp;
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                        sources.push(tempNode);
                    }
                }
            }
        }
        else {
            for (int k = 0; k < 8; k++) {
                if (slopes[k] != 0) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, node.row);
                    iCol = dem.getCol(k, node.col);
                    if (isOuterEdges(iRow, iCol, dem))
                        continue;
                    double temp = slopes[k] / fracSum;
                    accum.at(iRow, iCol) += accum.at(node.row, node.col) * temp;
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                        sources.push(tempNode);
                    }
                }
            }
        }
    }
}

void FlowAccumulationInnerImprove(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum) {
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }

    std::vector<std::vector<int>> dependencies(height, std::vector<int>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            accum.at(row, col) = 1;
            if (flowdirs.isNoData(row, col))
                continue;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dependencies[iRow][iCol]++;
                }
            }
        }
    }

    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (!dem.isNoData(row, col) && dependencies[row][col] == 0)  // dp=0, 加入队列，如果downstream只有一条，则单独计算
            {
                std::queue<Node> sources;
                uint8_t dir0 = flowdirs.at(row, col);
                int tempRow0 = row;
                int tempCol0 = col;
                int k0 = -1;
                bool visited0 = false;
                bool isPushed0 = true;
                while (BitCount2(dir0, k0) == 1) {
                    visited0 = true;
                    int iRow = dem.getRow(k0, tempRow0);
                    int iCol = dem.getCol(k0, tempCol0);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol) || isOuterEdges(iRow, iCol, dem)) {  // downstream cell 没有在这个范围内
                        isPushed0 = false;
                        break;
                    }
                    accum.at(iRow, iCol) += accum.at(tempRow0, tempCol0);
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        dependencies[iRow][iCol]--;
                        dir0 = flowdirs.at(iRow, iCol);
                        tempRow0 = iRow;
                        tempCol0 = iCol;
                        k0 = -1;
                    }
                    else {
                        isPushed0 = false;
                        break;
                    }
                }
                if (visited0) {
                    if (isPushed0) {
                        Node tempNode(tempRow0, tempCol0, dem.at(tempRow0, tempCol0));
                        sources.push(tempNode);
                    }
                }
                else {
                    Node node(row, col, dem.at(row, col));
                    sources.push(node);
                }
                while (!sources.empty()) {
                    Node tempNode = sources.front();
                    sources.pop();
                    int tempRow = tempNode.row;
                    int tempCol = tempNode.col;
                    dependencies[tempRow][tempCol]--;
                    if (dem.isNoData(tempRow, tempCol))
                        continue;
                    std::vector<double> slopes(8, 0);
                    double fracSum = 0;
                    uint8_t dir = flowdirs.at(tempRow, tempCol);
                    int k0 = -1;
                    bool visited = false;
                    bool isPushed = true;
                    while (BitCount2(dir, k0) == 1) {
                        visited = true;
                        int iRow = dem.getRow(k0, tempRow);
                        int iCol = dem.getCol(k0, tempCol);
                        if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol) || isOuterEdges(iRow, iCol, dem)) {  // downstream cell 没有在这个范围内
                            isPushed = false;
                            break;
                        }
                        accum.at(iRow, iCol) += accum.at(tempRow, tempCol);
                        dependencies[iRow][iCol]--;
                        if (dependencies[iRow][iCol] == 0) {
                            dependencies[iRow][iCol]--;
                            dir = flowdirs.at(iRow, iCol);
                            tempRow = iRow;
                            tempCol = iCol;
                            k0 = -1;
                        }
                        else {
                            isPushed = false;
                            break;
                        }
                    }
                    if (visited) {
                        if (isPushed) {
                            Node tempNode(tempRow, tempCol, dem.at(tempRow, tempCol));
                            sources.push(tempNode);
                        }
                        continue;
                    }
                    for (int k = 0; k < 8; k++) {

                        if (dir & (1 << k)) {
                            int iRow, iCol;
                            iRow = dem.getRow(k, tempRow);
                            iCol = dem.getCol(k, tempCol);
                            if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                                continue;
                            dSlope = (dem.at(tempRow, tempCol) - dem.at(iRow, iCol)) / cellSize;
                            if (k % 2 == 1) {
                                dSlope = dSlope / SQRT2;
                                slopes[k] = pow(dSlope, index[tempRow][tempCol]) * SQRT2 / 4;
                            }
                            else {
                                slopes[k] = pow(dSlope, index[tempRow][tempCol]) / 2;
                            }
                            fracSum += slopes[k];
                        }
                    }
                    if (fracSum == 0) {
                        int count = 0;
                        for (int k = 0; k < 8; k++) {
                            if (dir & (1 << k)) {
                                count++;
                            }
                        }
                        double temp = accum.at(tempRow, tempCol) / count;
                        for (int k = 0; k < 8; k++) {
                            if (dir & (1 << k)) {
                                int iRow, iCol;
                                iRow = dem.getRow(k, tempRow);
                                iCol = dem.getCol(k, tempCol);
                                if (isOuterEdges(iRow, iCol, dem))
                                    continue;
                                accum.at(iRow, iCol) += temp;
                                dependencies[iRow][iCol]--;
                                if (dependencies[iRow][iCol] == 0) {
                                    Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                                    sources.push(tempNode);
                                }
                            }
                        }
                    }
                    else {
                        for (int k = 0; k < 8; k++) {
                            if (slopes[k] == 0)
                                continue;
                            int iRow, iCol;
                            iRow = dem.getRow(k, tempRow);
                            iCol = dem.getCol(k, tempCol);
                            if (isOuterEdges(iRow, iCol, dem))
                                continue;
                            double temp = slopes[k] / fracSum;
                            accum.at(iRow, iCol) += accum.at(tempRow, tempCol) * temp;
                            dependencies[iRow][iCol]--;
                            if (dependencies[iRow][iCol] == 0) {
                                Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                                sources.push(tempNode);
                            }
                        }
                    }
                }
            }
        }
    }
}

bool is2RowCol(int row, int col, int i, int iRow, int iCol, uint8_t iDir) {
    std::vector<int> inverse = { 4, 5, 6, 7, 0, 1, 2, 3 };
    if (iDir & (1 << inverse[i]))
        return true;
    else
        return false;
}

double CheckNeighbor(int row, int col, Raster<float>& dem, Raster<uint8_t>& flowdirs, const double& cellSize, const std::vector<std::vector<double>>& index,
                     const std::vector<std::vector<double>>& fractionSum, Raster<double>& accum) {
    if (dem.isNoData(row, col))
        return 0;
    if (accum.isNoData(row, col)) {
        accum.at(row, col) = 1;
        for (int k = 0; k < 8; k++) {
            int iRow, iCol;
            iRow = dem.getRow(k, row);
            iCol = dem.getCol(k, col);
            if (isOuterEdges(iRow, iCol, dem))
                continue;
            if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                continue;
            double temp, dSlope;
            if (is2RowCol(row, col, k, iRow, iCol, flowdirs.at(iRow, iCol))) {
                if (fractionSum[iRow][iCol] == MAX_SINGLE) {
                    uint8_t dir = flowdirs.at(iRow, iCol);
                    int count = 0;
                    for (int k = 0; k < 8; k++) {
                        if (dir & (1 << k)) {
                            count++;
                        }
                    }
                    temp = 1.0 / count;
                }
                else {
                    dSlope = (dem.at(iRow, iCol) - dem.at(row, col)) / cellSize;
                    if (k % 2 == 1) {
                        dSlope = dSlope / SQRT2;
                        temp = pow(dSlope, index[iRow][iCol]) * SQRT2 / (4 * fractionSum[iRow][iCol]);
                    }
                    else {
                        temp = pow(dSlope, index[iRow][iCol]) / (2 * fractionSum[iRow][iCol]);
                    }
                }

                if (temp > 0) {
                    double value = accum.at(row, col) + temp * CheckNeighbor(iRow, iCol, dem, flowdirs, cellSize, index, fractionSum, accum);
                    accum.at(row, col) = value;
                }
            }
        }
    }
    return accum.at(row, col);
}

void FlowAccumulationInnerRecursion(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum) {
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }

    std::vector<std::vector<double>> fractionSum(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col) || flowdirs.isNoData(row, col))
                continue;
            double dSum = 0;
            uint8_t dir = flowdirs.at(row, col);
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, row);
                    iCol = dem.getCol(k, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (k % 2 == 1) {
                        dSlope = dSlope / SQRT2;
                        dSum = dSum + pow(dSlope, index[row][col]) * SQRT2 / 4;
                    }
                    else {
                        dSum = dSum + pow(dSlope, index[row][col]) / 2;
                    }
                }
            }
            fractionSum[row][col] = dSum == 0 ? MAX_SINGLE : dSum;
        }
    }

    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (!dem.isNoData(row, col))
                CheckNeighbor(row, col, dem, flowdirs, cellSize, index, fractionSum, accum);
        }
    }
}

void FlowAccumulationTradition(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c) {
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);

    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }

    std::vector<std::vector<int>> dependencies(height, std::vector<int>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            accum.at(row, col) = 1 + getOffset1(row, col, height, width, p2c);
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dependencies[iRow][iCol]++;
                }
            }
        }
    }
    std::queue<Node> sources;
    for (int row = 1; row < height - 1; row++)
        for (int col = 1; col < width - 1; col++)
            if (dependencies[row][col] == 0 && !dem.isNoData(row, col)) {
                Node node(row, col, dem.at(row, col));
                sources.push(node);
            }

    while (!sources.empty()) {
        Node node = sources.front();
        sources.pop();
        if (dem.isNoData(node.row, node.col))
            continue;
        std::vector<double> slopes(8, 0);
        double fracSum = 0;
        if (flowdirs.isNoData(node.row, node.col))
            continue;
        uint8_t dir = flowdirs.at(node.row, node.col);
        for (int k = 0; k < 8; k++) {
            if (dir & (1 << k)) {
                int iRow, iCol;
                iRow = dem.getRow(k, node.row);
                iCol = dem.getCol(k, node.col);
                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                    continue;
                dSlope = (dem.at(node.row, node.col) - dem.at(iRow, iCol)) / cellSize;
                if (k % 2 == 1) {
                    dSlope = dSlope / SQRT2;
                    slopes[k] = pow(dSlope, index[node.row][node.col]) * SQRT2 / 4;
                }
                else
                    slopes[k] = pow(dSlope, index[node.row][node.col]) / 2;
                fracSum += slopes[k];
            }
        }

        if (fracSum == 0) {
            int count = 0;
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    count++;
                }
            }
            double temp = accum.at(node.row, node.col) / count;
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, node.row);
                    iCol = dem.getCol(k, node.col);
                    if (isOuterEdges(iRow, iCol, dem) || dem.isNoData(iRow, iCol))
                        continue;
                    accum.at(iRow, iCol) += temp;
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                        sources.push(tempNode);
                    }
                }
            }
        }
        else {
            for (int k = 0; k < 8; k++) {
                if (slopes[k] != 0) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, node.row);
                    iCol = dem.getCol(k, node.col);
                    if (isOuterEdges(iRow, iCol, dem) || dem.isNoData(iRow, iCol))
                        continue;
                    double temp = slopes[k] / fracSum;
                    accum.at(iRow, iCol) += accum.at(node.row, node.col) * temp;
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                        sources.push(tempNode);
                    }
                }
            }
        }
    }
    int tempRow, tempCol;
    for (size_t i = 0; i < p2c->outer_accum.size(); i++) {
        serialToXYOuter(i, tempRow, tempCol, width, height);
        accum.at(tempRow, tempCol) = p2c->outer_accum[i];
    }
}

void FlowAccumulationImprove(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c) {
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;
    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }
    std::vector<std::vector<int>> dependencies(height, std::vector<int>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            accum.at(row, col) = 1 + getOffset1(row, col, height, width, p2c);
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dependencies[iRow][iCol]++;
                }
            }
        }
    }

    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (!dem.isNoData(row, col) && dependencies[row][col] == 0) {
                std::queue<Node> sources;
                uint8_t dir0 = flowdirs.at(row, col);
                int tempRow0 = row;
                int tempCol0 = col;
                int k0 = -1;
                bool visited0 = false;
                bool isPushed0 = true;
                while (BitCount2(dir0, k0) == 1) {
                    visited0 = true;
                    int iRow = dem.getRow(k0, tempRow0);
                    int iCol = dem.getCol(k0, tempCol0);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol) || isOuterEdges(iRow, iCol, dem)) {
                        isPushed0 = false;
                        break;
                    }
                    accum.at(iRow, iCol) += accum.at(tempRow0, tempCol0);
                    dependencies[iRow][iCol]--;
                    if (dependencies[iRow][iCol] == 0) {
                        dependencies[iRow][iCol]--;
                        dir0 = flowdirs.at(iRow, iCol);
                        tempRow0 = iRow;
                        tempCol0 = iCol;
                        k0 = -1;
                    }
                    else {
                        isPushed0 = false;
                        break;
                    }
                }
                if (visited0) {
                    if (isPushed0) {
                        Node tempNode(tempRow0, tempCol0, dem.at(tempRow0, tempCol0));
                        sources.push(tempNode);
                    }
                }
                else {
                    Node node(row, col, dem.at(row, col));
                    sources.push(node);
                }
                while (!sources.empty()) {
                    Node tempNode = sources.front();
                    sources.pop();
                    int tempRow = tempNode.row;
                    int tempCol = tempNode.col;
                    dependencies[tempRow][tempCol]--;
                    if (dem.isNoData(tempRow, tempCol))
                        continue;
                    std::vector<double> slopes(8, 0);
                    double fracSum = 0;
                    uint8_t dir = flowdirs.at(tempRow, tempCol);
                    int k0 = -1;
                    bool visited = false;
                    bool isPushed = true;
                    while (BitCount2(dir, k0) == 1) {
                        visited = true;
                        int iRow = dem.getRow(k0, tempRow);
                        int iCol = dem.getCol(k0, tempCol);
                        if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol) || isOuterEdges(iRow, iCol, dem)) {  // downstream cell 没有在这个范围内
                            isPushed = false;
                            break;
                        }
                        accum.at(iRow, iCol) += accum.at(tempRow, tempCol);
                        dependencies[iRow][iCol]--;
                        if (dependencies[iRow][iCol] == 0) {
                            dependencies[iRow][iCol]--;
                            dir = flowdirs.at(iRow, iCol);
                            tempRow = iRow;
                            tempCol = iCol;
                            k0 = -1;
                        }
                        else {
                            isPushed = false;
                            break;
                        }
                    }
                    if (visited) {
                        if (isPushed) {
                            Node tempNode(tempRow, tempCol, dem.at(tempRow, tempCol));
                            sources.push(tempNode);
                        }
                        continue;
                    }
                    for (int k = 0; k < 8; k++) {

                        if (dir & (1 << k)) {
                            int iRow, iCol;
                            iRow = dem.getRow(k, tempRow);
                            iCol = dem.getCol(k, tempCol);
                            if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                                continue;
                            dSlope = (dem.at(tempRow, tempCol) - dem.at(iRow, iCol)) / cellSize;
                            if (k % 2 == 1) {
                                dSlope = dSlope / SQRT2;
                                slopes[k] = pow(dSlope, index[tempRow][tempCol]) * SQRT2 / 4;
                            }
                            else {
                                slopes[k] = pow(dSlope, index[tempRow][tempCol]) / 2;
                            }
                            fracSum += slopes[k];
                        }
                    }
                    if (fracSum == 0) {
                        int count = 0;
                        for (int k = 0; k < 8; k++) {
                            if (dir & (1 << k)) {
                                count++;
                            }
                        }
                        double temp = accum.at(tempRow, tempCol) / count;
                        for (int k = 0; k < 8; k++) {
                            if (dir & (1 << k)) {
                                int iRow, iCol;
                                iRow = dem.getRow(k, tempRow);
                                iCol = dem.getCol(k, tempCol);
                                if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol) || isOuterEdges(iRow, iCol, dem))
                                    continue;
                                accum.at(iRow, iCol) += temp;
                                dependencies[iRow][iCol]--;
                                if (dependencies[iRow][iCol] == 0) {
                                    Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                                    sources.push(tempNode);
                                }
                            }
                        }
                    }
                    else {
                        for (int k = 0; k < 8; k++) {
                            if (slopes[k] == 0)
                                continue;
                            int iRow, iCol;
                            iRow = dem.getRow(k, tempRow);
                            iCol = dem.getCol(k, tempCol);
                            if (isOuterEdges(iRow, iCol, dem) || dem.isNoData(iRow, iCol))
                                continue;
                            double temp = slopes[k] / fracSum;
                            accum.at(iRow, iCol) += accum.at(tempRow, tempCol) * temp;
                            dependencies[iRow][iCol]--;
                            if (dependencies[iRow][iCol] == 0) {
                                Node tempNode(iRow, iCol, dem.at(iRow, iCol));
                                sources.push(tempNode);
                            }
                        }
                    }
                }
            }
        }
    }
    int tempRow, tempCol;
    for (size_t i = 0; i < p2c->outer_accum.size(); i++) {
        serialToXYOuter(i, tempRow, tempCol, width, height);
        accum.at(tempRow, tempCol) = p2c->outer_accum[i];
    }
}

double CheckNeighbor(int row, int col, Raster<float>& dem, Raster<uint8_t>& flowDir, const double& cellSize, const std::vector<std::vector<double>>& index,
                     const std::vector<std::vector<double>>& fractionSum, Raster<double>& accum, Producer2Consumer* p2c) {
    if (dem.isNoData(row, col))
        return 0;
    if (accum.isNoData(row, col)) {
        accum.at(row, col) = 1 + getOffset1(row, col, dem.getHeight(), dem.getWidth(), p2c);
        for (int k = 0; k < 8; k++) {
            int iRow, iCol;
            iRow = dem.getRow(k, row);
            iCol = dem.getCol(k, col);

            if (isOuterEdges(iRow, iCol, dem))
                continue;
            if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                continue;
            double temp = 0, dSlope = 0;
            if (is2RowCol(row, col, k, iRow, iCol, flowDir.at(iRow, iCol))) {
                if (fractionSum[iRow][iCol] == MAX_SINGLE) {  //平地
                    uint8_t dir = flowDir.at(iRow, iCol);
                    int count = 0;
                    for (int k = 0; k < 8; k++) {
                        if (dir & (1 << k)) {
                            count++;
                        }
                    }
                    temp = 1.0 / count;
                }
                else {
                    dSlope = (dem.at(iRow, iCol) - dem.at(row, col)) / cellSize;
                    if (k % 2 == 1) {
                        dSlope = dSlope / SQRT2;
                        temp = pow(dSlope, index[iRow][iCol]) * SQRT2 / (4 * fractionSum[iRow][iCol]);
                    }
                    else {
                        temp = pow(dSlope, index[iRow][iCol]) / (2 * fractionSum[iRow][iCol]);
                    }
                }

                if (temp > 0) {
                    double value = accum.at(row, col) + temp * CheckNeighbor(iRow, iCol, dem, flowDir, cellSize, index, fractionSum, accum, p2c);
                    accum.at(row, col) = value;
                }
            }
        }
    }
    return accum.at(row, col);
}

void FlowAccumulationRecursion(Raster<float>& dem, Raster<uint8_t>& flowdirs, Raster<double>& accum, Producer2Consumer* p2c) {
    double* geo = dem.getGeoTransformsPtr();
    double cellSize = abs(geo[1]);
    int height = flowdirs.getHeight();
    int width = flowdirs.getWidth();
    accum.init(height, width);
    accum.setAllValues(accum.NoDataValue);
    double dSlope;
    double a = 8.9;
    double b = 1.1;
    double slopeMax = 1;
    double slopeMin = 0;

    std::vector<std::vector<double>> index(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col))
                continue;
            double dMax = 0;
            uint8_t dirs = flowdirs.at(row, col);
            for (int i = 0; i < 8; i++) {
                if (dirs & (1 << i)) {
                    int iRow, iCol;
                    iRow = dem.getRow(i, row);
                    iCol = dem.getCol(i, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (i % 2 == 1)
                        dSlope = dSlope / SQRT2;
                    if (dSlope > dMax)
                        dMax = dSlope;
                }
            }
            if (dMax > slopeMax)
                index[row][col] = a + b;
            else if (dMax <= slopeMin)
                index[row][col] = b;
            else
                index[row][col] = a * dMax + b;
        }
    }

    std::vector<std::vector<double>> fractionSum(height, std::vector<double>(width, 0));
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (dem.isNoData(row, col) || flowdirs.isNoData(row, col))
                continue;
            double dSum = 0;
            uint8_t dir = flowdirs.at(row, col);
            for (int k = 0; k < 8; k++) {
                if (dir & (1 << k)) {
                    int iRow, iCol;
                    iRow = dem.getRow(k, row);
                    iCol = dem.getCol(k, col);
                    if (!dem.isInGrid(iRow, iCol) || dem.isNoData(iRow, iCol))
                        continue;
                    dSlope = (dem.at(row, col) - dem.at(iRow, iCol)) / cellSize;
                    if (k % 2 == 1) {
                        dSlope = dSlope / SQRT2;
                        dSum = dSum + pow(dSlope, index[row][col]) * SQRT2 / 4;
                    }
                    else {
                        dSum = dSum + pow(dSlope, index[row][col]) / 2;
                    }
                }
            }
            fractionSum[row][col] = dSum == 0 ? MAX_SINGLE : dSum;
        }
    }

    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            if (!dem.isNoData(row, col))
                CheckNeighbor(row, col, dem, flowdirs, cellSize, index, fractionSum, accum, p2c);
        }
    }

    int tempRow, tempCol;
    for (size_t i = 0; i < p2c->outer_accum.size(); i++) {
        serialToXYOuter(i, tempRow, tempCol, width, height);
        accum.at(tempRow, tempCol) = p2c->outer_accum[i];
    }
}
