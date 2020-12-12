#include "producer.h"
#include "perimeters.h"

#include <paradem/object_deleter.h>

#include <iomanip>
#include <iostream>
#include <queue>

const int dRow[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
const int dCol[8] = { 1, 1, 0, -1, -1, -1, 0, 1 };
const double Sqrt2 = 1.4142135623730950488016887242097;
const double a = 8.9;
const double b = 1.1;
const double slopeMax = 1;
const double slopeMin = 0;
const double NODATA = -9999;

int getRow(int dir, int row) {
    return row + dRow[dir];
}

int getCol(int dir, int col) {
    return (col + dCol[dir]);
}

bool isInGrid(int row, int col, int height, int width) {
    if (row < 0 || row >= height || col < 0 || col >= width)
        return false;
    return true;
}

atype getINode(const int k, const int row, const int col, const GridInfo& gridInfo, const TileInfo& tileInfo, const std::vector<TileInfo>& tileInfos) {
    atype node;
    node.gCol = -1;
    node.gRow = -1;
    node.s = -1;
    int width = tileInfo.width;
    int height = tileInfo.height;

    int iGridRow, iGridCol;
    int iRow = -1, iCol = -1;
    int n;

    if (row == 0 && col == 0) {
        switch (k) {
        case 3: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row + 1;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 4: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 5: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 6: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = col;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 7: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = col + 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        default:
            break;
        }
    }
    else if (row == 0 && col == width - 1) {
        switch (k) {
        case 0: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 1: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row + 1;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 5: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = col - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 6: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = col;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 7: {
            iGridRow = tileInfo.gridRow - 1;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = iTileInfo.height - 1;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        default:
            break;
        }
    }
    else if (row == height - 1 && col == 0) {
        switch (k) {
        case 1: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = col + 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 2: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 3: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 4: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 5: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol - 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row - 1;
            iCol = iTileInfo.width - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        default:
            break;
        }
    }
    else if (row == height - 1 && col == width - 1) {
        switch (k) {
        case 0: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 1: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 2: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = col;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 3: {
            iGridRow = tileInfo.gridRow + 1;
            iGridCol = tileInfo.gridCol;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = 0;
            iCol = col - 1;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        case 7: {
            iGridRow = tileInfo.gridRow;
            iGridCol = tileInfo.gridCol + 1;
            if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
                return node;
            TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
            iRow = row - 1;
            iCol = 0;
            n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
            node.gRow = iGridRow;
            node.gCol = iGridCol;
            node.s = n;
            break;
        }
        default:
            break;
        }
    }
    else if (row == 0) {
        iGridRow = tileInfo.gridRow - 1;
        iGridCol = tileInfo.gridCol;
        if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
            return node;
        TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
        iRow = iTileInfo.height - 1;
        switch (k) {
        case 5:
            iCol = col - 1;
            break;
        case 6:
            iCol = col;
            break;
        case 7:
            iCol = col + 1;
            break;
        default:
            break;
        }
        n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
        node.gRow = iGridRow;
        node.gCol = iGridCol;
        node.s = n;
    }
    else if (row == height - 1) {
        iGridRow = tileInfo.gridRow + 1;
        iGridCol = tileInfo.gridCol;
        if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
            return node;
        TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
        iRow = 0;
        switch (k) {
        case 1:
            iCol = col + 1;
            break;
        case 2:
            iCol = col;
            break;
        case 3:
            iCol = col - 1;
            break;
        default:
            break;
        }
        n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
        node.gRow = iGridRow;
        node.gCol = iGridCol;
        node.s = n;
    }
    else if (col == 0) {
        iGridRow = tileInfo.gridRow;
        iGridCol = tileInfo.gridCol - 1;
        if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
            return node;
        TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
        iCol = iTileInfo.width - 1;
        switch (k) {
        case 3:
            iRow = row + 1;
            break;
        case 4:
            iRow = row;
            break;
        case 5:
            iRow = row - 1;
            break;
        default:
            break;
        }
        n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
        node.gRow = iGridRow;
        node.gCol = iGridCol;
        node.s = n;
    }
    else if (col == width - 1) {
        iGridRow = tileInfo.gridRow;
        iGridCol = tileInfo.gridCol + 1;
        if (!isInGrid(iGridRow, iGridCol, gridInfo.gridHeight, gridInfo.gridWidth))
            return node;
        TileInfo iTileInfo = tileInfos[iGridRow * gridInfo.gridWidth + iGridCol];
        iCol = 0;
        switch (k) {
        case 7:
            iRow = row - 1;
            break;
        case 0:
            iRow = row;
            break;
        case 1:
            iRow = row + 1;
            break;
        default:
            break;
        }
        n = xyToSerialOuter(iRow, iCol, iTileInfo.width, iTileInfo.height);
        node.gRow = iGridRow;
        node.gCol = iGridCol;
        node.s = n;
    }

    return node;
}

void get_GlobalLink(Consumer2Producer& c2p, const GridInfo& gridInfo, const TileInfo& tileInfo, const std::vector<TileInfo>& tileInfos) {
    int width = tileInfo.width;
    int height = tileInfo.height;
    int mid = width * 2 + height * 2 - 4;
    c2p.outerGlobal_links.resize(mid);
    for (int i = 0; i < mid; i++) {
        elev_t elev = c2p.elevation[i];
        uint8_t dir = c2p.flowdirs[i];
        if (dir == 0 || elev == NODATA)
            continue;
        int row, col;
        serialToXYOuter(i, row, col, width, height);

        for (size_t k = 0; k < 8; k++) {
            if (dir & (1 << k)) {
                int iRow = getRow(k, row);
                int iCol = getCol(k, col);
                if (!isInGrid(iRow, iCol, height, width)) {
                    atype node = getINode(k, row, col, gridInfo, tileInfo, tileInfos);
                    if (!node.isNull())
                        c2p.outerGlobal_links[i].push_back(node);
                }
            }
        }
    }
}

void get_dependence(Consumer2Producer& c2p, Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) {
    auto& links = c2p.links;
    for (size_t i = 0; i < links.size(); i++) {
        for (size_t j = 0; j < links[i].size(); j++) {
            c2p.dependencies[links[i][j]]++;
        }
    }
    auto& globalLinks = c2p.outerGlobal_links;
    for (size_t i = 0; i < globalLinks.size(); i++) {
        for (size_t j = 0; j < globalLinks[i].size(); j++) {
            atype node = globalLinks[i][j];
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(node.gRow, node.gCol).get();
            c2p.dependencies[node.s]++;
        }
    }
}

void calculationFra(elev_t ele, double cellSize, const std::vector<elev_t>& link_elevs, const std::vector<elev_t>& globalLink_elevs, const std::vector<bool>& link_isDiagonal,
                    const std::vector<bool>& globalLink_isDiagonal, std::vector<double>& innerSlope, std::vector<double>& outerSlope) {
    double dMAX = 0;
    for (size_t m = 0; m < link_elevs.size(); m++) {
        innerSlope[m] = (ele - link_elevs[m]) / cellSize;
        if (link_isDiagonal[m])
            innerSlope[m] = innerSlope[m] / Sqrt2;
        dMAX = (dMAX > innerSlope[m] ? dMAX : innerSlope[m]);
    }
    for (size_t m = 0; m < globalLink_elevs.size(); m++) {
        outerSlope[m] = (ele - globalLink_elevs[m]) / cellSize;
        if (globalLink_isDiagonal[m])
            outerSlope[m] = outerSlope[m] / Sqrt2;
        dMAX = (dMAX > outerSlope[m] ? dMAX : outerSlope[m]);
    }
    double p;
    if (dMAX > slopeMax)
        p = a + b;
    else if (dMAX <= slopeMin)
        p = b;
    else
        p = a * dMAX + b;
    double dSum = 0;
    for (size_t j = 0; j < innerSlope.size(); j++) {
        if (link_isDiagonal[j])
            innerSlope[j] = pow(innerSlope[j], p) * Sqrt2 / 4;
        else
            innerSlope[j] = pow(innerSlope[j], p) / 2;
        dSum += innerSlope[j];
    }
    for (size_t j = 0; j < outerSlope.size(); j++) {
        if (globalLink_isDiagonal[j])
            outerSlope[j] = pow(outerSlope[j], p) * Sqrt2 / 4;
        else
            outerSlope[j] = pow(outerSlope[j], p) / 2;
        dSum += outerSlope[j];
    }
    if (dSum == 0) {
        double fraction = 1.0 / (link_elevs.size() + globalLink_elevs.size());
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] = fraction;
        for (size_t j = 0; j < outerSlope.size(); j++)
            outerSlope[j] = fraction;
    }
    else {
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] /= dSum;
        for (size_t j = 0; j < outerSlope.size(); j++)
            outerSlope[j] /= dSum;
    }
}

void getFraction(Consumer2Producer& c2p, const GridInfo& gridInfo, const TileInfo& tileInfo, const std::vector<TileInfo>& tileInfos, const double cellSize,
                 Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) {
    int mid = c2p.outerGlobal_links.size();

    for (int i = 0; i < mid; i++) {
        elev_t ele = c2p.elevation[i];
        int row, col;
        serialToXYOuter(i, row, col, tileInfo.width, tileInfo.height);
        auto& link = c2p.links[i];
        auto& globalLink = c2p.outerGlobal_links[i];
        std::vector<elev_t> link_elevs(link.size());
        std::vector<elev_t> globalLink_elevs(globalLink.size());
        for (size_t m = 0; m < link.size(); m++)
            link_elevs[m] = c2p.elevation[link[m]];
        for (size_t m = 0; m < globalLink.size(); m++) {
            auto& tempC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(globalLink[m].gRow, globalLink[m].gCol).get();
            globalLink_elevs[m] = tempC2p.elevation[globalLink[m].s];
        }
        std::vector<bool> link_isDiagonal(link.size(), true);
        std::vector<bool> globalLink_isDiagonal(globalLink.size(), true);
        int tempRow, tempCol;
        for (size_t m = 0; m < link.size(); m++) {
            serialToXYAll(link[m], tempRow, tempCol, tileInfo.width, tileInfo.height);
            if (row == tempRow || col == tempCol)
                link_isDiagonal[m] = false;
        }
        for (size_t m = 0; m < globalLink.size(); m++) {
            TileInfo tempTileInfo = tileInfos[globalLink[m].gRow * gridInfo.gridWidth + globalLink[m].gCol];
            serialToXYAll(globalLink[m].s, tempRow, tempCol, tempTileInfo.width, tempTileInfo.height);
            if (row == tempRow || col == tempCol)
                globalLink_isDiagonal[m] = false;
        }
        std::vector<double> innerSlope(link.size());
        std::vector<double> outerSlope(globalLink.size());
        calculationFra(ele, cellSize, link_elevs, globalLink_elevs, link_isDiagonal, globalLink_isDiagonal, innerSlope, outerSlope);
        c2p.fraction[i] = innerSlope;
        c2p.outerGlobal_fraction[i] = outerSlope;
    }
}

bool isDiagonal(const int row1, const int col1, const int row2, const int col2) {
    if (row1 == row2 || col1 == col2)
        return false;
    return true;
}

void calculationFra(const Node& node, const std::vector<Node>& toInner, const std::vector<GlobalPoint>& toOuter, const double cellSize, std::vector<double>& innerSlope,
                    std::vector<double>& outerSlope) {
    double elve = node.elevation;
    innerSlope.resize(toInner.size());
    std::vector<bool> innerFlag(toInner.size(), false);
    outerSlope.reserve(toOuter.size());
    std::vector<bool> outerFlag(toOuter.size(), false);
    double dMAX = 0;
    for (size_t j = 0; j < toInner.size(); j++) {
        innerSlope[j] = (elve - toInner[j].elevation) / cellSize;
        if (isDiagonal(node.row, node.col, toInner[j].row, toInner[j].col)) {
            innerSlope[j] = innerSlope[j] / Sqrt2;
            innerFlag[j] = true;
        }
        dMAX = (dMAX > innerSlope[j] ? dMAX : innerSlope[j]);
    }
    for (size_t j = 0; j < toOuter.size(); j++) {
        outerSlope[j] = (elve - toOuter[j].elevation) / cellSize;
        if (isDiagonal(node.row, node.col, toOuter[j].row, toOuter[j].col)) {
            outerSlope[j] = outerSlope[j] / Sqrt2;
            outerFlag[j] = true;
        }
        dMAX = (dMAX > outerSlope[j] ? dMAX : outerSlope[j]);
    }
    double p;
    if (dMAX > slopeMax)
        p = a + b;
    else if (dMAX <= slopeMin)
        p = b;
    else
        p = a * dMAX + b;
    double dSum = 0;
    for (size_t j = 0; j < innerSlope.size(); j++) {
        if (innerFlag[j])
            innerSlope[j] = pow(innerSlope[j], p) * Sqrt2 / 4;
        else
            innerSlope[j] = pow(innerSlope[j], p) / 2;
        dSum += innerSlope[j];
    }
    for (size_t j = 0; j < outerSlope.size(); j++) {
        if (outerFlag[j])
            outerSlope[j] = pow(outerSlope[j], p) * Sqrt2 / 4;
        else
            outerSlope[j] = pow(outerSlope[j], p) / 2;
        dSum += outerSlope[j];
    }
    if (dSum == 0) {
        double fraction = 1.0 / (toInner.size() + toOuter.size());
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] = fraction;
        for (size_t j = 0; j < outerSlope.size(); j++)
            outerSlope[j] = fraction;
    }
    else {
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] /= dSum;
        for (size_t j = 0; j < outerSlope.size(); j++)
            outerSlope[j] /= dSum;
    }
}

void calculationFra(const Node& node, const std::vector<Node>& toInner, const double cellSize, std::vector<double>& innerSlope) {
    double elve = node.elevation;
    innerSlope.resize(toInner.size());
    std::vector<bool> innerFlag(toInner.size(), false);
    double dMAX = 0;
    for (size_t j = 0; j < toInner.size(); j++) {
        innerSlope[j] = (elve - toInner[j].elevation) / cellSize;
        if (isDiagonal(node.row, node.col, toInner[j].row, toInner[j].col)) {
            innerSlope[j] = innerSlope[j] / Sqrt2;
            innerFlag[j] = true;
        }
        dMAX = (dMAX > innerSlope[j] ? dMAX : innerSlope[j]);
    }
    double p;
    if (dMAX > slopeMax)
        p = a + b;
    else if (dMAX <= slopeMin)
        p = b;
    else
        p = a * dMAX + b;
    double dSum = 0;
    for (size_t j = 0; j < innerSlope.size(); j++) {
        if (innerFlag[j])
            innerSlope[j] = pow(innerSlope[j], p) * Sqrt2 / 4;
        else
            innerSlope[j] = pow(innerSlope[j], p) / 2;
        dSum += innerSlope[j];
    }
    if (dSum == 0) {
        double fraction = 1.0 / toInner.size();
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] = fraction;
    }
    else {
        for (size_t j = 0; j < innerSlope.size(); j++)
            innerSlope[j] /= dSum;
    }
}

void Producer::process(const GridInfo& gridInfo, const std::vector<TileInfo>& tileInfos, Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) {
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    for (int gridRow = 0; gridRow < gridHeight; gridRow++) {
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();
            TileInfo tileDirInfo = tileInfos[gridRow * gridInfo.gridWidth + gridCol];
            get_GlobalLink(c2p, gridInfo, tileDirInfo, tileInfos);
        }
    }
    for (int gridRow = 0; gridRow < gridHeight; gridRow++)
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();
            c2p.dependencies.resize(c2p.elevation.size(), 0);
            c2p.outer_accum.resize(c2p.outerGlobal_links.size(), NODATA);
            for (size_t i = 0; i < c2p.outerGlobal_links.size(); i++) {
                if (c2p.elevation[i] != NODATA)
                    c2p.outer_accum[i] = 1;
            }
            c2p.inner_accum_setoff1.resize(c2p.inner_accum.size(), 0);
            c2p.inner_accum_setoff2.resize(c2p.inner_accum.size(), 0);
        }

    for (int gridRow = 0; gridRow < gridHeight; gridRow++) {
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();
            get_dependence(c2p, gridIConsumer2Producer);
        }
    }
    for (int gridRow = 0; gridRow < gridHeight; gridRow++) {
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();
            c2p.outerGlobal_fraction.resize(c2p.outerGlobal_links.size());
            TileInfo tileDirInfo = tileInfos[gridRow * gridInfo.gridWidth + gridCol];
            getFraction(c2p, gridInfo, tileDirInfo, tileInfos, gridInfo.cellSize, gridIConsumer2Producer);
        }
    }
    globalAccumOffset_improve(gridInfo, tileInfos, gridIConsumer2Producer);
    // globalAccumOffset_tradition(gridInfo, tileInfos, gridIConsumer2Producer);
}

bool isOuterEdges(int row, int col, int height, int width) {
    if (row == 0 || col == 0 || row == height - 1 || col == width - 1)
        return true;
    return false;
}

bool isInnerEdges(int row, int col, int height, int width) {
    if (row == 1 || col == 1 || row == height - 2 || col == width - 2)
        return true;
    return false;
}

int get_Mid(const TileInfo& tileInfo) {
    return tileInfo.height * 2 + tileInfo.width * 2 - 4;
}
int get_Mid(const int height, const int width) {
    return height * 2 + width * 2 - 4;
}

void Producer::globalAccumOffset_tradition(const GridInfo& gridInfo, const std::vector<TileInfo>& tileInfos, Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) {
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    std::queue<atype> q;
    for (int gridRow = 0; gridRow < gridHeight; gridRow++) {
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();
            for (size_t i = 0; i < c2p.dependencies.size(); i++) {
                if (c2p.dependencies[i] == 0) {
                    atype temp(gridRow, gridCol, i);
                    q.push(temp);
                }
            }
        }
    }
    while (!q.empty()) {
        atype gp = q.front();
        q.pop();
        auto& GPC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gp.gRow, gp.gCol).get();
        TileInfo GPTileInfo = tileInfos[gp.gRow * gridInfo.gridWidth + gp.gCol];
        double accum = 0, offset1 = 0, offset2 = 0;
        int mid = get_Mid(GPTileInfo);
        if (gp.s < mid) {
            accum = GPC2p.outer_accum[gp.s];
            if (accum == NODATA)
                continue;
            auto& link = GPC2p.links[gp.s];
            auto& frac = GPC2p.fraction[gp.s];
            auto& globalLink = GPC2p.outerGlobal_links[gp.s];
            auto& globalFrac = GPC2p.outerGlobal_fraction[gp.s];
            for (size_t k = 0; k < link.size(); k++) {
                double temp = accum * frac[k];
                if (link[k] < mid) {
                    GPC2p.outer_accum[link[k]] += temp;
                }
                else {
                    GPC2p.inner_accum_setoff1[link[k] - mid] += temp;
                }
                GPC2p.dependencies[link[k]]--;
                if (GPC2p.dependencies[link[k]] == 0) {
                    atype tempP(gp.gRow, gp.gCol, link[k]);
                    q.push(tempP);
                }
            }

            for (size_t k = 0; k < globalLink.size(); k++) {
                double temp = accum * globalFrac[k];
                auto& tempC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(globalLink[k].gRow, globalLink[k].gCol).get();
                tempC2p.outer_accum[globalLink[k].s] += temp;
                tempC2p.dependencies[globalLink[k].s]--;
                if (tempC2p.dependencies[globalLink[k].s] == 0) {
                    q.push(globalLink[k]);
                }
            }
        }
        else {
            accum = GPC2p.inner_accum[gp.s - mid];
            offset1 = GPC2p.inner_accum_setoff1[gp.s - mid];
            offset2 = GPC2p.inner_accum_setoff2[gp.s - mid];
            if (accum == NODATA)
                continue;
            auto& link = GPC2p.links[gp.s];
            auto& frac = GPC2p.fraction[gp.s];
            for (size_t k = 0; k < link.size(); k++) {
                if (link[k] < mid) {
                    double temp = (accum + offset1 + offset2) * frac[k];
                    GPC2p.outer_accum[link[k]] += temp;
                }
                else {
                    double temp = (offset1 + offset2) * frac[k];
                    GPC2p.inner_accum_setoff2[link[k] - mid] += temp;
                }
                GPC2p.dependencies[link[k]]--;
                if (GPC2p.dependencies[link[k]] == 0) {
                    atype tempP(gp.gRow, gp.gCol, link[k]);
                    q.push(tempP);
                }
            }
        }
    }
}

void Producer::globalAccumOffset_improve(const GridInfo& gridInfo, const std::vector<TileInfo>& tileInfos, Grid<std::shared_ptr<IConsumer2Producer>>& gridIConsumer2Producer) {
    int gridHeight = gridInfo.gridHeight;
    int gridWidth = gridInfo.gridWidth;
    std::queue<atype> q;
    for (int gridRow = 0; gridRow < gridHeight; gridRow++) {
        for (int gridCol = 0; gridCol < gridWidth; gridCol++) {
            auto& c2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gridRow, gridCol).get();

            for (size_t i = 0; i < c2p.dependencies.size(); i++) {
                if (c2p.dependencies[i] == 0) {
                    c2p.dependencies[i]--;
                    atype gp(gridRow, gridCol, i);
                    bool ispushed = false;
                    while (true) {
                        auto& GPC2p0 = *(Consumer2Producer*)gridIConsumer2Producer.at(gp.gRow, gp.gCol).get();
                        TileInfo GPTileInfo0 = tileInfos[gp.gRow * gridInfo.gridWidth + gp.gCol];
                        double accum0 = 0;
                        int mid0 = get_Mid(GPTileInfo0);
                        if (gp.s < mid0) {
                            accum0 = GPC2p0.outer_accum[gp.s];
                            if (accum0 == NODATA)
                                break;
                            auto& link0 = GPC2p0.links[gp.s];
                            auto& globalLink0 = GPC2p0.outerGlobal_links[gp.s];
                            if ((link0.size() == 1 && globalLink0.size() == 0) || (link0.size() == 0 && globalLink0.size() == 1)) {
                                if (link0.size() == 1) {
                                    if (link0[0] < mid0) {
                                        GPC2p0.outer_accum[link0[0]] += accum0;
                                    }
                                    else {
                                        GPC2p0.inner_accum_setoff1[link0[0] - mid0] += accum0;
                                    }
                                    GPC2p0.dependencies[link0[0]]--;
                                    if (GPC2p0.dependencies[link0[0]] == 0) {
                                        GPC2p0.dependencies[link0[0]]--;
                                        gp.s = link0[0];
                                    }
                                    else {
                                        break;
                                    }
                                }
                                else {
                                    auto& tempC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(globalLink0[0].gRow, globalLink0[0].gCol).get();
                                    tempC2p.outer_accum[globalLink0[0].s] += accum0;
                                    tempC2p.dependencies[globalLink0[0].s]--;
                                    if (tempC2p.dependencies[globalLink0[0].s] == 0) {
                                        tempC2p.dependencies[globalLink0[0].s]--;
                                        gp = globalLink0[0];
                                    }
                                    else {
                                        break;
                                    }
                                }
                            }
                            else {
                                ispushed = true;
                                break;
                            }
                        }
                        else {
                            accum0 = GPC2p0.inner_accum[gp.s - mid0];
                            int offset10 = GPC2p0.inner_accum_setoff1[gp.s - mid0];
                            int offset20 = GPC2p0.inner_accum_setoff2[gp.s - mid0];
                            if (accum0 == NODATA)
                                break;
                            auto& link0 = GPC2p0.links[gp.s];
                            if (link0.size() == 1) {
                                if (link0[0] < mid0) {
                                    double temp = (accum0 + offset10 + offset20);
                                    GPC2p0.outer_accum[link0[0]] += temp;
                                }
                                else {
                                    double temp = (offset10 + offset20);
                                    GPC2p0.inner_accum_setoff2[link0[0] - mid0] += temp;
                                }
                                GPC2p0.dependencies[link0[0]]--;
                                if (GPC2p0.dependencies[link0[0]] == 0) {
                                    GPC2p0.dependencies[link0[0]]--;
                                    atype tempP(gp.gRow, gp.gCol, link0[0]);
                                    gp = tempP;
                                }
                                else {
                                    break;
                                }
                            }
                            else {
                                ispushed = true;
                                break;
                            }
                        }
                    }

                    if (ispushed) {
                        q.push(gp);
                    }

                    while (!q.empty()) {
                        atype gp = q.front();
                        q.pop();
                        while (true) {

                            auto& GPC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(gp.gRow, gp.gCol).get();
                            TileInfo GPTileInfo = tileInfos[gp.gRow * gridInfo.gridWidth + gp.gCol];
                            double accum = 0, offset1 = 0, offset2 = 0;
                            int mid = get_Mid(GPTileInfo);
                            if (gp.s < mid) {
                                accum = GPC2p.outer_accum[gp.s];
                                if (accum == NODATA)
                                    break;
                                auto& link = GPC2p.links[gp.s];
                                auto& frac = GPC2p.fraction[gp.s];
                                auto& globalLink = GPC2p.outerGlobal_links[gp.s];
                                auto& globalFrac = GPC2p.outerGlobal_fraction[gp.s];
                                if ((link.size() == 0 && globalLink.size() == 1) || (link.size() == 1 && globalLink.size() == 0)) {
                                    if (link.size() == 1) {
                                        if (link[0] < mid) {
                                            GPC2p.outer_accum[link[0]] += accum;
                                        }
                                        else {
                                            GPC2p.inner_accum_setoff1[link[0] - mid] += accum;
                                        }
                                        GPC2p.dependencies[link[0]]--;
                                        if (GPC2p.dependencies[link[0]] == 0) {
                                            GPC2p.dependencies[link[0]]--;
                                            gp.s = link[0];
                                        }
                                        else {
                                            break;
                                        }
                                    }
                                    else {
                                        auto& tempC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(globalLink[0].gRow, globalLink[0].gCol).get();
                                        tempC2p.outer_accum[globalLink[0].s] += accum;
                                        tempC2p.dependencies[globalLink[0].s]--;
                                        if (tempC2p.dependencies[globalLink[0].s] == 0) {
                                            tempC2p.dependencies[globalLink[0].s]--;
                                            gp = globalLink[0];
                                        }
                                        else {
                                            break;
                                        }
                                    }
                                }
                                else {

                                    for (size_t k = 0; k < link.size(); k++) {
                                        double temp = accum * frac[k];
                                        if (link[k] < mid) {
                                            GPC2p.outer_accum[link[k]] += temp;
                                        }
                                        else {
                                            GPC2p.inner_accum_setoff1[link[k] - mid] += temp;
                                        }
                                        GPC2p.dependencies[link[k]]--;
                                        if (GPC2p.dependencies[link[k]] == 0) {
                                            GPC2p.dependencies[link[k]]--;
                                            atype tempP(gp.gRow, gp.gCol, link[k]);
                                            q.push(tempP);
                                        }
                                    }

                                    for (size_t k = 0; k < globalLink.size(); k++) {
                                        double temp = accum * globalFrac[k];
                                        auto& tempC2p = *(Consumer2Producer*)gridIConsumer2Producer.at(globalLink[k].gRow, globalLink[k].gCol).get();
                                        tempC2p.outer_accum[globalLink[k].s] += temp;
                                        tempC2p.dependencies[globalLink[k].s]--;
                                        if (tempC2p.dependencies[globalLink[k].s] == 0) {
                                            tempC2p.dependencies[globalLink[k].s]--;
                                            q.push(globalLink[k]);
                                        }
                                    }
                                    break;
                                }
                            }
                            else {
                                accum = GPC2p.inner_accum[gp.s - mid];
                                offset1 = GPC2p.inner_accum_setoff1[gp.s - mid];
                                offset2 = GPC2p.inner_accum_setoff2[gp.s - mid];
                                if (accum == NODATA)
                                    break;
                                auto& link = GPC2p.links[gp.s];
                                auto& frac = GPC2p.fraction[gp.s];
                                if (link.size() == 1) {
                                    if (link[0] < mid) {
                                        double temp = (accum + offset1 + offset2);
                                        GPC2p.outer_accum[link[0]] += temp;
                                    }
                                    else {
                                        double temp = (offset1 + offset2);
                                        GPC2p.inner_accum_setoff2[link[0] - mid] += temp;
                                    }
                                    GPC2p.dependencies[link[0]]--;
                                    if (GPC2p.dependencies[link[0]] == 0) {
                                        GPC2p.dependencies[link[0]]--;
                                        atype tempP(gp.gRow, gp.gCol, link[0]);
                                        gp.s = link[0];
                                    }
                                    else {
                                        break;
                                    }
                                }
                                else {

                                    for (size_t k = 0; k < link.size(); k++) {
                                        if (link[k] < mid) {
                                            double temp = (accum + offset1 + offset2) * frac[k];
                                            GPC2p.outer_accum[link[k]] += temp;
                                        }
                                        else {
                                            double temp = (offset1 + offset2) * frac[k];
                                            GPC2p.inner_accum_setoff2[link[k] - mid] += temp;
                                        }
                                        GPC2p.dependencies[link[k]]--;
                                        if (GPC2p.dependencies[link[k]] == 0) {
                                            GPC2p.dependencies[link[k]]--;
                                            atype tempP(gp.gRow, gp.gCol, link[k]);
                                            q.push(tempP);
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

std::shared_ptr<IProducer2Consumer> Producer::toConsumer(const IConsumer2Producer* ic2p) {
    Producer2Consumer* p2c = new Producer2Consumer();
    Consumer2Producer* c2p = (Consumer2Producer*)ic2p;
    p2c->inner_accum_setoff1 = c2p->inner_accum_setoff1;
    p2c->inner_accum_setoff2 = c2p->inner_accum_setoff2;
    p2c->outer_accum = c2p->outer_accum;
    return std::shared_ptr<IProducer2Consumer>(p2c, ObjectDeleter());
}

void Producer::free() {
    delete this;
}