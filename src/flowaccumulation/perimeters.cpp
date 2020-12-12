#include "perimeters.h"

#include <assert.h>

int xyToSerialOuter(const int row, const int col, const int width, const int height) {
    assert((col == 0 || col == width - 1 || row == 0 || row == height - 1) && col >= 0 && row >= 0 && col < width && row < height);

    if (row == 0)
        return col;

    if (col == width - 1)
        return (width - 1) + row;

    if (row == height - 1)
        return (width - 1) + (height) + col;

    return 2 * (width - 1) + (height - 1) + row;
}

void serialToXYOuter(const int serial, int& row, int& col, const int width, const int height) {
    if (serial < width) {
        col = serial;
        row = 0;
    }
    else if (serial < (width - 1) + height) {
        col = width - 1;
        row = serial - (width - 1);
    }
    else if (serial < 2 * (width - 1) + (height)) {
        col = serial - (width - 1) - (height - 1) - 1;
        row = height - 1;
    }
    else {
        col = 0;
        row = serial - 2 * (width - 1) - (height - 1);
    }

    assert((col == 0 || col == width - 1 || row == 0 || row == height - 1) && col >= 0 && row >= 0 && col < width && row < height);
}
int xyToSerialInner(const int row, const int col, const int width, const int height) {
    assert((col == 1 || col == width - 2 || row == 1 || row == height - 2) && col > 0 && row > 0 && col < width - 1 && row < height - 1);

    if (row == 1)
        return col - 1;

    if (col == width - 2)
        return (width - 3) + row - 1;

    if (row == height - 2)
        return (width - 3) + (height - 2) + (col - 1);

    return 2 * (width - 3) + (height - 3) + row - 1;
}

void serialToXYInner(const int serial, int& row, int& col, const int width, const int height) {
    if (serial < width - 2) {
        col = serial + 1;
        row = 1;
    }
    else if (serial < (width - 3) + height - 2) {
        col = width - 2;
        row = serial - (width - 4);
    }
    else if (serial < 2 * (width - 3) + (height - 2)) {
        col = serial - (width - 4) - (height - 3) - 1;
        row = height - 2;
    }
    else {
        col = 1;
        row = serial - 2 * (width - 4) - (height - 2);
    }
    assert((col == 1 || col == width - 2 || row == 1 || row == height - 2) && col > 0 && row > 0 && col < width - 1 && row < height - 1);
}

int xyToSerialAll(const int row, const int col, const int width, const int height) {
    assert(((col == 0 || col == width - 1 || row == 0 || row == height - 1) && col >= 0 && row >= 0 && col < width && row < height)
           || ((col == 1 || col == width - 2 || row == 1 || row == height - 2) && col > 0 && row > 0 && col < width - 1 && row < height - 1));
    if ((col == 0 || col == width - 1 || row == 0 || row == height - 1) && col >= 0 && row >= 0 && col < width && row < height) {
        return xyToSerialOuter(row, col, width, height);
    }
    else {
        return xyToSerialInner(row, col, width, height) + 2 * width + 2 * height - 4;
    }
}
void serialToXYAll(const int serial, int& row, int& col, const int width, const int height) {
    int mid = 2 * width + 2 * height - 4;
    if (serial < mid) {
        serialToXYOuter(serial, row, col, width, height);
    }
    else {
        serialToXYInner(serial - mid, row, col, width, height);
    }
}
