#pragma once
#include <paradem/raster.h>

#include <assert.h>
#include <vector>

int xyToSerialOuter(const int row, const int col, const int width, const int height);
void serialToXYOuter(const int serial, int& row, int& col, const int width, const int height);
int xyToSerialInner(const int row, const int col, const int width, const int height);
void serialToXYInner(const int serial, int& row, int& col, const int width, const int height);
void serialToXYAll(const int serial, int& row, int& col, const int width, const int height);
int xyToSerialAll(const int row, const int col, const int width, const int height);

template <class T> void GridInnerPerimToArray(Raster<T>& grid, std::vector<T>& vec) {
    assert(vec.size() == 0);  // Ensure receiving array is empty

    std::vector<T> vec2copy;

    vec2copy = grid.getRowData(1);  // Top
    vec.insert(vec.end(), vec2copy.begin() + 1, vec2copy.end() - 1);

    vec2copy = grid.getColData(grid.getWidth() - 2);  // Right
    vec.insert(vec.end(), vec2copy.begin() + 2, vec2copy.end() - 1);

    vec2copy = grid.getRowData(grid.getHeight() - 2);  // Bottom
    vec.insert(vec.end(), vec2copy.begin() + 1, vec2copy.end() - 2);

    vec2copy = grid.getColData(1);  // Left
    vec.insert(vec.end(), vec2copy.begin() + 2, vec2copy.end() - 2);
}

template <class T> void GridOuterPerimToArray(Raster<T>& grid, std::vector<T>& vec) {
    assert(vec.size() == 0);  // Ensure receiving array is empty

    std::vector<T> vec2copy;

    vec2copy = grid.getRowData(0);  // Top
    vec.insert(vec.end(), vec2copy.begin(), vec2copy.end());

    vec2copy = grid.getColData(grid.getWidth() - 1);  // Right
    vec.insert(vec.end(), vec2copy.begin() + 1, vec2copy.end());

    vec2copy = grid.getRowData(grid.getHeight() - 1);  // Bottom
    vec.insert(vec.end(), vec2copy.begin(), vec2copy.end() - 1);

    vec2copy = grid.getColData(0);  // Left
    vec.insert(vec.end(), vec2copy.begin() + 1, vec2copy.end() - 1);
}

template <class T> void GridAllPerimToArray(Raster<T>& grid, std::vector<T>& vec) {
    std::vector<T> outerVec;
    GridOuterPerimToArray(grid, outerVec);
    std::vector<T> innerVec;
    GridInnerPerimToArray(grid, innerVec);
    vec.insert(vec.end(), outerVec.begin(), outerVec.end());
    vec.insert(vec.end(), innerVec.begin(), innerVec.end());
}

template <class T> void GridInnerPerimToArray(std::vector<std::vector<T>>& grid, std::vector<T>& vec) {
    assert(vec.size() == 0);  // Ensure receiving array is empty

    int height = grid.size();
    int width = grid[0].size();

    std::vector<T> top_copy{ std::begin(grid[1]), std::end(grid[1]) };
    vec.insert(vec.end(), top_copy.begin() + 1, top_copy.end() - 1);

    std::vector<T> right_copy(height);
    for (int row = 0; row < height; row++)
        right_copy[row] = grid[row][width - 2];
    vec.insert(vec.end(), right_copy.begin() + 2, right_copy.end() - 1);

    std::vector<T> bot_copy{ std::begin(grid[grid.size() - 2]), std::end(grid[grid.size() - 2]) };
    vec.insert(vec.end(), bot_copy.begin() + 1, bot_copy.end() - 2);

    std::vector<T> left_copy(height);
    for (int row = 0; row < height; row++)
        left_copy[row] = grid[row][1];
    vec.insert(vec.end(), left_copy.begin() + 2, left_copy.end() - 2);
}

template <class T> void GridOuterPerimToArray(std::vector<std::vector<T>>& grid, std::vector<T>& vec) {
    assert(vec.size() == 0);

    int height = grid.size();
    int width = grid[0].size();

    std::vector<T> top_copy{ std::begin(grid[0]), std::end(grid[0]) };
    vec.insert(vec.end(), top_copy.begin(), top_copy.end());

    std::vector<T> right_copy(height);
    for (int row = 0; row < height; row++)
        right_copy[row] = grid[row][width - 1];
    vec.insert(vec.end(), right_copy.begin() + 1, right_copy.end());

    std::vector<T> bot_copy{ std::begin(grid[grid.size() - 1]), std::end(grid[grid.size() - 1]) };
    vec.insert(vec.end(), bot_copy.begin(), bot_copy.end() - 1);

    std::vector<T> left_copy(height);
    for (int row = 0; row < height; row++)
        left_copy[row] = grid[row][0];
    vec.insert(vec.end(), left_copy.begin() + 1, left_copy.end() - 1);
}

template <class T> void GridAllPerimToArray(std::vector<std::vector<T>>& grid, std::vector<T>& vec) {
    std::vector<T> outerVec;
    GridOuterPerimToArray(grid, outerVec);
    std::vector<T> innerVec;
    GridInnerPerimToArray(grid, innerVec);
    vec.insert(vec.end(), outerVec.begin(), outerVec.end());
    vec.insert(vec.end(), innerVec.begin(), innerVec.end());
}