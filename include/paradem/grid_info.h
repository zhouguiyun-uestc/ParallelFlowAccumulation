#ifndef PARADEM_GRIDINFO_H
#define PARADEM_GRIDINFO_H

#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <string>

class GridInfo {
private:
    friend class cereal::access;
    template <class Archive> void serialize(Archive& ar) {
        ar(tileHeight, tileWidth, gridHeight, gridWidth, grandHeight, grandWidth, inputFolder, outputFolder, cellSize);
    }

public:
    int tileHeight, tileWidth;                                  // tile dimensions
    int gridHeight, gridWidth;                                  // number of tiles along rows and cols
    int grandHeight, grandWidth;                                // total rows and cols in the grid
    std::string inputFolder, outputFolder, intermediateFolder;  // three folders
    double cellSize;

public:
    GridInfo() = default;
    ~GridInfo() = default;
};

#endif
