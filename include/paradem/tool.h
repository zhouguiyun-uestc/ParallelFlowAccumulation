#ifndef PARADEM_TOOL_H
#define PARADEM_TOOL_H

#include <paradem/grid_info.h>
#include <paradem/i_object_factory.h>
#include <paradem/tile_info.h>

#include <vector>

bool generateTiles(const char* filePath, int tileHeight, int tileWidth, const char* outputFolder, const std::string filename, const std::string dojob);
bool readGridInfo(const char* tileFolder, GridInfo& gridInfo);
void createTileInfoArray(GridInfo& gridInfo, std::vector<TileInfo>& tileInfos);
bool mergeTiles(GridInfo& gridInfo);
bool mergeDEM(GridInfo& gridInfo);
void createNewTif(const std::string path);

#endif