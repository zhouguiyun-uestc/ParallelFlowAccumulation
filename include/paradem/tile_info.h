#include<cereal/archives/binary.hpp>

#ifndef PARADEM_TILECOORD_H
#define PARADEM_TILECOORD_H

class TileInfo
{
private:
	friend class cereal::access;
	template<class Archive>
	void serialize(Archive &ar)
	{
		ar(gridRow, gridCol, height, width);
	}
public:
	int gridRow, gridCol; // tile coordinates in the Grid
	int height, width; //number of rows and cols in the tile, may be different from the standard tile sizes for border tiles
public:
	TileInfo()=default;
	~TileInfo()=default;
};

#endif

