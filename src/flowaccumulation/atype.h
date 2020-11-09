#pragma once

typedef int link_t;
class atype {
public:
    int gRow;  // X-coordinate of tile
    int gCol;  // Y-coordinate of tile
    link_t s;  // Serialized coordinate of cell's position on the tile's perimeter

    atype( int gRow, int gCol, link_t s ) {
        this->gRow = gRow;
        this->gCol = gCol;
        this->s    = s;
    }
    atype() {
        gRow = -1;
        gCol = -1;
        s    = -1;
    }

    bool atype::isNull() {
        if ( gRow == -1 && gCol == -1 && s == -1 )
            return true;
        return false;
    }
};