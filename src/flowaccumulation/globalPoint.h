#ifndef GLOBALPOINT_HEAD_H
#define GLOBALPOINT_HEAD_H

class GlobalPoint {
public:
    int row, col;
    int gridRow, gridCol;
    double elevation;
    double fraction;

    GlobalPoint() {
        gridRow   = -1;
        gridCol   = -1;
        row       = 0;
        col       = 0;
        elevation = 0;
        fraction  = 0;
    }

    GlobalPoint( int row, int col, int gridRow, int gridCol, double elevation, double fraction ) {
        this->row       = row;
        this->col       = col;
        this->gridRow   = gridRow;
        this->gridCol   = gridCol;
        this->elevation = elevation;
        this->fraction  = fraction;
    }

    GlobalPoint( int row, int col, int gridRow, int gridCol, double elevation ) {
        this->row       = row;
        this->col       = col;
        this->gridRow   = gridRow;
        this->gridCol   = gridCol;
        this->elevation = elevation;
        this->fraction  = 0;
    }

    GlobalPoint( int row, int col, int gridRow, int gridCol ) {
        this->row     = row;
        this->col     = col;
        this->gridRow = gridRow;
        this->gridCol = gridCol;
    }
};

#endif