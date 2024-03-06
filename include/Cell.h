/**
 * @file Cell.h
 * @brief The Cell class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#include "Surface.h"
#include "Point.h"
#include "LocalCoords.h"
#endif

class Universe;
class Surface;
class LocalCoords;

int cell_id();


/**
 * @enum cellType
 * @brief The type of cell.
*/
enum cellType {

    /** A cell filled by a material */
    MATERIAL,

    /** A cell filled by a universe */
    FILL
};

/**
 * @class Cell Cell.h "openmoc/src/host/Cell.h"
 * @brief Represents a cell inside of a universe.
 */
class Cell {

protected:

    /** A static counter for the number of cells in a simulation */
    static int _n;

    /** A static counter for the number of cell clones in a simulation */
    static int _clone_n;

    /** A monotonically increasing unique ID for each cell created */
    int _uid;

    /** A user-defined ID for each cell created */
    int _id;

    /** The type of cell (ie MATERIAL or FILL) */
    cellType _cell_type;

    /** The ID for the universe within which this cell resides */
    int _universe;

    /** Map of bounding surface ids (keys) to surface pointers (values). 
     *  Keys are +/- depending on the side of surface in which this cell
     *  resides. */
    std::map<int, Surface*> _surfaces;

public:
    Cell();
    Cell(int universe, int id=0);
    virtual ~Cell();
    int getUid() const;
    int getId() const;
    cellType getType() const;
    int getUniverse() const;
    int getNumSurfaces() const;
    std::map<int, Surface*> getSurfaces() const;

    /**
     * @brief Return the number of flat source regions in this cell. 
     * @details This method is used when the geometry recursively constructs 
     *          flat source regions.
     * @return the number of FSRs in this cell
     */
    virtual int getNumFSRs() =0;

    void setUniverse(int universe);
    void addSurface(int halfspace, Surface* surface);
    void setSurfacePointer(Surface* surface);

    bool cellContainsPoint(Point* point);
    bool cellContainsCoords(LocalCoords* coords);
    double minSurfaceDist(Point* point, double angle, 
                          Point* min_intersection);
    /**
     * @brief Convert this cellfill's attributes to a string format.
     * @return a character array of this cell's attributes
     */
    virtual std::string toString() =0;

    /**
     * @brief Prints a string representation of all of the cells's objects to
     *        the console.
     */
    virtual void printString() =0;
};


/**
 * @class CellBasic Cell.h "openmoc/src/host/Cell.h"
 * @brief Represents a cell filled with a material.
 */
class CellBasic: public Cell {

private: 

    /** A pointer to the material filling this cell */
    int _material;

    /** The number of rings sub-dividing this cell */
    int _num_rings;

    /** The number of sectors sub-dividing this cell */
    int _num_sectors;

    std::vector<CellBasic*> _rings;
    std::vector<CellBasic*> _sectors;
    std::vector<CellBasic*> _subcells;

    void ringify();
    void sectorize();

public:
    CellBasic(int universe, int material,
	      int rings=0, int sectors=0, int id=0);

    int getMaterial() const;
    int getNumRings();
    int getNumSectors();
    int getNumFSRs();

    void setNumRings(int num_rings);
    void setNumSectors(int num_sectors);
    CellBasic* clone();
    std::vector<CellBasic*> subdivideCell();

    std::string toString();
    void printString();
};


/**
 * @brief Represents a cell filled with a universe.
 * @param _universe_fill id for the universe filling this cell
 */
class CellFill: public Cell {

private:
    /** A pair of the ID and universe pointer filling this cell */
    std::pair<int, Universe*> _universe_fill;

public:
    CellFill(int universe, int universe_fill, int id=0);

    int getUniverseFillId() const;
    Universe* getUniverseFill() const;
    int getNumFSRs();

    void setUniverseFill(int universe_Fill);
    void setUniverseFillPointer(Universe* universe_fill);

    std::string toString();
    void printString();
};

#endif /* CELL_H_ */
