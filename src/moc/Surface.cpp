#include "Surface.h"

int Surface::_n = 0;

static int auto_id = DEFAULT_INIT_ID;

/**
 * @brief Returns an auto-generated unique surface ID.
 * @details This method is intended as a utility mehtod for user's writing
 *          OpenMOC input files. The method makes use of a static surface
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined surface IDs greater
 *          than or equal to 10000 is prohibited.
 */
int surface_id()
{
    int id = auto_id;
    auto_id++;
    return id;
}

/**
 * @brief Resets the auto-generated unique Surface ID counter to 1,000,000.
 */
void reset_surface_id()
{
    auto_id = DEFAULT_INIT_ID;
}

/**
 * @brief Maximize the auto-generated unique Surface ID counter.
 * @details This method updates the auto-generated unique Surface ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenMC compatibility module
 *          to ensure that the auto-generated Surface IDs do not
 *          collide with those created in OpenMC.
 * @param surface_id the id assigned to the auto-generated counter
 */
void maximize_surface_id(int surface_id)
{
    if (surface_id > auto_id)
        auto_id = surface_id;
}

/**
 * @brief Constructor assigns unique ID and user-defined ID for a surface.
 * @details Assigns a default boundary condition for this surface to
 *          BOUNDARY_NONE.
 * @param id an optional user-defined surface id
 */
Surface::Surface(const int id, const char *name)
{

    /* If the user did not define an optional ID, create one */
    if (id == 0)
        _id = surface_id();
    else if (id >= surface_id())
        log_printf(ERROR, "Unable to set the ID of a surface to %d since "
                          "surface IDs greater than or equal to 10000 is probibited "
                          "by OpenMOC.",
                   id);
    /* Use the user-defined ID */
    else
        _id = id;

    _uid = _n;
    _n++;

    _name = NULL;
    setName(name);
    _boundary_type = BOUNDARY_NONE;

    /* Initialize empty vectors of neighbor Cells for each halfspace */
    _neighbors[-1] = new std::vector<Cell *>();
    _neighbors[+1] = new std::vector<Cell *>();
}

/**
 * @brief Destructor.
 */
Surface::~Surface()
{
    if (_name != NULL)
        delete[] _name;

    if (!_neighbors.empty())
    {
        _neighbors[-1]->clear();
        _neighbors[+1]->clear();
        delete _neighbors[-1];
        delete _neighbors[+1];
        _neighbors.clear();
    }
}

/**
 * @brief Return the surface's unique ID.
 * @return the surface's unique ID
 */
int Surface::getUid() const
{
    return _uid;
}

/**
 * @brief Return the surface's user-defined id.
 * @return the surface's user-defined id
 */
int Surface::getId() const
{
    return _id;
}

/**
 * @brief Return the user-defined name of the Surface.
 * @return the Surface name
 */
char *Surface::getName() const
{
    return _name;
}

/**
 * @brief Returns the minimum coordinate in the axis direction of the
 *        surface.
 * @param axis The axis of interest (0 = x, 1 = y, 2 = z)
 * @param halfspace the halfspace to consider
 * @return the minimum coordinate in the axis direction
 */
double Surface::getMin(int axis, int halfspace)
{
    if (axis == 0)
        return getMinX(halfspace);
    else if (axis == 1)
        return getMinY(halfspace);
    else if (axis == 2)
        return getMinZ(halfspace);
    else
        log_printf(ERROR, "Could not retrieve minimum Surface coordinate since axis"
                          " is not recognized");
    return 0;
}

/**
 * @brief Returns the maximum coordinate in the axis direction of the
 *        surface.
 * @param axis The axis of interest (0 = x, 1 = y, 2 = z)
 * @param halfspace the halfspace to consider
 * @return the maximum coordinate in the axis direction
 */
double Surface::getMax(int axis, int halfspace)
{
    if (axis == 0)
        return getMaxX(halfspace);
    else if (axis == 1)
        return getMaxY(halfspace);
    else if (axis == 2)
        return getMaxZ(halfspace);
    else
        log_printf(ERROR, "Could not retrieve minimum Surface coordinate since axis"
                          " is not recognized");
    return 0;
}

/**
 * @brief Sets the name of the Surface.
 * @param name the Surface name string
 */
void Surface::setName(const char *name)
{
    int length = strlen(name);

    if (_name != NULL)
        delete[] _name;

    /* Initialize a character array for the Surface's name */
    _name = new char[length + 1];

    /* Copy the input character array Surface name to the class attribute name */
    for (int i = 0; i <= length; i++)
        _name[i] = name[i];
}
/**
 * @brief Return the type of surface (ie, XPLANE, ZCylinder, etc).
 * @return the surface type
 */
surfaceType Surface::getSurfaceType()
{
    return _surface_type;
}

/**
 * @brief Returns the type of boundary conditions for this surface (REFLECTIVE,
 *        VACUUM or BOUNDARY_NONE)
 * @return the type of boundary type for this surface
 */
boundaryType Surface::getBoundaryType()
{
    return _boundary_type;
}

/**
 * @brief Sets the boundary condition type (ie, VACUUM or REFLECTIVE) for this
 *        surface.
 * @param boundary_type the boundary condition forthis surface
 */
void Surface::setBoundaryType(boundaryType boundary_type)
{
    _boundary_type = boundary_type;
}

/**
 * @brief Adds a neighbor Cell to this Surface's collection of neighbors.
 * @param halfspace the +/-1 halfspace for the neighboring Cell
 * @param cell a pointer to the neighboring Cell
 */
void Surface::addNeighborCell(int halfspace, Cell *cell)
{

    if (halfspace != -1 && halfspace != +1)
        log_printf(ERROR, "Unable to add neighbor Cell %d to Surface %d since the "
                          "halfspace %d is not -1 or 1",
                   cell->getId(), _id, halfspace);

    /* Get pointer to vector of neighbor Cells for this halfspace */
    std::vector<Cell *> *neighbors = _neighbors[halfspace];

    /* Add the neighbor Cell if the collection does not already contain it */
    if (std::find(neighbors->begin(), neighbors->end(), cell) == neighbors->end())
        neighbors->push_back(cell);

    /* Update Cells with the neighbor Cells on the opposite Surface halfspace */
    std::vector<Cell *>::iterator iter1;
    std::vector<Cell *>::iterator iter2;
    for (iter1 = _neighbors[-1]->begin();
         iter1 != _neighbors[-1]->end(); ++iter1)
    {
        for (iter2 = _neighbors[1]->begin();
             iter2 != _neighbors[1]->end(); ++iter2)
            (*iter1)->addNeighborCell(*iter2);
    }

    for (iter1 = _neighbors[1]->begin();
         iter1 != _neighbors[1]->end(); ++iter1)
    {
        for (iter2 = _neighbors[-1]->begin();
             iter2 != _neighbors[-1]->end(); ++iter2)
            (*iter1)->addNeighborCell(*iter2);
    }
}

/**
 * @brief Return true or false if a point is on or off of a surface.
 * @param point pointer to the point of interest
 * @return on (true) or off (false) the surface
 */
bool Surface::isPointOnSurface(Point *point)
{

    /* Uses a threshold to determine whether the point is on the surface */
    if (abs(evaluate(point)) < ON_SURFACE_THRESH)
        return true;
    else
        return false;
}

/**
 * @brief Return true or false if a localcoord is on or off of a surface.
 * @param coord pointer to the localcoord of interest
 * @return on (true) or off (false) the surface
 */
bool Surface::isCoordOnSurface(LocalCoords *coord)
{
    return isPointOnSurface(coord->getPoint());
}

/**
 * @brief Finds the minimum distance to a Surface.
 * @details Finds the miniumum distance to a Surface from a LocalCoords
 *          with a trajectory defined by an angle to this Surface. If the
 *          trajectory will not intersect the Surface, returns INFINITY.
 * @param coords a pointer to a localcoords object
 * @return the minimum distance to the Surface
 */
double Surface::getMinDistance(LocalCoords *coords)
{

    Point *point = coords->getPoint();
    double phi = coords->getPhi();
    double polar = coords->getPolar();

    /* Point array for intersections with this Surface */
    Point intersections[2];

    /* Find the intersection Point(s) */
    int num_inters = this->intersection(point, phi, polar, intersections);
    double distance = INFINITY;

    /* If there is one intersection Point */
    if (num_inters == 1)
        distance = intersections[0].distanceToPoint(point);

    /* If there are two intersection Points */
    else if (num_inters == 2)
    {
        double dist1 = intersections[0].distanceToPoint(point);
        double dist2 = intersections[1].distanceToPoint(point);

        /* Determine which intersection Point is nearest */
        if (dist1 < dist2)
            distance = dist1;
        else
            distance = dist2;
    }

    return distance;
}
// /**
//  * @brief Constructor.
//  * @param id the surface id
//  * @param A the first coefficient in \f$ A * x + B * y + C = 0 \f$
//  * @param B the second coefficient in \f$ A * x + B * y + C = 0 \f$
//  * @param C the third coefficient in \f$ A * x + B * y + C = 0 \f$
//  */
// Plane::Plane(const double A, const double B, const double C, const int id) : Surface(id)
// {
//     _surface_type = PLANE;
//     _A = A;
//     _B = B;
//     _C = C;
// }

/**
 * @brief Prints a string representation of all of the Surface's objects to
 *        the console.
 */
void Surface::printString() {
  log_printf(RESULT, toString().c_str());
}
/**
 * @brief Constructor.
 * @param A the first coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param B the second coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param C the third coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param D the fourth coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param id the optional Surface ID
 * @param name the optional name of the Surface
 */
Plane::Plane(const double A, const double B,
             const double C, const double D, const int id, const char* name):
  Surface(id, name) {

  _surface_type = PLANE;
  _A = A;
  _B = B;
  _C = C;
  _D = D;
}

// /**
//  * @brief Returns the minimum x value of -INFINITY on this surface.
//  * @return the minimum x value of -INFINITY
//  */
// double Plane::getXMin()
// {
//     log_printf(ERROR, "Plane::getXMin() not yet implemented");
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the maximum x value of INFINITY on this surface.
//  * @return the maximum x value of INFINITY
//  */
// double Plane::getXMax()
// {
//     log_printf(ERROR, "Plane::getXMax() not yet implemented");
//     return std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the minimum y value of -INFINITY on this surface.
//  * @return the minimum y value of -INFINITY
//  */
// double Plane::getYMin()
// {
//     log_printf(ERROR, "Plane::getYMin not yet implemented");
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the maximum y value of INFINITY on this surface.
//  * @return the maximum y value of INFINITY
//  */
// double Plane::getYMax()
// {
//     log_printf(ERROR, "Plane::getYMax not yet implemented");
//     return std::numeric_limits<double>::infinity();
// }

/**
 * @brief Returns the minimum x value of -INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum x value of -INFINITY
 */
double Plane::getMinX(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum x value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum x value of INFINITY
 */
double Plane::getMaxX(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum y value of -INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum y value of -INFINITY
 */
double Plane::getMinY(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum y value of INFINITY
 */
double Plane::getMaxY(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum z value of -INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum z value of -INFINITY
 */
double Plane::getMinZ(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum z value of INFINITY
 */
double Plane::getMaxZ(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the A coefficient multiplying x in the surface equation
 * @return the value for the A coefficient
 */
double Plane::getA() {
  return _A;
}


/**
 * @brief Returns the B coefficient multiplying y in the surface equation
 * @return the value for the B coefficient
 */
double Plane::getB() {
  return _B;
}


/**
 * @brief Returns the C coefficient multiplying z in the surface equation
 * @return the value for the C coefficient
 */
double Plane::getC() {
  return _C;
}


/**
 * @brief Returns the D constant coefficient
 * @return the value for the D coefficient
 */
double Plane::getD() {
  return _D;
}

/**
 * @brief Finds the intersection point with this plane from a given point and
 *        trajectory defined by an angle.
 * @param point pointer to the point of interest
 * @param angle the angle defining the trajectory in radians
 * @param points pointer to a point to store the intersection point
 * @return the number of intersection points (0 or 1)
 */
inline int Plane::intersection(Point *point, double angle, Point *points)
{

    double x0 = point->getX();
    double y0 = point->getY();

    int num = 0;         /* number of intersections */
    double xcurr, ycurr; /* coordinates of current intersection point */

    /* The track is vertical */
    if ((fabs(angle - (M_PI / 2))) < 1.0e-10)
    {

        /* The plane is also vertical => no intersections */
        if (_B == 0)
            return 0;

        /* The plane is not vertical */
        else
        {
            xcurr = x0;
            ycurr = (-_A * x0 - _C) / _B;
            points->setCoords(xcurr, ycurr);

            /* Check that point is in same direction as angle */
            if (angle < M_PI && ycurr > y0)
                num++;
            else if (angle > M_PI && ycurr < y0)
                num++;
            return num;
        }
    }

    /* If the track isn't vertical */
    else
    {
        double m = sin(angle) / cos(angle);

        /* The plane and track are parallel, no intersections */
        if (fabs(-_A / _B - m) < 1e-11 && _B != 0)
            return 0;

        else
        {
            xcurr = -(_B * (y0 - m * x0) + _C) / (_A + _B * m);
            ycurr = y0 + m * (xcurr - x0);
            points->setCoords(xcurr, ycurr);

            if (angle < M_PI && ycurr > y0)
                num++;
            else if (angle > M_PI && ycurr < y0)
                num++;

            return num;
        }
    }
}

/**
 * @brief Converts this plane's attributes to a character array.
 * @details The character array returned conatins the type of plane (ie,
 *          PLANE) and the A, B, and C coefficients in the
 *          quadratic surface equation.
 * @return a character array of this plane's attributes
 */
std::string Plane::toString()
{
    std::stringstream string;

    string << "Surface id = " << _id << ", type = PLANE "
           << ", A = "
        //    << _A << ", B = " << _B << ", C = " << _C;
         << ", A = " << _A << ", B = " << _B << ", C = " << _C << ", D = " << _D;

    return string.str();
}

// /**
//  * @brief Prints a string representation of all of the plane's objects to
//  *        the console.
//  */
// void Plane::printString()
// {
//     log_printf(RESULT, toString().c_str());
// }

/**
 * @brief Constructor for a plane perpendicular to the x-axis.
 * @param id the user-defined surface id
 * @param x the location of the plane along the x-axis
 */
// XPlane::XPlane(const double x, const int id) : Plane(0, 1, -x, id)
XPlane::XPlane(const double x, const int id, const char* name):
  Plane(1, 0, 0, -x, id, name) {
    _surface_type = XPLANE;
    _x = x;
}

/**
 * @brief Set the location of this xplane on the x-axis.
 * @param x the location of the xplane on the x-axis
 */
void XPlane::setX(const double x)
{
    _x = x;
    _D = -x;
}

/**
 * @brief Returns the location of the xplane on the x-axis.
 * @return the location of the xplane on the x-axis
 */
double XPlane::getX()
{
    return _x;
}

/**
 * @brief Returns the minimum x value on the xplane.
 * @return the minimum x value
 */
double XPlane::getXMin()
{
    // return _x;
    if (halfspace == +1)
        return _x;
    else
        return -std::numeric_limits<double>::infinity();

}

/**
 * @brief Returns the maximum x value on the xplane.
 * @return the maximum x value
 */
double XPlane::getXMax()
{
    // return _x;
  if (halfspace == -1)
    return _x;
  else
    return std::numeric_limits<double>::infinity();
}

// /**
//  * @brief Returns the minimum y value of -INFINITY on the xplane.
//  * @return the minimum y value of -INFINITY
//  */
// double XPlane::getYMin()
// {
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * Returns the maximum y value of INFINITY on this xplane.
//  * @return the maximum y value of INFINITY
//  */
// double XPlane::getYMax()
// {
//     return std::numeric_limits<double>::infinity();
// }

/**
 * @brief Converts this xplane's attributes to a character array.
 * @details The character array returned conatins the type of plane (ie,
 *          XPLANE) and the A, B, and C coefficients in the
 *          quadratic surface equation and the location of the plane on
 *          the x-axis.
 * @return a character array of this xplane's attributes
 */
std::string XPlane::toString()
{
    std::stringstream string;

    // string << "Surface id = " << _id << ", type = XPLANE "
    //        << ", A = "
    //        << _A << ", B = " << _B << ", C = " << _C << ", x = " << _x;

      string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = XPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", x = " << _x;
    return string.str();
}

// /**
//  * @brief Constructor for a plane perpendicular to the y-axis.
//  * @param id the surface id
//  * @param y the location of the plane along the y-axis
//  */
// YPlane::YPlane(const double y, const int id) : Plane(1, 0, -y, id)
// {
/**
 * @brief Constructor for a Plane perpendicular to the y-axis.
 * @param y the location of the Plane along the y-axis
 * @param id the optional Surface id
 * @param name the optional Surface name
 */
YPlane::YPlane(const double y, const int id, const char* name):
  Plane(0, 1, 0, -y, id, name) {    
    _surface_type = YPLANE;
    _y = y;
}

/**
 * @brief Set the location of this yplane on the y-axis.
 * @param y the location of the yplane on the y-axis
 */
void YPlane::setY(const double y)
{
    _y = y;
      _D = -y;

}

/**
 * @brief Returns the location of the yplane on the y-axis.
 * @return the location of the yplane on the y-axis
 */
double YPlane::getY()
{
    return _y;
}

// /**
//  * @brief Returns the minimum x value of -INFINITY on this yplane.
//  * @return the minimum x value of -INFINITY
//  */
// double YPlane::getXMin()
// {
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the maximum x value of INFINITY on this yplane.
//  * @return the maximum x value of INFINITY
//  */
// double YPlane::getXMax()
// {
//     return std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the minimum y value on this yplane.
//  * @return the minimum y value
//  */
// double YPlane::getYMin()
// {
//     return _y;
// }

// /**
//  * @brief Returns the maximum y value on this yplane.
//  * @return the maximum y value
//  */
// double YPlane::getYMax()
// {
//     return _y;
// }



/**
 * @brief Returns the minimum y value for one of this YPlane's halfspaces.
 * @param halfspace the halfspace of the YPlane to consider
 * @return the minimum y value
 */
double YPlane::getMinY(int halfspace) {
  if (halfspace == +1)
    return _y;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value for one of this YPlane's halfspaces.
 * @param halfspace the halfspace of the YPlane to consider
 * @return the maximum y value
 */
double YPlane::getMaxY(int halfspace) {
  if (halfspace == -1)
    return _y;
  else
    return std::numeric_limits<double>::infinity();
}

/**
 * @brief Converts this yplane's attributes to a character array
 * @details The character array returned conatins the type of plane (ie,
 *          YPLANE) and the A, B, and C coefficients in the
 *          quadratic surface equation and the location of the plane on
 *          the y-axis.
 * @return a character array of this yplane's attributes
 */
std::string YPlane::toString()
{
    std::stringstream string;

    // string << "Surface id = " << _id << ", type = YPLANE "
    //        << ", A = "
    //        << _A << ", B = " << _B << ", C = " << _C << ", y = " << _y;
    string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = YPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", y = " << _y;
    return string.str();
}

// /**
//  * @brief Prints a string representation of all of the yplane's objects to
//  *        the console.
//  */
// void YPlane::printString()
// {
//     log_printf(RESULT, toString().c_str());
// }

/**
 * @brief Constructor for a plane perpendicular to the z-axis.
 * @param z the location of the Plane along the z-axis
 * @param id the optional Surface ID
 * @param name the optional Surface name * /*/
// ZPlane::ZPlane(const double z, const int id) : Plane(0, 0, -z, id)
// {
ZPlane::ZPlane(const double z, const int id, const char* name):
  Plane(0, 0, 1, -z, id, name) {
    _surface_type = ZPLANE;
    _z = z;
}

/**
 * @brief Set the location of this zplane on the z-axis.
 * @param z the location of the zplane on the z-axis
 */
void ZPlane::setZ(const double z)
{
    _z = z;
    _D = -z;
}

/**
 * @brief Returns the location of the zplane on the z-axis.
 * @return the location of the zplane on the z-axis
 */
double ZPlane::getZ()
{
    return _z;
}

// /**
//  * @brief Returns the minimum x value of -INFINITY on this zplane.
//  * @return the minimum x value of -INFINITY
//  */
// double ZPlane::getXMin()
// {
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the maximum x value of INFINITY on this zplane.
//  * @return the maximum x value of INFINITY
//  */
// double ZPlane::getXMax()
// {
//     return std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the minimum y value of -INFINITY on this zplane.
//  * @return the minimum y value of -INFINITY
//  */
// double ZPlane::getYMin()
// {
//     return -std::numeric_limits<double>::infinity();
// }

// /**
//  * @brief Returns the maximum y value on this zplane.
//  * @return the maximum y value
//  */
// double ZPlane::getYMax()
// {
//     return std::numeric_limits<double>::infinity();
// }


/**
 * @brief Returns the minimum z value for one of this ZPlane's halfspaces.
 * @param halfspace the halfspace of the ZPlane to consider
 * @return the minimum z value
 */
double ZPlane::getMinZ(int halfspace) {
  if (halfspace == +1)
    return _z;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value for one of this ZPlane's halfspaces.
 * @param halfspace the halfspace of the ZPlane to consider
 * @return the maximum z value
 */
double ZPlane::getMaxZ(int halfspace) {
  if (halfspace == -1)
    return _z;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Converts this zplane's attributes to a character array.
 * @details The character array returned conatins the type of plane (ie,
 *          ZPLANE) and the A, B, and C coefficients in the
 *          quadratic surface equation and the location of the plane along
 *          the z-axis.
 * @return a character array of this zplane's attributes
 */
std::string ZPlane::toString()
{
    std::stringstream string;

    // string << "Surface id = " << _id << ", type = ZPLANE "
    //        << ", A = "
    //        << _A << ", B = " << _B << ", C = " << _C << ", z = " << _z;


  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = ZPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", z = " << _z;
    return string.str();
}

// /**
//  * @brief Prints a string representation of all of the zplane's objects to
//  *        the console.
//  */
// void ZPlane::printString()
// {
//     log_printf(RESULT, toString().c_str());
// }

/**
 * @brief constructor.
 * @param id the surface id
 * @param x the x-coordinte of the ZCylinder center
 * @param y the y-coordinate of the ZCylinder center
 * @param radius the radius of the ZCylinder
 */
// Circle::Circle(const double x, const double y,
//                const double radius, const int id) : Surface(id)
// {
ZCylinder::ZCylinder(const double x, const double y,
               const double radius, const int id, const char* name):
  Surface(id, name) {
    _surface_type = ZCYLINDER;
    _A = 1.;
    _B = 1.;
    _C = -2. * x;
    _D = -2. * y;
    _E = x * x + y * y - radius * radius;
    _radius = radius;
    _center.setX(x);
    _center.setY(y);
}

/**
 * @brief Return the x-coordinate of the ZCylinder's center point.
 * @return the x-coordinate of the ZCylinder center
 */
double ZCylinder::getX0()
{
    return _center.getX();
}

/**
 * @brief Return the y-coordinate of the circle's center point.
 * @return the y-coordinate of the circle center
 */
double ZCylinder::getY0()
{
    return _center.getY();
}

/**
 * @brief Returns the minimum x value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum x value
 */
double ZCylinder::getMinX(int halfspace) {
  if (halfspace == -1)
    return _center.getX() - _radius;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum x value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum x value
 */
double ZCylinder::getMaxX(int halfspace) {
  if (halfspace == -1)
    return _center.getX() + _radius;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum y value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum y value
 */
double ZCylinder::getMinY(int halfspace) {
  if (halfspace == -1)
    return _center.getY() - _radius;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum y value
 */
double ZCylinder::getMaxY(int halfspace) {
  if (halfspace == -1)
    return _center.getY() + _radius;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum z value of -INFINITY.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum z value of -INFINITY
 */
double ZCylinder::getMinZ(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value of INFINITY.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum z value of INFINITY
 */
double ZCylinder::getMaxZ(int halfspace) {
  return std::numeric_limits<double>::infinity();
}

/**
 * @brief Finds the intersection Point with this zcylinder from a given Point
 *        and trajectory defined by an azim/polar angles (0, 1, or 2 points).
 * @param point pointer to the Point of interest
 * @param azim the azimuthal angle defining the trajectory in radians
 * @param polar the polar angle defining the trajectory in radians
 * @param points pointer to a an array of Points to store intersection Points,
 *        the array is not sorted by how close the points are to point
 * @return the number of intersection Points (0, 1 or 2)
 */
int ZCylinder::intersection(Point* point, double azim, double polar, Point* points) {

  double x0 = point->getX();
  double y0 = point->getY();
  double z0 = point->getZ();
  double xcurr = 0;
  double ycurr = 0;
  double zcurr = 0;

  /* Number of intersection Points */
  int num = 0;

  /* Vertical tracks only intersect Z cylinders at infinity */
  if (fabs(polar) < FLT_EPSILON || fabs(polar - M_PI) < FLT_EPSILON)
    return 0;

  /* If the track is vertical in y */
  if ((fabs(azim - M_PI_2)) < 1.0e-10 || (fabs(azim - 3.0 * M_PI_2)) < 1.0e-10) {

    /* Solve for where the line x = x0 and the Surface F(x,y) intersect
     * Find the y where F(x0, y) = 0
     * Substitute x0 into F(x,y) and rearrange to put in
     * the form of the quadratic formula: ay^2 + by + c = 0 */
    double a = 1.0;
    double b = _D;
    double c = _A * x0 * x0 + _C * x0 + _E;

    double discr = b*b - 4*c;

    /* There are no intersections */
    if (discr <= -ON_SURFACE_THRESH)
      return 0;

    /* There is one intersection (ie on the Surface) */
    else if (fabs(discr) < ON_SURFACE_THRESH) {
      double xcurr = x0;
      double ycurr = -b / 2;
      zcurr = z0 + fabs(ycurr - y0) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Check that point is in same direction as angle */
      if (azim < M_PI && ycurr > y0)
        num++;
      else if (azim > M_PI && ycurr < y0)
        num++;
    }

    /* There are two intersections */
    else {

      /* Compute first intersection point */
      xcurr = x0;
      ycurr = (-b + sqrt(discr)) / (2 * a);
      zcurr = z0 + fabs(ycurr - y0) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);
      if (azim < M_PI && ycurr > y0)
        num++;
      else if (azim > M_PI && ycurr < y0)
        num++;

      /* Compute second intersection point */
      ycurr = (-b - sqrt(discr)) / (2 * a);
      zcurr = z0 + fabs(ycurr - y0) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);
      if (azim < M_PI && ycurr > y0)
        num++;
      else if (azim > M_PI && ycurr < y0)
        num++;
    }
  }

  /* If the track isn't vertical */
  else {
    /* Solve for where the line y-y0 = m*(x-x0) and the Surface F(x,y)
     * intersect. Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
     * Substitute the point-slope formula for y into F(x,y) and
     * rearrange to put in the form of the quadratic formula:
     * ax^2 + bx + c = 0
     */
    double m = tan(azim);
    double q = y0 - m * x0;
    double a = 1 + m * m;
    double b = 2 * m * q + _C + _D * m;
    double c = q * q + _D * q + _E;

    double discr = b*b - 4*a*c;

    /* Boolean value describing whether the track is traveling to the right */
    bool right = azim < M_PI / 2. || azim > 3. * M_PI / 2.;

    /* There are no intersections */
    if (discr <= -ON_SURFACE_THRESH)
      return 0;

    /* There is one intersection (ie on the Surface) */
    else if (fabs(discr) < ON_SURFACE_THRESH) {
      xcurr = -b / (2*a);
      ycurr = y0 + m * (xcurr - x0);
      double interior = pow(ycurr - y0, 2.0) + pow(xcurr - x0, 2.0);
      zcurr = z0 + sqrt(interior) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;
    }

    /* There are two intersections */
    else {

      /* Determine first point of intersection */
      xcurr = (-b + sqrt(discr)) / (2*a);
      ycurr = y0 + m * (xcurr - x0);
      double interior = pow(ycurr - y0, 2.0) + pow(xcurr - x0, 2.0);
      zcurr = z0 + sqrt(interior) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      /* Determine second point of intersection */
      xcurr = (-b - sqrt(discr)) / (2*a);
      ycurr = y0 + m * (xcurr - x0);
      interior = pow(ycurr - y0, 2.0) + pow(xcurr - x0, 2.0);
      zcurr = z0 + sqrt(interior) * tan(M_PI_2 - polar);
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;
    }
  }
  return num;
}


// /**
//  * @brief Finds the intersection point with this circle from a given point and
//  *        trajectory defined by an angle (0, 1, or 2 points).
//  * @param point pointer to the point of interest
//  * @param angle the angle defining the trajectory in radians
//  * @param points pointer to a an array of points to store intersection points
//  * @return the number of intersection points (0 or 1)
//  */
// int ZCylinder::intersection(Point *point, double angle, Point *points)
// {

//     double x0 = point->getX();
//     double y0 = point->getY();
//     double xcurr, ycurr;
//     int num = 0; /* Number of intersection points */
//     double a, b, c, q, discr;

//     /* If the track is vertical */
//     if ((fabs(angle - (M_PI / 2))) < 1.0e-10)
//     {

//         /* Solve for where the line x = x0 and the surface F(x,y) intersect
//          * Find the y where F(x0, y) = 0
//          * Substitute x0 into F(x,y) and rearrange to put in
//          * the form of the quadratic formula: ay^2 + by + c = 0
//          */
//         a = _B * _B;
//         b = _D;
//         c = _A * x0 * x0 + _C * x0 + _E;

//         discr = b * b - 4 * a * c;

//         /* There are no intersections */
//         if (discr < 0)
//             return 0;

//         /* There is one intersection (ie on the surface) */
//         else if (discr == 0)
//         {
//             xcurr = x0;
//             ycurr = -b / (2 * a);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//                 num++;
//             else if (angle > M_PI && ycurr < y0)
//                 num++;
//             return num;
//         }

//         /* There are two intersections */
//         else
//         {
//             xcurr = x0;
//             ycurr = (-b + sqrt(discr)) / (2 * a);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//                 num++;
//             else if (angle > M_PI && ycurr < y0)
//                 num++;

//             xcurr = x0;
//             ycurr = (-b - sqrt(discr)) / (2 * a);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//                 num++;
//             else if (angle > M_PI && ycurr < y0)
//                 num++;
//             return num;
//         }
//     }

//     /* If the track isn't vertical */
//     else
//     {
//         /* Solve for where the line y-y0 = m*(x-x0) and the surface F(x,y)
//          * intersect. Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
//          * Substitute the point-slope formula for y into F(x,y) and
//          * rearrange to put in the form of the quadratic formula:
//          * ax^2 + bx + c = 0
//          */
//         double m = sin(angle) / cos(angle);
//         q = y0 - m * x0;
//         a = _A + _B * _B * m * m;
//         b = 2 * _B * m * q + _C + _D * m;
//         c = _B * q * q + _D * q + _E;

//         discr = b * b - 4 * a * c;

//         /* There are no intersections */
//         if (discr < 0)
//             return 0;

//         /* There is one intersection (ie on the surface) */
//         else if (discr == 0)
//         {
//             xcurr = -b / (2 * a);
//             ycurr = y0 + m * (points[0].getX() - x0);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//                 num++;
//             else if (angle > M_PI && ycurr < y0)
//                 num++;
//             return num;
//         }

//         /* There are two intersections */
//         else
//         {
//             xcurr = (-b + sqrt(discr)) / (2 * a);
//             ycurr = y0 + m * (xcurr - x0);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//             {
//                 num++;
//             }
//             else if (angle > M_PI && ycurr < y0)
//             {
//                 num++;
//             }

//             xcurr = (-b - sqrt(discr)) / (2 * a);
//             ycurr = y0 + m * (xcurr - x0);
//             points[num].setCoords(xcurr, ycurr);
//             if (angle < M_PI && ycurr > y0)
//             {
//                 num++;
//             }
//             else if (angle > M_PI && ycurr < y0)
//             {
//                 num++;
//             }

//             return num;
//         }
//     }
// }

// /**
//  * @brief Converts this circle's attributes to a character array.
//  * @details The character array returned conatins the type of plane (ie,
//  *          CIRCLE) and the A, B, C, D and E coefficients in the
//  *          quadratic surface equation.
//  * @return a character array of this circle's attributes
//  */
// std::string Circle::toString()
// {
//     std::stringstream string;

//     string << "Surface id = " << _id << ", type = CIRCLE "
//            << ", A = "
//            << _A << ", B = " << _B << ", C = " << _C << ", D = " << _D
//            << ", E = " << _E << ", x0 = " << _center.getX() << ", y0 = "
//            << _center.getY() << ", radius = " << _radius;

//     return string.str();
// }

/**
 * @brief Converts this ZCylinder's attributes to a character array.
 * @details The character array returned contains the type of Plane (ie,
 *          ZCYLINDER) and the A, B, C, D and E coefficients in the
 *          quadratic Surface equation.
 * @return a character array of this ZCylinder's attributes
 */
std::string ZCylinder::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = ZCYLINDER "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D << ", E = " << _E
         << ", x0 = " << _center.getX()
         << ", y0 = " << _center.getY()
         << ", radius = " << _radius;

    return string.str();
}

// /**
//  * @brief Prints a string representation of all of the circle's objects to
//  *        the console.
//  */
// void Circle::printString()
// {
//     log_printf(RESULT, toString().c_str());
// }

// /**
//  * @brief Returns the minimum x value on this circle.
//  * @return the minimum y value
//  */
// double Circle::getXMin()
// {
//     return _center.getX() - _radius;
// }

// /**
//  * @brief Returns the maximum x value on this circle.
//  * @return the maximum x value
//  */
// double Circle::getXMax()
// {
//     return _center.getX() + _radius;
// }

// /**
//  * @brief Returns the minimum y value on this circle.
//  * @return the minimum y value
//  */
// double Circle::getYMin()
// {
//     return _center.getY() - _radius;
// }

// /**
//  * @brief Returns ths maximum y value on this circle.
//  * @return the maximum y value
//  */
// double Circle::getYMax()
// {
//     return _center.getY() + _radius;
// }

