#ifndef ROVER_SENSOR_H
#define ROVER_SENSOR_H

#include <set>
#include <memory>
#include "location.h"

class Sensor {
public:
    virtual bool is_safe(coordinate_t x, coordinate_t y) = 0;
};

using sensor_ptr = std::unique_ptr<Sensor>;
using sensors_collection_t = std::set<sensor_ptr>;

#endif // ROVER_SENSOR_H
