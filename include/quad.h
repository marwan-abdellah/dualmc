#ifndef QUAD_H
#define QUAD_H

// c includes
#include <cstdint>

namespace dualmc
{

/// Quad indices structure
struct Quad
{
    /// Non-initializing constructor
    Quad();

    /// Initializing constructor
    Quad(int32_t i0, int32_t i1,int32_t i2, int32_t i3);

    // quad indices
    int32_t i0,i1,i2,i3;

};

}

#endif // QUAD_H

