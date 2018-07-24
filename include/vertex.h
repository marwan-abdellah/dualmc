#ifndef VERTEX_H
#define VERTEX_H

// c includes
#include <cstdint>

namespace dualmc
{

/// Vertex structure for dual points
struct Vertex
{
    /// Non-initializing constructor
    Vertex();

    /// Initializing constructor
    Vertex(float x, float y, float z);

    /// Initializing constructor
    Vertex(Vertex const & v);

    // Components
    float x,y,z;
};

}

#endif // VERTEX_H

