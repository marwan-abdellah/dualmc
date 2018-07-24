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
    Vertex()
    {
        /// EMPTY
    }

    /// Initializing constructor
    Vertex( float x, float y, float z )
        : x( x ),
          y( y ),
          z( z )
    {
        /// EMPTY
    }

    /// Initializing constructor
    Vertex( Vertex const & v )
        : x( v.x ),
          y( v.y ),
          z( v.z )
    {
        /// EMPTY
    }

    // Components
    float x,y,z;
};

}

#endif // VERTEX_H

