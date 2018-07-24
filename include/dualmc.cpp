#include "dualmc.h"

namespace dualmc
{

/**
 * @brief DualMC::index
 * @param x
 * @param y
 * @param z
 * @return
 */
int32_t DualMC::_index( const int32_t x, const int32_t y, const int32_t z ) const
{
    return x + _volumeDimensions[0] * (y + _volumeDimensions[1] * z);
}


/**
 * @brief DualMC::getCellCode
 * @param x
 * @param y
 * @param z
 * @param isoValue
 * @return
 */
int DualMC::_getCellCode( const int32_t x, const int32_t y, const int32_t z,
                         const uint8_t isoValue ) const
{
    // Determine for each cube corner if it is outside or inside
    int code = 0;

    if( _volumeGrid[_index( x, y, z )] >= isoValue )
        code |= 1;
    if( _volumeGrid[_index( x + 1, y, z )] >= isoValue )
        code |= 2;
    if( _volumeGrid[_index( x, y + 1, z )] >= isoValue )
        code |= 4;
    if( _volumeGrid[_index( x + 1, y + 1, z )] >= isoValue )
        code |= 8;
    if( _volumeGrid[_index( x, y, z + 1 )] >= isoValue )
        code |= 16;
    if( _volumeGrid[_index( x + 1, y, z + 1 )] >= isoValue )
        code |= 32;
    if( _volumeGrid[_index( x, y + 1, z + 1 )] >= isoValue )
        code |= 64;
    if( _volumeGrid[_index( x + 1, y + 1, z + 1 )] >= isoValue )
        code |= 128;

    return code;
}

/**
 * @brief DualMC::getDualPointCode
 * @param x
 * @param y
 * @param z
 * @param isoValue
 * @param edge
 * @return
 */
int DualMC::_getDualPointCode( const int32_t x, const int32_t y, const int32_t z,
                              const uint8_t isoValue,
                              const DMC_EDGE_CODE edge) const
{
    // Get the code of the cube that corresponds to the given XYZ voxel
    int cubeCode = _getCellCode( x, y, z, isoValue );

    // Is manifold dual marching cubes desired?
    if(generateManifold)
    {
        // The Manifold Dual Marching Cubes approach from Rephael Wenger as
        // described in chapter 3.3.5 of his book "Isosurfaces: Geometry,
        // Topology, and Algorithms" is implemente here.
        // If a problematic C16 or C19 configuration shares the ambiguous face
        // with another C16 or C19 configuration we simply invert the cube code
        // before looking up dual points.
        // Doing this for these pairs ensures manifold meshes.
        // But this removes the dualism to marching cubes.

        // Check if we have a potentially problematic configuration
        const uint8_t direction = problematicConfigs[uint8_t( cubeCode )];

        // If the direction code is in {0,...,5} we have a C16 or C19 configuration.
        if( direction != 255 )
        {
            // We have to check the neighboring cube, which shares the ambiguous
            // face. For this we decode the direction. This could also be done
            // with another lookup table.
            // Copy current cube coordinates into an array.
            int32_t neighborCoords[] = {x,y,z};

            // Get the dimension of the non-zero coordinate axis
            unsigned int const component = direction >> 1;

            // Get the sign of the direction
            int32_t delta = (direction & 1) == 1 ? 1 : -1;

            // Modify the correspong cube coordinate
            neighborCoords[component] += delta;

            // Have we left the volume in this direction?
            if( neighborCoords[component] >= 0 &&
                neighborCoords[component] < ( _volumeDimensions[component] - 1 ))
            {
                // Get the cube configuration of the relevant neighbor
                int neighborCubeCode = _getCellCode( neighborCoords[0],
                                                    neighborCoords[1],
                                                    neighborCoords[2],
                                                    isoValue );

                // Look up the neighbor configuration ambiguous face direction.
                // If the direction is valid we have a C16 or C19 neighbor.
                // As C16 and C19 have exactly one ambiguous face this face is
                // guaranteed to be shared for the pair.
                if( problematicConfigs[uint8_t( neighborCubeCode )] != 255 )
                {
                    // Replace the cube configuration with its inverse.
                    cubeCode ^= 0xff;
                }
            }
        }
    }

    for( int i = 0; i < 4; ++i )
    {
        if( dualPointsList[ cubeCode ][ i ] & edge )
        {
            return dualPointsList[ cubeCode ][ i ];
        }
    }

    return 0;
}


/**
 * @brief DualMC::calculateDualPoint
 * @param x
 * @param y
 * @param z
 * @param isoValue
 * @param pointCode
 * @param v
 */
void DualMC::_calculateDualPoint( const int32_t x,
                                  const int32_t y,
                                  const int32_t z,
                                  const uint8_t isoValue,
                                  const int pointCode,
                                  Vertex & v ) const
{
    // Initialize the point with lower voxel coordinates
    v.x = x;
    v.y = y;
    v.z = z;

    // Compute the dual point as the mean of the face vertices belonging to the
    // original marching cubes face
    Vertex p;
    p.x = 0; p.y = 0; p.z = 0;

    int points = 0;

    // Sum edge intersection vertices using the point code
    if( pointCode & EDGE0 )
    {
        p.x += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y, z )]) /
                (( float ) _volumeGrid[_index( x + 1, y, z )] -
                ( float ) _volumeGrid[_index( x, y, z )]);
        points++;
    }

    if( pointCode & EDGE1 )
    {
        p.x += 1.0f;
        p.z += (( float ) isoValue - ( float ) _volumeGrid[_index(x + 1, y, z )]) /
                (( float ) _volumeGrid[_index( x + 1, y, z + 1 )] -
                ( float ) _volumeGrid[_index(x + 1, y, z )]);
        points++;
    }

    if( pointCode & EDGE2 )
    {
        p.x += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y, z + 1 )]) /
                (( float ) _volumeGrid[_index( x + 1, y, z + 1 )] -
                ( float ) _volumeGrid[_index( x, y, z + 1 )]);
        p.z += 1.0f;
        points++;
    }

    if( pointCode & EDGE3 )
    {
        p.z += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y, z )] ) /
                (( float ) _volumeGrid[_index( x, y, z + 1)] -
                ( float ) _volumeGrid[_index( x, y, z )]);
        points++;
    }

    if( pointCode & EDGE4 )
    {
        p.x += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y + 1, z )]) /
                (( float ) _volumeGrid[_index( x + 1, y + 1, z )] -
                ( float ) _volumeGrid[_index( x, y + 1, z )]);
        p.y += 1.0f;
        points++;
    }

    if( pointCode & EDGE5 )
    {
        p.x += 1.0f;
        p.z += (( float ) isoValue - ( float ) _volumeGrid[_index( x + 1, y + 1, z )]) /
                (( float ) _volumeGrid[_index( x + 1, y + 1, z + 1 )] -
                ( float ) _volumeGrid[_index( x + 1, y + 1, z )]);
        p.y += 1.0f;
        points++;
    }

    if( pointCode & EDGE6 )
    {
        p.x += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y + 1, z + 1 )]) /
                (( float ) _volumeGrid[_index( x + 1, y + 1, z + 1 )] -
                ( float ) _volumeGrid[_index( x, y + 1, z + 1 )]);
        p.z += 1.0f;
        p.y += 1.0f;
        points++;
    }

    if( pointCode & EDGE7 )
    {
        p.z += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y + 1 , z )]) /
                (( float ) _volumeGrid[_index( x, y + 1, z + 1 )] -
                ( float ) _volumeGrid[_index( x, y + 1 , z )]);
        p.y += 1.0f;
        points++;
    }

    if( pointCode & EDGE8 )
    {
        p.y += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y, z )]) /
                (( float ) _volumeGrid[_index(x, y + 1, z )] -
                ( float ) _volumeGrid[_index( x, y, z )]);
        points++;
    }

    if( pointCode & EDGE9 )
    {
        p.x += 1.0f;
        p.y += (( float ) isoValue - ( float ) _volumeGrid[_index( x + 1, y, z )]) /
                (( float ) _volumeGrid[_index(x + 1, y + 1, z )] -
                ( float ) _volumeGrid[_index( x + 1, y, z )]);
        points++;
    }

    if( pointCode & EDGE10 )
    {
        p.x += 1.0f;
        p.y += (( float ) isoValue - ( float ) _volumeGrid[_index(x + 1, y, z + 1 )]) /
                (( float ) _volumeGrid[_index( x + 1, y + 1, z + 1 )] -
                ( float ) _volumeGrid[_index( x + 1, y, z + 1 )]);
        p.z += 1.0f;
        points++;
    }

    if( pointCode & EDGE11 )
    {
        p.z += 1.0f;
        p.y += (( float ) isoValue - ( float ) _volumeGrid[_index( x, y, z + 1 )]) /
                (( float ) _volumeGrid[_index( x, y + 1, z + 1 )] -
                ( float ) _volumeGrid[_index( x, y, z + 1 )]);
        points++;
    }

    // Divide by number of accumulated points
    float invPoints = 1.0f / ( float ) points;
    p.x *= invPoints; p.y *= invPoints; p.z *= invPoints;

    // Offset point by voxel coordinates
    v.x += p.x;
    v.y += p.y;
    v.z += p.z;
}

/**
 * @brief DualMC::getSharedDualPointIndex
 * @param x
 * @param y
 * @param z
 * @param isoValue
 * @param edge
 * @param vertices
 * @return
 */
int32_t DualMC::_getSharedDualPointIndex( const int32_t x,
                                          const int32_t y,
                                          const int32_t z,
                                          const uint8_t isoValue,
                                          const DMC_EDGE_CODE edge,
                                          std::vector<Vertex> & vertices )
{
    // Create a key for the dual point from its linearized cell ID and point code
    DualPointKey key;
    key.linearizedCellID = _index(x,y,z);
    key.pointCode = _getDualPointCode(x,y,z,isoValue,edge);

    // have we already computed the dual point?
    auto iterator = pointToIndex.find( key );
    if( iterator != pointToIndex.end())
    {
        // Just return the dual point index
        return iterator->second;
    }
    else
    {
        // Create new vertex and vertex id
        int32_t newVertexId = vertices.size();
        vertices.emplace_back();
        _calculateDualPoint( x, y, z, isoValue, key.pointCode, vertices.back());

        // Insert vertex ID into map and also return it
        pointToIndex[key] = newVertexId;
        return newVertexId;
    }
}

/**
 * @brief DualMC::build
 * @param data
 * @param x
 * @param y
 * @param z
 * @param isoValue
 * @param generateManifold
 * @param generateSoup
 * @param vertices
 * @param quads
 */
void DualMC::build(const uint8_t* data,
                   const int32_t x, const int32_t y, const int32_t z,
                   const uint8_t isoValue,
                   const bool generateManifold,
                   const bool generateSoup,
                   std::vector<Vertex> & vertices,
                   std::vector<Quad> & quads)
{

    // Set members
    this->_volumeDimensions[0] = x;
    this->_volumeDimensions[1] = y;
    this->_volumeDimensions[2] = z;
    this->_volumeGrid = data;
    this->generateManifold = generateManifold;

    // Clear vertices and quad indices
    vertices.clear();
    quads.clear();

    // Generate quad soup or shared vertices quad list
    if( generateSoup )
    {
        _buildQuadSoup( isoValue, vertices, quads );
    }
    else
    {
        _buildSharedVerticesQuads( isoValue, vertices, quads );
    }
}

/**
 * @brief DualMC::buildSharedVerticesQuads
 * @param isoValue
 * @param vertices
 * @param quads
 */
void DualMC::_buildSharedVerticesQuads( const uint8_t isoValue,
                                       std::vector<Vertex> & vertices,
                                       std::vector<Quad> & quads)
{
    // TODO: Why the volume dimensions are reduced by two ?!!
    int32_t const reducedX = _volumeDimensions[0] - 2;
    int32_t const reducedY = _volumeDimensions[1] - 2;
    int32_t const reducedZ = _volumeDimensions[2] - 2;

    int32_t i0, i1, i2, i3;

    pointToIndex.clear();

    // Iterate voxels
    for( int32_t z = 0; z < reducedZ; ++z )
    {
        for( int32_t y = 0; y < reducedY; ++y )
        {
            for( int32_t x = 0; x < reducedX; ++x )
            {
                // Construct quads for X edge
                if( z > 0 && y > 0 )
                {
                    bool const entering = _volumeGrid[_index( x, y, z)] < isoValue &&
                                          _volumeGrid[_index( x + 1, y, z)] >= isoValue;
                    bool const exiting  = _volumeGrid[_index( x, y, z)] >= isoValue &&
                                          _volumeGrid[_index( x + 1, y, z)] < isoValue;
                    if( entering || exiting )
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex( x, y, z,
                                                      isoValue, EDGE0, vertices );
                        i1 = _getSharedDualPointIndex( x, y, z - 1,
                                                      isoValue, EDGE2, vertices );
                        i2 = _getSharedDualPointIndex( x, y - 1, z - 1,
                                                      isoValue, EDGE6, vertices );
                        i3 = _getSharedDualPointIndex( x, y - 1, z,
                                                      isoValue, EDGE4, vertices );

                        if( entering )
                        {
                            quads.emplace_back( i0, i1, i2, i3);
                        }
                        else
                        {
                            quads.emplace_back( i0, i3, i2, i1 );
                        }
                    }
                }

                // Construct quads for y edge
                if( z > 0 && x > 0 )
                {
                    bool const entering = _volumeGrid[_index( x, y, z )] < isoValue &&
                                          _volumeGrid[_index( x, y + 1, z )] >= isoValue;
                    bool const exiting  = _volumeGrid[_index( x, y, z )] >= isoValue &&
                                          _volumeGrid[_index( x, y + 1, z )] < isoValue;

                    if( entering || exiting )
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex( x, y, z,
                                                      isoValue, EDGE8, vertices );
                        i1 = _getSharedDualPointIndex( x, y, z - 1,
                                                      isoValue, EDGE11, vertices );
                        i2 = _getSharedDualPointIndex( x - 1, y, z - 1,
                                                      isoValue, EDGE10, vertices );
                        i3 = _getSharedDualPointIndex( x - 1, y, z,
                                                      isoValue, EDGE9, vertices );

                        if( exiting )
                        {
                            quads.emplace_back( i0, i1, i2, i3 );
                        }
                        else
                        {
                            quads.emplace_back( i0, i3, i2, i1 );
                        }
                    }
                }

                // Construct quads for z edge
                if( x > 0 && y > 0 )
                {
                    bool const entering = _volumeGrid[_index( x, y, z )] < isoValue &&
                                          _volumeGrid[_index( x, y, z + 1 )] >= isoValue;
                    bool const exiting  = _volumeGrid[_index( x, y, z )] >= isoValue &&
                                          _volumeGrid[_index( x, y, z + 1 )] < isoValue;
                    if( entering || exiting )
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex( x, y, z,
                                                      isoValue, EDGE3, vertices);
                        i1 = _getSharedDualPointIndex( x - 1, y, z,
                                                      isoValue, EDGE1, vertices);
                        i2 = _getSharedDualPointIndex( x - 1, y - 1, z,
                                                      isoValue, EDGE5, vertices);
                        i3 = _getSharedDualPointIndex( x, y - 1, z,
                                                      isoValue, EDGE7, vertices );

                        if( exiting )
                        {
                            quads.emplace_back( i0, i1, i2, i3 );
                        }
                        else
                        {
                            quads.emplace_back( i0, i3, i2, i1 );
                        }
                    }
                }
            }
        }
    }
}


void DualMC::_buildQuadSoup(uint8_t const isoValue,
    std::vector<Vertex> & vertices,
    std::vector<Quad> & quads
    ) {

    int32_t const reducedX = _volumeDimensions[0] - 2;
    int32_t const reducedY = _volumeDimensions[1] - 2;
    int32_t const reducedZ = _volumeDimensions[2] - 2;

    Vertex vertex0;
    Vertex vertex1;
    Vertex vertex2;
    Vertex vertex3;
    int pointCode;

    // iterate voxels
    for(int32_t z = 0; z < reducedZ; ++z)
        for(int32_t y = 0; y < reducedY; ++y)
            for(int32_t x = 0; x < reducedX; ++x) {
                // construct quad for x edge
                if(z > 0 && y > 0) {
                    // is edge intersected?
                    bool const entering = _volumeGrid[_index(x,y,z)] < isoValue && _volumeGrid[_index(x+1,y,z)] >= isoValue;
                    bool const exiting  = _volumeGrid[_index(x,y,z)] >= isoValue && _volumeGrid[_index(x+1,y,z)] < isoValue;
                    if(entering || exiting){
                        // generate quad
                        pointCode = _getDualPointCode(x,y,z,isoValue,EDGE0);
                        _calculateDualPoint(x,y,z,isoValue,pointCode, vertex0);

                        pointCode = _getDualPointCode(x,y,z-1,isoValue,EDGE2);
                        _calculateDualPoint(x,y,z-1,isoValue,pointCode, vertex1);

                        pointCode = _getDualPointCode(x,y-1,z-1,isoValue,EDGE6);
                        _calculateDualPoint(x,y-1,z-1,isoValue,pointCode, vertex2);

                        pointCode = _getDualPointCode(x,y-1,z,isoValue,EDGE4);
                        _calculateDualPoint(x,y-1,z,isoValue,pointCode, vertex3);

                        if(entering) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }

                // construct quad for y edge
                if(z > 0 && x > 0) {
                    // is edge intersected?
                    bool const entering = _volumeGrid[_index(x,y,z)] < isoValue && _volumeGrid[_index(x,y+1,z)] >= isoValue;
                    bool const exiting  = _volumeGrid[_index(x,y,z)] >= isoValue && _volumeGrid[_index(x,y+1,z)] < isoValue;
                    if(entering || exiting){
                        // generate quad
                        pointCode = _getDualPointCode(x,y,z,isoValue,EDGE8);
                        _calculateDualPoint(x,y,z,isoValue,pointCode, vertex0);

                        pointCode = _getDualPointCode(x,y,z-1,isoValue,EDGE11);
                        _calculateDualPoint(x,y,z-1,isoValue,pointCode, vertex1);

                        pointCode = _getDualPointCode(x-1,y,z-1,isoValue,EDGE10);
                        _calculateDualPoint(x-1,y,z-1,isoValue,pointCode, vertex2);

                        pointCode = _getDualPointCode(x-1,y,z,isoValue,EDGE9);
                        _calculateDualPoint(x-1,y,z,isoValue,pointCode, vertex3);

                        if(exiting) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }

                // construct quad for z edge
                if(x > 0 && y > 0) {
                    // is edge intersected?
                    bool const entering = _volumeGrid[_index(x,y,z)] < isoValue && _volumeGrid[_index(x,y,z+1)] >= isoValue;
                    bool const exiting  = _volumeGrid[_index(x,y,z)] >= isoValue && _volumeGrid[_index(x,y,z+1)] < isoValue;
                    if(entering || exiting){
                        // generate quad
                        pointCode = _getDualPointCode(x,y,z,isoValue,EDGE3);
                        _calculateDualPoint(x,y,z,isoValue,pointCode, vertex0);

                        pointCode = _getDualPointCode(x-1,y,z,isoValue,EDGE1);
                        _calculateDualPoint(x-1,y,z,isoValue,pointCode, vertex1);

                        pointCode = _getDualPointCode(x-1,y-1,z,isoValue,EDGE5);
                        _calculateDualPoint(x-1,y-1,z,isoValue,pointCode, vertex2);

                        pointCode = _getDualPointCode(x,y-1,z,isoValue,EDGE7);
                        _calculateDualPoint(x,y-1,z,isoValue,pointCode, vertex3);

                        if(exiting) {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex1);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex3);
                        } else {
                            vertices.emplace_back(vertex0);
                            vertices.emplace_back(vertex3);
                            vertices.emplace_back(vertex2);
                            vertices.emplace_back(vertex1);
                        }
                    }
                }
            }

    // generate triangle soup quads
    size_t const numQuads = vertices.size() / 4;
    quads.reserve(numQuads);
    for (size_t i = 0; i < numQuads; ++i) {
        quads.emplace_back(i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3);
    }
}

/**
 * @brief DualMC::DualPointKey::operator ==
 * @param other
 * @return
 */
bool DualMC::DualPointKey::operator==( typename DualMC::DualPointKey const & other ) const
{
    return linearizedCellID == other.linearizedCellID && pointCode == other.pointCode;
}



}
