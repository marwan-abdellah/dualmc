#ifndef DUALMC_H_INCLUDED
#define DUALMC_H_INCLUDED

// C includes
#include <cstdint>

// STL includes
#include <unordered_map>
#include <vector>

#include "quad.h"
#include "vertex.h"
#include "tables.h"

namespace dualmc
{
    
/**
 * @brief The DualMC class
 * Class which implements the dual marching cubes algorithm from Gregory M. Nielson.
 * Faces and vertices of the standard marching cubes algorithm correspond to
 * vertices and faces in the dual algorithm. As a vertex in standard marching cubes
 * usually is shared by 4 faces, the dual mesh is entirely made from quadrangles.
 * Unfortunately, under rare circumstances the original algorithm can create
 * non-manifold meshes. See the remarks of the original paper on this.
 * The class optionally can guarantee manifold meshes by taking the Manifold
 * Dual Marching Cubes approach from Rephael Wenger as described in
 * chapter 3.3.5 of his book "Isosurfaces: Geometry, Topology, and Algorithms".
 */
class DualMC
{
public:

    /**
     * @brief build
     * Extracts the iso surface for a given volume and iso value.
     * Output is a list of vertices and a list of indices, which connect
     * vertices to quads.
     * The quad mesh either uses shared vertex indices or is a quad soup if
     * desired.
     * @param volumeGrid
     * @param x
     * @param y
     * @param z
     * @param isoValue
     * @param generateManifold
     * @param generateSoup
     * @param vertices
     * @param quads
     */
    void build( const uint8_t* volumeGrid,
                int32_t const x, int32_t const y, int32_t const z,
                uint8_t const isoValue,
                bool const _generateManifold, bool const generateSoup,
                std::vector<Vertex> & vertices, std::vector<Quad> & quads );

private:

    /**
     * @brief _buildSharedVerticesQuads
     * Extract quad mesh with shared vertex indices.
     * @param iso
     * @param vertices
     * @param quads
     */
    void _buildSharedVerticesQuads( const uint8_t iso,
                                   std::vector<Vertex> & vertices,
                                   std::vector<Quad> & quads );

    /**
     * @brief _buildQuadSoup
     * Extract quad soup.
     * @param isoValue
     * @param vertices
     * @param quads
     */
    void _buildQuadSoup( const uint8_t isoValue,
                        std::vector<Vertex> & vertices,
                        std::vector<Quad> & quads );

private:

    /**
     * @brief _getCellCode
     * Get the 8-bit in-out mask for the voxel corners of the cell cube at
     * ( x, y, z ) and the given iso value.
     * @param x
     * @param y
     * @param z
     * @param isoValue
     * @return
     */
    int _getCellCode( const int32_t x, const int32_t y, const int32_t z,
                      const uint8_t isoValue ) const;

    /**
     * @brief _getDualPointCode
     * Get the 12-bit dual point code mask, which encodes the traditional
     * marching cube vertices of the traditional marching cubes face which
     * corresponds to the dual point.
     * This is also where the manifold dual marching cubes algorithm is
     * implemented.
     * @param x
     * @param y
     * @param z
     * @param isoValue
     * @param edge
     * @return
     */
    int _getDualPointCode( const int32_t x, const int32_t y, const int32_t z,
                           const uint8_t isoValue,
                           const DMC_EDGE_CODE edge ) const;

    /**
     * @brief _calculateDualPoint
     * Given a dual point code and iso value, compute the dual point.
     * @param x
     * @param y
     * @param z
     * @param isoValue
     * @param pointCode
     * @param v
     */
    void _calculateDualPoint( const int32_t x, const int32_t y, const int32_t z,
                              uint8_t const isoValue, int const pointCode,
                              Vertex &v ) const;

    /**
     * @brief _getSharedDualPointIndex
     * Get the shared index of a dual point which is uniquly identified by its
     * cell cube index and a cube edge. The dual point is computed, if it has
     * not been computed before.
     * @param x
     * @param y
     * @param cz
     * @param isoValue
     * @param edge
     * @param vertices
     * @return
     */
    int32_t _getSharedDualPointIndex( const int32_t x,
                                      const int32_t y,
                                      const int32_t cz,
                                      const uint8_t isoValue,
                                      const DMC_EDGE_CODE edge,
                                      std::vector<Vertex> & vertices );

    /**
     * @brief _index
     * Compute a linearized cell cube index.
     * @param x
     * @param y
     * @param z
     * @return
     */
    int32_t _index( const int32_t x, const int32_t y, const int32_t z ) const;
    
private:

    /**
     * @brief _volumeDimensions
     * Volume dimensions.
     */
    int32_t _volumeDimensions[3];

    /**
     * @brief _volumeGrid
     * The input volume grid.
     */
    uint8_t const * _volumeGrid;
    

    /**
     * @brief _generateManifold
     * Store whether the manifold dual marching cubes algorithm should be applied.
     */
    bool _generateManifold;

    /**
     * @brief The DualPointKey struct
     * Dual point key structure for hashing of shared vertices
     */
    struct DualPointKey
    {
        // A dual point can be uniquely identified by ite linearized volume cell
        // id and point code
        int32_t linearizedCellID;
        int pointCode;

        /// Equal operator for unordered map
        bool operator==( DualPointKey const & other ) const;
    };

    /**
     * @brief The DualPointKeyHash struct
     * Functor for dual point key hash generation
     */
    struct DualPointKeyHash
    {
        size_t operator()( DualPointKey const & k ) const
        {
            return size_t( k.linearizedCellID ) | ( size_t( k.pointCode ) << 32u );
        }
    };

    /**
     * @brief pointToIndex
     * Hash map for shared vertex index computations
     */
    std::unordered_map< DualPointKey, int32_t, DualPointKeyHash > pointToIndex;
};
}
#endif // DUALMC_H_INCLUDED
