#ifndef EDGES_H
#define EDGES_H

/**
 * @brief The DMC_EDGE_CODE enum
 * Enum with edge codes for a '12-bit' voxel edge mask to indicate
 * grid edges which intersect the ISO surface of classic marching cubes.
 */
enum DMC_EDGE_CODE
{
    EDGE0 = 1,
    EDGE1 = 1 << 1,
    EDGE2 = 1 << 2,
    EDGE3 = 1 << 3,
    EDGE4 = 1 << 4,
    EDGE5 = 1 << 5,
    EDGE6 = 1 << 6,
    EDGE7 = 1 << 7,
    EDGE8 = 1 << 8,
    EDGE9 = 1 << 9,
    EDGE10 = 1 << 10,
    EDGE11 = 1 << 11,
    FORCE_32BIT = 0xffffffff
};

#endif // EDGES_H

