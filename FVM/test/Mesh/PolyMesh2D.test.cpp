#include "TestUtils.h"
#include <Mesh/2D/Unstructured/PolyMesh2D.h>
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <sstream>

class TestPolyMesh2D : public testing::Test, protected PolyMesh2D
{
protected:
    TestPolyMesh2D()
        : PolyMesh2D(std::stringstream(
#include "poiseuille_unifrom_11x15.msh"
              ))
    {
    }

    static constexpr Index  nx = 11, ny = 15;
    static constexpr Scalar dx = 1.0 / nx, dy = 1.0 / ny;
    static constexpr Index  totalCells = nx * ny;
    static constexpr Index  totalFaces = 2 * nx * ny + nx + ny;

    static constexpr Scalar tolerance = 1e-5;
};

TEST_F(TestPolyMesh2D, CellCentroidArrayTest)
{
    ASSERT_EQ(m_cellCentroids.size(), totalCells);

    for (Index x = 0; x < nx; x++)
    {
        for (Index y = 0; y < ny; y++)
        {
            Index  cellIdx = y * nx + x;
            Vector expectedCentroid = {x * dx + dx / 2, y * dy + dy / 2, 0};
            Vector centroid = getCellCentroid(cellIdx);

            EXPECT_LT((centroid - expectedCentroid).norm(), tolerance) << "cell idx is " << cellIdx;
        }
    }
}

TEST_F(TestPolyMesh2D, CellVolumesArrayTest)
{
    ASSERT_EQ(m_cellVolumes.size(), totalCells);

    Scalar expectedCellVolume = dx * dy;
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Scalar cellVolume = getCellVolume(cellIdx);

        EXPECT_NEAR(cellVolume, expectedCellVolume, tolerance) << "cell idx is " << cellIdx;
    }
}

TEST_F(TestPolyMesh2D, CellFacesArrayTest)
{
    ASSERT_EQ(m_cellFaces.size(), totalCells);

    for (Index x = 0; x < nx; x++)
    {
        for (Index y = 0; y < ny; y++)
        {
            Index cellIdx = y * nx + x;
            Index f = y * (2 * nx + 1) + x;

            List<Index> expectedFaces = {f, f + nx, f + nx + 1, f + 2 * nx + 1};
            List<Index> faces = getCellFaces(cellIdx);
            std::sort(faces.begin(), faces.end());

            EXPECT_EQ(faces, expectedFaces) << "cell idx is " << cellIdx;
        }
    }
}

TEST_F(TestPolyMesh2D, FaceVectorArrayTest)
{
    ASSERT_EQ(m_faceVectors.size(), totalFaces);

    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Vector faceVector = getFaceVector(faceIdx);
        Vector expectedFaceVector = (faceIdx % (2 * nx + 1) < nx ? Vector{0, dx, 0} : Vector{dy, 0, 0});
        if (faceIdx < nx || faceIdx % (2 * nx + 1) == nx)
            expectedFaceVector *= -1;

        EXPECT_LT((faceVector - expectedFaceVector).norm(), tolerance) << "face idx is " << faceIdx;
    }
}

TEST_F(TestPolyMesh2D, FaceCentroidArrayTest)
{
    ASSERT_EQ(m_faceCentroids.size(), totalFaces);

    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Index y = faceIdx / (2 * nx + 1);
        Index x = faceIdx % (2 * nx + 1);
        if (x >= nx)
            x -= nx;
        Vector expectedCentroid = {x * dx, y * dy, 0};
        expectedCentroid += (faceIdx % (2 * nx + 1) < nx ? Vector{dx / 2, 0, 0} : Vector{0, dy / 2, 0});
        Vector centroid = getFaceCentroid(faceIdx);

        EXPECT_LT((centroid - expectedCentroid).norm(), tolerance) << "face idx is " << faceIdx;
    }
}

TEST_F(TestPolyMesh2D, FaceNeighborsArrayTest)
{
    ASSERT_EQ(m_faceNeighbors.size(), totalFaces);

    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Index y = faceIdx / (2 * nx + 1);
        Index x = faceIdx % (2 * nx + 1);
        if (x >= nx)
            x -= nx;
        Array<Index, 2> expectedNeighbours = {-1, nx * y + x};
        expectedNeighbours[0] = (faceIdx % (2 * nx + 1) < nx ? nx * (y - 1) + x : nx * y + x - 1);
        if (faceIdx < nx)
            expectedNeighbours = {x, -1};
        else if (faceIdx >= totalFaces - nx)
            expectedNeighbours = {nx * (y - 1) + x, -1};
        else if (faceIdx % (2 * nx + 1) == nx)
            expectedNeighbours = {nx * y + x, -1};
        else if (faceIdx % (2 * nx + 1) == 2 * nx)
            expectedNeighbours = {nx * y + x - 1, -1};

        auto neighbours = getFaceNeighbors(faceIdx);

        EXPECT_EQ(neighbours, expectedNeighbours) << "face idx is " << faceIdx;
    }
}

TEST_F(TestPolyMesh2D, BoundariesTest)
{
    const Boundaries leftBoundaries = {{{0.01, 0, 0}, BoundaryConditionType::FIXED_VALUE},
                                       {0, BoundaryConditionType::FIXED_GRADIENT}};
    const Boundaries rightBoundaries = {{{0, 0, 0}, BoundaryConditionType::FIXED_GRADIENT},
                                        {0, BoundaryConditionType::FIXED_VALUE}};
    const Boundaries topBoundaries = {{{0, 0, 0}, BoundaryConditionType::FIXED_VALUE},
                                      {0, BoundaryConditionType::FIXED_GRADIENT}};
    const Boundaries bottomBoundaries = {{{0, 0, 0}, BoundaryConditionType::FIXED_VALUE},
                                         {0, BoundaryConditionType::FIXED_GRADIENT}};

    ASSERT_EQ(m_boundariesMap.size(), 2 * nx + 2 * ny);

    // top
    for (Index faceIdx = 0; faceIdx < nx; faceIdx++)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, topBoundaries);

    // bottom
    for (Index faceIdx = totalFaces - nx; faceIdx < totalFaces; faceIdx++)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, bottomBoundaries);

    // left
    for (Index faceIdx = nx; faceIdx < totalFaces; faceIdx += 2 * nx + 1)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, leftBoundaries);

    // right
    for (Index faceIdx = 2 * nx; faceIdx < totalFaces; faceIdx += 2 * nx + 1)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, rightBoundaries);
}