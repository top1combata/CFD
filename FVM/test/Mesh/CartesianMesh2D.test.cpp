#include "TestUtils.h"
#include <Mesh/2D/Structured/CartesianMesh2D.h>
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>

using CartesianMeshTestParam = std::tuple<Index, Index, Scalar, Scalar>;

class TestCartesinaMesh2D : public testing::TestWithParam<CartesianMeshTestParam>, protected CartesianMesh2D
{
protected:
    TestCartesinaMesh2D()
        : CartesianMesh2D(std::get<0>(TestWithParam::GetParam()), std::get<1>(TestWithParam::GetParam()),
                          std::get<2>(TestWithParam::GetParam()), std::get<3>(TestWithParam::GetParam())),
          nx(std::get<0>(TestWithParam::GetParam())), ny(std::get<1>(TestWithParam::GetParam())), totalCells(nx * ny),
          totalFaces(2 * nx * ny + nx + ny), dx(std::get<2>(TestWithParam::GetParam()) / nx),
          dy(std::get<3>(TestWithParam::GetParam()) / ny)
    {
    }

    const Index             nx, ny, totalCells, totalFaces;
    const Scalar            dx, dy;
    static constexpr Scalar tolerance = 1e-5;
};

TEST_P(TestCartesinaMesh2D, ConstructorCheck)
{
    EXPECT_EQ(m_nx, nx);
    EXPECT_EQ(m_ny, ny);
    EXPECT_EQ(m_boundariesMap.size(), 0);
}

TEST_P(TestCartesinaMesh2D, CellCentroidArrayTest)
{
    ASSERT_EQ(m_cellCentroids.size(), totalCells);

    for (Index x = 0; x < nx; x++)
    {
        for (Index y = 0; y < ny; y++)
        {
            Index  cellIdx = y * nx + x;
            Vector expectedCentroid = {x * dx + dx / 2, y * dy + dy / 2, 0};
            Vector centroid = getCellCentroid(cellIdx);

            EXPECT_LT((centroid - expectedCentroid).norm(), tolerance) << std::format("cell idx is {}", cellIdx);
        }
    }
}

TEST_P(TestCartesinaMesh2D, CellVolumesArrayTest)
{
    ASSERT_EQ(m_cellVolumes.size(), totalCells);

    Scalar expectedCellVolume = dx * dy;
    for (Index cellIdx = 0; cellIdx < totalCells; cellIdx++)
    {
        Scalar cellVolume = getCellVolume(cellIdx);

        EXPECT_NEAR(cellVolume, expectedCellVolume, tolerance) << std::format("cell idx is {}", cellIdx);
    }
}

TEST_P(TestCartesinaMesh2D, CellFacesArrayTest)
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

            EXPECT_EQ(faces, expectedFaces) << std::format("cell idx is {}", cellIdx);
        }
    }
}

TEST_P(TestCartesinaMesh2D, FaceVectorArrayTest)
{
    ASSERT_EQ(m_faceVectors.size(), totalFaces);

    for (Index faceIdx = 0; faceIdx < totalFaces; faceIdx++)
    {
        Vector faceVector = getFaceVector(faceIdx);
        Vector expectedFaceVector = (faceIdx % (2 * nx + 1) < nx ? Vector{0, dx, 0} : Vector{dy, 0, 0});
        if (faceIdx < nx || faceIdx % (2 * nx + 1) == nx)
            expectedFaceVector *= -1;

        EXPECT_LT((faceVector - expectedFaceVector).norm(), tolerance) << std::format("face idx is {}", faceIdx);
    }
}

TEST_P(TestCartesinaMesh2D, FaceCentroidArrayTest)
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

        EXPECT_LT((centroid - expectedCentroid).norm(), tolerance) << std::format("face idx is {}", faceIdx);
    }
}

TEST_P(TestCartesinaMesh2D, FaceNeighborsArrayTest)
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

        EXPECT_EQ(neighbours, expectedNeighbours) << std::format("face idx is {}", faceIdx);
    }
}

TEST_P(TestCartesinaMesh2D, LeftBoundariesTest)
{
    using enum BoundaryConditionType;
    Boundaries boundaries = {{{1, 2, 3}, FIXED_GRADIENT}, {-3, FIXED_VALUE}};

    setLeftBoundary(boundaries);
    EXPECT_EQ(m_boundariesMap.size(), ny);

    for (Index faceIdx = nx; faceIdx < totalFaces; faceIdx += 2 * nx + 1)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, boundaries);
}

TEST_P(TestCartesinaMesh2D, RightBoundariesTest)
{
    using enum BoundaryConditionType;
    Boundaries boundaries = {{{1, 2, 3}, FIXED_GRADIENT}, {-3, FIXED_VALUE}};

    setRightBoundary(boundaries);
    EXPECT_EQ(m_boundariesMap.size(), ny);

    for (Index faceIdx = 2 * nx; faceIdx < totalFaces; faceIdx += 2 * nx + 1)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, boundaries);
}

TEST_P(TestCartesinaMesh2D, TopBoundariesTest)
{
    using enum BoundaryConditionType;
    Boundaries boundaries = {{{1, 2, 3}, FIXED_GRADIENT}, {-3, FIXED_VALUE}};

    setTopBoundary(boundaries);
    EXPECT_EQ(m_boundariesMap.size(), nx);

    for (Index faceIdx = 0; faceIdx < nx; faceIdx++)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, boundaries);
}

TEST_P(TestCartesinaMesh2D, BottomBoundariesTest)
{
    using enum BoundaryConditionType;
    Boundaries boundaries = {{{1, 2, 3}, FIXED_GRADIENT}, {-3, FIXED_VALUE}};

    setBottomBoundary(boundaries);
    EXPECT_EQ(m_boundariesMap.size(), nx);

    for (Index faceIdx = totalFaces - nx; faceIdx < totalFaces; faceIdx++)
        EXPECT_MAP_VALUE(m_boundariesMap, faceIdx, boundaries);
}

INSTANTIATE_TEST_SUITE_P(TestParamTest, TestCartesinaMesh2D,
                         testing::Values(CartesianMeshTestParam{2, 2, 1.0, 2.3}, CartesianMeshTestParam{4, 9, 1.0, 1.0},
                                         CartesianMeshTestParam{21, 21, 10.0, 10.0},
                                         CartesianMeshTestParam{16, 6, 0.4, 0.2}));