#include "TestUtils.h"
#include <Discretization/LinearCombination.h>
#include <algorithm>
#include <gtest/gtest.h>

#include <unordered_set>

TEST(TestLinearCombination, ValidConstructScalar)
{
    int                       i = 3;
    double                    d = 5;
    unsigned short            us = 17;
    LinearCombination<Scalar> lci(i), lcd(d), lcus(us);

    EXPECT_EQ(lci.bias, i);
    EXPECT_EQ(lci.terms.size(), 0);

    EXPECT_EQ(lcd.bias, d);
    EXPECT_EQ(lcd.terms.size(), 0);

    EXPECT_EQ(lcus.bias, us);
    EXPECT_EQ(lcus.terms.size(), 0);

    LinearCombination<Scalar> lc({{3.5, 0}, {7, 2}, {1.5, 1}});
    EXPECT_EQ(lc.bias, 0);
    sortIndices(lc);
    EXPECT_EQ(lc.terms, List<Term<Scalar>>({{3.5, 0}, {1.5, 1}, {7, 2}}));
}

TEST(TestLinearCombination, ValidConstructVector)
{
    {
        Vector v1 = {1, 2, 3};
        Tensor t{{1, 2, 3}, {2, 4, 6}, {3, 6, 9}};

        LinearCombination<Vector, Vector> lc(v1 * v1.transpose());
        EXPECT_EQ(lc.terms.size(), 0);
        EXPECT_EQ(lc.bias, t);
    }

    {
        LinearCombination<Scalar, Vector> lc = {{{1, -1, -1}, 3}, {{-2, 7, 10}, 1}};
        EXPECT_EQ(lc.bias, Vector(0, 0, 0));
        sortIndices(lc);
        EXPECT_EQ(lc.terms, List<Term<Vector>>({{{-2, 7, 10}, 1}, {{1, -1, -1}, 3}}));
    }
}

TEST(TestLinearCombination, AddAssign)
{
    LinearCombination<Scalar> lc1 = {{3.5, 0}, {2, 1}, {1.5, 2}};
    lc1 += Scalar(10);
    ASSERT_EQ(lc1.bias, 10);
    lc1 += -3.5;
    ASSERT_EQ(lc1.bias, 6.5);
    sortIndices(lc1);
    ASSERT_EQ(lc1.terms, List<Term<Scalar>>({{3.5, 0}, {2, 1}, {1.5, 2}}));

    LinearCombination<Scalar> lc = lc1;

    // full overlap
    lc += LinearCombination<Scalar>({{2.5, 0}, {4, 2}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{6, 0}, {2, 1}, {5.5, 2}}));

    // partial overlap
    lc += LinearCombination<Scalar>({{1.5, 1}, {2.5, 0}, {3, 4}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{8.5, 0}, {3.5, 1}, {5.5, 2}, {3, 4}}));

    // no overlap
    lc += LinearCombination<Scalar>({{1.25, 5}, {5, 3}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{8.5, 0}, {3.5, 1}, {5.5, 2}, {5, 3}, {3, 4}, {1.25, 5}}));
}

TEST(TestLinearCombination, SubtractAssign)
{
    LinearCombination<Scalar> lc1 = {{3.5, 0}, {2, 1}, {1.5, 2}};
    lc1 -= 10;
    ASSERT_EQ(lc1.bias, -10);
    lc1 -= -3.5;
    ASSERT_EQ(lc1.bias, -6.5);
    sortIndices(lc1);
    ASSERT_EQ(lc1.terms, List<Term<Scalar>>({{3.5, 0}, {2, 1}, {1.5, 2}}));

    LinearCombination<Scalar> lc = lc1;

    // full overlap
    lc -= LinearCombination<Scalar>({{2.5, 0}, {4, 2}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{1, 0}, {2, 1}, {-2.5, 2}}));

    // partial overlap
    lc -= LinearCombination<Scalar>({{1.5, 1}, {2.5, 0}, {3, 4}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{-1.5, 0}, {0.5, 1}, {-2.5, 2}, {-3, 4}}));

    // no overlap
    lc -= LinearCombination<Scalar>({{1.25, 5}, {5, 3}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{-1.5, 0}, {0.5, 1}, {-2.5, 2}, {-5, 3}, {-3, 4}, {-1.25, 5}}));
}

TEST(TestLinearCombination, nonMemberOperators)
{
    auto randomCombination = []()
    {
        constexpr Index sz = 4;

        constexpr Index coeffMin = -30;
        constexpr Index coeffMax = 30;

        constexpr Index idxMin = 0;
        constexpr Index idxMax = 10;

        LinearCombination<Scalar, Scalar> res;
        res.bias = randomIndex(coeffMin, coeffMax);
        res.terms.resize(sz);

        std::unordered_set<Index> usedIndices;
        for (Index i = 0; i < sz; i++)
        {
            Scalar coeff = randomIndex(coeffMin, coeffMax);
            Index  idx;
            while (usedIndices.count(idx = randomIndex(idxMin, idxMax)))
                ;
            usedIndices.insert(idx);
            res.terms[i] = {coeff, idx};
        }
        return res;
    };

    {
        auto lc1 = randomCombination();
        auto lc2 = randomCombination();
        auto lc3 = lc1;
        EXPECT_EQ(lc3 += lc2, lc1 + lc2);
    }
    {
        auto lc1 = randomCombination();
        auto lc2 = lc1;

        Scalar s = randomScalar();
        lc2 += s;

        EXPECT_EQ(lc1 + s, lc2);
        EXPECT_EQ(s + lc1, lc2);
    }
    {
        auto lc1 = randomCombination();
        auto lc2 = randomCombination();
        auto lc3 = lc1;
        EXPECT_EQ(lc3 -= lc2, lc1 - lc2);
    }
    {
        auto   lc1 = randomCombination();
        auto   lc2 = lc1;
        Scalar s = std::rand();
        lc2 -= s;
        EXPECT_EQ(lc1 - s, lc2);
        EXPECT_EQ(s - lc1, -lc2);
    }
    {
        auto   lc1 = randomCombination();
        auto   lc2 = lc1;
        Scalar s = std::rand();
        lc2 *= s;
        EXPECT_EQ(lc1 * s, lc2);
        EXPECT_EQ(s * lc1, lc2);
    }
}

TEST(TestLinearCombination, MulDivAssign)
{
    LinearCombination<Scalar> lc = {{3.5, 0}, {2, 1}, {1.5, 2}};
    // implicit cast
    lc += 10;
    ASSERT_EQ(lc.bias, 10);
    sortIndices(lc);
    ASSERT_EQ(lc.terms, List<Term<Scalar>>({{3.5, 0}, {2, 1}, {1.5, 2}}));

    // *=(Scalar)
    auto lc1 = lc;
    lc1 *= -2;
    EXPECT_EQ(lc1.bias, -20);
    sortIndices(lc1);
    EXPECT_EQ(lc1.terms, List<Term<Scalar>>({{-7, 0}, {-4, 1}, {-3, 2}}));

    // /=(Scalar)
    auto lc2 = lc;
    lc2 /= 2;
    EXPECT_EQ(lc2.bias, 5);
    sortIndices(lc2);
    EXPECT_EQ(lc2.terms, List<Term<Scalar>>({{1.75, 0}, {1, 1}, {.75, 2}}));

    auto lc3 = lc * Scalar(3.5);
    EXPECT_EQ(lc3.bias, 35);
    sortIndices(lc3);
    EXPECT_EQ(lc3.terms, List<Term<Scalar>>({{12.25, 0}, {7, 1}, {5.25, 2}}));
}

TEST(TestLinearCombination, MulScalarByVector)
{
    LinearCombination<Scalar, Scalar> lc = {{0.5, 2}, {1.5, 1}, {3, 3}};
    lc -= 25;
    auto vlc = lc * Vector{-1, 3, 7};
    EXPECT_EQ(vlc.bias, Vector(25, -75, -175));
    sortIndices(vlc);
    EXPECT_EQ(vlc.terms, List<Term<Vector>>({{{-1.5, 4.5, 10.5}, 1}, {{-0.5, 1.5, 3.5}, 2}, {{-3, 9, 21}, 3}}));

    auto vlc2 = Vector{-1, 3, 7} * lc;
    EXPECT_EQ(vlc2.bias, Vector(25, -75, -175));
    sortIndices(vlc2);
    EXPECT_EQ(vlc.terms, List<Term<Vector>>({{{-1.5, 4.5, 10.5}, 1}, {{-0.5, 1.5, 3.5}, 2}, {{-3, 9, 21}, 3}}));
}

TEST(TestLinearCombination, MulVectorByVector)
{
    {
        Scalar c1 = randomScalar(), c2 = randomScalar(), c3 = randomScalar();
        Vector bias = randomVector();

        LinearCombination<Vector, Scalar> lc = {{c3, 3}, {c2, 2}, {c1, 1}};
        lc += bias;

        Vector multiplier = randomVector();

        auto result = lc * multiplier;
        sortIndices(result);
        List<Term<Vector>> expectedTerms = {{c1 * multiplier, 1}, {c2 * multiplier, 2}, {c3 * multiplier, 3}};
        EXPECT_EQ(result.bias, lc.bias * multiplier.transpose());
        EXPECT_EQ(result.terms, expectedTerms);

        result = multiplier * lc;
        sortIndices(result);
        EXPECT_EQ(result.bias, multiplier * lc.bias.transpose());
        EXPECT_EQ(result.terms, expectedTerms);
    }
    {
        Vector v1 = randomVector(), v2 = randomVector(), v3 = randomVector();
        Vector bias = randomVector();

        LinearCombination<Scalar, Vector> lc = {{v2, 2}, {v3, 3}, {v1, 1}};
        lc += bias;

        Vector multiplier = randomVector();
        auto   transposed = multiplier.transpose();

        auto result = lc * multiplier;
        sortIndices(result);
        List<Term<Tensor>> expectedTerms = {{v1 * transposed, 1}, {v2 * transposed, 2}, {v3 * transposed, 3}};
        EXPECT_EQ(result.bias, bias * transposed);
        EXPECT_EQ(result.terms, expectedTerms);

        result = multiplier * lc;
        sortIndices(result);
        for (auto& [coeff, _] : expectedTerms)
            coeff.transposeInPlace();
        EXPECT_EQ(result.bias, multiplier * bias.transpose());
        EXPECT_EQ(result.terms, expectedTerms);
    }
}

TEST(TestLinearCombination, EvaluateOnScalarField)
{
    Field<Scalar> field(6);

    LinearCombination<Scalar, Scalar> lc1 = {{0.5, 0}, {-0.75, 2}, {2, 3}};

    lc1 += -3;
    field << 2, 1, -1, 3, 2, 1;
    EXPECT_EQ(lc1.evaluate(field), 4.75);

    lc1 += 5.25;
    field << 1, -1, 1, 0.5, 1.5, -2;
    EXPECT_EQ(lc1.evaluate(field), 3);

    LinearCombination<Scalar, Vector> lc2 = {{{1, 1, 1}, 1}, {{-2.5, 0.5, -1}, 3}, {{0, 0.25, 2}, 4}};

    lc2 += {-1, 1, 1};
    field << 1, 1, 1, -1, -1, -1;
    EXPECT_EQ(lc2.evaluate(field), Vector(2.5, 1.25, 1));

    lc2 -= {0.5, 0.25, -1};
    field << -2, -3, 1.5, 1.25, 1.75, -0.5;
    EXPECT_EQ(lc2.evaluate(field), Vector(-7.625, -1.1875, 1.25));
}

TEST(TestLinearCombination, EvaluateOnVectorField)
{
    Vector v1(1, -1, 2), v2(-0.5, 3, -1.75), v3(4, -5, 8), v4(-2, -2, 8), v5(-4, -10, -7);
    Tensor t{{1, 1, 1}, {-1, -2, -3}, {2, 5, -3}};

    LinearCombination<Vector, Vector> lc = {{v1, 1}, {v2, 3}, {v3, 0}};
    lc += t;

    Field<Vector> field(5);
    field << v1, v2, v3, v4, v5;

    Tensor ans = t + v2 * v1.transpose() + v1 * v3.transpose() + v4 * v2.transpose();

    EXPECT_EQ(lc.evaluate(field), ans);
}

TEST(TestLinearCombination, DotProductOnScalarVariable)
{
    Vector                            v1(1, -1, 2), v2(-0.5, 3, -1.75), v3(4, -5, 8), v4(-2, -2, 8), v5(-4, -10, -7);
    LinearCombination<Scalar, Vector> lc = {{v1, 1}, {v3, 0}, {v5, 3}};
    lc += v5;

    auto slc = lc.dot(v2);
    EXPECT_EQ(slc.bias, v2.dot(v5));
    sortIndices(slc);
    EXPECT_EQ(slc.terms, List<Term<Scalar>>({{v2.dot(v3), 0}, {v1.dot(v2), 1}, {v2.dot(v5), 3}}));

    slc = lc.dot(v4);
    EXPECT_EQ(slc.bias, v4.dot(v5));
    sortIndices(slc);
    EXPECT_EQ(slc.terms, List<Term<Scalar>>({{v4.dot(v3), 0}, {v1.dot(v4), 1}, {v4.dot(v5), 3}}));
}

TEST(TestLinearCombination, DotProductOnVectorVariable)
{
    Vector                            v1(-3, 3, 2), v2(-2, 3.5, 0.25), v3(-4, -5, -3), v4(1, 4, 8), v5(1, -1, 1);
    LinearCombination<Vector, Vector> lc = {{v2, 0}, {v4, 2}, {v5, 4}};
    Tensor                            t{{-1, 2, 2}, {4, 1, 8}, {-3, 2, 7}};
    lc += t;

    auto slc = lc.dot(v1);
    EXPECT_EQ(slc.bias, t * v1);
    sortIndices(slc);
    EXPECT_EQ(slc.terms, List<Term<Scalar>>({{v1.dot(v2), 0}, {v1.dot(v4), 2}, {v1.dot(v5), 4}}));

    slc = lc.dot(v3);
    EXPECT_EQ(slc.bias, t * v3);
    sortIndices(slc);
    EXPECT_EQ(slc.terms, List<Term<Scalar>>({{v2.dot(v3), 0}, {v3.dot(v4), 2}, {v3.dot(v5), 4}}));
}