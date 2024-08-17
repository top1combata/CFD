#include <gtest/gtest.h>
#include <Discretization/LinearCombination/LinearCombination.h>
#include <algorithm>
#include <random>


template<class T> requires std::convertible_to<T,Scalar>
static Scalar norm(T value)
{
    return std::abs(value);
}

template<class T> requires requires(T value) {{value.norm()};}
static Scalar norm(T const& value)
{
    return value.norm();
}

template<class T>
static bool operator==(Term<T> const& t1, Term<T> const& t2)
{
    constexpr Scalar tolerance = 1e-5;
    return t1.idx == t2.idx && norm(t1.coeff - t2.coeff) < tolerance;
}

template<class U, class V>
static bool operator==(LinearCombination<U,V> const& lc1, LinearCombination<U,V> const& lc2)
{
    constexpr Scalar tolerance = 1e-5;
    return norm(lc1.bias - lc2.bias) < tolerance && lc1.terms == lc2.terms;
}

template<class T>
static std::ostream& operator<<(std::ostream& os, Term<T> const& term)
{
    return os << "{" << term.coeff << "," << term.idx << "}";
}

template<class U, class V>
void sortIndices(LinearCombination<U,V>& lc)
{
    std::sort(lc.terms.begin(), lc.terms.end(), [](auto t1, auto t2){return t1.idx < t2.idx;});
}

TEST(TestLinearCombination, ValidConstructScalar)
{
    int i = 3;
    double d = 5;
    unsigned short us = 17;
    LinearCombination<Scalar> lci(i), lcd(d), lcus(us);

    EXPECT_EQ(lci.bias, i);
    EXPECT_EQ(lci.terms.size(), 0);

    EXPECT_EQ(lcd.bias, d);
    EXPECT_EQ(lcd.terms.size(), 0);

    EXPECT_EQ(lcus.bias, us);
    EXPECT_EQ(lcus.terms.size(), 0);

    LinearCombination<Scalar> lc({{3.5, 0}, {7,2}, {1.5, 1}});
    EXPECT_EQ(lc.bias, 0);
    sortIndices(lc);
    EXPECT_EQ(lc.terms, std::vector<Term<Scalar>>({{3.5, 0}, {1.5, 1}, {7,2}}));
}

TEST(TestLinearCombination, ValidConstructVector)
{
    {
    Vector v1 = {1,2,3};
    Tensor t{{1,2,3},{2,4,6},{3,6,9}};

    LinearCombination<Vector, Vector> lc(v1*v1.transpose());
    EXPECT_EQ(lc.terms.size(), 0);
    EXPECT_EQ(lc.bias, t);
    }

    {
    LinearCombination<Scalar,Vector> lc = {{{1,-1,-1},3}, {{-2,7,10}, 1}};
    EXPECT_EQ(lc.bias, Vector(0,0,0));
    sortIndices(lc);
    EXPECT_EQ(lc.terms, std::vector<Term<Vector>>({ {{-2,7,10}, 1}, {{1,-1,-1},3}}));
    }
}

TEST(TestLinearCombination, AddAssign)
{
    LinearCombination<Scalar> lc1 = {{3.5, 0}, {2,1}, {1.5,2}};
    lc1 += Scalar(10);
    ASSERT_EQ(lc1.bias, 10);
    lc1 += -3.5;
    ASSERT_EQ(lc1.bias, 6.5);
    sortIndices(lc1);
    ASSERT_EQ(lc1.terms, std::vector<Term<Scalar>>({{3.5,0},{2,1},{1.5,2}}));

    LinearCombination<Scalar> lc = lc1;

    // full overlap
    lc += LinearCombination<Scalar>({{2.5, 0},{4,2}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{6,0},{2,1},{5.5,2}}));

    // partial overlap
    lc += LinearCombination<Scalar>({{1.5, 1},{2.5,0},{3,4}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{8.5,0},{3.5,1},{5.5,2},{3,4}}));

    // no overlap
    lc += LinearCombination<Scalar>({{1.25,5},{5,3}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, 6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{8.5,0},{3.5,1},{5.5,2},{5,3},{3,4},{1.25,5}}));
}


TEST(TestLinearCombination, SubtractAssign)
{
    LinearCombination<Scalar> lc1 = {{3.5, 0}, {2,1}, {1.5,2}};
    lc1 -= 10;
    ASSERT_EQ(lc1.bias, -10);
    lc1 -= -3.5;
    ASSERT_EQ(lc1.bias, -6.5);
    sortIndices(lc1);
    ASSERT_EQ(lc1.terms, std::vector<Term<Scalar>>({{3.5,0},{2,1},{1.5,2}}));

    LinearCombination<Scalar> lc = lc1;

    // full overlap
    lc -= LinearCombination<Scalar>({{2.5, 0},{4,2}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{1,0},{2,1},{-2.5,2}}));

    // partial overlap
    lc -= LinearCombination<Scalar>({{1.5, 1},{2.5,0},{3,4}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{-1.5,0},{0.5,1},{-2.5,2},{-3,4}}));

    // no overlap
    lc -= LinearCombination<Scalar>({{1.25,5},{5,3}});
    sortIndices(lc);
    ASSERT_EQ(lc.bias, -6.5);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{-1.5,0},{0.5,1},{-2.5,2},{-5,3},{-3,4},{-1.25,5}}));
}


TEST(TestLinearCombination, nonMemberOperators)
{
    auto randomCombination = []()
    {
        Index sz = 4;
        std::mt19937 generator{(std::size_t)std::rand()};
        std::uniform_int_distribution<int> distribution(-30, 30);
        std::uniform_int_distribution<int> idxDistribution(0, 10);

        LinearCombination<Scalar, Scalar> res;
        res.bias = distribution(generator);
        res.terms.resize(sz);
        for (Index idx = 0; idx < sz; idx++)
            res.terms[idx] = {(Scalar)distribution(generator), idxDistribution(generator)};
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
        Scalar s = std::rand();
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
        auto lc1 = randomCombination();
        auto lc2 = lc1;
        Scalar s = std::rand();
        lc2 -= s;
        EXPECT_EQ(lc1 - s, lc2);
        EXPECT_EQ(s - lc1, -lc2);
    }
    {
        auto lc1 = randomCombination();
        auto lc2 = lc1;
        Scalar s = std::rand();
        lc2 *= s;
        EXPECT_EQ(lc1 * s, lc2);
        EXPECT_EQ(s * lc1, lc2);
    }
}


TEST(TestLinearCombination, MulDivAssign)
{
    LinearCombination<Scalar> lc = {{3.5, 0}, {2,1}, {1.5,2}};
    // implicit cast
    lc += 10;
    ASSERT_EQ(lc.bias, 10);
    sortIndices(lc);
    ASSERT_EQ(lc.terms, std::vector<Term<Scalar>>({{3.5,0},{2,1},{1.5,2}}));

    // *=(Scalar)
    auto lc1 = lc;
    lc1 *= -2;
    EXPECT_EQ(lc1.bias, -20);
    sortIndices(lc1);
    EXPECT_EQ(lc1.terms, std::vector<Term<Scalar>>({{-7,0},{-4,1},{-3,2}}));

    // /=(Scalar)
    auto lc2 = lc;
    lc2 /= 2;
    EXPECT_EQ(lc2.bias, 5);
    sortIndices(lc2);
    EXPECT_EQ(lc2.terms, std::vector<Term<Scalar>>({{1.75,0},{1,1},{.75,2}}));

    auto lc3 = lc*3.5;
    EXPECT_EQ(lc3.bias, 35);
    sortIndices(lc3);
    EXPECT_EQ(lc3.terms, std::vector<Term<Scalar>>({{12.25,0},{7,1},{5.25,2}}));
}


TEST(TestLinearCombination, MulScalarByVector)
{
    LinearCombination<Scalar, Scalar> lc = {{0.5, 2}, {1.5, 1}, {3, 3}};
    lc -= 25;
    auto vlc = lc * Vector{-1,3,7};
    EXPECT_EQ(vlc.bias, Vector(25, -75, -175));
    sortIndices(vlc);
    EXPECT_EQ(vlc.terms, std::vector<Term<Vector>>({ {{-1.5,4.5,10.5}, 1}, {{-0.5,1.5,3.5}, 2}, {{-3,9,21}, 3} }));

    auto vlc2 = Vector{-1,3,7} * lc;
    EXPECT_EQ(vlc2.bias, Vector(25, -75, -175));
    sortIndices(vlc2);
    EXPECT_EQ(vlc.terms, std::vector<Term<Vector>>({ {{-1.5,4.5,10.5}, 1}, {{-0.5,1.5,3.5}, 2}, {{-3,9,21}, 3} }));
}


TEST(TestLinearCombination, MulVectorByVector)
{
    LinearCombination<Vector, Scalar> lc = {{0.5, 2}, {1.5, 1}, {3, 3}};
    lc += {1,-2,3};
    auto vlc = lc * Vector{-1,3,7};
    EXPECT_EQ(vlc.bias, Tensor({{-1,3,7},{2,-6,-14},{-3,9,21}}));
    sortIndices(vlc);
    EXPECT_EQ(vlc.terms, std::vector<Term<Vector>>({ {{-1.5,4.5,10.5}, 1}, {{-0.5,1.5,3.5}, 2}, {{-3,9,21}, 3} }));
}


TEST(TestLinearCombination, EvaluateOnScalarField)
{
    Field<Scalar> field(6);

    LinearCombination<Scalar, Scalar> lc1 = {{0.5, 0}, {-0.75,2},{2,3}};

    lc1 += -3;
    field << 2,1,-1,3,2,1;
    EXPECT_EQ(lc1.evaluate(field), 4.75);

    lc1 += 5.25;
    field << 1,-1,1,0.5,1.5,-2;
    EXPECT_EQ(lc1.evaluate(field), 3);

    LinearCombination<Scalar, Vector> lc2 = {{{1,1,1},1}, {{-2.5,0.5,-1},3}, {{0,0.25,2},4}};

    lc2 += {-1,1,1};
    field << 1,1,1,-1,-1,-1;
    EXPECT_EQ(lc2.evaluate(field), Vector(2.5, 1.25, 1));

    lc2 -= {0.5,0.25,-1};
    field << -2,-3,1.5,1.25,1.75,-0.5;
    EXPECT_EQ(lc2.evaluate(field), Vector(-7.625, -1.1875, 1.25));
}


TEST(TestLinearCombination, EvaluateOnVectorField)
{
    Vector v1(1,-1,2), v2(-0.5,3,-1.75), v3(4,-5,8), v4(-2,-2,8), v5(-4,-10,-7);
    Tensor t{{1,1,1},{-1,-2,-3},{2,5,-3}};

    LinearCombination<Vector, Vector> lc = {{v1,1},{v2,3},{v3,0}};
    lc += t;

    Field<Vector> field(5);
    field << v1,v2,v3,v4,v5;

    Tensor ans = t + v2*v1.transpose() + v1*v3.transpose() + v4*v2.transpose();

    EXPECT_EQ(lc.evaluate(field), ans);
}


TEST(TestLinearCombination, DotProductOnScalarVariable)
{
    Vector v1(1,-1,2), v2(-0.5,3,-1.75), v3(4,-5,8), v4(-2,-2,8), v5(-4,-10,-7);
    LinearCombination<Scalar, Vector> lc = {{v1,1},{v3,0},{v5,3}};
    lc += v5;

    auto slc = lc.dot(v2);
    EXPECT_EQ(slc.bias, v2.dot(v5));
    sortIndices(slc);
    EXPECT_EQ(slc.terms, std::vector<Term<Scalar>>({{v2.dot(v3),0}, {v1.dot(v2),1}, {v2.dot(v5),3}}));

    slc = lc.dot(v4);
    EXPECT_EQ(slc.bias, v4.dot(v5));
    sortIndices(slc);
    EXPECT_EQ(slc.terms, std::vector<Term<Scalar>>({{v4.dot(v3),0}, {v1.dot(v4),1}, {v4.dot(v5),3}}));
}


TEST(TestLinearCombination, DotProductOnVectorVariable)
{
    Vector v1(-3,3,2), v2(-2,3.5,0.25), v3(-4,-5,-3), v4(1,4,8), v5(1,-1,1);
    LinearCombination<Vector, Vector> lc = {{v2,0},{v4,2},{v5,4}};
    Tensor t{{-1,2,2},{4,1,8},{-3,2,7}};
    lc += t;

    auto slc = lc.dot(v1);
    EXPECT_EQ(slc.bias, t*v1);
    sortIndices(slc);
    EXPECT_EQ(slc.terms, std::vector<Term<Scalar>>({{v1.dot(v2),0}, {v1.dot(v4),2}, {v1.dot(v5),4}}));

    slc = lc.dot(v3);
    EXPECT_EQ(slc.bias, t*v3);
    sortIndices(slc);
    EXPECT_EQ(slc.terms, std::vector<Term<Scalar>>({{v2.dot(v3),0}, {v3.dot(v4),2}, {v3.dot(v5),4}}));
}