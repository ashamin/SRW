#include "gtest/gtest.h"
#include "../areas/sampleArea2d.h"
#include "../precond/SSORpar.h"

#include <iostream>

#define abs_error 10e-6

class AreaTest : public ::testing::Test{
protected:
    virtual void SetUp(){
        I = J = 30;
        test_area = new sampleArea2d(I, J);
    }

    virtual void TearDown(){
        delete test_area;
    }

    int I, J;
    sampleArea2d* test_area;
};

/*
TODO:
	2) add zeros to all diags
	3) tests module stucture
    4) operation dx_l[i] = dx_u[i] = value - may be undefined result
    5) ASSERT_DOUBLE_EQ => ASSERT_NEAR 
*/

TEST_F (AreaTest, MainDiag) { 
    for (int i = 0; i < test_area->getN(); i++){
        ASSERT_NEAR (-3364, test_area->getAp()[i], abs_error);
    }
}

TEST_F (AreaTest, NorthDiag){
    for (int i = 0; i < test_area->getN()-1; i++){
        ASSERT_NEAR (841, test_area->getAn()[i], abs_error);
    }
    ASSERT_NEAR (0, test_area->getAn()[test_area->getN()-1], abs_error);
}


TEST_F (AreaTest, SouthDiag){
    ASSERT_DOUBLE_EQ (0, test_area->getAs()[0]);
    for (int i = 1; i < test_area->getN(); i++){
        ASSERT_DOUBLE_EQ (841, test_area->getAs()[i]);
    }
}

TEST_F (AreaTest, EastDiag){
    for (int i = 0; i < test_area->getN(); i++){
        ASSERT_DOUBLE_EQ (841, test_area->getAe()[i]);
    }
}

TEST_F (AreaTest, WestDiag){
    for (int i = 0; i < test_area->getN(); i++){
        ASSERT_DOUBLE_EQ (841, test_area->getAw()[i]);
    }
}


class PrecondTest : public ::testing::Test{
protected:
    virtual void SetUp(){
        area = new sampleArea2d(30, 30);
        n = area->getN();
        precond = new SSORpar(.4, n, area->h1);
    }

    virtual void TearDown(){
        delete area;
        delete precond;
    }

    int n;
    sampleArea2d* area;
    SSORpar* precond;
};

TEST_F (PrecondTest, DX_D){
    for (int i = 0; i<n; i++)
        ASSERT_DOUBLE_EQ (-8410, precond->dx_d[i]);
}

TEST_F (PrecondTest, DX_L){
    ASSERT_DOUBLE_EQ (0, precond->dx_l[0]);
    for (int i = 1; i<n; i++)
        ASSERT_DOUBLE_EQ (1682, precond->dx_l[i]);
}

TEST_F (PrecondTest, DX_U){
    for (int i = 0; i<n-1; i++)
        ASSERT_DOUBLE_EQ (1682, precond->dx_u[i]);
    ASSERT_DOUBLE_EQ (0, precond->dx_u[n-1]);
}

TEST_F (PrecondTest, DY_D){
    for (int i = 0; i<n; i++)
        ASSERT_NEAR (0.4, precond->dy_d[i], abs_error);
}

TEST_F (PrecondTest, DY_L){
    ASSERT_DOUBLE_EQ (0, precond->dy_l[0]);
    for (int i = 1; i<n; i++)
        ASSERT_NEAR (-0.08, precond->dy_l[i], abs_error);
}

TEST_F (PrecondTest, DY_U){
    for (int i = 0; i<n-1; i++)
        ASSERT_NEAR (-0.08, precond->dy_u[i], abs_error);
    ASSERT_DOUBLE_EQ (0, precond->dy_u[n-1]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}