#include "gtest/gtest.h"
#include "../areas/sampleArea2d.h"

#include <iostream>

/*
TODO:
	1) make test fixture
	2) find out where zeros to n, s diags are added and why there are no zeros on e, w diags
	3) tests module stucture
*/

TEST (AreaSuite, MainDiag) { 
	int I = 30, J = 30;
	sampleArea2d* test_area = new sampleArea2d(I, J);
    for (int i = 0; i < test_area->getN(); i++){
    	ASSERT_DOUBLE_EQ (-3364, test_area->getAp()[i]);
    }
    delete test_area;
}

TEST (AreaSuite, NorthDiag){
	int I = 30, J = 30;
	sampleArea2d* test_area = new sampleArea2d(I, J);
    for (int i = 0; i < test_area->getN(); i++){
    	ASSERT_DOUBLE_EQ (841, test_area->getAn()[i]);
    }
    delete test_area;
}


TEST (AreaSuite, SouthDiag){
	int I = 30, J = 30;
	sampleArea2d* test_area = new sampleArea2d(I, J);
    for (int i = 0; i < test_area->getN(); i++){
    	ASSERT_DOUBLE_EQ (841, test_area->getAs()[i]);
    }
    delete test_area;
}

TEST (AreaSuite, EastDiag){
	int I = 30, J = 30;
	sampleArea2d* test_area = new sampleArea2d(I, J);
    for (int i = 0; i < test_area->getN(); i++){
    	ASSERT_DOUBLE_EQ (841, test_area->getAe()[i]);
    }
    delete test_area;
}

TEST (AreaSuite, WestDiag){
	int I = 30, J = 30;
	sampleArea2d* test_area = new sampleArea2d(I, J);
    for (int i = 0; i < test_area->getN(); i++){
    	ASSERT_DOUBLE_EQ (841, test_area->getAw()[i]);
    }
    delete test_area;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}