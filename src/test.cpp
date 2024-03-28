#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

void Filling(S21Matrix& other, const double* value) {
  int n = 0;
  for (int i = 0; i < other.GetRows(); i++) {
    for (int j = 0; j < other.GetCols(); j++) {
      other(i, j) = value[n++];
    }
  }
}

TEST(Matrix, Constructor) {
  S21Matrix one;
  one(0, 0) = 50.0000008;
  S21Matrix two(1, 1);
  two(0, 0) = 50.0000008;
  ASSERT_EQ(one, two);
}

TEST(Matrix, Copy) {
  S21Matrix one(2, 2);
  double two[] = {1.1, 2.2, 3.3, 4.4};
  Filling(one, two);
  S21Matrix copy(one);
  ASSERT_EQ(one, copy);
}

TEST(Matrix, Move) {
  S21Matrix one(2, 2);
  double two[] = {1.1, 2.2, 3.3, 4.4};
  Filling(one, two);
  S21Matrix copy(one);
  S21Matrix move(std::move(one));
  EXPECT_EQ(move, copy);
  ASSERT_FALSE(move == one);
}

TEST(Matrix, Operators) {
  S21Matrix one(3, 3);
  double two[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Filling(one, two);
  S21Matrix sum(3, 3);
  S21Matrix vrem(3, 3);
  sum = one + one;
  sum -= one;
  sum += one;
  vrem = one;
  vrem *= 2.0;
  EXPECT_EQ(vrem, sum);
  vrem = one;
  vrem * 2.0;
  EXPECT_EQ(vrem, sum);
  sum = sum - one;
  EXPECT_EQ(sum, one);
  sum = one * one;
  vrem = one;
  vrem *= one;
  EXPECT_EQ(vrem, sum);
}

TEST(Matrix, Inverse) {
  S21Matrix R(3, 3);
  double R_1[] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
  Filling(R, R_1);
  S21Matrix R2(3, 3);
  S21Matrix A(3, 3);
  double A_1[] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  Filling(A, A_1);
  R2 = A.InverseMatrix();
  ASSERT_TRUE(R2 == R);
}

TEST(Matrix, SetGetAdd) {
  S21Matrix one(3, 3);
  EXPECT_EQ(one.GetRows(), 3);
  EXPECT_EQ(one.GetCols(), 3);
  one.SetRows(1);
  EXPECT_EQ(one.GetRows(), 1);
  one.SetCols(10);
  EXPECT_EQ(one.GetCols(), 10);
  S21Matrix two(20, 20);
  two.SetRows(30);
  EXPECT_EQ(two.GetRows(), 30);
  two.SetCols(30);
  EXPECT_EQ(two.GetCols(), 30);
}

TEST(Matrix, Errors) {
  EXPECT_THROW(S21Matrix error(0, 0), std::invalid_argument);
  EXPECT_THROW(S21Matrix error(2, 0), std::invalid_argument);
  S21Matrix error_1(2, 4);
  S21Matrix error_2(5, 3);
  EXPECT_THROW(error_1 += error_2, std::invalid_argument);
  EXPECT_THROW(error_1 -= error_2, std::invalid_argument);
  EXPECT_THROW(error_1 *= error_2, std::invalid_argument);
  EXPECT_THROW(error_1.CalcComplements(), std::invalid_argument);
  EXPECT_THROW(error_1.Determinant(), std::invalid_argument);
  EXPECT_THROW(error_1.InverseMatrix(), std::invalid_argument);
  EXPECT_THROW(error_2(-2, 20), std::invalid_argument);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}