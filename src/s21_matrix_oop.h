#ifndef _S21_MATRIX_OOP_H_
#define _S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>

using namespace std;

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  // FUNCTIONS
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  S21Matrix Minor(int i_rows, int i_cols);
  double Determinant();
  S21Matrix InverseMatrix();

  // CONSTRUCTORS
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);

  // DESTRUCTOR
  ~S21Matrix();

  // OPERATORS
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix operator=(const S21Matrix& other);
  bool operator==(const S21Matrix& other) const;
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(const double num);
  double& operator()(int i, int j);

  // TOOLS
  void Print();
  void Alloc(double*** matrix, int rows, int cols);
  void Del(double** matrix);

  // ACESSORS
  int GetRows();
  int GetCols();

  // MUTATORS
  void SetRows(int rows);
  void SetCols(int cols);
};

#endif