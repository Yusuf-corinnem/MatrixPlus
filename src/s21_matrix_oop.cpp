#include "s21_matrix_oop.h"

// FUNCTIONS
bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
          return false;
        }
      }
    }
  } else
    return false;
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ && cols_ != other.cols_)
    throw invalid_argument("Incorrect matrix");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ && cols_ != other.cols_)
    throw invalid_argument("Incorrect matrix");
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) throw invalid_argument("Incorrect matrix");
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (cols_ != rows_) throw invalid_argument("Incorrect matrix");
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[i][j] = pow(-1, (i + j)) * Minor(i, j).Determinant();
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (cols_ != rows_) throw invalid_argument("Incorrect matrix");
  double result = 0;
  if (rows_ == 1)
    result = matrix_[0][0];
  else if (rows_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    int sign = 1;
    for (int i = 0; i < cols_; i++) {
      result += matrix_[0][i] * Minor(0, i).Determinant() * sign;
      sign = -sign;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) throw invalid_argument("Incorrect matrix");
  double det = Determinant();
  if (fabs(det) < 1e-7) throw invalid_argument("Incorrect matrix");

  S21Matrix inverse(rows_, cols_);
  inverse = CalcComplements().Transpose();
  inverse.MulNumber(1 / det);

  return inverse;
}

S21Matrix S21Matrix::Minor(int i_rows, int i_cols) {
  if (i_rows < 0 || i_rows > 10 || i_rows >= rows_ || i_cols < 0 ||
      i_cols > 10 || i_cols >= cols_)
    throw invalid_argument("Incorrect matrix");
  int rows = rows_ - 1, cols = cols_ - 1;
  if (rows == 0) rows = 1;
  if (cols == 0) cols = 1;
  S21Matrix minor(rows, cols);
  int minor_rows = 0;
  for (int i = 0; i < rows_; i++) {
    if (i == i_rows) continue;
    int minor_cols = 0;
    for (int j = 0; j < cols_; j++) {
      if (j == i_cols) continue;
      minor.matrix_[minor_rows][minor_cols] = matrix_[i][j];
      minor_cols++;
    }
    minor_rows++;
  }
  return minor;
}

// CONSTRUCTORS
S21Matrix::S21Matrix() : rows_(1), cols_(1) { Alloc(&matrix_, rows_, cols_); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  Alloc(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

// DESTRUCTOR
S21Matrix::~S21Matrix() {
  if (matrix_) {
    Del(matrix_);
  }
}

// OPERATORS
double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || j < 0) throw invalid_argument("invalid size");
  return matrix_[i][j];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result = *this;
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result = *this;
  MulNumber(num);
  return result;
}

S21Matrix S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix(other.rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return S21Matrix();
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

void S21Matrix::operator+=(const S21Matrix& other) { return SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { return SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { return MulMatrix(other); }

void S21Matrix::operator*=(const double num) { return MulNumber(num); }

// ACESSORS
int S21Matrix::GetRows() { return rows_; }
int S21Matrix::GetCols() { return cols_; }

// MUTATORS
void S21Matrix::SetRows(int rows) {
  double** new_matrix;
  Alloc(&new_matrix, rows, cols_);
  if (rows > rows_) {
    for (int i = rows_; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        new_matrix[i][j] = 0;
      }
    }
  } else if (rows < rows_) {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        new_matrix[i][j] = matrix_[i][j];
      }
    }
  }
  Del(matrix_);
  rows_ = rows;
  matrix_ = new_matrix;
}
void S21Matrix::SetCols(int cols) {
  double** new_matrix;
  Alloc(&new_matrix, rows_, cols);
  if (cols > cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = cols_; j < cols; j++) {
        new_matrix[i][j] = 0;
      }
    }
  } else if (cols < cols_) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        new_matrix[i][j] = matrix_[i][j];
      }
    }
  }
  Del(matrix_);
  cols_ = cols;
  matrix_ = new_matrix;
}

// TOOLS
void S21Matrix::Print() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      if (i > 0) cout << endl;
      for (int j = 0; j < cols_; j++) {
        cout << matrix_[i][j];
      }
    }
    cout << endl;
  }
}

void S21Matrix::Alloc(double*** matrix, int rows, int cols) {
  if (rows < 1 || cols < 1) throw invalid_argument("Incorrect matrix");
  *matrix = new double*[rows]();
  for (int i = 0; i < rows; i++) {
    (*matrix)[i] = new double[cols]();
  }
}

void S21Matrix::Del(double** matrix) {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;
}