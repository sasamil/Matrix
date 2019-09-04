// Class name:   Matrix
// Authors:      NN(design 75%, implementation 15%)
//              Sasa Milenkovic(design 25%, implementation 85%)
// Description:  This class performs most of common matrix operations

#ifndef _lmjjsMatrix
#define _lmjjsMatrix

#include <iostream>
#include <fstream>
#include <string.h>

#define _error 0

using namespace std;
//---------------------------------------------------------------------------
inline void error(const char *str)
{
  cout << str << endl;
  exit(_error);
}

//---------------------------------------------------------------------------
class Matrix
{
public:
  Matrix();
  Matrix(int m, int n);
  Matrix(const Matrix &);
  ~Matrix();

  void release();
  const int getRows() const
  {
    return rows;
  };
  const int getCols() const
  {
    return cols;
  };

  Matrix &operator=(const Matrix &);
  friend int operator==(const Matrix &, const Matrix &);
  friend int operator!=(const Matrix &, const Matrix &);
  const double &operator()(int i, int j) const
  {
    return mat[i - 1][j - 1];
  };
  double &operator()(int i, int j)
  {
    return mat[i - 1][j - 1];
  };

  friend Matrix E(int);
  friend Matrix operator+(const Matrix &);
  friend Matrix operator-(const Matrix &);
  friend Matrix operator+(const Matrix &, const Matrix &);
  friend Matrix operator-(const Matrix &, const Matrix &);
  friend Matrix operator*(const Matrix &, const Matrix &);
  friend Matrix operator*(const Matrix &, const double &);
  friend Matrix operator*(const double &, const Matrix &);
  friend Matrix operator/(const Matrix &, const double &);

  friend double innerproduct(const Matrix &, const Matrix &);
  friend double trag(const Matrix &);
  friend Matrix trans(const Matrix &);
  friend Matrix comp(const Matrix &, int i, int j);
  friend double det(const Matrix &);
  friend double cofact(const Matrix &, int i, int j);
  friend Matrix adj(const Matrix &);
  friend Matrix inv(const Matrix &);
  friend Matrix invh(const Matrix &);

  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);
  Matrix &operator*=(const Matrix &);
  Matrix &operator*=(const double &);

  friend istream &operator>>(istream &, Matrix &);
  friend ostream &operator<<(ostream &, const Matrix &);

  friend Matrix gaussj(const Matrix &);
  friend Matrix cholesky(const Matrix &);
  friend Matrix lowerL(const Matrix &);
  friend Matrix Lyb(const Matrix &, const Matrix &);
  friend Matrix Ltxy(const Matrix &, const Matrix &);
  friend pair<Matrix, Matrix> ludecomp(const Matrix &);
  friend Matrix factor(const Matrix &A);
  friend Matrix factor2(const Matrix &A);
  friend Matrix solveUC(const Matrix &, const Matrix &);
  friend Matrix solveH(const Matrix &, const Matrix &);
  friend Matrix solveLU(const Matrix &, const Matrix &);
  Matrix submatrix(int srow, int erow, int scol, int ecol) const;

protected:
  void copy(const Matrix &);
  void allocate();

private:
  int cols, rows;
  double **mat;
};

#endif
