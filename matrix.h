//Class name:   Matrix
//Authors:      NN(design 75%, implementation 15%)
//              Sasha Milenkovic(design 25%, implementation 85%)
//Date:         15.06.2002.
//Description:  This class performs most of common matrix operations
//
//Changes:      11.10.2002. Choleskey's decomposition added
//
// 8128402, 8989890 ms, 8982637, 8929549, 9155736, 627071, 9174708, 9195674, 9195446!
// 8406698, 7942164, 9147768, 8246641, 7930088, 9063498, 9063498, 8909241, 8882657, 9155863
// 8671116, 9335589, 9176783, 8284305
//
// 8490361, 8465071, 8305473,709291, 8742276, 8726472, 8882657
// 483100, 330697, 334045, 315672, 309501, 308656, 715761
// 8505881, 622905, 8559210, 8929550, 8553784, 8216692
// 

#ifndef _lmjjsMatrix
#define _lmjjsMatrix

//#include <stdio>
#include <iostream>
#include <fstream>
#include <string.h>


#define _error 0

using namespace std;
//---------------------------------------------------------------------------
inline void error(const char* str)
{
  cout << str << endl;
  //::MessageBox(NULL, str, "Matrix Operation Error", 0x10010);
  //MessageDlg(str, mtError, TMsgDlgButtons()<<mbOK, 0);
  exit(_error);
}

//---------------------------------------------------------------------------
class Matrix{
public:
  Matrix();
  Matrix(int m, int n); // matrica dimenzija mxn
  Matrix(const Matrix&);
  ~Matrix();

  void release();
  const int getRows() const {return rows;};
  const int getCols() const {return cols;};

  Matrix& operator= (const Matrix&);
  friend int operator==(const Matrix&, const Matrix&);
  friend int operator!=(const Matrix&, const Matrix&);
  const double& operator() (int i, int j) const {return mat[i-1][j-1];};
  double& operator()(int i, int j) {return mat[i-1][j-1];};

  friend Matrix E(int);
  friend Matrix operator+ (const Matrix&);
  friend Matrix operator- (const Matrix&);
  friend Matrix operator+ (const Matrix&, const Matrix&);
  friend Matrix operator- (const Matrix&, const Matrix&);
  friend Matrix operator* (const Matrix&, const Matrix&);
  friend Matrix operator* (const Matrix&, const double&);
  friend Matrix operator* (const double&, const Matrix&);
  friend Matrix operator/ (const Matrix&, const double&);

  friend double innerproduct (const Matrix&, const Matrix&);
  friend double trag(const Matrix&);
  friend Matrix trans(const Matrix&);
  friend Matrix comp(const Matrix&, int i, int j);
  friend double det(const Matrix&);
  friend double cofact(const Matrix&, int i, int j);
  friend Matrix adj(const Matrix&);
  friend Matrix inv(const Matrix&);
  friend Matrix invh(const Matrix&);

  Matrix& operator+= (const Matrix&);
  Matrix& operator-= (const Matrix&);
  Matrix& operator*= (const Matrix&);
  Matrix& operator*= (const double&);

  friend istream& operator>> (istream&, Matrix&);
  friend ostream& operator<< (ostream&, const Matrix&);

//friend Matrix gaussj(double**, int, double**, int);
  friend Matrix gaussj(const Matrix&); //##sm
  friend Matrix cholesky(const Matrix&); //##sm
  friend Matrix lowerL(const Matrix&); //##sm
  friend Matrix Lyb(const Matrix&, const Matrix&); //##sm
  friend Matrix Ltxy(const Matrix&, const Matrix&); //##sm
  friend pair<Matrix, Matrix> ludecomp(const Matrix&); //##sm
  friend Matrix factor(const Matrix& A); //###sm
  friend Matrix factor2(const Matrix& A); //###sm
  friend Matrix solveUC(const Matrix&, const Matrix&); //##sm
  friend Matrix solveH(const Matrix&, const Matrix&); //##sm
  friend Matrix solveLU(const Matrix&, const Matrix&); //##sm
  /*friend void fKolokacija(const Matrix&, const Matrix&,
                          const Matrix&, const Matrix&,
                          const Matrix&); //##sm*/
  Matrix submatrix(int srow, int erow, int scol, int ecol) const;//##sm

protected:
  void allocate();
  void copy(const Matrix&);

private:
  int cols, rows;
  double **mat;
};

#endif
