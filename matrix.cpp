//Class name:   Matrix
//Authors:      NN(design 75%, implementation 15%)
//              Sasa Milenkovic(design 25%, implementation 85%)
//Description:  This class performs most of common matrix operations

#include <stdlib.h>
#include <math.h>

#include "matrix.h"

//---------------------------------------------------------------------------
inline void swap(double& a, double& b)
{
  double temp=a;
  a=b;
  b=temp;
}

//---------------------------------------------------------------------------
void Matrix::allocate()
{
  mat = new double*[rows];
  for(int i=0; i<rows; i++)
  {
    mat[i] = new double[cols];
    for(int j=0; j<cols; j++)
      mat[i][j] = 0.0;
  }
}
//---------------------------------------------------------------------------
void Matrix::release()
{
  for(int i=0; i<rows; i++)
  {
    delete[] mat[i];
  }
  delete[] mat;
  rows = cols = 1;
}
//---------------------------------------------------------------------------
void Matrix::copy(const Matrix& right)
{
  if(right.cols != cols || right.rows != rows)
    error("Matrix::copy:\nBroj redova/kolona nije isti!");

  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      mat[i][j] = right.mat[i][j];
}
//---------------------------------------------------------------------------
Matrix::Matrix() : cols(1), rows(1)
{
  allocate();
}
//---------------------------------------------------------------------------
Matrix::Matrix(int m, int n) : cols(n), rows(m)
{
  allocate();
}
//---------------------------------------------------------------------------
Matrix::Matrix(const Matrix& right) : cols(right.cols), rows(right.rows)
{
  allocate();
  copy(right);
}
//---------------------------------------------------------------------------
Matrix::~Matrix()
{
  release();
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator= (const Matrix& right)
{
  if(&right == this)
    return *this;

  release();
  rows=right.rows;
  cols=right.cols;
  allocate();
  copy(right);
  return *this;
}
//---------------------------------------------------------------------------
Matrix Matrix::submatrix(int srow, int erow, int scol, int ecol) const
{

  if(srow>erow || scol>ecol || srow<1 || erow>rows ||scol<1 || ecol>cols)
    error("Matrix::submatrix()\nNeodgovarajuci argumenti!");

  Matrix m(erow-srow+1,ecol-scol+1);
  for(int r=srow-1, i=0; r<erow; r++, i++)
    for(int c=scol-1, j=0; c<ecol; c++, j++)
      m.mat[i][j]=this->mat[r][c];

  return m;
}
//---------------------------------------------------------------------------
int operator== (const Matrix& left, const Matrix& right)
{
  if(right.cols != left.cols || right.rows != left.rows)
    return 0;

  for(int i=0; i<left.rows; i++)
    for(int j=0; j<left.cols; j++)
      if(left.mat[i][j] != right.mat[i][j]) return 0;

  return 1;
}
//---------------------------------------------------------------------------
int operator!= (const Matrix& left, const Matrix& right)
{
  return !(left == right);
}
//---------------------------------------------------------------------------
Matrix operator+ (const Matrix& m) {return m;}
//---------------------------------------------------------------------------
Matrix operator- (const Matrix& m) {return -1.0*m;}
//---------------------------------------------------------------------------
Matrix operator+ (const Matrix& left, const Matrix& right)
{
  Matrix m = left;
  return m += right;
}
//---------------------------------------------------------------------------
Matrix operator- (const Matrix& left, const Matrix& right)
{
  Matrix m = left;
  return m -= right;
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator+= (const Matrix& right)
{
  if(right.cols != cols || right.rows != rows)
    error("operator +=\nBroj redova/kolona nije isti!");

  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      mat[i][j] += right.mat[i][j];

  return *this;
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator-= (const Matrix& right)
{
  if(right.cols != cols || right.rows != rows)
    error("operator -=\nBroj redova/kolona nije isti!");

  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
    	mat[i][j] -= right.mat[i][j];

  return *this;
}
//---------------------------------------------------------------------------
Matrix operator* (const Matrix& left, const Matrix& right)
{
  if(left.cols != right.rows)
    error("operator *\nA.cols != B.rows!");

  Matrix m(left.rows, right.cols);
  for(int i=0; i<left.rows; i++)
    for(int j=0; j<right.cols; j++)
      for(int k=0; k<left.cols; k++)
        m.mat[i][j] += left.mat[i][k]*right.mat[k][j];

  return m;
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator*= (const double& t)
{
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      mat[i][j] *= t;

  return *this;
}
//---------------------------------------------------------------------------
Matrix operator* (const Matrix& left, const double& t)
{
  Matrix m = left;
  return m *= t;
}
//---------------------------------------------------------------------------
Matrix operator/ (const Matrix& left, const double& t)
{
  Matrix m = left;
  return m *= 1/t;
}
//---------------------------------------------------------------------------
Matrix operator* (const double& t, const Matrix& right)
{
  return right*t;
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator*= (const Matrix& right)
{
  return *this = *this*right;
}
//---------------------------------------------------------------------------
istream& operator>> (istream& is, Matrix& m)
{
  for(int i=0; i<m.rows; i++)
    for(int j=0; j<m.cols; j++)
    {
      if(is.eof())
        error("Matrix:: istream >> \n Ulazni podaci ne odgovaraju velicini matrice.");
      is >> m.mat[i][j];
    }

  is.seekg(0, ios::beg);
  return is;
}
//---------------------------------------------------------------------------
ostream& operator<< (ostream& os, const Matrix& m)
{
  os << "{\n";
  os.precision(4);
  os.setf(ios::fixed);
  for(int i=0; i<m.rows; i++)
  {
    os << '{';
    for(int j=0; j<m.cols; j++)
      os << m.mat[i][j] << ((j<m.cols-1) ? "    " : "}");

    os << ((i<m.rows-1)?"\n":"\n}\n");
  }
  return os;
}
//---------------------------------------------------------------------------
Matrix E(int n)
{
  Matrix m(n, n);
  for(int i=0; i<n; i++)
    m.mat[i][i] = 1.0;

  return m;
}
//---------------------------------------------------------------------------
double innerproduct(const Matrix& left, const Matrix& right)
{
  if(left.cols != right.rows)
    error("innerproduct()\nA.cols != B.rows!");
  if(left.rows != 1 || right.cols != 1)
    error("innerproduct()\n Both A and B should be vectors!");

  double d = 0.0;
  for(int k=0; k<left.cols; k++)
    d += left.mat[0][k]*right.mat[k][0];

  return d;
}
//---------------------------------------------------------------------------
double trag(const Matrix& m)
{
  if(m.cols != m.rows)
    error("F-ja trag()\nNije kvadratna matrica!");

  double tr=0.0;
  for(int j=0; j<m.rows; j++)
  {
    tr += m.mat[j][j];
  }
  return tr;
}
//---------------------------------------------------------------------------
Matrix trans(const Matrix& m)
{
  Matrix mtemp(m.cols, m.rows);
  for(int i=0; i<m.rows; i++)
    for(int j=0; j<m.cols; j++)
      mtemp.mat[j][i] = m.mat[i][j];

  return mtemp;
}

//---------------------------------------------------------------------------
//complement
Matrix comp(const Matrix& m, int r, int c)
{
  if(m.cols<1 || m.rows<1)
    error("F-ja comp()\n A.cols < 1 ili A.rows <1 !");

  Matrix mtemp(m.rows-1, m.cols-1);
  for(int i=0,x=0; i<m.rows; i++)
  {
    if(i == r-1) continue;
    for(int j=0,y=0; j<m.cols; j++)
    {
      if(j == c-1) continue;
      mtemp.mat[x][y] = m.mat[i][j];
      y++;
    }
    x++;
  }
  
  return mtemp;
}
//---------------------------------------------------------------------------
//determinant
double det(const Matrix& m)
{
  if(m.cols != m.rows)
    error("F-ja det()\nNije kvadratna matrica!");

  if(m.rows == 1) return m.mat[0][0];
  double temp = 0.0;//double(0);
  for(int j=0,c=1; j<m.cols; j++,c*=-1)
    temp += (c*m.mat[0][j]*det(comp(m, 1, j+1)));

  return temp;
}
//---------------------------------------------------------------------------
// cofactor
double cofact(const Matrix& m, int r, int c)
{
  if(m.cols != m.rows)
    error("F-ja cofact()\nNije kvadratna matrica!");

  if((r+c)%2 != 0)
    return (-det(comp(m, r, c)));
  else
    return (det(comp(m, r, c)));
}
//---------------------------------------------------------------------------
// adjunct matrix
Matrix adj(const Matrix& m)
{
  if(m.cols != m.rows)
    error("F-ja adj()\nNije kvadratna matrica!");

  Matrix mtemp(m.rows, m.cols);
  for(int i=0; i<m.rows; i++)
    for(int j=0; j<m.cols; j++)
      mtemp.mat[i][j] = cofact(m, i+1, j+1);

  return trans(mtemp);
}

//---------------------------------------------------------------------------
Matrix inv(const Matrix& m)
{
  return gaussj(m);
}
//---------------------------------------------------------------------------
Matrix invh(const Matrix& m)
{
  return cholesky(m);
}
//---------------------------------------------------------------------------
//Gauss-Jordan
Matrix gaussj(const Matrix& aa)
{
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv;

  if (aa.rows != aa.cols)
    error("F-ja gaussj()\nNije kvadratna matrica!");

  int n=aa.rows;

  int* indexc = new int[n];
  int* indexr = new int[n];
  int* ipiv = new int[n];
  Matrix inv(n, n);
  inv.copy(aa);

  for(j=0; j<n; ++j)
   ipiv[j]=0;

  for(i=0; i<n; ++i)
  {
    big=0.0;
    for(j=0; j<n; ++j)
      if(ipiv[j]!=1)
        for(k=0; k<n; ++k)
        {
          if(ipiv[k]==0)
          {
            if(fabs(inv.mat[j][k])>=big)
            {
              big=fabs(inv.mat[j][k]);
              irow=j;
              icol=k;
            }
          }
          else if(ipiv[k] > 1)
            error("F-ja gaussj()\nipiv[k] > 1"); //error1
        }

    ++(ipiv[icol]);

    if(irow != icol)
    {
      for(l=0; l<n; ++l)
       swap(inv.mat[irow][l], inv.mat[icol][l]);
    }

    indexr[i]=irow;
    indexc[i]=icol;
    if(inv.mat[icol][icol]==0)
      error("F-ja gaussj()\ninv.mat[icol][icol]==0"); //error2

    pivinv=1.0/inv.mat[icol][icol];
    inv.mat[icol][icol]=1.0;

    for(l=0; l<n; ++l)
      inv.mat[icol][l] *= pivinv;

    for(ll=0; ll<n; ++ll)
      if(ll != icol)
      {
        dum=inv.mat[ll][icol];
        inv.mat[ll][icol]=0.0;
        for(l=0; l<n; ++l)
          inv.mat[ll][l] -= inv.mat[icol][l]*dum;
      }
  }

  for(l=n-1; l>=0; --l)
  {
    if(indexr[l] != indexc[l])
      for(k=0; k<n; ++k)
        swap(inv.mat[k][indexr[l]], inv.mat[k][indexc[l]]);
  }

  delete[] ipiv;
  delete[] indexc;
  delete[] indexr;

  return inv;
}
//---------------------------------------------------------------------------
//Cholesky
//aa symetric and positive definite
Matrix lowerL(const Matrix& aa)
{
 if (aa.rows != aa.cols)
    error("F-ja lowerL()\nNije kvadratna matrica!");

 int i, j, k;
 int n=aa.rows;
 double sum;
 Matrix l(n, n); //left lower matrix of aa

 for(j=1; j<=n; j++)
   for(k=j; k<=n; k++)
   {
     for(i=1, sum=0.0; i<j; i++)   sum += l(j,i)*l(k,i);

     if(j==k)
       l(j,j) = sqrt(aa(j,j)-sum);
     else
       l(k,j) = (aa(j,k)-sum)/l(j,j);
   } 


 return l;
}

//---------------------------------------------------------------------------
// L*y=b
// L(n,n) left diagonal
Matrix Lyb(const Matrix& L, const Matrix& b)
{
 if (L.rows != L.cols)
    error("F-ja Lyb()\nL nije kvadratna matrica!");
 if (L.rows != b.rows)
    error("F-ja Lyb()\nL i b nemaju isti broj redova!");
 if (b.cols != 1)
    error("F-ja Lyb()\nb nije vektor!");

 int n=L.rows;
 double sum;
 Matrix y(n, 1); 

 for(int k=0; k<n; k++)
 {
   sum=0.0;
   for(int i=0; i<=k-1; i++)
     sum += L.mat[k][i]*y.mat[i][0];

   y.mat[k][0] = (b.mat[k][0]-sum)/L.mat[k][k];
 }

 return y;

}
//---------------------------------------------------------------------------
// L*x=y
// L(n,n) upper diagonal
Matrix Ltxy(const Matrix& L, const Matrix& y)
{
 if (L.rows != L.cols)
    error("F-ja Ltxy()\nL nije kvadratna matrica!");
 if (L.rows != y.rows)
    error("F-ja Ltxy()\nL i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 int n=L.rows;
 double sum;
 Matrix x(n, 1); 

 for(int k=n-1; k>=0; k--)
 {
   sum=0.0;
   for(int i=n-1; i>=k+1; i--)
     sum += L.mat[k][i]*x.mat[i][0];

   x.mat[k][0] = (y.mat[k][0]-sum)/L.mat[k][k];
 }

 return x;
}

//---------------------------------------------------------------------------
pair<Matrix, Matrix> ludecomp(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja ludecomp()\nA nije kvadratna matrica!");

 int k;
 int n=A.rows;
 double sum;
 Matrix L = E(n);  
 Matrix U(n, n);

 double eps=1e-12; 

 for(int it=1; it<=n; it++)
 {
   for(int j=it; j<=n; j++)
   {
     for(sum=0.0, k=1; k<it; k++)   sum += L(it,k)*U(k,j);
     U(it,j) = A(it,j) - sum;
     if(it==j && fabs(U(it,j)) < eps)
     {
       cout << "linearna zavisnost:  " << it << ". red/kolona" << endl;
       error("F-ja ludecomp()\nA nepotpun rang!");
     }
   }

   for(int i=it+1; i<=n; i++)
   {
     for(sum=0.0, k=1; k<it; k++)   sum += L(i,k)*U(k,it);
     L(i,it) = (A(i,it) - sum)/U(it,it);
   }
 }

 return make_pair(L, U);
}

//---------------------------------------------------------------------------
Matrix factor(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja factor()\nA nije kvadratna matrica!");

 int i, j, k;
 int n=A.rows;
 double sum, sum2;
 double* p = new double[n];
 Matrix L(n,n);

 for(i=1; i<=n; i++)
 {
   for(j=1, memset(p,0,(i-1)*sizeof(double)); j<i; j++)
   {
     for(k=1, sum=0.0; k<j; k++)   sum += p[k]*L(j,k);
     p[j] = A(i,j) - sum;
   }

   for(j=1; j<i; j++)
   {
     L(i,j) = (p[j]-A(i,j))/L(j,j);
   }

   for(j=1, sum2=0.0; j<i; j++)
   {
     L(i,j) += A(i,j)/L(j,j);
     sum2 += p[j] * L(i,j);
   }

   L(i,i) = A(i,i) - sum2;
 }

 delete[] p;
 return L;
}

//---------------------------------------------------------------------------
Matrix factor2(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja factor()\nA nije kvadratna matrica!");

 int i, j, k;
 int n=A.rows;
 double sum, sum2;
 double* p = new double[n]; 
 Matrix L(n,n);

 for(i=1; i<=n; i++)
 {
   for(j=1, sum2=0.0, memset(p,0,(i-1)*sizeof(double)); j<i; j++)
   {
     for(k=1, sum=0.0; k<j; k++)  sum += p[k]*L(j,k);

     p[j] = A(i,j) - sum;
     L(i,j) = p[j] / L(j,j;
     sum2 += p[j] * L(i,j);
   }

   L(i,i) = A(i,i) - sum2;
 }

 delete[] p;
 return L;
}

//---------------------------------------------------------------------------
Matrix solveUC(const Matrix& A, const Matrix& y)
{
 if (A.rows != A.cols)
    error("F-ja solve()\nA nije kvadratna matrica!");
 if (A.rows != y.rows)
    error("F-ja solve()\nA i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 int i, k;
 int n=A.rows;
 double sum;

 Matrix p = factor(A);

 Matrix y2(n, 1); 
 for(i=1; i<=n; i++)
 {
   for(k=1, sum=0.0; k<i; k++)   sum += p(i,k)*y2(k,1); 
   y2(i,1) = (y(i,1) - sum);
 }

 Matrix x(n, 1); 
 for(i=n; i>=1; i--)
 {
   for(k=n, sum=0.0; k>i; k--)   sum += p(k,i)*x(k,1); 
   x(i,1) = y2(i,1)/p(i,i) - sum;
 }

 return x;
}

//---------------------------------------------------------------------------
Matrix solveH(const Matrix& A, const Matrix& y)
{
 if (A.rows != A.cols)
    error("F-ja solve()\nA nije kvadratna matrica!");
 if (A.rows != y.rows)
    error("F-ja solve()\nA i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 Matrix L = lowerL(A);  
 Matrix y2 = Lyb(L,y); 

 return Ltxy(trans(L),y2); 
}

//---------------------------------------------------------------------------
Matrix solveLU(const Matrix& A, const Matrix& y)
{
 if (A.rows != A.cols)
    error("F-ja solve()\nA nije kvadratna matrica!");
 if (A.rows != y.rows)
    error("F-ja solve()\nA i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 pair<Matrix, Matrix> p = ludecomp(A);  
 Matrix y2 = Lyb(p.first,y); 

 return Ltxy(p.second,y2); 
}

//---------------------------------------------------------------------------
Matrix cholesky(const Matrix& aa)
{
  if (aa.rows != aa.cols)
    error("F-ja cholesky()\nNije kvadratna matrica!");

  int n=aa.rows;
  Matrix inv(n,n);

  Matrix L=lowerL(aa);
  Matrix C(n,n); 
  Matrix E1=E(n);

  Matrix b(n,1), y(n,1), x(n,1);
  for(int i=1; i<=n; i++)
  {
    b=E1.submatrix(1,n,i,i);
    y=Lyb(L, b);
    for(int j=0; j<n; j++)
      C.mat[j][i-1]=y.mat[j][0];
  }

  for(int i=1; i<=n; i++)
  {
    b=C.submatrix(1,n,i,i);
    x=Ltxy(trans(L), b);
    for(int j=0; j<n; j++)
      inv.mat[j][i-1]=x.mat[j][0];
  }

  return inv;
}

