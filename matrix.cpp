#include <stdlib.h>
#include <math.h>
//#include <string.h> //mozda zbog memset

#include "matrix.h"

//---------------------------------------------------------------------------
inline void swap(double& a, double& b)
{
  double temp=a;
  a=b;
  b=temp;
}
/*
//---------------------------------------------------------------------------
inline const int Matrix::getRows() const
{
  return rows;
}
//---------------------------------------------------------------------------
inline const int Matrix::getCols() const
{
  return cols;
}
*/
//---------------------------------------------------------------------------
void Matrix::allocate()
{
  mat = new double*[rows];
  for(int i=0; i<rows; i++)
  {
    mat[i] = new double[cols];
    for(int j=0; j<cols; j++)
      mat[i][j] = 0.0;//double(0);
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

  //if(right.cols != cols || right.rows != rows)
  //error();

  release();
  rows=right.rows;
  cols=right.cols;
  allocate();
  copy(right);
  return *this;
}
/*
//---------------------------------------------------------------------------
//vraca clanove matrice ali indeks pocinje od 1 (ne od 0)
//A(1,1) je u stvari A[0][0]
inline const double& Matrix::operator() (int i, int j) const
{
  if(i<1 || i>rows || j<1 || j>cols)
  {
    cout << i << "  " << j << endl;
    error("Matrix::operator()\nNeodgovarajuci argumenti!");
  }

  return mat[i-1][j-1];
}

//---------------------------------------------------------------------------
//vraca clanove matrice ali indeks pocinje od 1 (ne od 0)
//A(1,1) je u stvari A[0][0]
inline double& Matrix::operator() (int i, int j)
{
  if(i<1 || i>rows || j<1 || j>cols)
  {
    cout << i << "  " << j << endl;
    error("Matrix::operator()\nNeodgovarajuci argumenti!");
  }
  return mat[i-1][j-1];
}
*/
//---------------------------------------------------------------------------
//vraca submatricu
//srow - startrow
//erow - endrow
//...itd, itd...
Matrix Matrix::submatrix(int srow, int erow, int scol, int ecol) const
{

  if(srow>erow || scol>ecol || srow<1 || erow>rows ||scol<1 || ecol>cols)
    error("Matrix::submatrix()\nNeodgovarajuci argumenti!");

  Matrix m(erow-srow+1,ecol-scol+1);
  for(int r=srow-1, i=0; r<erow; r++, i++)
  {
    for(int c=scol-1, j=0; c<ecol; c++, j++)
    {
      m.mat[i][j]=this->mat[r][c];
    }
  }

  return m;
}
//---------------------------------------------------------------------------
int operator== (const Matrix& left, const Matrix& right)
{
  if(right.cols != left.cols || right.rows != left.rows)
    //error("operator ==\nBroj redova/kolona nije isti!");
    return 0;

  for(int i=0; i<left.rows; i++)
  {
    for(int j=0; j<left.cols; j++)
    {
      if(left.mat[i][j] != right.mat[i][j]) return 0;
    }
  }
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
  {
    for(int j=0; j<cols; j++)
    {
      mat[i][j] += right.mat[i][j];
    }
  }
  return *this;
}
//---------------------------------------------------------------------------
Matrix& Matrix::operator-= (const Matrix& right)
{
  if(right.cols != cols || right.rows != rows)
    error("operator -=\nBroj redova/kolona nije isti!");

  for(int i=0; i<rows; i++)
  {
    for(int j=0; j<cols; j++)
    {

			mat[i][j] -= right.mat[i][j];
    }
  }
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
  {
    for(int j=0; j<cols; j++)
    {
      mat[i][j] *= t;
    }
  }
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
//cita matricu iz inputstream-a
istream& operator>> (istream& is, Matrix& m)
{
  for(int i=0; i<m.rows; i++)
  {
    for(int j=0; j<m.cols; j++)
    {
//fajl manji od matrice
      if(is.eof())
        error("Matrix:: istream >> \n Ulazni podaci ne odgovaraju velicini matrice.");
      is >> m.mat[i][j];
    }
  }
/* UKINUTO!!
//matrica manja od fajla
  double d;
  is >> d;
  if(!is.eof())
    error("Matrix:: istream >> \n Ulazni podaci ne odgovaraju velicini matrice.");
*/
//go again to the start of file
  is.seekg(0, ios::beg);
  return is;
}
//---------------------------------------------------------------------------
//upisuje matricu u outputstream
ostream& operator<< (ostream& os, const Matrix& m)
{
  os << "{\n";
  os.precision(4);
  os.setf(ios::fixed);
  for(int i=0; i<m.rows; i++)
  {
    os << '{';
    for(int j=0; j<m.cols; j++)
    {
      os << m.mat[i][j] << ((j<m.cols-1) ? "    " : "}");
    }
    os << ((i<m.rows-1)?"\n":"\n}\n");
  }
  return os;
}
//---------------------------------------------------------------------------
//pravi jedinacnu matricu reda n
Matrix E(int n)
{
  Matrix m(n, n);
  for(int i=0; i<n; i++)
  {
    m.mat[i][i] = 1.0;
  }

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
//trag matrice
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
//transpozicija
Matrix trans(const Matrix& m)
{
  Matrix mtemp(m.cols, m.rows);
  for(int i=0; i<m.rows; i++)
  {
    for(int j=0; j<m.cols; j++)
    {
      mtemp.mat[j][i] = m.mat[i][j];
    }
  }
  return mtemp;
}
//---------------------------------------------------------------------------
//komplement
//izbacuje red - r, kolonu - c
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
//determinanta
double det(const Matrix& m)
{
  if(m.cols != m.rows)
    error("F-ja det()\nNije kvadratna matrica!");

  if(m.rows == 1) return m.mat[0][0];
  double temp = 0.0;//double(0);
  for(int j=0,c=1; j<m.cols; j++,c*=-1)
  {
    temp += (c*m.mat[0][j]*det(comp(m, 1, j+1)));
  }
  return temp;
}
//---------------------------------------------------------------------------
//kofaktor (r,c) clana matrice m
double cofact(const Matrix& m, int r, int c)
{
  if(m.cols != m.rows)
    error("F-ja cofact()\nNije kvadratna matrica!");

  if((r+c)%2 != 0)
  {
    return (-det(comp(m, r, c)));
  }
  else
  {
    return (det(comp(m, r, c)));
  }
}
//---------------------------------------------------------------------------
//adjungovana matrica
Matrix adj(const Matrix& m)
{
  if(m.cols != m.rows)
    error("F-ja adj()\nNije kvadratna matrica!");

  Matrix mtemp(m.rows, m.cols);
  for(int i=0; i<m.rows; i++)
  {
    for(int j=0; j<m.cols; j++)
    {
      mtemp.mat[i][j] = cofact(m, i+1, j+1);
    }
  }
  return trans(mtemp);
}
//---------------------------------------------------------------------------
//UKINUTO !!
//inverzna matrica preko adjungovane matrice
//Matrix inv(const Matrix& m)
//{
//	if(m.cols != m.rows) error();
//	double temp = det(m);
//	if(temp == 0.0) error();
//	return (1.0/temp)*adj(m);
//}
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
//Gauss-Jordan-ova Eliminacija
//(sa pivotizacijom)
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
          }//ifipiv[k]=0
          else if(ipiv[k] > 1)
            error("F-ja gaussj()\nipiv[k] > 1"); //error1
        }//k loop

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
      }//ll loop
  }//i loop

  for(l=n-1; l>=0; --l)
  {
    if(indexr[l] != indexc[l])
      for(k=0; k<n; ++k)
        swap(inv.mat[k][indexr[l]], inv.mat[k][indexc[l]]);
  }//l loop

  delete[] ipiv;
  delete[] indexc;
  delete[] indexr;

  return inv;
}
//---------------------------------------------------------------------------
//Cholesky-jeva dekompozicija
//aa mora biti simetricna i pozitivno definitna
Matrix lowerL(const Matrix& aa)
{
 if (aa.rows != aa.cols)
    error("F-ja lowerL()\nNije kvadratna matrica!");

 //double MINVAL = 0.00001;
 //double MINVALKV = MINVAL*MINVAL;

 int i, j, k;
 int n=aa.rows;
 double sum;
 Matrix l(n, n); //left lower matrix of aa

 for(j=1; j<=n; j++)
 {
   for(k=j; k<=n; k++)
   {
     for(i=1, sum=0.0; i<j; i++)   sum += l(j,i)*l(k,i);

     if(j==k)
       l(j,j) = sqrt(aa(j,j)-sum);
     //if(j==k)
       //l.mat[j][j] = aa.mat[j][j]>sum+MINVALKV ? sqrt(aa.mat[j][j]-sum): MINVAL;
     else
       l(k,j) = (aa(j,k)-sum)/l(j,j);
   } //k loop
 }//j loop


 return l;
}
//---------------------------------------------------------------------------
//resavanje linearne jednacine - L*y=b
//L(n,n) je leva dijagonalna matrica; y(n) i b(n) su vektori
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
 Matrix y(n, 1); //resenje jednacine

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
//Resavanje linearne jednacine - L*x=y
//L(n,n) je gornja dijagonalna matrica; x(n) i y(n) su vektori
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
 Matrix x(n, 1); //resenje jednacine

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
//procedura za racunanje LU dekompozicije kvadratne matrice
//(za sada bez pivotiranja)
pair<Matrix, Matrix> ludecomp(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja ludecomp()\nA nije kvadratna matrica!");

 int k;
 int n=A.rows;
 double sum;
 Matrix L = E(n);  //L(i,i)=1.0;
 Matrix U(n, n);

 double eps=0.0000000001; // kriterijum za nepotpun rang - eps=10^(-10)

 //idemo po dijagonali i racunamo red desno i kolonu ispod
 for(int it=1; it<=n; it++)
 {
   //izracunaj red "it"; sve elemente U(it,j), (j=it,...,n)
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

   //izracunaj kolonu "it"; sve elemente L(i,it), (i=it+1,...,n)
   for(int i=it+1; i<=n; i++)
   {
     for(sum=0.0, k=1; k<it; k++)   sum += L(i,k)*U(k,it);
     L(i,it) = (A(i,it) - sum)/U(it,it);
   }
 }

 return make_pair(L, U);
}

//---------------------------------------------------------------------------
//brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
//(za sada bez pivotiranja)
Matrix factor(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja factor()\nA nije kvadratna matrica!");

 int i, j, k;
 int n=A.rows;
 double sum, sum2;
 double* p = new double[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
 Matrix L(n,n);

 for(i=1; i<=n; i++)
 {
   //racunanje cesto koriscenih proizvoda - p
   for(j=1, memset(p,0,(i-1)*sizeof(double)); j<i; j++)
   {
     for(k=1, sum=0.0; k<j; k++)   /*p[j] -=*/ sum += p[k]*L(j,k);
     p[j] = A(i,j) - sum;
     //p[j] += A(i,j);
   }

   //uz pomoc izracunatog "p" racunamo poddijagonalne clanove (L(i,j),j<i) u i-tom redu...
   for(j=1; j<i; j++)
   {
     L(i,j) = (p[j]-A(i,j))/L(j,j);
   }

   for(j=1, sum2=0.0; j<i; j++)
   {
     L(i,j) += A(i,j)/L(j,j);
     sum2 += p[j] * L(i,j);
   }

   //... a zatim i dijagonalni clan L(i,i) = d(i)
   L(i,i) = A(i,i) - sum2;
 }

 delete[] p;
 return L;
}

//---------------------------------------------------------------------------
//brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
//(za sada bez pivotiranja)
Matrix factor2(const Matrix& A)
{
 if (A.rows != A.cols)
    error("F-ja factor()\nA nije kvadratna matrica!");

 int i, j, k;
 int n=A.rows;
 double sum, sum2;
 double* p = new double[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
 Matrix L(n,n);

 for(i=1; i<=n; i++)
 {
   //racunanje cesto koriscenih proizvoda - p
   for(j=1, sum2=0.0, memset(p,0,(i-1)*sizeof(double)); j<i; j++)
   {
     for(k=1, sum=0.0; k<j; k++)  sum += p[k]*L(j,k);

     p[j] = A(i,j) - sum;
     //cout << "p(" << i << "," << j << ")=" << p[j] << endl;

    //uz pomoc izracunatog "p" racunamo poddijagonalne clanove (L(i,j),j<i) u i-tom redu...
     L(i,j) = (p[j]/*-A(i,j)*/) / L(j,j);
     //cout << "L(" << i << "," << j << ")=" << L(i,j) << endl;
     sum2 += p[j] * L(i,j);
   }

   //... a zatim i dijagonalni clan L(i,i) = d(i)
   L(i,i) = A(i,i) - sum2;
   //cout << " L(" << i << "," << i << ")=" << L(i,i) << endl << endl;
 }

 delete[] p;
 return L;
}

//---------------------------------------------------------------------------
//procedura za resavanje sistema jednacina
//preko brze uopstene holeskijeve dekompozicije (UC)
//(A mora da bude simetricna matrica)
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

 //uopstena holeskijeva dekompozicija
 //u pair<Matrix, Matrix> p = factor(A);
 Matrix p = factor(A);

 Matrix y2(n, 1); //promenjeni vektor y - resenje sistema L*y2 = y
 for(i=1; i<=n; i++)
 {
   for(k=1, sum=0.0; k<i; k++)   sum += p(i,k)*y2(k,1); //Dk * Lik * Yk
   y2(i,1) = (y(i,1) - sum);
 }

 Matrix x(n, 1); //vektor x - resenje sistema D*trans(L)*x = y2 (x se upisuje u y2)
 for(i=n; i>=1; i--)
 {
   for(k=n, sum=0.0; k>i; k--)   sum += p(k,i)*x(k,1); //Lki * Xk
   x(i,1) = y2(i,1)/p(i,i) - sum;
   //(y2(i,1) /= p(i,i)) -= sum; //moze i ovo
   /*
   y2(i,1) /= p(i,i);
   y2(i,1) -= sum;*/
 }

 return x;
}

//---------------------------------------------------------------------------
//procedura za resavanje sistema jednacina
//preko stroge holeskijeve dekompozicije
//(A mora da bude simetricna i pozitivno definitna matrica)
Matrix solveH(const Matrix& A, const Matrix& y)
{
 if (A.rows != A.cols)
    error("F-ja solve()\nA nije kvadratna matrica!");
 if (A.rows != y.rows)
    error("F-ja solve()\nA i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 Matrix L = lowerL(A);  //stroga Holeskijeva dekompozicija
 Matrix y2 = Lyb(L,y); //promenjeni vektor y - resenje sistema L*y2 = y

 return Ltxy(trans(L),y2); //return - resenje sistema trans(L)*x = y2
}

//---------------------------------------------------------------------------
//procedura za resavanje sistema jednacina
//preko LU holeskijeve dekompozicije
//(A je PROIZVOLJNA matrica)
Matrix solveLU(const Matrix& A, const Matrix& y)
{
 if (A.rows != A.cols)
    error("F-ja solve()\nA nije kvadratna matrica!");
 if (A.rows != y.rows)
    error("F-ja solve()\nA i b nemaju isti broj redova!");
 if (y.cols != 1)
    error("F-ja solve()\ny nije vektor!");

 pair<Matrix, Matrix> p = ludecomp(A);  //LU dekompozicija
 Matrix y2 = Lyb(p.first,y); //promenjeni vektor y - resenje sistema L*y2 = y

 return Ltxy(p.second,y2); //return - resenje sistema trans(L)*x = y2
}

//---------------------------------------------------------------------------
//procedura za racunanje inverzije matrice
//preko holeskijeve dekompozicije
Matrix cholesky(const Matrix& aa)
{
  if (aa.rows != aa.cols)
    error("F-ja cholesky()\nNije kvadratna matrica!");

  int n=aa.rows;
  Matrix inv(n,n);

  Matrix L=lowerL(aa);
  Matrix C(n,n); //C<=>trans(L)*inv(aa)
  Matrix E1=E(n);

  Matrix b(n,1), y(n,1), x(n,1);
//resavanje sistema L*C=E
  for(int i=1; i<=n; i++)
  {
    b=E1.submatrix(1,n,i,i);
    y=Lyb(L, b);
    for(int j=0; j<n; j++)
      C.mat[j][i-1]=y.mat[j][0];
  }

//resavanje sistema trans(L)*inv(aa)=C
  for(int i=1; i<=n; i++)
  {
    b=C.submatrix(1,n,i,i);
    x=Ltxy(trans(L), b);
    for(int j=0; j<n; j++)
      inv.mat[j][i-1]=x.mat[j][0];
  }

  return inv;
}

