/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   smilenkovic@rgz.gov.rs                                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
#ifndef _spMatrix
#define _spMatrix

#include <vector>
#include <map>
#include <list>
#include <cmath>
#include <algorithm>
#include <iterator> 

#include "matrix.h"

using namespace std;

//////////////////////////   SparseMatrix   /////////////////////////////////
template<typename T = double>
class SparseMatrix
{
  typedef T value_type;
  typedef map<unsigned int, T, less<unsigned int> > row_type;
  typedef map<unsigned int, T, less<unsigned int> > col_type;

  typedef typename vector<row_type*>::iterator rowsp_iterator;
  typedef typename vector<row_type>::iterator rows_iterator;
  typedef typename vector<col_type>::iterator cols_iterator;
  typedef typename vector<row_type>::const_iterator rows_const_iterator;
  typedef typename vector<col_type>::const_iterator cols_const_iterator;
  typedef typename vector<col_type>::const_reverse_iterator cols_const_riterator;

  typedef typename row_type::iterator onerow_iterator;
  typedef typename col_type::iterator onecol_iterator;
  typedef typename row_type::const_iterator onerow_const_iterator;
  typedef typename col_type::const_iterator onecol_const_iterator;
  typedef typename col_type::const_reverse_iterator onecol_const_riterator;

//---------------------------------------------------------------------------
  friend SparseMatrix<T> operator* (const SparseMatrix<T>& left, const SparseMatrix<T>& right)
  {
    //if(left.cols() != right.rows())
      //error("operator *\nA.cols != B.rows!");

    unsigned int i, j;

    typename SparseMatrix<T>::rows_const_iterator itrleft;
    typename SparseMatrix<T>::onerow_const_iterator itr1, itr2;

    unsigned int r = left.rows();
    unsigned int c = right.cols();
    SparseMatrix<T> m(r,c);
    list<pair<unsigned int, T> > poss;
    typename list<pair<unsigned int, T> >::iterator it;

    for( itrleft=left.rows_.begin(), i=1;
         itrleft!=left.rows_.end();
         ++itrleft, ++i)
    {
      poss.clear();
      poss.push_back(make_pair(c, left.nullVal)); //trick

      for( itr1=itrleft->begin(), j=itr1->first;
           itr1!=itrleft->end();
           ++itr1, j=itr1->first )
      	{
        for(itr2=right.rows_[j].begin(), it=poss.begin();
            itr2!=right.rows_[j].end();
            ++itr2)
        {
          while(itr2->first > it->first)
            ++it;

          if(c == it->first)
          {
            poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
            while(++itr2 != right.rows_[j].end())
            {
              poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
            }

            break; 
          }

          else if(itr2->first < it->first)
          {
            it = poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
          }

          else
          {
            it->second += itr1->second * itr2->second;
          }
        }
      }

      for(it=poss.begin(); it!=poss.end(); ++it)
        m.push_fast(i, it->first+1, it->second);
    }

    return m;
  }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i, j;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...

  SparseMatrix<T> L(n,n);

  for(it=A.rows_.begin(), i=0, memset(d,0,n*sizeof(T));
      it!=A.rows_.end(); ++it, ++i)
   {
     if( it->empty() )  continue;  

     for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
         itr!=it->end() && itr->first <= i; ++itr)
       a[itr->first] = itr->second;

     for(j=it->begin()->first, memset(p,0,i*sizeof(T)); j<i; ++j)
     {
       for(itr=L.rows_[j].begin(); itr!=L.rows_[j].end() && itr->first < j; ++itr)
         p[j] -= itr->second * p[itr->first];

       p[j] += a[j];

       if(p[j]==A.nullVal)  continue;

       l[j] = p[j]/d[j];
       d[i] -= p[j] * l[j];
       L.push_fast(i+1, j+1, l[j]);
     }

     d[i] += a[i];
     L.push_fast(i+1, i+1, d[i]);
   }

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor2(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i, j;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   for(it=A.rows_.begin(), i=0, memset(d,0,n*sizeof(T));
      it!=A.rows_.end(); ++it, ++i)
   {
     if( it->empty() )  continue;  

     for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
         itr!=it->end() && itr->first <= i; ++itr)
       a[itr->first] = itr->second;

     for(j=it->begin()->first, memset(p,0,i*sizeof(T)); j<i; ++j)
     {
       p[j] += a[j];

       for(itc=++L.cols_[j].find(j);
           itc!=L.cols_[j].end() && itc->first<i; ++itc)
       {
         p[itc->first] -= p[j] * itc->second;
       }
     }

     for(j=0; j<i; ++j)
     {
       if(p[j]==A.nullVal)  continue;

       l[j] = p[j]/d[j];
       d[i] -= p[j] * l[j];
       L.push_fast(i+1, j+1, l[j]);
     }

     d[i] += a[i];
     L.push_fast(i+1, i+1, d[i]);
   }

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor3(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i, j, k;
   unsigned int j2; 
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...

   SparseMatrix<T> L(n,n);

   for(it=A.rows_.begin(), i=0, memset(d,0,n*sizeof(T));
      it!=A.rows_.end(); ++it, ++i)
   {
     if( !it->empty() )  
     {
       for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
           itr!=it->end() && itr->first <= i; ++itr)
         a[itr->first] = itr->second;

       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first;
       lt.clear();  lt.push_back(j);
       lit0=lt.begin();
       while(j<i)
       {
         if(j==j2)
         {
           p[j] += a[j];
           j2 = (++itr)->first; 
       }

         for(itc=++L.cols_[j].find(j), k=itc->first, lit = lit0;
             itc!=L.cols_[j].end() && k<i;
             ++itc, k=itc->first)
         {
           p[k] -= p[j] * itc->second;

           if(lit==lt.end())
             lt.push_back(k);
           else
           {
             while( ++lit!=lt.end() && k > *lit);
             if(lit==lt.end() || k<*lit)
               lit=lt.insert(lit, k);
           }
         }

         if( ++lit0 != lt.end() && *lit0<=j2) 
           j = *lit0;
         else 
         {
           j=j2;
           if(j<i)  
             lit0 = lt.insert(lit0, j);
         }
       }

       for(lit=lt.begin(), j=*lit; lit!=lt.end() && j<i; ++lit, j=*lit)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);
       }

       d[i] += a[i];
       L.push_fast(i+1, i+1, d[i]);
     }
   }

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor4(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2; 
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...

   SparseMatrix<T> L(n,n);

   it=A.rows_.begin();
   memset(d,0,n*sizeof(T));
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   
     {
       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first;
       lt.clear();
       lit0 = lt.insert(lt.end(), n); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           p[j] += itr->second;

           if(++itr != it->end())
             j2 = itr->first;
           else
             j2=i;

           if( j != *lit0 )
            lit0 = lt.insert(lit0, j);
         }

         lit = lit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           k = itc->first;
           p[k] -= p[j] * itc->second;

           while(k > *(++lit));

           if(n == *lit) 
           {
             lit=lt.insert(lit, k);
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               p[k] -= p[j] * itc->second;
               lit=lt.insert(++lit, k);
             }
             break; 
           }

           else if(k < *lit)
             lit=lt.insert(lit, k);

         }

         if( *(++lit0) <= j2 ) 
           j = *lit0;
         else 
           j = j2;
       }

       lit = lt.begin();  j = *lit;
       while(j<i)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);

         j=*(++lit);
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] p;
         delete[] d;
         delete[] l;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }
     ++it; ++i;
   }

   delete[] p;
   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor45(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2; 
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...

   SparseMatrix<T> L(n,n);

   it=A.rows_.begin();
   memset(d,0,n*sizeof(T));
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )  
     {
       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first; 
       lt.clear();
       lit0 = lt.insert(lt.end(), n); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           p[j] += itr->second;

           if(++itr != it->end())
             j2 = itr->first;
           else
             j2=i;

           if( j != *lit0 )
            lit0 = lt.insert(lit0, j);
         }

         lit = lit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           k = itc->first;
           p[k] -= p[j] * itc->second;

           while(k > *(++lit));

           if(n == *lit) 
           {
             lit=lt.insert(lit, k);
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               p[k] -= p[j] * itc->second;
               lit=lt.insert(++lit, k);
             }
             break; 
           }

           else if(k < *lit)
             lit=lt.insert(lit, k);

         }

         if( *(++lit0) <= j2 ) 
           j = *lit0;
         else 
           j = j2;
       }

       lit = lt.begin();  j = *lit;
       while(j<i)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);

         j=*(++lit);
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] p;
         delete[] d;
         delete[] l;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }
     ++it; ++i;
   }

   delete[] p;
   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor5(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<pair<unsigned int, T> > lt;   //f
   typename list<pair<unsigned int, T> >::iterator lit0;
   typename list<pair<unsigned int, T> >::iterator lit;

   unsigned int n=A.rows();
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...

   SparseMatrix<T> L(n,n);

   memset(d,0,n*sizeof(T));
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )  
     {
       itr=it->begin();
       j = j2 = itr->first;
       lt.clear();
       lit0 = lt.insert(lt.end(), make_pair(n, A.nullVal)); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           if( j != lit0->first )
             lit0 = lt.insert(lit0, make_pair(j, itr->second));
           else
             lit0->second += itr->second;  

           if(++itr != it->end())
             j2 = itr->first;
           else
             j2=i;
         }

         lit = lit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           k = itc->first;

           while(k > (++lit)->first);

           if(n == lit->first) 
           {
             lit=lt.insert(lit, make_pair(k, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               lit=lt.insert(++lit, make_pair(k, -lit0->second*itc->second));
             }
             break; 
           }

	   else if(k < lit->first)
             lit=lt.insert(lit, make_pair(k, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }

         if( (++lit0)->first <= j2 ) 
           j = lit0->first;
         else 
           j = j2;
       }

       lit = lt.begin();  j = lit->first;
       while(j<i)
       {
         l[j] = lit->second/d[j];
         d[i] -= lit->second * l[j];
         L.push_fast(i+1, j+1, l[j]);

         j=(++lit)->first;
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         delete[] l;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }
     ++it; ++i;
   }

   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  // fast decomposition
  friend SparseMatrix<T> factor6(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j;
   unsigned int j2; //for double traversing
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<pair<unsigned int, T> > lt;   //f
   typename list<pair<unsigned int, T> >::iterator lit0, lit;

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; 

   SparseMatrix<T> L(n,n);

   memset(d,0,n*sizeof(T));

   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   
     {
       itr=it->begin();
       j = j2 = itr->first; 
       lt.clear();
       lit0 = lt.insert(lt.end(), make_pair(n, A.nullVal)); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           if( j != lit0->first )
             lit0 = lt.insert(lit0, make_pair(j, itr->second));
           else
             lit0->second += itr->second;  

           if(++itr != it->end())
             j2 = itr->first;
           else
             j2=i;
         }

         lit = lit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           while(itc->first > (++lit)->first);

           if(n == lit->first) 
           {
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               lit=lt.insert(++lit, make_pair(itc->first, -lit0->second*itc->second));

             break; 
           }

	   else if(itc->first < lit->first)
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }

         lj = lit0->second/d[j];
         L.push_fast(i+1, j+1, lj);
         d[i] -= lit0->second * lj;

         if( (++lit0)->first <= j2 ) 
           j = lit0->first;
         else 
           j = j2;
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }

     ++it; ++i;
   }

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor7(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j;
   unsigned int j2; 
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   map<unsigned int, T> mp;  
   typename map<unsigned int, T>::iterator mit0, mit;

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; 

   SparseMatrix<T> L(n,n);

   memset(d,0,n*sizeof(T));
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   
     {
       itr=it->begin();
       j = j2 = itr->first; 
       mp.clear();
       mit0 = mp.insert(make_pair(n, A.nullVal)).first; 
       while(j<i)
       {
         if(j==j2)
         {
           if( j != mit0->first )
             mit0 = mp.insert(mit0, make_pair(j, itr->second));
           else
             mit0->second += itr->second; 

           if(++itr != it->end())
             j2 = itr->first;
           else
             j2=i;
         }

         mit = mit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           while(itc->first > (++mit)->first);

           if(n == mit->first) 
           {
             mp.insert(make_pair(itc->first, -mit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               mp.insert(make_pair(itc->first, -mit0->second*itc->second));

             break; 
           }

	   else if(itc->first < mit->first)
             mit = mp.insert(mit, make_pair(itc->first, -mit0->second*itc->second));
           else
             mit->second -= mit0->second * itc->second;
         }

         lj = mit0->second/d[j];
         L.push_fast(i+1, j+1, lj);
         d[i] -= mit0->second * lj;

         if( (++mit0)->first <= j2 ) 
           j = mit0->first;
         else 
           j = j2;
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }

     ++it; ++i;
   }

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor8(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   map<unsigned int, T> mp;   //f
   typename map<unsigned int, T>::iterator mit; 

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; 

   SparseMatrix<T> L(n,n);

   memset(d,0,n*sizeof(T));
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )  
     {
       mp.clear();
       itr=it->begin();
       while(itr!=it->end() && itr->first<i)
       {
         mp.insert(make_pair(itr->first, itr->second));
         ++itr;
       }

       mit=mp.begin();
       while(mit!=mp.end())
       {
         itc = L.cols_[mit->first].begin();
         while(++itc != L.cols_[mit->first].end())
         {  
           mp[itc->first] -= mit->second * itc->second;
         }

         lj = mit->second/d[mit->first];
         L.push_fast(i+1, mit->first+1, lj);
         d[i] -= mit->second * lj;

         ++mit;
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }

     ++it; ++i;
   }

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> factor9(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<pair<unsigned int, T> > lt;  
   typename list<pair<unsigned int, T> >::iterator lit0, lit;

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; 

   SparseMatrix<T> L(n,n);

   memset(d,0,n*sizeof(T));
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   
     {
       lt.clear();
       itr=it->begin();
       while(itr!=it->end() && itr->first<i)
       {
         lt.push_back(make_pair(itr->first, itr->second)); 
         ++itr;
       }

       lt.push_back(make_pair(n, T()));
       lit0 = lt.begin();
       while( lit0->first < i )
       {
         int j = lit0->first;   
         T pj = lit0->second;    

         lit = lit0;
         itc = L.cols_[lit0->first].begin();
         while(++itc != L.cols_[lit0->first].end())
         {
           while(itc->first > (++lit)->first);

           if(n == lit->first) 
           {
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               lit=lt.insert(++lit, make_pair(itc->first, -lit0->second*itc->second));

             break; 
           }

	   else if(itc->first < lit->first)
             lit = lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }

         lj = lit0->second/d[lit0->first];
         L.push_fast(i+1, lit0->first+1, lj);
         d[i] -= lit0->second * lj;

         ++lit0;
       }

       if( itr->first == i )
         d[i] += itr->second;

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }

     ++it; ++i;
   }

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> ltl(SparseMatrix<T>& A)
  {
   //if (A.rows() != A.rows())
      //error("F-ja lowerL()\nNije kvadratna matrica!");

   unsigned int i, j, k;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr1;
   typename SparseMatrix<T>::onecol_const_iterator itr2;

   unsigned int n=A.rows();
   T* sum = new T[n]; 
   T* a = new T[n]; 
   T dd;
   SparseMatrix<T> L(n, n); 

   for(it=A.rows_.begin(), j=1, memset(sum,0,n*sizeof(T));
       it!=A.rows_.end(); ++it, ++j, memset(sum, 0, n*sizeof(T)))
   {
     for(itr1=it->begin(), memset(a,0,n*sizeof(T)); 
         itr1!=it->end(); ++itr1)
       a[itr1->first] = itr1->second;

     for(itr1=L.rows_[j-1].begin(), i=itr1->first;
         itr1!=L.rows_[j-1].end(); ++itr1, i=itr1->first)
     {
       for(itr2=L.cols_[i].find(j-1); itr2!=L.cols_[i].end(); ++itr2)
         sum[itr2->first] += itr1->second * itr2->second;
     }

     for(k=j; k<=n; k++)
     {
       if(j==k)
       {
         dd = sqrt(a[j-1]-sum[k-1]);
         L.push_fast( j,j, dd );
       }
       else
         L.push_fast( k,j, (a[k-1]-sum[k-1])/dd );
     } 
   }

 return L;
 }

//---------------------------------------------------------------------------
  friend ostream& operator<< (ostream& os, const SparseMatrix<T>& sm)
  {
    if(sm.rows_.empty())
      return os << "{}" << endl;

    os << "{" << endl;
    os.precision(4);
    os.setf(ios::fixed);
    for(unsigned int i=1; i<=sm.rows(); ++i)
    {
      os << '{';
      for(unsigned int j=1; j<=sm.cols(); ++j)
      {
        os << sm(i,j);
        os << ( j<sm.cols() ? "    " : "}" ) ;
      }
      os << endl;
      if( i==sm.rows() )
        os << "}" << endl;
    }

    return os;
  };

//---------------------------------------------------------------------------
private:
  T nullVal;
  vector<row_type> rows_;
  vector<col_type> cols_;

public:
//---------------------------------------------------------------------------
  //destructor
  ~SparseMatrix() { clear(); }

//---------------------------------------------------------------------------
  //default constructor
  SparseMatrix() : nullVal()
  {
    clear();
  };

//---------------------------------------------------------------------------
  //copy constructor
  SparseMatrix(const SparseMatrix<T>& sm2) : nullVal()
  {
    copy(sm2);
  };

//---------------------------------------------------------------------------
  //creates SparseMatrix of mxn dimension
  SparseMatrix(unsigned int m, unsigned int n) : nullVal()
  {
    //fill row-map
    rows_.clear();
    for(unsigned int i=0; i<m; ++i)
    {
      map<unsigned int, T> temprow;
      temprow.clear();
      rows_.push_back(temprow);
    }

    cols_.clear();
    for(unsigned int i=0; i<n; ++i)
    {
      map<unsigned int, T> tempcol;
      tempcol.clear();
      cols_.push_back(tempcol);
    }
  };

//---------------------------------------------------------------------------
  void clear()
  {
    for(rows_iterator itr=rows_.begin(); itr!=rows_.end(); ++itr)
      itr->clear();

    for(cols_iterator itc=cols_.begin(); itc!=cols_.end(); ++itc)
      itc->clear();

    rows_.clear();
    cols_.clear();

  }

//---------------------------------------------------------------------------
  unsigned int rows() const  { return rows_.size();};

//---------------------------------------------------------------------------
  unsigned int cols() const  { return cols_.size();};

//---------------------------------------------------------------------------
  bool copy(const SparseMatrix<T>& sm2)
  {
    clear();

    for(rows_const_iterator itr=sm2.rows_.begin(); itr!=sm2.rows_.end(); ++itr)
      rows_.push_back(*itr);

    for(cols_const_iterator itc=sm2.cols_.begin(); itc!=sm2.cols_.end(); ++itc)
      cols_.push_back(*itc);

    return true;
};

//---------------------------------------------------------------------------
  SparseMatrix<T>& operator= (const SparseMatrix<T>& sm2)
  {
    if( this != &sm2 )
      copy(sm2);

    return *this;
  };

//---------------------------------------------------------------------------
  friend SparseMatrix<T> trans(const SparseMatrix<T>& sm2)
  {
    unsigned int r = sm2.rows();
    unsigned int c = sm2.cols();
    SparseMatrix<T> t(r,c);
    t.cols_ = sm2.rows_;
    t.rows_ = sm2.cols_;

    return t;
  };

//---------------------------------------------------------------------------
  pair<int, double> findN()
  {
   int n=0;
   double p=0.0;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;

   for(it=rows_.begin(); it!=rows_.end(); ++it)
     n += it->size();

   p = static_cast<double>(n)/rows_.size()/cols_.size();

   return make_pair(n, p);
  }

//---------------------------------------------------------------------------
  friend Matrix solveUC(const SparseMatrix<T>& A, const Matrix& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();

    if (n != c)
      error("F-ja solveUC()\nA nije kvadratna matrica!");
    if (n != y.getRows())
      error("F-ja solveUC()\nA i b nemaju isti broj redova!");
    if (y.getCols() != 1)
      error("F-ja solveUC()\nb nije vektor!");
    //*/

    Matrix y2(n, 1); 
    SparseMatrix<T> p = factor6(A);
    for(i=1; i<=n; i++)
      for(k=1, y2(i,1) = y(i,1); k<i; k++)
        y2(i,1) -= p(i,k)*y2(k,1); 

    Matrix x(n, 1); 
    for(i=n; i>=1; i--)
      for(k=n, x(i,1) = y2(i,1)/p(i,i); k>i; k--)
        x(i,1) -= p(k,i)*x(k,1); 

    return x;
  };

//---------------------------------------------------------------------------
  friend vector<double> solveUC2(const SparseMatrix<T>& A, const vector<double>& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();

    if (n != c)
      error("F-ja solveUC2()\nA nije kvadratna matrica!");
    if (n != y.size())
      error("F-ja solveUC()\nA i y nemaju isti broj redova!");

    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;
    typename SparseMatrix<T>::cols_const_iterator it2, itn, it0;
    typename SparseMatrix<T>::onecol_const_iterator itc, itcend, itcbeg;

    vector<double> y2(n); 
    SparseMatrix<T> p = factor6(A);
    for(it=p.rows_.begin(), i=0; it!=p.rows_.end(); ++it, ++i)
      for(itr=it->begin(), k=itr->first, y2[i] = y[i]; itr!=--it->end(); ++itr, k=itr->first)
        y2[i] -= itr->second * y2[k];

    vector<double> x(n); 
    itn = --p.cols_.end();
    it0 = --p.cols_.begin();

    for(it2=itn, i=n; it2!=it0; --it2, --i)
    {
      itcend = --it2->end();
      itcbeg = it2->begin();
      for(itc=itcend, k=itc->first, x[i-1] = y2[i-1]/p(i,i); itc!=itcbeg; --itc, k=itc->first)
        x[i-1] -= itc->second * x[k];
    }

    return x;
  };

//---------------------------------------------------------------------------
  friend vector<double> solveUC3(const SparseMatrix<T>& A, const vector<double>& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();
    
    if (n != c)
      error("F-ja solveUC2()\nA nije kvadratna matrica!");
    if (n != y.size())
      error("F-ja solveUC()\nA i y nemaju isti broj redova!");

    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;
    typename SparseMatrix<T>::cols_const_riterator it2;
    typename SparseMatrix<T>::onecol_const_riterator itc;

    vector<double> y2(n); 
    SparseMatrix<T> p = factor6(A);
    for(it=p.rows_.begin(), i=0; it!=p.rows_.end(); ++it, ++i)
      for(itr=it->begin(), k=itr->first, y2[i] = y[i]; itr!=--it->end(); ++itr, k=itr->first)
        y2[i] -= itr->second * y2[k];

    vector<double> x(n); 
    for(it2=p.cols_.rbegin(), i=n-1; it2!=p.cols_.rend(); ++it2, --i)
	  for(itc=it2->rbegin(), k=itc->first, x[i] = y2[i]/p(i+1,i+1); itc!=--it2->rend(); ++itc, k=itc->first)
        x[i] -= itc->second * x[k];

    return x;
  };

//---------------------------------------------------------------------------
  friend SparseMatrix<double>& operator*=(SparseMatrix<double>& A, double x)
  {
    typename SparseMatrix<T>::rows_iterator it;
    typename SparseMatrix<T>::onerow_iterator itr;

    for(it=A.rows_.begin(); it!=A.rows_.end(); ++it)
      for(itr=it->begin(); itr!=it->end(); ++itr)
        itr->second *= x;

    return A;
  };

//---------------------------------------------------------------------------
  friend vector<double> operator*(const SparseMatrix<double>& A, const vector<double>& x)
  {
    if(A.cols() != x.size())
      error("operator *\nA.cols != x.size()!");

    unsigned int i, k;
    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;
    vector<double>::const_iterator itx, itm;

    vector<double> m(A.rows());
    for(it=A.rows_.begin(), i=0; it!=A.rows_.end(); ++it, ++i)
      for(itr=it->begin(), k=itr->first; itr!=it->end(); ++itr, k=itr->first)
        m[i] += itr->second * x[k];

    return m;
  };

//---------------------------------------------------------------------------
  //opt
  friend Matrix operator*(const SparseMatrix<T>& A, const Matrix& x)
  {
    if(A.cols() != x.getRows())
      error("operator *\nA.cols != x.rows!");

    unsigned int i=1, k=1;
    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;

    Matrix m(A.rows(), x.getCols());
    for(it=A.rows_.begin(); it!=A.rows_.end(); ++it, ++i)
      for(unsigned int j=1; j<=x.getCols(); j++)
        for(itr=it->begin(), k=1+itr->first; itr!=it->end(); ++itr, k=1+itr->first)
          m(i,j) += itr->second * x(k,j);

    return m;
  };

//---------------------------------------------------------------------------
  T& operator() (unsigned int i, unsigned int j)
  {
    //if( i>rows() || j>cols() )
      //return nullVal;

    onerow_iterator it = rows_[i-1].find(j-1);
    return it != rows_[i-1].end() ? it->second : nullVal;
  };

//---------------------------------------------------------------------------
  const T operator() (unsigned int i, unsigned int j) const //would return ref to temp
  {
    if( i>rows() || j>cols() )
      return nullVal;

    onerow_const_iterator it = rows_[i-1].find(j-1);
    return it != rows_[i-1].end() ? it->second : nullVal;
  };

//---------------------------------------------------------------------------
  //push provides only write-access
  void push(unsigned int i, unsigned int j, const T& val)
  {
    if( i<=rows() && j<=cols() )
    {
      onerow_iterator it = rows_[i-1].find(j-1);
      if( it != rows_[i-1].end() )
      {
        if( val != nullVal )
          cols_[j-1][i-1] = it->second = val;
        else
        {
          rows_[i-1].erase(j-1);
          cols_[j-1].erase(i-1);
        }
      }
      else if( val != nullVal )
        cols_[j-1].insert( make_pair(i-1,
                           rows_[i-1].insert(make_pair(j-1,val)).first->second) );
    }
  };

//---------------------------------------------------------------------------
  //it is to be used only when there isn't (i,j) element (100% sure)
  //and i, j are within rows(), cols()
  void push_fast(unsigned int i, unsigned int j, const T& val)
  {
    if(val != nullVal)
      cols_[j-1].insert( make_pair(i-1,
                         rows_[i-1].insert(make_pair(j-1,val)).first->second) );
  };

//---------------------------------------------------------------------------


};
//////////////////////////   SparseMatrix   /////////////////////////////////

#endif
