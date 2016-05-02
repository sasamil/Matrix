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
//#include <iostream> //u //f samo za testiranja mozda kasnije za copy zajebancije
#include <iterator> //u //f samo za testiranja mozda kasnije za copy zajebancije

#include "matrix.h"

using namespace std;

//////////////////////////   SparseMatrix   /////////////////////////////////
template<typename T = double>
class SparseMatrix
{
  typedef T value_type;
  typedef map<unsigned int, T, less<unsigned int> > row_type;
  typedef map<unsigned int, T, less<unsigned int> > col_type;

  typedef typename vector<row_type*>::iterator rowsp_iterator; //f 07.01.2007.
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
  //mnozenje matrica. std::map se koristi kao posrednik kod umetanja
  //(nije najbrza ali je odmah posle najbrze varijante)
  friend SparseMatrix<T> operator- (const SparseMatrix<T>& left, const SparseMatrix<T>& right)
  {
    //if(left.cols() != right.rows())
      //error("operator *\nA.cols != B.rows!");

    unsigned int i, j;

    typename SparseMatrix<T>::rows_const_iterator itrleft;
    typename SparseMatrix<T>::onerow_const_iterator itr1, itr2;

    unsigned int r = left.rows();
    unsigned int c = right.cols();
    SparseMatrix<T> m(r,c);
    map<unsigned int, T> skeys;
    typename map<unsigned int, T>::iterator its;

    for( itrleft=left.rows_.begin(), i=1;
         itrleft!=left.rows_.end();
         ++itrleft, ++i, skeys.clear() )
    {
      //if( itrleft->empty() )  continue;
      for( itr1=itrleft->begin(), j=itr1->first;
           itr1!=itrleft->end();
           ++itr1, j=itr1->first )
      {
        //if( right.rows_[j].empty() )  continue;  //u
        for(itr2=right.rows_[j].begin(); itr2!=right.rows_[j].end(); ++itr2)
        {
          its = skeys.find(itr2->first);
          if(its!=skeys.end())
            its->second += itr1->second * itr2->second;
          else
            skeys.insert(make_pair(itr2->first, itr1->second * itr2->second));
        }
      }//for every left-row member - iterate corresponding right row

      //fill up the product-row
      //for(unsigned int k=1; k<=c; ++k)
        //m.push_fast(i, k, vals[k-1]);
      for(its=skeys.begin(); its!=skeys.end(); ++its)
        m.push_fast(i, its->first+1, its->second);
    }//for every non-empty row in left matrix

    return m;
  }

//---------------------------------------------------------------------------
  //mnozenje matrica. std::list (koji se kasnije sortira) se koristi kao posrednik kod umetanja
  friend SparseMatrix<T> operator+ (const SparseMatrix<T>& left, const SparseMatrix<T>& right)
  {
    //if(left.cols() != right.rows())
      //error("operator *\nA.cols != B.rows!");

    unsigned int i, j;

    typename SparseMatrix<T>::rows_const_iterator itrleft;
    typename SparseMatrix<T>::onerow_const_iterator itr1, itr2;

    unsigned int r = left.rows();
    unsigned int c = right.cols();
    SparseMatrix<T> m(r,c);
    T* vals = new T[c];
    list<unsigned int> poss;
    typename list<unsigned int>::const_iterator its;

    for( itrleft=left.rows_.begin(), i=1, memset(vals, 0, c*sizeof(T));
         itrleft!=left.rows_.end();
         ++itrleft, ++i, memset(vals, 0, c*sizeof(T)), poss.clear() )
    {
      //if( itrleft->empty() )  continue;
      for( itr1=itrleft->begin(), j=itr1->first;
           itr1!=itrleft->end();
           ++itr1, j=itr1->first )
      {
        //if( right.rows_[j].empty() )  continue;  //u
        for(itr2=right.rows_[j].begin(); itr2!=right.rows_[j].end(); ++itr2)
        {
          vals[itr2->first] += itr1->second * itr2->second;
          poss.push_back(itr2->first);
        }
      }//for every left-row member - iterate corresponding right row

      //fill up the product-row
      //for(unsigned int k=1; k<=c; ++k)
        //m.push_fast(i, k, vals[k-1]);
      poss.sort();
      poss.unique();
      for(its=poss.begin(); its!=poss.end(); ++its)
        m.push_fast(i, *its+1, vals[*its]);
    }//for every non-empty row in left matrix

    delete[] vals;
    return m;
  }

//---------------------------------------------------------------------------
  //mnozenje matrica. std::list (sortiranje na licu mesta) se koristi kao posrednik kod umetanja
  //ovo je najbrza metoda
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

    //double fpp = 0;

    for( itrleft=left.rows_.begin(), i=1;
         itrleft!=left.rows_.end();
         ++itrleft, ++i)
    {
      poss.clear();
      poss.push_back(make_pair(c, left.nullVal)); //trick

      //if( itrleft->empty() )  continue;
      for( itr1=itrleft->begin(), j=itr1->first;
           itr1!=itrleft->end();
           ++itr1, j=itr1->first )
      {
        //if( right.rows_[j].empty() )  continue;  //u
        for(itr2=right.rows_[j].begin(), it=poss.begin();
            itr2!=right.rows_[j].end();
            ++itr2)
        {
          while(itr2->first > it->first)
            ++it;

          if(c == it->first)
          {
            poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
            //++fpp;
            while(++itr2 != right.rows_[j].end())
            {
              poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
              //++fpp;
            }

            break; //go out of for(itr2) loop
          }

          //insert if first time on this position
					//sdf  sdf   jfjjhd
          else if(itr2->first < it->first)
          {
            it = poss.insert(it, make_pair(itr2->first, itr1->second*itr2->second) );
            //++fpp;
          }

          //change if we've already been here (itr2->first == it->first)
          else
          {
            it->second += itr1->second * itr2->second;
            //++fpp;
          }
        }
      }//for every left-row member - iterate corresponding right row

      //fill up the product-row
      for(it=poss.begin(); it!=poss.end(); ++it)
        m.push_fast(i, it->first+1, it->second);
    }//for every non-empty row in left matrix

    //std::cout << "fpp = " << fpp << "   #####################" << std::endl;

    return m;
  }

//---------------------------------------------------------------------------
  friend SparseMatrix<T> operator/ (const SparseMatrix<T>& left, const SparseMatrix<T>& right)
  {
    //if(left.cols() != right.rows())
      //error("operator *\nA.cols != B.rows!");

    unsigned int i, j;

    typename SparseMatrix<T>::rows_const_iterator itrleft;
    typename SparseMatrix<T>::onerow_const_iterator itr1, itr2;

    unsigned int r = left.rows();
    unsigned int c = right.cols();
    SparseMatrix<T> m(r,c);
    T* vals = new T[c];

    for( itrleft=left.rows_.begin(), i=1, memset(vals, 0, c*sizeof(T));
         itrleft!=left.rows_.end();
         ++itrleft, ++i, memset(vals, 0, c*sizeof(T)) )
    {
      //if( itrleft->empty() )  continue;
      for( itr1=itrleft->begin(), j=itr1->first;
           itr1!=itrleft->end();
           ++itr1, j=itr1->first )
      {
        //if( right.rows_[j].empty() )  continue;  //u
        for(itr2=right.rows_[j].begin(); itr2!=right.rows_[j].end(); ++itr2)
          vals[itr2->first] += itr1->second * itr2->second;
      }//for every left-row member - iterate corresponding right row

      //fill up the product-row
      for(unsigned int k=1; k<=c; ++k)
        m.push_fast(i, k, vals[k-1]);
    }//for every non-empty row in left matrix

    delete[] vals;
    return m;
  }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //PRVOBITNA VARIJANTA - KAO U KLASI Matrix
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
   //Matrix L(n,n);

  for(it=A.rows_.begin(), i=0, memset(d,0,n*sizeof(T));
      it!=A.rows_.end(); ++it, ++i)
   {
     if( it->empty() )  continue;   //f neki text

     //izvlacenje i-tog reda matrice A
     for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
         itr!=it->end() && itr->first <= i; ++itr)
       a[itr->first] = itr->second;

     //racunanje cesto koriscenih proizvoda - p
     for(j=it->begin()->first, memset(p,0,i*sizeof(T)); j<i; ++j)
     {
       for(itr=L.rows_[j].begin(); itr!=L.rows_[j].end() && itr->first < j; ++itr)
         p[j] -= itr->second * p[itr->first];

       p[j] += a[j];
       //u cout << "p(" << i+1 << "," << j+1 << ")=" << p[j] << endl;

       if(p[j]==A.nullVal)  continue;

       //uz pomoc izracunatog "p" racunamo poddijagonalne clanove (L(i,j),j<i) u i-tom redu...
       l[j] = p[j]/d[j];
       d[i] -= p[j] * l[j];
       L.push_fast(i+1, j+1, l[j]);
       //cout << "L(" << i+1 << "," << j+1 << ")=" << l[j] << endl; //u
     }

     //... a zatim i dijagonalni clan L(i,i) = d(i)
     d[i] += a[i];//A(i,i);
     L.push_fast(i+1, i+1, d[i]);
     //cout << " L(" << i+1 << "," << i+1 << ")=" << L(i,i) << endl << endl; //u
   }

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //SLICNO KAO U KLASI Matrix; REDOSLED RACUNANJA/ITERIRANJA POBOLJSAN
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
     if( it->empty() )  continue;   //f neki text

     //izvlacenje i-tog reda matrice A
     for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
         itr!=it->end() && itr->first <= i; ++itr)
       a[itr->first] = itr->second;

     //racunanje cesto koriscenih proizvoda - p
     //j = it->begin()->first;
     //memset(p,0,i*sizeof(T));
     //for(itr=it->begin(), j=itr->first, memset(p,0,i*sizeof(T));
       //  j<i; j=(++itr)->first)
     for(j=it->begin()->first, memset(p,0,i*sizeof(T)); j<i; ++j)
     {
       p[j] += a[j];

       //for(k=j+1; k<i; ++k)
       for(itc=++L.cols_[j].find(j);
           itc!=L.cols_[j].end() && itc->first<i; ++itc)
       {
         p[itc->first] -= p[j] * itc->second;
         //p[k] -= p[j] * L(k+1,j+1);
       }
       /*
       if(p[j]==A.nullVal)  continue;

       l[j] = p[j]/d[j];
       d[i] -= p[j] * l[j];
       L.push_fast(i+1, j+1, l[j]);
       //cout << "L(" << i+1 << "," << j+1 << ")=" << l[j] << endl; //u
       */
     }

     for(j=0; j<i; ++j)
     {
       if(p[j]==A.nullVal)  continue;

       l[j] = p[j]/d[j];
       d[i] -= p[j] * l[j];
       L.push_fast(i+1, j+1, l[j]);
     }

     //... a zatim i dijagonalni clan L(i,i) = d(i)
     d[i] += a[i];//A(i,i);
     L.push_fast(i+1, i+1, d[i]);
     //cout << " L(" << i+1 << "," << i+1 << ")=" << L(i,i) << endl << endl; //u
   }

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //UBACENA JE LISTA; SVAKI RED MATRICE L (ZAJEDNO SA FILLIN-ovima)
  //SE PUNI IZ LISTE
  friend SparseMatrix<T> factor3(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i, j, k;
   unsigned int j2; //for double traversing
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   //f
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...
   //list<unsigned int> lt;

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //iteriranje iz reda u red
   for(it=A.rows_.begin(), i=0, memset(d,0,n*sizeof(T));
      it!=A.rows_.end(); ++it, ++i)
   {
     if( !it->empty() )   //f neki text
     {
       //izvlacenje i-tog reda matrice A
       for(itr=it->begin(), memset(a,0,(i+1)*sizeof(T));
           itr!=it->end() && itr->first <= i; ++itr)
         a[itr->first] = itr->second;

       //racunanje cesto koriscenih proizvoda - p
       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       lt.clear();  lt.push_back(j);
       lit0=lt.begin();
       //lit0 = /*lt.end(); //*/lt.insert(lt.end(), n); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           p[j] += a[j];
           j2 = (++itr)->first; //f nema provere za itr->end() podrazumeva se da je dijagonala puna
           //lit0 = lt.insert(lit0, j);
         }

         //f find(j) ide u begin() u dobro dezajniranoj simetricnoj matrici

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
               //mora lit=lt.insert jer bi sledeci ++lit zapravo NEOPRAVDANO preskocio ispitivanje trenutnog lit
           }
         }//for all subdiagonal items in column j

         //find out the next item
         if( ++lit0 != lt.end() && *lit0<=j2) //next item in list for next iteration
           j = *lit0;
         else //A(i,j2) is overjumped OR we have just reached the diagonal
         {
           j=j2;
           if(j<i)  //definately, A(i,j2) is overjumped! So, it becomes the next list item
             lit0 = lt.insert(lit0, j);
         }
       }//while j<i

       //cout << "(" << i << "," << j << ")   ";
       for(lit=lt.begin(), j=*lit; lit!=lt.end() && j<i; ++lit, j=*lit)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);
         //cout << j << " ";
       }
       //cout << endl;

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       d[i] += a[i];//A(i,i);
       L.push_fast(i+1, i+1, d[i]);
     }//if not empty
   }//for every row

   delete[] p;
   delete[] d;
   delete[] l;
   delete[] a;
   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //UBACENA JE LISTA; SVAKI RED MATRICE L (ZAJEDNO SA FILLIN-ovima)
  //SE PUNI IZ LISTE; POBOLJSANO ITERIRANJE PO LISTI;
  //umesto itc->find(j) ide itc.begin(); SVUDA WHILE-LOOP
  friend SparseMatrix<T> factor4(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2; //for double traversing
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   //f
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   //T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...
   //list<unsigned int> lt;

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //iteriranje iz reda u red
   it=A.rows_.begin();
   memset(d,0,n*sizeof(T));
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       lt.clear();
       lit0 = lt.insert(lt.end(), n); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           /*p[j] += a[j];*/
           p[j] += itr->second;

           if(++itr != it->end())//mora provera; pravi probleme ako je a(k>=i)=0
             j2 = itr->first;
           else
             j2=i;

           if( j != *lit0 )
            lit0 = lt.insert(lit0, j);
         }

         lit = lit0;
         //itc = ++L.cols_[j].find(j);
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           k = itc->first;
           p[k] -= p[j] * itc->second;

           while(k > *(++lit));

           if(n == *lit) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             lit=lt.insert(lit, k);
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               p[k] -= p[j] * itc->second;
               lit=lt.insert(++lit, k);
             }
             break; //exit while(++itc != L.cols_[j].end()) loop
           }

           else if(k < *lit)
             lit=lt.insert(lit, k);//umecemo novi ako je < a ako je == onda ne

         }//for all subdiagonal items in column j

         //find out the next item
         if( *(++lit0) <= j2 ) //next item in list for next iteration
           j = *lit0;
         else //A(i,j2) is overjumped OR we have just reached the diagonal
           j = j2;
       }//while j<i

       //u cout << "(" << i << "," << j << ")   ";
       lit = lt.begin();  j = *lit;
       while(j<i)
       //for(lit=lt.begin(), j=*lit; j<i; ++lit, j=*lit)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);

         //u cout << j << " ";
         j=*(++lit);
       }
       //u cout << endl;

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

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
     }//if not empty
     ++it; ++i;
   }//while iterating every row

   delete[] p;
   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //UBACENA JE LISTA; SVAKI RED MATRICE L (ZAJEDNO SA FILLIN-ovima)
  //SE PUNI IZ LISTE; POBOLJSANO ITERIRANJE PO LISTI;
  //umesto itc->find(j) ide itc.begin(); SVUDA WHILE-LOOP
  friend SparseMatrix<T> factor45(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2; //for double traversing
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<unsigned int> lt;   //f
   list<unsigned int>::iterator lit0, lit;

   unsigned int n=A.rows();
   T* p = new T[n]; // p - proizvod l(ij)*d(j) ; da bi se izbeglo visestruko izracunavanje
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   T* l = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   //T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...
   //list<unsigned int> lt;

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //iteriranje iz reda u red
   it=A.rows_.begin();
   memset(d,0,n*sizeof(T));
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       memset(p,0,i*sizeof(T));
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       lt.clear();
       lit0 = lt.insert(lt.end(), n); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           /*p[j] += a[j];*/
           p[j] += itr->second;

           if(++itr != it->end())//mora provera; pravi probleme ako je a(k>=i)=0
             j2 = itr->first;
           else
             j2=i;

           if( j != *lit0 )
            lit0 = lt.insert(lit0, j);
         }

         lit = lit0;
         //itc = ++L.cols_[j].find(j);
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           k = itc->first;
           p[k] -= p[j] * itc->second;

           while(k > *(++lit));

           if(n == *lit) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             lit=lt.insert(lit, k);
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               p[k] -= p[j] * itc->second;
               lit=lt.insert(++lit, k);
             }
             break; //exit for loop
           }

           else if(k < *lit)
             lit=lt.insert(lit, k);//umecemo novi ako je < a ako je == onda ne

         }//for all subdiagonal items in column j

         //find out the next item
         if( *(++lit0) <= j2 ) //next item in list for next iteration
           j = *lit0;
         else //A(i,j2) is overjumped OR we have just reached the diagonal
           j = j2;
       }//while j<i

       //u cout << "(" << i << "," << j << ")   ";
       lit = lt.begin();  j = *lit;
       while(j<i)
       //for(lit=lt.begin(), j=*lit; j<i; ++lit, j=*lit)
       {
         l[j] = p[j]/d[j];
         d[i] -= p[j] * l[j];
         L.push_fast(i+1, j+1, l[j]);

         //u cout << j << " ";
         j=*(++lit);
       }
       //u cout << endl;

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

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
     }//if not empty
     ++it; ++i;
   }//while iterating every row

   delete[] p;
   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //ISTO KAO factor4; IZBACEN memset
  friend SparseMatrix<T> factor5(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j, k;
   unsigned int j2; //for double traversing
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
   //Matrix L(n,n);

   //svi dijagonalni elementi za pocetak nula
   memset(d,0,n*sizeof(T));
   //iteriranje iz reda u red
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       lt.clear();
       lit0 = lt.insert(lt.end(), make_pair(n, A.nullVal)); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           if( j != lit0->first )
             lit0 = lt.insert(lit0, make_pair(j, itr->second));
           else
             lit0->second += itr->second;  //p[j] = itr->second;

           if(++itr != it->end())//mora provera; pravi probleme ako je a(k>=i)=0
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

           if(n == lit->first) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             lit=lt.insert(lit, make_pair(k, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
             {
               k = itc->first;
               lit=lt.insert(++lit, make_pair(k, -lit0->second*itc->second));
             }
             break; //exit outer "while(++itc!=L.cols_[j].end())" loop
           }

           //umecemo novi ako je < a ako je == onda ne; nego samo
	   //korigujemo staru vrednost i idemo u sledecu iteraciju
	   else if(k < lit->first)
             lit=lt.insert(lit, make_pair(k, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }//while(++itc != L.cols_[j].end()) i.e. for all subdiagonal items in column j

         //find out the next item in list
         if( (++lit0)->first <= j2 ) //next item in list for next iteration is fill-in
           j = lit0->first;
         else //next item is A(i,j2) OR we have just reached the diagonal
           j = j2;
       }//while j<i i.e for all elements in row before the diagonal element

       lit = lt.begin();  j = lit->first;
       while(j<i)
       //for(lit=lt.begin(), j=*lit; j<i; ++lit, j=*lit)
       {
         l[j] = lit->second/d[j];
         d[i] -= lit->second * l[j];
         L.push_fast(i+1, j+1, l[j]);

         j=(++lit)->first;
       }

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         delete[] l;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }//if not empty
     ++it; ++i;
   }//while iterating every row

   delete[] d;
   delete[] l;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //ISTO KAO factor5; uz neka poboljsanja;
  //izbacen k, izbacena petlja za l[i],d[i], izbacen niz l[]
  //NAJBRZA VARIJANTA!!!  (13.09.2006.)
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
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //svi dijagonalni elementi za pocetak nula
   memset(d,0,n*sizeof(T));

   //iteriranje iz reda u red
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //ako red nije prazan
     {
       //racunanje cesto koriscenih proizvoda - p
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       lt.clear();
       lit0 = lt.insert(lt.end(), make_pair(n, A.nullVal)); //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           if( j != lit0->first )
             lit0 = lt.insert(lit0, make_pair(j, itr->second));
           else
             lit0->second += itr->second;  //p[j] = itr->second;

           if(++itr != it->end())//mora provera; pravi probleme ako je a(k>=i)=0
             j2 = itr->first;
           else
             j2=i;
         }

         lit = lit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           while(itc->first > (++lit)->first);

           if(n == lit->first) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               lit=lt.insert(++lit, make_pair(itc->first, -lit0->second*itc->second));

             break; //exit outer "while(++itc!=L.cols_[j].end())" loop
           }

       //umecemo novi ako je < a ako je == onda ne; nego samo
	   //korigujemo staru vrednost i idemo u sledecu iteraciju
	   else if(itc->first < lit->first)
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }//while(++itc != L.cols_[j].end()) i.e. for all subdiagonal items in column j

         lj = lit0->second/d[j];
         L.push_fast(i+1, j+1, lj);
         d[i] -= lit0->second * lj;

         //find out the next item in list
         if( (++lit0)->first <= j2 ) //next item in list for next iteration is fill-in
           j = lit0->first;
         else //next item is A(i,j2) OR we have just reached the diagonal
           j = j2;
       }//while j<i i.e for all elements in row before the diagonal element

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }//if not empty

     ++it; ++i;
   }//while iterating every row

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //UMESTO LISTE - UPOTREBA sdt::map<unsigned, T>
  //ali se map upotrebljava kao list, razlika je samu u brzini iteriranja
  friend SparseMatrix<T> factor7(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0, j;
   unsigned int j2; //for double traversing
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   map<unsigned int, T> mp;   //f
   typename map<unsigned int, T>::iterator mit0, mit;

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //svi dijagonalni elementi za pocetak nula
   memset(d,0,n*sizeof(T));
   //iteriranje iz reda u red
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       itr=it->begin();
       j = j2 = itr->first; //j2 - sledeci A tj. A(i,j2)
       mp.clear();
       mit0 = mp.insert(make_pair(n, A.nullVal)).first; //trick!!
       while(j<i)
       {
         if(j==j2)
         {
           if( j != mit0->first )
             mit0 = mp.insert(mit0, make_pair(j, itr->second));
           else
             mit0->second += itr->second;  //p[j] = itr->second;

           if(++itr != it->end())//mora provera; pravi probleme ako je a(k>=i)=0
             j2 = itr->first;
           else
             j2=i;
         }

         mit = mit0;
         itc = L.cols_[j].begin();
         while(++itc != L.cols_[j].end())
         {
           while(itc->first > (++mit)->first);

           if(n == mit->first) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             mp.insert(make_pair(itc->first, -mit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               mp.insert(make_pair(itc->first, -mit0->second*itc->second));

             break; //exit outer "while(++itc!=L.cols_[j].end())" loop
           }

           //umecemo novi ako je < a ako je == onda ne; nego samo
	   //korigujemo staru vrednost i idemo u sledecu iteraciju
	   else if(itc->first < mit->first)
             mit = mp.insert(mit, make_pair(itc->first, -mit0->second*itc->second));
           else
             mit->second -= mit0->second * itc->second;
         }//while(++itc != L.cols_[j].end()) i.e. for all subdiagonal items in column j

         lj = mit0->second/d[j];
         L.push_fast(i+1, j+1, lj);
         d[i] -= mit0->second * lj;

         //find out the next item in list
         if( (++mit0)->first <= j2 ) //next item in list for next iteration is fill-in
           j = mit0->first;
         else //next item is A(i,j2) OR we have just reached the diagonal
           j = j2;
       }//while j<i i.e for all elements in row before the diagonal element

       //... a zatim racunamo i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }//if not empty

     ++it; ++i;
   }//while iterating every row

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //UMESTO LISTE - UPOTREBA sdt::map<unsigned, T>
  //ali se map upotrebljava standardno kao map
  friend SparseMatrix<T> factor8(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   map<unsigned int, T> mp;   //f
   typename map<unsigned int, T>::iterator mit; //, mit2; //1. varijanta

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje
   //pair<typename map<unsigned int, T>::iterator, bool> pr; // 2. varijanta

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //svi dijagonalni elementi za pocetak nula
   memset(d,0,n*sizeof(T));
   //iteriranje iz reda u red
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       mp.clear();
       itr=it->begin();
       while(itr!=it->end() && itr->first<i)
       {
         mp.insert(make_pair(itr->first, itr->second)); //ubaci osnovne elemente
         ++itr;
       }

       //ubaci fill-ins-ove
       mit=mp.begin();
       //while(j<i)
       while(mit!=mp.end())
       //while(mit->first < i)
       {
         itc = L.cols_[mit->first].begin();
         while(++itc != L.cols_[mit->first].end())
         {  /*  1. varijanta
           mit2 = mp.find(itc->first);

           //korigujemo staru vrednost i idemo u sledecu iteraciju
           if(mit2 != mp.end())
             mit2->second -= mit->second * itc->second;

           //umecemo novi ako je nije pronadjen odgovarajuci key
           else
             mp.insert(make_pair(itc->first, -mit->second*itc->second));
           */

           /*  2. varijanta
           pr = mp.insert(make_pair(itc->first, -mit->second*itc->second));
           if(!pr.second)
             pr.first->second -= mit->second * itc->second;
           */

           // 3. varijanta
           mp[itc->first] -= mit->second * itc->second;

           //++itc;
         }//while(itc != L.cols_[mit->first].end()) i.e. for all subdiagonal items in column j

         //racunamo l[j], d[i] i inkrementiramo mit za sledecu operaciju
         lj = mit->second/d[mit->first];
         L.push_fast(i+1, mit->first+1, lj);
         d[i] -= mit->second * lj;
         //u cout << i << "," << mit->first/*j*/ << "," << lj << endl; //u

         ++mit;
       }//while j<i i.e for all elements in row before the diagonal element

       //... a zatim racunamo i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }//if not empty

     ++it; ++i;
   }//while iterating every row

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  //brza UCt dekompozicija simetricne kvadratne matrice (Doolitle)
  //(za sada bez pivotiranja)
  //SLICNO KAO factor6;
  //RAZIKA: prvo idu copy elementi p tek naknadno fillins-ovi
  friend SparseMatrix<T> factor9(const SparseMatrix<T>& A)
  {
   //f if (A.rows != A.cols)
      //f error("F-ja factor()\nA nije kvadratna matrica!");

   unsigned int i=0;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr;
   typename SparseMatrix<T>::onecol_const_iterator itc;
   list<pair<unsigned int, T> > lt;   //f
   typename list<pair<unsigned int, T> >::iterator lit0, lit;

   unsigned int n=A.rows();
   T lj;
   T* d = new T[n]; // d - dijagonalni clan d(j,j) ; da bi se izbeglo visestruko pretrazivanje

   SparseMatrix<T> L(n,n);
   //Matrix L(n,n);

   //svi dijagonalni elementi za pocetak nula
   memset(d,0,n*sizeof(T));
   //iteriranje iz reda u red
   it=A.rows_.begin();
   while( it!=A.rows_.end() )
   {
     if( !it->empty() )   //f neki text
     {
       //racunanje cesto koriscenih proizvoda - p
       lt.clear();
       itr=it->begin();
       while(itr!=it->end() && itr->first<i)
       {
         lt.push_back(make_pair(itr->first, itr->second)); //ubaci osnovne elemente
         ++itr;
       }

       lt.push_back(make_pair(n, T())); //Trick !!!
       lit0 = lt.begin();
       while( lit0->first < i )
       {
         int j = lit0->first;   //u
         T pj = lit0->second;      //u

         lit = lit0;
         itc = L.cols_[lit0->first].begin();
         while(++itc != L.cols_[lit0->first].end())
         {
           while(itc->first > (++lit)->first);

           if(n == lit->first) //dosli smo do kraja i sada samo dodajemo na kraj
           {
             lit=lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
             while(++itc!=L.cols_[j].end())
               lit=lt.insert(++lit, make_pair(itc->first, -lit0->second*itc->second));

             break; //exit outer "while(++itc!=L.cols_[j].end())" loop
           }

           //umecemo novi ako je < a ako je == onda ne; nego samo
	   //korigujemo staru vrednost i idemo u sledecu iteraciju
	   else if(itc->first < lit->first)
             lit = lt.insert(lit, make_pair(itc->first, -lit0->second*itc->second));
           else
             lit->second -= lit0->second * itc->second;
         }//while(++itc != L.cols_[j].end()) i.e. for all subdiagonal items in column j

         lj = lit0->second/d[lit0->first];
         L.push_fast(i+1, lit0->first+1, lj);
         d[i] -= lit0->second * lj;

         ++lit0;
       }//while j<i i.e for all elements in row before the diagonal element

       //... a zatim i dijagonalni clan L(i,i) = d(i)
       if( itr->first == i )
         d[i] += itr->second;
         //d[i] += a[i];//A(i,i);

       if(d[i]==0.0)
       {
         cout << "Singularity!" << endl;
         L.clear();
         delete[] d;
         return L;
       }

       L.push_fast(i+1, i+1, d[i]);
     }//if not empty

     ++it; ++i;
   }//while iterating every row

   delete[] d;

   return L;
 }

//---------------------------------------------------------------------------
  //brza Holeskijeva dekompozicija simetricne pozitivno definitne
  //kvadratne matrice (za sada bez pivotiranja)
  friend SparseMatrix<T> ltl(/*const //f */SparseMatrix<T>& A)
  {
   //if (A.rows() != A.rows())
      //error("F-ja lowerL()\nNije kvadratna matrica!");

   //double MINVAL = 0.00001;
   //double MINVALKV = MINVAL*MINVAL;

   unsigned int i, j, k;
   typename SparseMatrix<T>::rows_const_iterator it;
   typename SparseMatrix<T>::onerow_const_iterator itr1;
   typename SparseMatrix<T>::onecol_const_iterator itr2;

   unsigned int n=A.rows();
   T* sum = new T[n]; // l - poddijagonalni clanovi (L(i,j),j<i) u i-tom redu...
   T* a = new T[n]; // a - poddijagonalni clanovi (A(i,j),j<=i) u i-tom redu...
   T dd;
   SparseMatrix<T> L(n, n); //left lower matrix of aa

   for(it=A.rows_.begin(), j=1, memset(sum,0,n*sizeof(T));//f n
       it!=A.rows_.end(); ++it, ++j, memset(sum, 0, n*sizeof(T)))
   {
     //izvlacenje i-tog reda matrice A
     for(itr1=it->begin(), memset(a,0,n*sizeof(T)); //n
         itr1!=it->end(); ++itr1)
       a[itr1->first] = itr1->second;

     //racunanje privreminih suma
     for(itr1=L.rows_[j-1].begin(), i=itr1->first;
         itr1!=L.rows_[j-1].end(); ++itr1, i=itr1->first)
     {
       for(itr2=L.cols_[i].find(j-1); itr2!=L.cols_[i].end(); ++itr2)
         sum[itr2->first] += itr1->second * itr2->second;
     }//for every left-row member - iterate corresponding right row

     for(k=j; k<=n; k++)
     {
       if(j==k)
       {
         dd = sqrt(a[j-1]/*A(j,j)*/-sum[k-1]);
         L.push_fast( j,j, dd );
       }
       //if(j==k)
         //l.mat[j][j] = aa.mat[j][j]>sum+MINVALKV ? sqrt(aa.mat[j][j]-sum): MINVAL;
       else
         L.push_fast( k,j, (a[k-1]/*A(j,k)*/-sum[k-1])/dd/*L(j,j)*/ );
     } //k loop
   }//j loop

 //delete[] sum;
 return L;
 }

//---------------------------------------------------------------------------
  //friend istream& operator>> (istream&, SparseMatrix&);
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
        //if( sm(i,j) >= sm.nullVal )  os << " ";
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
  //const T nullValC;
  vector<row_type> rows_;
  vector<col_type> cols_;

public:
//---------------------------------------------------------------------------
  //destructor
  ~SparseMatrix() { clear(); }

//---------------------------------------------------------------------------
  //default constructor
  SparseMatrix() : nullVal()//, nullValC()
  {
    clear();
  };

//---------------------------------------------------------------------------
  //copy constructor
  SparseMatrix(const SparseMatrix<T>& sm2) : nullVal()//, nullValC()
  {
    copy(sm2);
  };

//---------------------------------------------------------------------------
  //creates SparseMatrix of mxn dimension
  SparseMatrix(unsigned int m, unsigned int n) : nullVal()//, nullValC()
  {
    //fill row-map
    rows_.clear();
    for(unsigned int i=0; i<m; ++i)
    {
      map<unsigned int, T> temprow;
      temprow.clear();
      rows_.push_back(temprow);
    }

    //fill column-map
    cols_.clear();
    for(unsigned int i=0; i<n; ++i)
    {
      map<unsigned int, T> tempcol;
      tempcol.clear();
      cols_.push_back(tempcol);
    }
  };

//---------------------------------------------------------------------------
  //clear the entire content
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
  //returns the number of rows
  unsigned int rows() const  { return rows_.size();};

//---------------------------------------------------------------------------
  //returns the number of columns
  unsigned int cols() const  { return cols_.size();};

//---------------------------------------------------------------------------
  //f postoje stl::copy zajebancije
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
  //opt
  SparseMatrix<T>& operator= (const SparseMatrix<T>& sm2)
  {
    if( this != &sm2 )
      copy(sm2);

    return *this;
  };

//---------------------------------------------------------------------------
  //opt
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
  //nalazi broj clanova razlicitih od nule kao i "stepen gustine"
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
  //opt
  friend Matrix solveUC(const SparseMatrix<T>& A, const Matrix& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();
    ///* //f
    if (n != c)
      error("F-ja solveUC()\nA nije kvadratna matrica!");
    if (n != y.getRows())
      error("F-ja solveUC()\nA i b nemaju isti broj redova!");
    if (y.getCols() != 1)
      error("F-ja solveUC()\nb nije vektor!");
    //*/

    Matrix y2(n, 1); //promenjeni vektor y - resenje sistema L*y2 = y
    SparseMatrix<T> p = factor6(A);
    for(i=1; i<=n; i++)
      for(k=1, y2(i,1) = y(i,1); k<i; k++)
        y2(i,1) -= p(i,k)*y2(k,1); //Dk * Lik * Yk

    Matrix x(n, 1); //vektor x - resenje sistema D*trans(L)*x = y2 (x se upisuje u y2)
    for(i=n; i>=1; i--)
      for(k=n, x(i,1) = y2(i,1)/p(i,i); k>i; k--)
        x(i,1) -= p(k,i)*x(k,1); //Lki * Xk

    return x;
  };

//---------------------------------------------------------------------------
  //opt
  //Pazi: moze i da se specificira tip (umesto <T> ide <double>)
  friend vector<double> solveUC2(const SparseMatrix<T>& A, const vector<double>& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();
    ///* //f
    if (n != c)
      error("F-ja solveUC2()\nA nije kvadratna matrica!");
    if (n != y.size())
      error("F-ja solveUC()\nA i y nemaju isti broj redova!");
    //*/

    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;
    typename SparseMatrix<T>::cols_const_iterator it2, itn, it0;
    typename SparseMatrix<T>::onecol_const_iterator itc, itcend, itcbeg;

    vector<double> y2(n); //promenjeni vektor y - resenje sistema L*y2 = y
    SparseMatrix<T> p = factor6(A);
    for(it=p.rows_.begin(), i=0; it!=p.rows_.end(); ++it, ++i)
      for(itr=it->begin(), k=itr->first, y2[i] = y[i]; itr!=--it->end(); ++itr, k=itr->first)
        y2[i] -= itr->second * y2[k];

    vector<double> x(n); //vektor x - resenje sistema D*trans(L)*x = y2 (x se upisuje u y2)
    itn = --p.cols_.end();
    it0 = --p.cols_.begin();
    //for(it2=p.cols_.rbegin(), i=n; it2!=p.cols_.rend(); --it2, --i)
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
  //opt
  //Pazi: moze i da se specificira tip (umesto <T> ide <double>)
  friend vector<double> solveUC3(const SparseMatrix<T>& A, const vector<double>& y)
  {
    unsigned int i, k;
    unsigned int n = A.rows();
    unsigned int c = A.cols();
    ///* //f
    if (n != c)
      error("F-ja solveUC2()\nA nije kvadratna matrica!");
    if (n != y.size())
      error("F-ja solveUC()\nA i y nemaju isti broj redova!");
    //*/

    typename SparseMatrix<T>::rows_const_iterator it;
    typename SparseMatrix<T>::onerow_const_iterator itr;
    //typename SparseMatrix<T>::cols_const_iterator it2, itn, it0;
    typename SparseMatrix<T>::cols_const_riterator it2;
    //typename SparseMatrix<T>::onecol_const_iterator itc, itcend, itcbeg;
    typename SparseMatrix<T>::onecol_const_riterator itc;

    vector<double> y2(n); //promenjeni vektor y - resenje sistema L*y2 = y
    SparseMatrix<T> p = factor6(A);
    for(it=p.rows_.begin(), i=0; it!=p.rows_.end(); ++it, ++i)
      for(itr=it->begin(), k=itr->first, y2[i] = y[i]; itr!=--it->end(); ++itr, k=itr->first)
        y2[i] -= itr->second * y2[k];

    vector<double> x(n); //vektor x - resenje sistema D*trans(L)*x = y2 (x se upisuje u y2)
    for(it2=p.cols_.rbegin(), i=n-1; it2!=p.cols_.rend(); ++it2, --i)
	  for(itc=it2->rbegin(), k=itc->first, x[i] = y2[i]/p(i+1,i+1); itc!=--it2->rend(); ++itc, k=itc->first)
        x[i] -= itc->second * x[k];

    return x;
  };

//---------------------------------------------------------------------------
  //opt
  //Pazi: moze i da se specificira tip (umesto <T> ide <double>)
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
  //opt
  //Pazi: moze i da se specificira tip (umesto <T> ide <double>)
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
  //f ovo moze dodatno da se optimizije //opt
  //f verovatno moze expresion template za npr. a(i,j)=3.89;
  //operator() provides only read-access, at the moment
  T& operator() (unsigned int i, unsigned int j)
  {
    //if( i>rows() || j>cols() )
      //return nullVal;

    onerow_iterator it = rows_[i-1].find(j-1);
    return it != rows_[i-1].end() ? it->second : nullVal;
    //onecol_iterator it = cols_[j-1].find(i-1);
    //return it != cols_[j-1].end() ? it->second : nullVal;
  };

//---------------------------------------------------------------------------
  ///*//f //u
  //operator() provides only read-access, at the moment
  const T operator() (unsigned int i, unsigned int j) const //would return ref to temp
  {
    if( i>rows() || j>cols() )
      return nullVal;

    onerow_const_iterator it = rows_[i-1].find(j-1);
    return it != rows_[i-1].end() ? it->second : nullVal;
    //onecol_iterator it = cols_[j-1].find(i-1);
    //return it != cols_[j-1].end() ? it->second : nullVal;
  };
  //*/

//---------------------------------------------------------------------------
  //push provides only write-access, at the moment
 //first (upper left) element is (1,1) not (0,0)
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
        //rows_[i-1].insert( pair<unsigned int, T>(j-1, val) );
        //cols_[j-1].insert( pair<unsigned int, T&>(i-1, val) );

      //cout << i << "," << j << "  " << val << "  " << cols_[j-1][i-1] << endl;
    }
  };

//---------------------------------------------------------------------------
  //push provides only write-access, at the moment
  //it is to be used only when there isn't (i,j) element (100% sure)
  //and i, j are within rows(), cols()
 //first (upper left) element is (1,1) not (0,0)
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
