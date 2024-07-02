#ifndef LU_DECOMP_H
#define LU_DECOMP_H


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>


#include "util.h"

template <class number>
class LU_Decomposition
{

public:
    LU_Decomposition(){}
    void ring_order(int order[][2],int place[][2],int v_len,int dimension);
    int  lu_decomp(number **a,int *place,int num);
    void solve_lu(number **a,number *x,number *b,int *place,int num);
    void diag3(number a[][3],number x[],number b[],int size,int num);
    int inverse(number **a,number **res,int num);
};




template <class number>
void LU_Decomposition<number>::ring_order(int order[][2],int place[][2],int v_len,int dimension)

/*   ***************************************************
         Subroutine to calculate the grey codes
         for a Hyper-Cube.  The program also
         finds the links to the neighbors.

         Return Values:
         order[i][0] : Number of the previous node
         order[i][1] : Number of the next node
         place[i][0] : Node number that is in the
         i'th position.
         place[i][1] : Position of node i in the ring

     ***************************************************

         */
{
  int i,count,here;
  int m,*bits;

  bits = new int[v_len];

  bits[0] = 0;
  bits[1] = 1;
  m = 2;
  for (i=1;i<dimension;m=1<<++i)
        for (count=0;count<m;++count)
            bits[m+count] = m | bits[m-1-count];

  m = 1 << dimension;
  for (i=0;i<m;++i)
        {
      here = bits[i];
      place[i][0] = here;
      place[here][1] = i;
      count = (m+i-1) % m;
      order[here][0] = bits[count];
      count = (i+1) % m;
      order[here][1] = bits[count];
        }

  delete bits;
}






template <class number>
int  LU_Decomposition<number>::lu_decomp(number **a,int *place,int num)
{

    /*
*******************************************************
Subroutine to perform an L/U decomposition
in place on the matrix A[][]
*******************************************************
*/


  number *work;
  int i,j,k;
  int tmp,pivot;


  ArrayUtils<number> arrays;
  work = arrays.onetensor(num);
  if(work==NULL) {
    perror("Error - Could not allocate memory in lu_decomp\n");
    exit(2);
  }


  for (i=0;i<num;++i)
  {
    place[i] = i;                   // initialize the permutation vector
    work[i] = fabs(a[i][0]);        // determine which entry in each row is "largest"
    for (j=1;j<num;++j)
      if (work[i] < fabs(a[i][j]))
                work[i] = fabs(a[i][j]);
  }

  for (i=0;i<num;++i)
  {
    pivot = i;                                                     // Determine if we need to pivot
    for (j=i+1;j<num;++j)                                          // see which row has biggest ratio
            if (fabs(a[place[j]][i])/work[place[j]] >
                    fabs(a[place[pivot]][i])/work[place[pivot]])
                pivot = j;                                         // This ratio is more bigger than
                                                                   // the current one.

    tmp = place[i];
    place[i] = place[pivot];
    place[pivot] = tmp;

    if (fabs(a[place[i]][i]) < 1.0E-9)
            {
                arrays.delonetensor(work);
                return(0);
            }

    for (j=i+1;j<num;++j)
            {
                if (fabs(a[place[j]][i]) > 1.0E-12)
                {
                    a[place[j]][i] = a[place[j]][i]/a[place[i]][i];
                    for (k=i+1;k<num;++k)
                        a[place[j]][k] -= a[place[j]][i]*a[place[i]][k];
                }
            }

  }

  //for(i=0;i<num;++i)
  //      std::cout << i << " " << place[i] << std::endl;

  arrays.delonetensor(work);
  return(1);
}






template <class number>
void LU_Decomposition<number>::solve_lu(number **a,number *x,number *b,int *place,int num)
/*
    Subroutine to solve the linear equation Ax=b for x.
    A is a matrix that has been decomposed into its L/U
    decomposition.  x and b are vectors where b is given
    and x is returned.  place is the index vector returned
    by the L/U decomposition.  size is the size of the
    matrix as it was declared and num is the number of
    values used in this application.
*/

{
    number *c;
    int i,j;

    ArrayUtils<number> arrays;
    c = arrays.onetensor(num);

    c[0] = b[place[0]];
    for (i=1;i<num;++i)
    {
        c[i] = b[place[i]];
        for (j=0;j<i;++j)
            c[i] -= a[place[i]][j]*c[j];
    }

    x[num-1] = c[num-1]/a[place[num-1]][num-1];
    for (i=num-2;i>=0;--i)
    {
        x[i] = c[i];
        for (j=i+1;j<num;++j)
            x[i] -= x[j]*a[place[i]][j];
        x[i] = x[i]/a[place[i]][i];
    }

    arrays.delonetensor(c);

}



template <class number>
void LU_Decomposition<number>::diag3(number a[][3],number x[],number b[],int,int num)
/*
**********************************************
Subroutine to solve the linear system Ax=b.
A is a tridiagonal matrix, and x and b are
vectors.  In this routine A is destroyed.
**********************************************
*/

{
    int i;
    number fact;

    for (i=0;i<num-1;++i)
    {
            fact = a[i+1][0]/a[i][1];
            a[i+1][1] -= fact*a[i][2];
            b[i+1] -= fact*b[i];
    }

    x[num-1] = b[num-1]/a[num-1][1];
    for (i=num-2;i>=0;--i)
    x[i] = (b[i]-x[i+1]*a[i][2])/a[i][1];

}




template <class number>
int LU_Decomposition<number>::inverse(number **a,number **res,int num)

/*
***************************************************
Subroutine to invert the matrix A.  The matrix is
inverted by using an L/U decomposition in place
on A and solving the necessary systems.  A is
returned as its L/U decomposition and the result
is returned in RES.
***************************************************
*/

{
  number *wrk1,*wrk2;
  int i,j;
  int *place;

  ArrayUtils<number> arrays;
  wrk1 = arrays.onetensor(num);
  wrk2 = arrays.onetensor(num);
  place = new int[num];
  if (place==nullptr) {
    printf("Could not allocate memory in inverse\n");
    exit(2);
  }

  if (!lu_decomp(a,place,num))
        {
            free(wrk1);
            free(wrk2);
            return(0);
        }

  for (i=0;i<num;++i)
        wrk1[i] = 0.0;

  wrk1[0] = 1.0;
  solve_lu(a,wrk2,wrk1,place,num);
  for (i=0;i<num;++i)
        res[i][0] = wrk2[place[i]];

  for (j=1;j<num;++j)
        {
      wrk1[j-1] = 0.0;
      wrk1[j] = 1.0;
      solve_lu(a,wrk2,wrk1,place,num);
      for (i=0;i<num;++i)
                res[i][j] = wrk2[place[i]];
        }

  arrays.delonetensor(wrk1);
  arrays.delonetensor(wrk2);
  delete place;
  return(1);

}



#endif // LU_DECOMP_H
