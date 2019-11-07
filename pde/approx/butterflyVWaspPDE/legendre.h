#ifndef LEGENDRE_H
#define LEGENDRE_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "util.h"



template <class number>
class Legendre
{

public:
    Legendre();

    // Take the exponential of x (raised to n)
    static number pythag(number a, number b);



  /*  QL Algortihm to find the eigen vectors and eigen
      values of a real symmetric tridiagonal matrix.
      Found on page 480 of Numerical Recipes in C, 2nd ed.

            d[ ] = diagonal elements
            e[ ] = subdiagonal elements (e[1] arbitrary)

  */
    static void tqli(number *d,number *e,int n,int size,number **z);



  /* Routine to find the Legendre-Gauus-Lobatto quadrature.
     Uses Golub's method, SIAM Review, Vol. 15, No.2 April, 73,
     Page 318.
  */
    static void leg_quad(number *x,number *w,int n);


  /* Routine to find the Legendre-Gauus quadrature.
     Uses Golub's method, SIAM Review, Vol. 15, No.2 April, 73,
     Page 318.
  */
    static void leg_quad_Gauss(number *x,number *w,int n);


    // Values of the Legendre polynomials at the point defined by *x
    // Returns a matrix with all the polynomials up to a given
    // degree:
    //      interp[i][j] = L_i(x_j)
    //
    static void leg_val(number **interp,number *x,int n,int rows);

    // Values of the Legendre polynomials at the point defined by *x
    // for the degree given by degree.
    // Returns only the polynomials to the indicated degree:
    //    interp[i] = L_num(x_i)
    //
    static void leg_val(number *interp ,number *x,int n,int degree);



    // Values of the derivatives of the Legendre polynomials at the
    // specified points (x) for the given degree (degree)
    //      interp[i] = L'_degree(x_i)
    //         i = 0...num
    static void leg_der_val(number *interp,number *x,int num,int degree);


    // Values of the derivatives of the Legendre polynomials at the
    // specified points (x) for the given degree (degree)
    //      interp[i] = L''_degree(x_i)
    //         i = 0...num
    static void leg_2_der_val(number *interp,number *x,int num,int degree);


    // Legendre (Gauss-Lobatto) collocation derivative matrix
    static void leg_der(number **d1,number **lval,number x[],int n,int size);

    // Legendre (Gauss) collocation derivative matrix
    static void leg_gauss_der(number **d1,number *x,int n);

    // Routine to calculate the entries in the stiffness
    // matrix for a Legendre-collocation type of method.
    static void stiffLeg(number **stiff,number *w,number **D1,int N);

    // Values of the Lagrange interpolates (Gauss) at the endpoints
    static void gauss_collocation_end_values(number *left,number *right,
                                                                                     number *x,int n);

    // Petrov-Galerkin mass and stiffness matrices
    static void petrov(number **mass,number **stiff,number *x,number *w,
                                         number **lval,int n,int size);


    // Legendre-Tau stiffness and mass matrices
    static void tauMatrices(number **stiff,number **mass,number **der,
                                                    int numuse);

    /* ******************************************
         Subroutine to initialize the stiffnes,
         mass, and derivative matrices for the
         galerkin method.
         ********************************************* */
    static void galerkinMatrices(number **stiff,number **mass,number **der,
                                                             int numuse);


    /* ******************************************
         Subroutine to initialize the stiffnes,
         mass, and derivative matrices for the
         galerkin method. (infinite domain)
         ********************************************* */
    static void galerkinInfinite(number **rightStiff,number **leftStiff,
                                                             number **mass,number **der,
                                                             int num);

};




static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b)>0.0 ? fabs(a) : -fabs(a))







// Take the expnential of x (raised to n)
template <class number>
number power(number x,int n) {
  number y;
  int i;

  for (i=0,y=1.0;i<abs(n);++i)
    y *= x;

  if (n<0)
    return(1.0/y);
  else
    return(y);

}



// Find the length of the side of a triang;e
// Taken from numerical recipes in C
template <class number>
number Legendre<number>::pythag(number a,number b) {
  number absa,absb;

  absa = fabs(a);
  absb = fabs(b);
  if (absa>absb) {
    /* s = absb/absa; */
    return( absa*sqrt(1.0+SQR(absb/absa)));
  }
  else {
    /* s = absa/absb; */
    return( absb*sqrt(1.0+SQR(absa/absb)));
  }
}


template <class number>
void Legendre<number>::tqli(number *d,number *e,int n,int,number **z) {
  /*  QL Algortihm to find the eigen vectors and eigen
      values of a real symmetric tridiagonal matrix.
      Found on page 480 of Numerical Recipes in C, 2nd ed.

            d[ ] = diagonal elements
            e[ ] = subdiagonal elements (e[1] arbitrary)

  */

  int m,l,iter,i,k;
  number s,r,p,g,f,dd,c,b;

  for(i=1;i<=n;++i) {
    for(k=1;k<=n;++k)
      z[i][k] = 0.0;
    z[i][i] = 1.0;
  }

  for(i=2;i<=n;++i) e[i-1] = e[i];
  e[n] = 0.0;
  for(l=1;l<=n;++l) {
    iter = 0;
    do {

      for(m=l;m<=n-1;m++) {
                dd = fabs(d[m]) + fabs(d[m+1]);
                if (static_cast<number>(fabs(e[m]+dd)) == dd) break;
      }

      if(m != l) {
                if (iter++==300) {
                    perror("To Many iterations in tqli\n");
                    exit(3);
                }
                g = (d[l+1]-d[l])/(2.0*e[l]);
                r = pythag(g,1.0);
                g = d[m] - d[l] + e[l]/(g+SIGN(r,g));
                s=(c=1.0);
                p = 0.0;

                for(i=m-1;i>=l;--i) {
                    f = s*e[i];
                    b = c*e[i];
                    e[i+1] = (r=pythag(f,g));
                    if(fabs(r)==0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f/r;
                    c = g/r;
                    g = d[i+1] - p;
                    r = (d[i]-g)*s+2.0*c*b;
                    d[i+1] = g + (p=s*r);
                    g = c*r-b;

                    for(k=1;k<=n;++k) {
                        f = z[k][i+1];
                        z[k][i+1] = s*z[k][i]+c*f;
                        z[k][i]   = c*z[k][i]-s*f;
                    }
                }

                if(r == 0.0 && i) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
      }

    } while(m!=l);
  }
}


template <class number>
void Legendre<number>::leg_quad(number *x,number *w,int n) {
  /* Routine to find the Legendre-Gauus-Lobatto quadrature.
     Uses Golub's method, SIAM Review, Vol. 15, No.2 April, 73,
     Page 318.
  */

  number *d,*e,**v;
  int i,j,k;
  number norm,p;

  d = ArrayUtils<number>::onetensor(n+2);
  e = ArrayUtils<number>::onetensor(n+2);
  v = ArrayUtils<number>::twotensor(n+2,n+2);


  d[1] = 0.0;
  e[1] = 1.0/sqrt(3);
  for(i=2;i<=n;++i)  {
    norm = static_cast<number>(i-1);
    e[i] = sqrt((2*norm-1)/(2*norm+1))*norm/(2*norm-1);
    d[i] = 0.0;
  }
  d[n+1] = 0.0;
  norm = static_cast<number>(n);
  e[n+1] = sqrt(norm/(2*norm-1));

  tqli(d,e,n+1,n+2,v);

  for(i=0;i<=n;++i) {
    norm = 0.0;
    for(j=1;j<=n+1;++j)
      norm += v[j][i+1]*v[j][i+1];
    w[i] = v[1][i+1];
    w[i] = w[i]*w[i]/norm*2.0;
    x[i] = d[i+1];
  }

  ArrayUtils<number>::delonetensor(d);
  ArrayUtils<number>::delonetensor(e);
  ArrayUtils<number>::deltwotensor(v);

  for(i=0;i<n;++i) {
    p = x[k=i];
    for(j=i+1;j<=n;++j)
      if(x[j]>=p) p=x[k=j];
    if (k!=i) {
      x[k] = x[i];
      x[i] = p;
      p = w[k];
      w[k] = w[i];
      w[i] = p;
    }
  }

}


template <class number>
void Legendre<number>::leg_quad_Gauss(number *x,number *w,int n) {
  /* Routine to find the Legendre-Gauus quadrature.
     Uses Golub's method, SIAM Review, Vol. 15, No.2 April, 73,
     Page 318.
  */

  number *d,*e,**v;
  int i,j,k;
  number norm,p;

  d = ArrayUtils<number>::onetensor(n+2);
  e = ArrayUtils<number>::onetensor(n+2);
  v = ArrayUtils<number>::twotensor(n+2,n+2);


  d[1] = 0.0;
  e[1] = 1.0/sqrt(3.0);
  for(i=2;i<=n+1;++i)  {
    norm = static_cast<number>(i-1);
    e[i] = sqrt((2*norm-1)/(2*norm+1))*norm/(2*norm-1);
    d[i] = 0.0;
  }

  tqli(d,e,n+1,n+2,v);

  for(i=0;i<=n;++i) {
    norm = 0.0;
    for(j=1;j<=n+1;++j)
      norm += v[j][i+1]*v[j][i+1];
    w[i] = v[1][i+1];
    w[i] = w[i]*w[i]/norm*2.0;
    x[i] = d[i+1];
  }

  ArrayUtils<number>::delonetensor(d);
  ArrayUtils<number>::delonetensor(e);
  ArrayUtils<number>::deltwotensor(v);

  for(i=0;i<n;++i) {
    p = x[k=i];
    for(j=i+1;j<=n;++j)
      if(x[j]>=p) p=x[k=j];
    if (k!=i) {
      x[k] = x[i];
      x[i] = p;
      p = w[k];
      w[k] = w[i];
      w[i] = p;
    }
  }

}



// Values of the Legendre polynomials at the point defined by *x
// Returns a matrix with all the polynomials up to a given
// degree:
//      interp[i][j] = L_i(x_j)
//
template <class number>
void Legendre<number>::leg_val(number **interp,number *x,int n,int rows) {
  number x1,x2;
  int i,j,k;
  number l;


  for(i=0;i<=n;++i) {
    interp[0][i] = 1.0;
    interp[1][i] = x[i];
  }

  for(k=2;k<=rows;++k) {
    l = static_cast<number>(k);
    x1 = (2*l-1)/l;
    x2 = (l-1)/l;
    for(j=0;j<=n;++j)
      interp[k][j] = (x1*x[j]*interp[k-1][j] -
                                            x2*interp[k-2][j]);
  }



}


// Values of the Legendre polynomials at the point defined by *x
// for the degree given by degree.
// Returns only the polynomials to the indicated degree:
//    interp[i] = L_num(x_i)
//
template <class number>
void Legendre<number>::leg_val(number *interp,number *x,int num,int degree) {
  number *l1 = NULL,
    *l2 = NULL;

  number *tmp;

  number x1,x2;
  int i,j,k;
  number l;

  l1 = ArrayUtils<number>::onetensor(num+1);
  l2 = ArrayUtils<number>::onetensor(num+1);

  for(i=0;i<=num;++i) {
    l1[i] = 1.0;
    l2[i] = x[i];
  }

  for(k=2;k<=degree;++k) {
    l = static_cast<number>(k);
    x1 = (2.0*l-1.0)/l;
    x2 = (l-1.0)/l;
    for(j=0;j<=num;++j)
      l1[j] = x1*x[j]*l2[j] - x2*l1[j];

    tmp = l2;
    l2 = l1;
    l1 = tmp;

  }

  for(j=0;j<=num;++j)
    interp[j] = l2[j];

  ArrayUtils<number>::delonetensor(l1);
  ArrayUtils<number>::delonetensor(l2);

}



// Values of the derivatives of the Legendre polynomials at the
// specified points (x) for the given degree (degree)
//      interp[i] = L'_degree(x_i)
//         i = 0...num
template <class number>
void Legendre<number>::leg_der_val(number *interp,
                                                                     number *x,
                                                                     int num,
                                                                     int degree) {

  number *l1 = NULL,
    *l2 = NULL,
    *l3 = NULL,
    *l4 = NULL;

  number *tmp;

  number x1,x2,x3;
  int i,j,k;
  number l;

  l1 = ArrayUtils<number>::onetensor(num+1);
  l2 = ArrayUtils<number>::onetensor(num+1);
  l3 = ArrayUtils<number>::onetensor(num+1);
  l4 = ArrayUtils<number>::onetensor(num+1);


  for(i=0;i<=num;++i) {
    l1[i] = 1.0;
    l2[i] = x[i];
    l3[i] = 0.0;
    l4[i] = 1.0;
  }

  for(k=2;k<=degree;++k) {
    l = static_cast<number>(k);
    x3 = (2.0*l-1.0);
    x1 = x3/l;
    x2 = (l-1.0)/l;
    for(j=0;j<=num;++j) {
      l1[j] = x1*x[j]*l2[j] - x2*l1[j];
      l3[j] = l3[j] + x3*l2[j];
    }
    tmp = l2;
    l2 = l1;
    l1 = tmp;

    tmp = l4;
    l4 = l3;
    l3 = tmp;

  }

  for(j=0;j<=num;++j)
    interp[j] = l4[j];

  ArrayUtils<number>::delonetensor(l1);
  ArrayUtils<number>::delonetensor(l2);
  ArrayUtils<number>::delonetensor(l3);
  ArrayUtils<number>::delonetensor(l4);

}


// Values of the derivatives of the Legendre polynomials at the
// specified points (x) for the given degree (degree)
//      interp[i] = L''_degree(x_i)
//         i = 0...num
template <class number>
void Legendre<number>::leg_2_der_val(number *interp,
                                                                         number *x,
                                                                         int num,
                                                                         int degree) {

  number *l1 = NULL,
    *l2 = NULL,
    *l3 = NULL,
    *l4 = NULL,
    *l5 = NULL,
    *l6 = NULL;

  number *tmp;

  number x1,x2,x3;
  int i,j,k;
  number l;

  l1 = ArrayUtils<number>::onetensor(num+1);
  l2 = ArrayUtils<number>::onetensor(num+1);
  l3 = ArrayUtils<number>::onetensor(num+1);
  l4 = ArrayUtils<number>::onetensor(num+1);
  l5 = ArrayUtils<number>::onetensor(num+1);
  l6 = ArrayUtils<number>::onetensor(num+1);


  for(i=0;i<=num;++i) {
    l1[i] = 1.0;
    l2[i] = x[i];
    l3[i] = 0.0;
    l4[i] = 1.0;
    l5[i] = 0.0;
    l6[i] = 0.0;
  }

  for(k=2;k<=degree;++k) {
    l = static_cast<number>(k);
    x3 = (2.0*l-1.0);
    x1 = x3/l;
    x2 = (l-1.0)/l;
    for(j=0;j<=num;++j) {
      l1[j] = x1*x[j]*l2[j] - x2*l1[j];
      l3[j] += x3*l2[j];
      l5[j] += x3*l4[j];
    }
    tmp = l2;
    l2 = l1;
    l1 = tmp;

    tmp = l4;
    l4 = l3;
    l3 = tmp;

    tmp = l6;
    l6 = l5;
    l5 = tmp;

  }

  for(j=0;j<=num;++j)
    interp[j] = l6[j];

  ArrayUtils<number>::delonetensor(l1);
  ArrayUtils<number>::delonetensor(l2);
  ArrayUtils<number>::delonetensor(l3);
  ArrayUtils<number>::delonetensor(l4);
  ArrayUtils<number>::delonetensor(l5);
  ArrayUtils<number>::delonetensor(l6);

}


// Legendre (Gauss-Lobatto) collocation derivative matrix
template <class number>
void Legendre<number>::leg_der(number **d1,number **lval,number *x,int n,int) {
  number nx;
  int l,j;

  nx = static_cast<number>(n);
  for(l=0;l<=n;++l)
    for(j=0;j<=n;++j) {
      if (l!=j)
                d1[l][j] = lval[n][l]/(lval[n][j]*(x[l]-x[j]));
      else if((l==0)&&(j==0))
                d1[0][0] = 0.25*nx*(nx+1.0);
      else if ((l==n)&&(j==n))
                d1[n][n] = -0.25*nx*(nx+1.0);
      else
                d1[l][j] = 0.0;
    }
}

// Legendre (Gauss) collocation derivative matrix
template <class number>
void Legendre<number>::leg_gauss_der(number **d1,number *x,int n) {

  int i,j;
  number *p1;
  number *p2;

  p1 = ArrayUtils<number>::onetensor(n+1);
  p2 = ArrayUtils<number>::onetensor(n+1);

  leg_2_der_val(p2,x,n,n+1);
  leg_der_val(p1,x,n,n+1);

  for(i=0;i<=n;++i) {
    for(j=0;j<i;++j)
      d1[i][j] = p1[i]/((x[i]-x[j])*p1[j]);
    d1[i][i] = p2[i]*0.5/p1[j];
    for(j=i+1;j<=n;++j)
      d1[i][j] = p1[i]/((x[i]-x[j])*p1[j]);
  }

  ArrayUtils<number>::delonetensor(p1);
  ArrayUtils<number>::delonetensor(p2);

}

// Routine to calculate the entries in the stiffness
// matrix for a Legendre-collocation type of method.
template <class number>
void Legendre<number>::stiffLeg(number **stiff,number *w,number **D1,int N){

  int i,j,m;

  for (m=0;m<=N;++m)
    for(i=0;i<=N;++i) {

      stiff[m][i] = 0.0;
      for(j=0;j<=N;++j)
        stiff[m][i] -= D1[j][m]*D1[j][i]*w[j];

        }


}


// Values of the Lagrange interpolates (Gauss) at the endpoints
template <class number>
void Legendre<number>::gauss_collocation_end_values(number *left,number *right,
                                                                                                        number *x,int n){
  int i;
  number *p1;
  p1 = ArrayUtils<number>::onetensor(n+1);
  leg_der_val(p1,x,n,n+1);

  for(i=0;i<=n;++i) {
    right[i] = 1.0/( ( 1-x[i])*p1[i]);
    left[i]  = 1.0/( (-1-x[i])*p1[i]) * ((n%2==1) ? 1.0 : -1.0);
    // cout << right[i] << "  " << left[i] << endl;
  }

  ArrayUtils<number>::delonetensor(p1);

}



// Petrov-Galerkin mass and stiffness matrices
template <class number>
void Legendre<number>::petrov(number **mass,number **stiff,number *x,number *w,
                                                            number **lval,int n,int) {
  int i,j,k;
  number *gamma;
  number *diag,*off;
  number x1,x2,x3;

  gamma = ArrayUtils<number>::onetensor(n+1);
  diag  = ArrayUtils<number>::onetensor(n+1);
  off   = ArrayUtils<number>::onetensor(n+1);

  gamma[n] = 0.5*static_cast<number>(n);
  for(i=0;i<n;++i)
    gamma[i] = ((static_cast<number>(2*i))+1)*0.5;

  diag[0] = 4.0/3;
  diag[1] = 4.0/15.0;
  off[0] = -4.0/15.0;
  off[1] = -4.0/35;
  for(i=2;i<=n;++i) {
    x1 = static_cast<number>(i);
    x2 = 1.0/(2.0*x1+1.0);
    x3 = x1 + 1.0;
    diag[i] = x2*2.0*(1.0-x3*x3*x2/(2.0*x1+3.0)-x1*x1*x2/(2.0*x1-1.0));
    off[i-2] = - 2.0*x1*x2*(x1-1.0)/((2.0*x1-3.0)*(2.0*x1-1.0));
  }

  for(k=0;k<=n;++k) {
    mass[0][k] = w[k]*(gamma[0] + gamma[1]/3.0*x[k]);
    mass[n][k] = w[k]*(gamma[0] - gamma[1]/3.0*x[k]);
  }

  for(k=0;k<=n;++k) {
    mass[1][k] = w[k]*
      (gamma[0]*diag[0]*lval[0][k] +
       gamma[2]*off[0]*lval[2][k]);
    mass[2][k] = w[k]*
      (gamma[1]*diag[1]*lval[1][k] +
       gamma[3]*off[1]*lval[3][k]);
  }

  for(j=2;j<=n-2;++j)
    for(k=0;k<=n;++k) {
      mass[j+1][k] = w[k]*
                (gamma[j]*diag[j]*lval[j][k] +
                 gamma[j-2]*off[j-2]*lval[j-2][k] +
                 gamma[j+2]*off[j]*lval[j+2][k]);
    }

  diag[0] = 0.0;
  diag[1] = 0.0;
  for(i=2;i<=n;++i) {
    x1 = static_cast<number>(i);
    diag[i] = 	2.0*x1*(1.0-x1)/(2.0*x1+1.0);
  }

  for(i=0;i<=n;++i) {
    stiff[0][i] = 0.0;
    stiff[n][i] = 0.0;
  }
  stiff[0][0] = -0.5;
  stiff[0][n] = 0.5;
  stiff[n][0] = 0.5;
  stiff[n][n] = -0.5;

  for(j=0;j<=n-2;++j)
    for(k=0;k<=n;++k) {
      x2 = gamma[j]*lval[j][k]*diag[j];
      for(i=j+2;i<=n;i+=2)
                x2 += 4.0*gamma[i]*lval[i][k];
      stiff[j+1][k] = w[k]*x2;
    }

  ArrayUtils<number>::delonetensor(gamma);
  ArrayUtils<number>::delonetensor(diag);
  ArrayUtils<number>::delonetensor(off);

}



// Legendre-Tau stiffness and mass matrices
template <class number>
void Legendre<number>::tauMatrices(number **stiff,number **mass,number **der,int num) {

    /*
        number len,number stretchm,
        number stretchd1,number stretchd2 ) {
    */

  // number mulm,muls;
  number tmp0,tmp1;

  int i,j;


  for(i=0;i<=num;++i)
    for(j=0;j<=num;++j) {
      mass[i][j] = 0.0;
      stiff[i][j] = 0.0;
      der[i][j] = 0.0;
        }

  // mulm = len*0.5;
  // muls = 1.0/mulm;
  // stretchm = mulm;
  // stretchd1 = muls;
  // stretchd2 = muls*muls;
  for(i=0;i<num-1;++i) {
    tmp0 = static_cast<number>(i);
    tmp1 = tmp0/(2.0*tmp0+1.0);
    mass[i][i]  = ( 1.0/(tmp0+0.5)
                                        - 2.0/(2.0*tmp0+3.0)*(tmp0+1.0)/(2.0*tmp0+1.0)
                                        *(tmp0+1.0)/(2.0*tmp0+1.0)
                                        - 2.0/(2.0*tmp0-1.0)*tmp1*tmp1 );
    stiff[i][i] = 2.0*tmp1*(1.0-tmp0);
    }

  for(j=0;j<num-1;++j)
    for(i=j+2;i<=num;i+=2)
      stiff[j][i] = 4.0;

  /* row (num-1) is found by integrating against (1+x)/2
     row (num)   is found by integrating against (1-x)/2 */
  for(i=1;i<=num;i+=2) {
    stiff[num-1][i] = -1.0;
    stiff[num][i] = 1.0;
  }

  for(i=0;i<num-1;++i) {
    tmp0 = static_cast<number>(i);
    tmp1 = static_cast<number>(i+2);
    mass[i][i+2] = -tmp1/(2.0*tmp1+1.0)*(tmp0+1.0)
      /(2.0*tmp0+1.0)*2.0/(2.0*tmp1-1.0);
  }

  for(i=2;i<num-1;++i) {
    tmp0 = static_cast<number>(i);
    tmp1 = static_cast<number>(i-2);
    mass[i][i-2] = -(tmp1+1.0)/(2.0*tmp1+1.0)*tmp0
      /(2.0*tmp0+1.0)*2.0/(2.0*tmp1+3.0);
  }

  mass[num-1][0] = 1.0;
  mass[num-1][1] = 1.0/3.0;
  mass[num][0] = 1.0;
  mass[num][1] = -1.0/3.0;

  for(i=0;i<=num-2;++i) {
    tmp0 = static_cast<number>(i);
    der[i][i+1] = 2.0*(tmp0+1.0)*(tmp0+2.0)/
      ( (2.0*tmp0+1.0)*(2.0*tmp0+3.0));
  }

  for(i=1;i<=num-2;++i) {
    tmp0 = static_cast<number>(i);
    der[i][i-1] = -2.0*tmp0*(tmp0-1.0)/
      ( (2.0*tmp0-1.0)*(2.0*tmp0+1.0) );
  }

  der[num-1][0] = 0.0;
  der[num][0] = 0.0;
  tmp0 = -1.0;
  for(i=1;i<=num;++i) {
    der[num-1][i] = 1.0;
    der[num][i] = (tmp0 *= -1.0);
  }

}





/* ******************************************
   Subroutine to initialize the stiffnes,
   mass, and derivative matrices for the
   galerkin method.
     ********************************************* */
template <class number>
void Legendre<number>::galerkinMatrices(number **stiff,number **mass,number **der,int num) {

  int i,j,k;


  for(i=0;i<=num;++i)
    for(j=0;j<=num;++j) {
      mass[i][j] = 0.0;
      stiff[i][j] = 0.0;
      der[i][j] = 0.0;
        }


  stiff[0][0] = -1.0/2.0;
  stiff[1][0] = 1.0/2.0;
  stiff[0][1] = 1.0/2.0;
  stiff[1][1] = -1.0/2.0;


  for(i=2;i<=num;++i) {
    stiff[i][i] = 2.0*(-2.0*static_cast<number>(i)+1);
  }


  mass[0][0] = 2.0/3.0;
  mass[1][0] = 1.0/3.0;
  mass[2][0] = -1.0;
  mass[3][0] = -1.0/3.0;

  mass[0][1] = 1.0/3.0;
  mass[1][1] = 2.0/3.0;
  mass[2][1] = -1.0;
  mass[3][1] = 1.0/3.0;

  mass[0][2] = -1.0;
  mass[1][2] = -1.0;
  mass[2][2] = (2.0/5.0 + 2.0);
  mass[4][2] = -2.0/5.0;

  mass[0][3] = -1.0/3.0;
  mass[1][3] = 1.0/3.0;
  mass[3][3] = 2.0/7.0+2.0/3.0;
  mass[5][3] = -2.0/7.0;

  for(j=4;j<=num;++j) {
    k = j - 2;
    mass[k][j] = mass[k][j] - 2.0/(2.0*static_cast<number>(k)+1);
    if(j<num-1) {
      k = j + 2;
      mass[k][j] = mass[k][j] - 2.0/(2.0*static_cast<number>(j)+1);
    }
    mass[j][j] = mass[j][j] +
      (2.0/(2.0*static_cast<number>(j)+1)+2.0/(2.0*static_cast<number>(j)-3));
  }


  der[0][0] = 0.5;
  der[0][1] = -0.5;
  der[0][2] = 1.0;

  der[1][0] = 0.5;
  der[1][1] = -0.5;
  der[1][2] = -1.0;

  der[2][0] = -1.0;
  der[2][1] = 1.0;
  der[2][3] = 2.0;

  for(j=3;j<=num;++j) {
    der[j][j-1] = -2.0;
    if(j<num)
      der[j][j+1] = 2.0;
  }



}






/* ******************************************
   Subroutine to initialize the stiffnes,
   mass, and derivative matrices for the
   galerkin method. (infinite domain)
     ********************************************* */
template <class number>
void Legendre<number>::galerkinInfinite(number **rightStiff,number **leftStiff,
                                                                                number **mass,number **,
                                                                                int num)
{


  number Lparam = 2.0;
  number ix;
  number lp2;
  //int N1;
  int i;

  //N1 = num + 1;



  mass[0][0] = 2./3.;
  mass[0][1] = 1./3.;
  mass[0][2] = -1.0;
  mass[0][3] = -1./3.;

  mass[1][0] = 1./3.;
  mass[1][1] = 2./3.;
  mass[1][2] = -1.0;
  mass[1][3] = 1./3.;

  mass[2][0] = -1.0;
  mass[2][1] = -1.0;
  mass[2][2] = (2./5. + 2.);
  mass[2][4] = -2./5.;

  mass[3][0] = -1./3.;
  mass[3][1] = 1./3.;
  mass[3][3] = 2./7.+2./3.;
  mass[3][5] = -2./7.;

  for(i=4;i<=num;++i) {
    ix = static_cast<number>(i);
    mass[i][i] = 2.*(1./(2.*ix+1.)+1./(2.*ix-3.));
    mass[i][i-2] = -2./(2.*ix-3.);

    if(ix<num-1) {
      mass[i][i+2] = -2./(2.*ix+1.);
    }

  }



  lp2 = Lparam*Lparam;

  rightStiff[0][0] = -1./5./lp2;
  rightStiff[0][1] = 1./5./lp2;
  rightStiff[0][2] = 6./5./lp2;
  rightStiff[0][3] = -10./7./lp2;
  rightStiff[0][4] = 4./5./lp2;
  rightStiff[0][5] = -6./35./lp2;

  rightStiff[1][0] = 6./5./lp2;
  rightStiff[1][1] = -6./5./lp2;
  rightStiff[1][2] = -24./5./lp2;
  rightStiff[1][3] = 24./7./lp2;
  rightStiff[1][4] = -6./5./lp2;
  rightStiff[1][5] = 6./35./lp2;

  rightStiff[2][0] = 4./5./lp2;
  rightStiff[3][0] = -4./7./lp2;
  rightStiff[4][0] = 1./5./lp2;
  rightStiff[5][0] = -1./35./lp2;

  rightStiff[2][1] = -4./5./lp2;
  rightStiff[3][1] = 4./7./lp2;
  rightStiff[4][1] = -1./5./lp2;
  rightStiff[5][1] = 1./35./lp2;

  rightStiff[2][2] = -192./35./lp2;
  rightStiff[2][3] = 8./lp2;
  rightStiff[2][4] = -92./15./lp2;
  rightStiff[2][5] = 12./5./lp2;
  rightStiff[2][6] = -8./21./lp2;

  rightStiff[3][2] = 36./7./lp2;
  rightStiff[3][3] = -220./21./lp2;
  rightStiff[3][4] = 12./lp2;
  rightStiff[3][5] = -612./77./lp2;
  rightStiff[3][6] = 20./7./lp2;
  rightStiff[3][7] = -100./231./lp2;

  rightStiff[4][2] = -14./5./lp2;
  rightStiff[5][2] = 6./7./lp2;
  rightStiff[6][2] = -4./35./lp2;


  rightStiff[4][3] = 26./3./lp2;
  rightStiff[5][3] = -346./77./lp2;
  rightStiff[6][3] = 4./3./lp2;
  rightStiff[7][3] = -40./231./lp2;


  for(i=4;i<=num;++i) {

    ix = static_cast<number>(i);

    rightStiff[i][i] =
      -(35.*ix*ix*ix*ix-70.*ix*ix*ix-130.*1.*ix*ix+165.*ix+126.)/
      (2.*ix+3.)/(2.*ix+1.)/(2.*ix-5.)/(2.*ix-3.)*(2.*ix-1.)/lp2;

    if(i>4) {
      rightStiff[i][i-1] = (7.*ix*ix-7.*ix-6.)/(2.*ix+1.)/lp2;
      if(i>5)  {
                rightStiff[i][i-2] =
                    -(ix-2.)*(14.*ix*ix*ix-70.*ix*ix+59.*ix+51.)/
                    (2.*ix-7.)/(2.*ix+1.)/(2.*ix-3.)/lp2;
                if(i>6)  {
                    rightStiff[i][i-3] = (ix-2.)*(ix-3.)/(2.*ix-3.)/lp2;
                    if(i>7)
                        rightStiff[i][i-4] = -(ix-2.)*(ix-4.)*(ix-3.)*(ix-3.)
                            /lp2/(2.*ix-7.)/(2.*ix-5.)/(2.*ix-3.)/2.;
                }
      }
    }

    if(i<num)  {
      rightStiff[i][i+1] = (7.*ix*ix-7.*ix-6.)/(2.*ix-3.)/lp2;
      if(i<num-1)  {
                rightStiff[i][i+2] = -(ix+1.)*
                    (14.*ix*ix*ix+28.*ix*ix-39.*ix-54.)/
                    (2.*ix+5.)/(2.*ix-3.)/(2.*ix+1.)/lp2;
                if(i<num-2)  {
                    rightStiff[i][i+3] = (ix+2.)*(ix+1.)/(2.*ix+1.)/lp2;
                    if(i<num-3)
                        rightStiff[i][i+4] = -(ix+3.)*(ix+1.)*(ix+2.)*(ix+2.)/
                            (2.*ix+1.)/(2.*ix+3.)/(2.*ix+5.)/lp2/2.;
                }
      }
    }

  }





  leftStiff[0][0] = -6./5./lp2;
  leftStiff[0][1] = 6./5./lp2;
  leftStiff[0][2] = -24./5./lp2;
  leftStiff[0][3] = -24./7./lp2;
  leftStiff[0][4] = -6./5./lp2;
  leftStiff[0][5] = -6./35./lp2;

  leftStiff[1][0] = 1./5./lp2;
  leftStiff[1][1] = -1./5./lp2;
  leftStiff[1][2] = 6./5./lp2;
  leftStiff[1][3] = 10./7./lp2;
  leftStiff[1][4] = 4./5./lp2;
  leftStiff[1][5] = 6./35./lp2;

  leftStiff[2][0] = -4./5./lp2;
  leftStiff[3][0] = -4./7./lp2;
  leftStiff[4][0] = -1./5./lp2;
  leftStiff[5][0] = -1./35./lp2;

  leftStiff[2][1] = 4./5./lp2;
  leftStiff[3][1] = 4./7./lp2;
  leftStiff[4][1] = 1./5./lp2;
  leftStiff[5][1] = 1./35./lp2;

  leftStiff[2][2] = -192./35./lp2;
  leftStiff[2][3] = -8./lp2;
  leftStiff[2][4] = -92./15./lp2;
  leftStiff[2][5] = -12./5./lp2;
  leftStiff[2][6] = -8./21./lp2;

  leftStiff[3][2] = -36./7./lp2;
  leftStiff[3][3] = -220./21./lp2;
  leftStiff[3][4] = -12./lp2;
  leftStiff[3][5] = -612./77./lp2;
  leftStiff[3][6] = -20./7./lp2;
  leftStiff[3][7] = -100./231./lp2;

  leftStiff[4][2] = -14./5./lp2;
  leftStiff[5][2] = -6./7./lp2;
  leftStiff[6][2] = -4./35./lp2;


  leftStiff[4][3] = -26./3./lp2;
  leftStiff[5][3] = -346./77./lp2;
  leftStiff[6][3] = -4./3./lp2;
  leftStiff[7][3] = -40./231./lp2;


  for(i=4;i<=num;++i) {

    ix = static_cast<number>(i);

    leftStiff[i][i] =
      -(35.*ix*ix*ix*ix-70.*ix*ix*ix-130.*ix*ix+165.*ix+126.)/
      (2.*ix-5.)/(2.*ix-3.)/(2.*ix+3.)/(2.*ix+1.)*(2.*ix-1.)/lp2;

    if(i>4) {
      leftStiff[i][i-1] = -(7.*ix*ix-7.*ix-6.)/(2.*ix+1.)/lp2;

      if(i>5)  {
                leftStiff[i][i-2] = -(ix-2.)*
                    (14.*ix*ix*ix-70.*ix*ix+59.*ix+51.)/
                    (2.*ix-7.)/(2.*ix+1.)/(2.*ix-3.)/lp2;

                if(i>6)  {
                    leftStiff[i][i-3] = -(ix-2.)*(ix-3.)/(2.*ix-3.)/lp2;

                    if(i>7)
                        leftStiff[i][i-4] = -(ix-2.)*(ix-4.)*(ix-3.)*(ix-3.)/
                            lp2/(2.*ix-7.)/(2.*ix-5.)/(2.*ix-3.)/2.;

                }
      }
    }

    if(i<num)  {
      leftStiff[i][i+1] = -(7.*ix*ix-7.*ix-6.)/(2.*ix-3.)/lp2;

      if(i<num-1)  {
                leftStiff[i][i+2] = -(ix+1.)*
                    (14.*ix*ix*ix+28.*ix*ix-39.*ix-54.)/
                    (2.*ix+5.)/(2.*ix-3.)/(2.*ix+1.)/lp2;

                if(i<num-2)  {
                    leftStiff[i][i+3] = -(ix+2.)*(ix+1.)/(2.*ix+1.)/lp2;

                    if(i<num-3)
                        leftStiff[i][i+4] = -(ix+3.)*(ix+1.)*(ix+2.)*(ix+2.)/
                            (2.*ix+1.)/(2.*ix+3.)/(2.*ix+5.)/lp2/2.;

                }
      }
    }

  }







}




#endif // LEGENDRE_H
