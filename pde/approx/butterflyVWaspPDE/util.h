#ifndef UTIL_H
#define UTIL_H

#include<stdio.h>
#include<stdlib.h>

template <class number>
class ArrayUtils
{

public:
    ArrayUtils(){}
    static number *****fivetensor(int n1,int n2,int n3,int n4,int n5);
    static number ****fourtensor(int n1,int n2,int n3,int n4);
    static number ***threetensor(int n1,int n2,int n3);
    static number **twotensor(int n1,int n2);
    static number *onetensor(int n1);

    static void delfivetensor(number *****u);
    static void delfourtensor(number ****u);
    static void delthreetensor(number ***u);
    static void deltwotensor(number **u);
    static void delonetensor(number *u);

};



using namespace std;


template <class number>
number *****ArrayUtils<number>::fivetensor(int n1,int n2,int n3,int n4,int n5) {
  number *****u;
  int s,i,j,k,l;

  u = new number****[n1];

  if(u==NULL) {
    printf("Error - fivetensor. Could not allocate memory\n");
    exit(2);
  }

  u[0] = new number***[n1*n2];
  if(u[0]==NULL) {
    printf("Error - fivetensor. Could not allocate memory for *** vector.\n");
    exit(2);
   }

  u[0][0] = new number**[n1*n2*n3];
  if(u[0][0]==NULL) {
    printf("Error - fivetensor. Could not allocate memory for ** vector.\n");
    exit(2);
  }

  u[0][0][0] = new number*[n1*n2*n3*n4];
  if(u[0][0][0]==NULL) {
    printf("Error - fivetensor. Could not allocate memory for * vector.\n");
    exit(2);
  }

  u[0][0][0][0] = new number[n1*n2*n3*n4*n5];
  if(u[0][0][0][0]==NULL) {
    printf("Error - fivetensor. Could not allocate memory for * vector.\n");
    exit(2);
  }


  for(s=0;s<n1;++s) {
      u[s] = u[0] + s*n2;
      for(i=0;i<n2;++i) {
          u[s][i] = u[0][0] + n3*n2*s + n3*i;
          for(j=0;j<n3;++j) {
              u[s][i][j] = u[0][0][0] + n4*n3*n2*s + n4*n3*i + n4*j;
              for(k=0;k<n4;++k) {
                  u[s][i][j][k] = u[0][0][0][0] + n5*n4*n3*n2*s + n5*n4*n3*i + n5*n4*j + n5*k;
                  for(l=0;l<n5;++l)
                      u[s][i][j][k][l] = 0.0;
              }
          }
      }
  }

  return(u);

}


template <class number>
number ****ArrayUtils<number>::fourtensor(int n1,int n2,int n3,int n4) {
  number ****u;
  int s,i,j,k;

  u = new number***[n1];

  if(u==NULL) {
    printf("Error - fourtensor. Could not allocate memory\n");
    exit(2);
  }

  u[0] = new number**[n1*n2];
  if(u[0]==NULL) {
    printf("Error - fourtensor. Could not allocate memory for *** vector.\n");
    exit(2);
   }

  u[0][0] = new number*[n1*n2*n3];
  if(u[0][0]==NULL) {
    printf("Error - fourtensor. Could not allocate memory for ** vector.\n");
    exit(2);
  }

  u[0][0][0] = new number[n1*n2*n3*n4];
  if(u[0][0][0]==NULL) {
    printf("Error - fourtensor. Could not allocate memory for * vector.\n");
    exit(2);
  }

  for(s=0;s<n1;++s) {
      u[s] = u[0] + s*n2;
      for(i=0;i<n2;++i) {
      u[s][i] = u[0][0] + n3*n2*s + n3*i;
      for(j=0;j<n3;++j) {
          u[s][i][j] = u[0][0][0] + n4*n3*n2*s + n4*n3*i + n4*j;
          for(k=0;k<n4;++k)
          u[s][i][j][k] = 0.0;
      }
      }
  }

  return(u);

}


template <class number>
number ***ArrayUtils<number>::threetensor(int n1,int n2,int n3) {
  number ***u;
  int s,i,j;

  u = new number**[n1];

  if(u==NULL) {
    printf("Error - threetensor. Could not allocate memory.\n");
    exit(2);
  }

  u[0] = new number*[n1*n2];
  if(u[0]==NULL) {
    printf("Error - threetensor. Could not allocate memory for ** vector.\n");
    exit(2);
  }

  u[0][0] = new number[n1*n2*n3];
  if(u[0][0]==NULL) {
    printf("Error - threetensor. Could not allocate memory for * vector.\n");
    exit(2);
  }

  for(s=0;s<n1;++s) {
      u[s] = u[0] + s*n2;
      for(i=0;i<n2;++i) {
      u[s][i] = u[0][0] + n3*n2*s + n3*i;
      for(j=0;j<n3;++j)
          u[s][i][j] = 0.0;
      }
  }

  return(u);

}



template <class number>
number **ArrayUtils<number>::twotensor(int n1,int n2){
  number **u;
  int s,i;

  u = new number*[n1];


  if(u==NULL) {
    printf("Error - twotensor. Could not allocate memory.\n");
    exit(2);
  }

  u[0] = new number[n1*n2];
  if(u[0]==NULL) {
    printf("Error - twotensor. Could not allocate memory for vector.\n");
    exit(2);
  }

  for(s=0;s<n1;++s) {
    u[s] = u[0] + s*n2;

    for(i=0;i<n2;++i)
      u[s][i] = 0.0;

  }

  return(u);

}



template <class number>
number *ArrayUtils<number>::onetensor(int n1) {
  number *u;
  int i;

  u = new number[n1];

  if (u==NULL) {
    printf("Error - onetensor. Could not allocate memory.\n");
    exit(2);
  }

  for(i=0;i<n1;++i)
    u[i] = 0.0;

  return(u);

}


template <class number>
void ArrayUtils<number>::delfivetensor(number *****u) {

    if(u==NULL)
    return;

    delete u[0][0][0][0];
    delete u[0][0][0];
    delete u[0][0];
    delete u[0];
    delete [] u;

}


template <class number>
void ArrayUtils<number>::delfourtensor(number ****u) {

    if(u==NULL)
    return;

    delete u[0][0][0];
    delete u[0][0];
    delete u[0];
    delete [] u;

}



template <class number>
void ArrayUtils<number>::delthreetensor(number ***u) {

    if(u==NULL)
    return;

  delete u[0][0];
  delete u[0];
  delete [] u;

}



template <class number>
void ArrayUtils<number>::deltwotensor(number **u){


    if(u==NULL)
    return;

  delete u[0];
  delete [] u;


}


template <class number>
void ArrayUtils<number>::delonetensor(number *u) {

    if(u==NULL)
    return;

  delete [] u;

}



#endif // UTIL_H
