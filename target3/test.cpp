#include <stdlib.h>
#include <stdio.h>
int main () {
  int a[2][3];

  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;
  a[1][0] = 4;
  a[1][1] = 5;
  a[1][2] = 6;

  int *p = a[1];
  int *p2 = a[1];
  int (*t) = NULL;

  *a[1]=*a[0];
  *a[0]=*p;
  
  
  //a[1] = p;

  //p=p2;
  //p2=t;

  printf("%i\n", a[0]);
  
  for (int i=0; i<2; i++)
    for (int j=0; j<3; j++)
      printf("%i\n", a[i][j]);
  return 0;
      
  
}


