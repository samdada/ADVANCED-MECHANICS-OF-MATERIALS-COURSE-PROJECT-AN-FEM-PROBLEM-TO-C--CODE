
#include <stdio.h>
#include<math.h>
void thetareturn(void);
double matrixreturn(void);
double inversematrix(void);
double dispcal(void);
double reactioncal(void);
double stressandstrains(void);
int n,t;
 double l[25] ,m[25],X[25],Y[25],le[25],E,A;
 int ni[25] ,nj[25];
 double C1,C2,k[25][25],stress[25][25];


 double K[25][25];
 double inv[25][25];
 double disp[25][25],tdisp[25][25],reactions[25][25];
int main()
{
    double d;
    thetareturn();
    matrixreturn();
    inversematrix();
    dispcal();
    reactioncal();
    stressandstrains();
  return 0;
}
void thetareturn()
{
    printf("enter the no. of nodes\n");
    scanf("%d",&n);

    printf("enter the no. of elements\n");
  scanf("%d",&t);
  printf("ln element connectivity table\n");
  for(int i=0; i<t; i++)
  {
  printf("enter the nodes of element %d\n",i+1);
  scanf("%d\t%d",&ni[i],&nj[i]);
  }
  printf("element node i node j\n");
  for(int i=0; i<t; i++)
  {
  printf("%d\t %d\t %d\n",i+1,ni[i],nj[i]);
  }
  printf("\n enter the coordinates of nodes\n");
  for(int i=0;i<n;i++)
  {
  printf("enter the X & Y coordinates of %d\n",i+1);
  scanf("%lf\t%lf",&X[i],&Y[i]);
  }

  printf("Enter the modulus of elasticity in N/mm2\n");
  scanf("%lf",&E);
  printf("Enter the value of area mm2\n");
  scanf("%lf",&A);
  for(int i=0;i<t;i++)
  {
  le[i]=sqrt((pow((X[nj[i]-1]-X[ni[i]-1]),2)+(pow((Y[nj[i]-1]-Y[ni[i]-1]),2))));
  l[i]=(X[nj[i]-1]-X[ni[i]-1])/le[i];
  m[i]=(Y[nj[i]-1]-Y[ni[i]-1])/le[i];
  }
}
 double matrixreturn(void)
 {
     double l1,l2,m1,m2;
     double a,b,c,d,e,f;


     l1=l[0];
     l2=l[1];
     m1=m[0];
     m2=m[1];
     a=l1*l1;
     b=m1*m1;
     c=m1*l1;
     d=l2*l2;
     e=m2*m2;
     f=m2*l2;
     C1=((A*E)/le[0]);
     C2=((A*E)/le[1]);
     double K1[6][6]={

         {a,c,-a,-c,0,0},
         {c,b,-c,-b,0,0},
         {-a,-c,a,c,0,0},
         {-c,-b,c,b,0,0},
         {0,0,0,0,0,0},
         {0,0,0,0,0,0}
     };
    double K2[6][6]=
     {
         {0,0,0,0,0,0},
         {0,0,0,0,0,0},
         {0,0,d,f,-d,-f},
         {0,0,f,e,-f,-e},
         {0,0,-d,-f,d,f},
         {0,0,-f,-e,f,e}

     };
       printf("\n\nThe 1st element stiffness matrix is\n");
       for(int i=0;i<4;i++)
      {
        for(int j=0;j<4;j++)
      {
          printf("%lf\t",C1*K1[i][j]);
      }
      printf("\n");
      }
       printf("\n\n The 2st element stiffness matrix is\n");
       for(int i=2;i<6;i++)
      {
        for(int j=2;j<6;j++)
      {
          printf("%lf\t",C2*K2[i][j]);
      }
      printf("\n");
      }
       printf("\n\n");



     for(int i=0;i<6;i++)
      {

        for(int j=0;j<6;j++)
         {
              K[i][j]=C1*K1[i][j]+C2*K2[i][j];
         }
      }
      printf("The global stiffness matrix is\n");
      for(int i=0;i<6;i++)
      {
        for(int j=0;j<6;j++)
      {
          printf("%lf\t",K[i][j]);
      }
      printf("\n");
      }
      printf("\n\n");
      printf("The reduced stiffness matrix is\n");
      for(int i=2;i<=3;i++)
      {
          for(int j=2;j<=3;j++)
          {
              printf("%f\t",K[i][j]);
          }
          printf("\n");
      }
      for(int i=0;i<2;i++)
      {
          for(int j=0;j<2;j++)
          {
              k[i][j]=K[i+2][j+2];

          }
      }
 }
double inversematrix(void)
{
    //determinent of reduced stiffness matrix
    double d=(k[0][0]*k[1][1])-(k[0][1]*k[1][0]);

    //inverse of matrix
     double b[25][25],fac[25][25];
   fac[0][0]=k[1][1];
   fac[0][1]=-1*k[1][0];
   fac[1][0]=-1*k[0][1];
   fac[1][1]=k[0][0];


      //for transpose
      for(int i=0;i<2;i++)
      {
          for(int j=0;j<2;j++)
          {
              inv[i][j]=fac[j][i];
          }
      }
      printf("\n\n The inverse matrix is\n\n ");
       for(int i=0;i<=1;i++)
      {
          for(int j=0;j<=1;j++)
          { inv[i][j]=inv[i][j]/d;
              printf("%lf\t",inv[i][j]);
          }
          printf("\n");
      }
}
double dispcal(void)
{
    double f[25];
    printf("\n Enter the value of reduced force vector\n");
    for(int i=0;i<2;i++)
      {
              scanf("%lf",&f[i]);
      }

    //matrix multiplication
    for(int k=0;k<1;k++)
    {
            for(int i=0;i<2;i++)
            { disp[i][k]=0;
                for(int j=0;j<2;j++)
                {
                    disp[i][k]=disp[i][k]+inv[i][j]*f[j];
                }
         }
    }
     printf("\n\n The displacements are\n\n ");
       for(int i=0;i<=1;i++)
      {
          for(int j=0;j<1;j++)
          {
              printf("%lf\t",disp[i][j]);
          }
          printf("\n");
      }
}
double reactioncal()
{
    tdisp[0][0]=0;
    tdisp[1][0]=0;
    tdisp[2][0]=disp[0][0];
    tdisp[3][0]=disp[1][0];
    tdisp[4][0]=0;
    tdisp[5][0]=0;
    printf("\n\n the total displacements are\n");

       for(int k=0;k<1;k++)
    {
            for(int i=0;i<6;i++)
            { reactions[i][k]=0;
                for(int j=0;j<6;j++)
                {
                    reactions[i][k]=reactions[i][k]+K[i][j]*tdisp[j][k];

                }
         }
    }

      printf("\n\n The reactions are\n\n ");
       for(int i=0;i<6;i++)
      {
          for(int j=0;j<1;j++)
          {
              printf("R%d = %lf",i+1,reactions[i][j]);
          }

          printf("\n");
      }

}
double stressandstrains()
{


      double cs[2][4]={{-l[0],-m[0],l[0],m[0]},
                         {-l[1],-m[1],l[1],m[1]} };

    for(int k=0;k<1;k++)
    {
            for(int i=0;i<2;i++)
            { stress[i][k]=0;
                for(int j=0;j<4;j++)
                {
                    stress[i][k]=stress[i][k]+(E/le[i])*(cs[i][j]*tdisp[j+2*i][k]);
                }
            }
    }
    for(int r=0;r<t;r++)
    {
    printf("\n The stress in element %d is %lf",r+1,stress[r][0]);
    }

    for(int i=0;i<t;i++)
    {
        for(int j=0;j<1;j++)
        {
            printf("\n The strain in element %d is %lf",i+1,stress[i][j]/E);
        }
    }
}
