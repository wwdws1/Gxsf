//Powered by Walter
//在与此文件同目录下新建“a1.txt”和“a2.txt”，以便承接结果

#include<fstream>
#include<iostream>

using namespace std;

static double x0[7] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};
static double y0[7];

void fx();
void setRes(double a[]);
void cal(double a[], double b[]);
double lag(double *add);

int main()
{
   ofstream outfile;
   
   double x1[800];
   double y1[800];

   fx();
   setRes(x1);
   cal(x1,y1);

   outfile.open("a1.txt");
   for(int i = 0; i < 800; i++)
   {
      outfile << x1[i] << ",";
   }
   outfile << endl;
   for(int i = 0; i < 800; i++)
   {
      outfile << y1[i] << ",";
   }
   outfile << endl;
   outfile.close();

   outfile.open("a2.txt");
   for (int i = 0; i < 7; i++)
   {
       outfile << x0[i] << ",";
   }
   outfile << endl;
   for (int i = 0; i < 7; i++)
   {
       outfile << y0[i] << ",";
   }
   outfile << endl;
   outfile.close();

   return 0;
}

void fx()
{
   for(int i = 0; i < 7; i++)
   {
      y0[i] = 30 - x0[i] + 15*x0[i]*x0[i] - 6*x0[i]*x0[i]*x0[i] - 0.5*x0[i]*x0[i]*x0[i]*x0[i] +5*x0[i]*x0[i]*x0[i]*x0[i]*x0[i] - 1.5*x0[i]*x0[i]*x0[i]*x0[i]*x0[i]*x0[i];
   }

   return;
};

void setRes(double a[])
{
   for(int i = 0; i < 800; i++)
   {
      a[i] = 0.01 * i;
   }
   return;
};

void cal(double a[], double b[])
{
   for(int i = 0; i < 800; i++)
   {
      b[i] = lag(&a[i]);
   }
   return;
};

double lag(double *add)
{
   double sg = 0.0;
   double pro1, pro2;
   for(int i = 0; i < 7; i++)
   {
      pro1 = 1.0;
      pro2 = 1.0;
      for(int j = 0; j < 7; j++)
      {
         if(j == i)
         continue;
         pro1 *= (*add - x0[j]);
         pro2 *= (x0[i] - x0[j]);
      }
      sg += y0[i] * pro1 / pro2;
   }
   return sg;
};