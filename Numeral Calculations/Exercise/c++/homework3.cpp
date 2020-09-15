//Powered by Walter
//原始数据存放于“data1.txt”，输出文件“ax.txt”仅用于Matlab作带误差线的折线图用
//定稿日期：2019年10月27日
#include <fstream>
#include <iostream>
#include <math.h>

#define delta 0.01

using namespace std;

static double x[30], y[30], deltay[30];

void getInitialData(void);
double Pxa(double a[3], double i);
double Qa(double a[3]);
void Grad(double a[3]);
void output(double a[3]);
void exchange(double a[3], double b[3]);

int main()
{
    double a[3] = {-6.0, 0.0, 0.0};

    getInitialData();

    Grad(a);

    output(a);

    system("pause");

    return 0;
}

void getInitialData(void)
{
    ifstream infile;
    ofstream outfile;
    double temp[90];
    int j = 0;
    infile.open("data1.txt");
    for (int i = 0; i < 90; i++)
    {
        infile >> temp[i];
        j = i / 3;
        if (i % 3 == 0)
        {
            x[j] = temp[i];
        }
        else if (i % 3 == 1)
        {
            y[j] = temp[i];
        }
        else
        {
            deltay[j] = temp[i];
        }
    }
    infile.close();

    outfile.open("ax.txt");
    for (int i = 0; i < 30; i++)
    {
        outfile << x[i] << ",";
    }
    outfile << endl;
    for (int i = 0; i < 30; i++)
    {
        outfile << y[i] << ",";
    }
    outfile << endl;
    for (int i = 0; i < 30; i++)
    {
        outfile << deltay[i] << ",";
    }
    outfile.close();

    return;
}

double Pxa(double a[3], int i)
{
    double res;

    res = (1 + a[0] * x[i] + a[1] * pow(x[i], 2)) * exp(-1 * a[2] * pow(x[i], 2));

    return res;
}

double Qa(double a[3])
{
    double res = 0;

    for (int i = 0; i < 30; i++)
    {
        res += pow((y[i] - Pxa(a, i)), 2) / pow(deltay[i], 2);
    }

    return res;
}

void Grad(double a[3])
{
    double res0, res1;
    double a0[3], b[3];
    res0 = Qa(a);
    exchange(a, a0);
    for (b[0] = a[0]; b[0] < 6; b[0] += delta)
    {
        for (b[1] = a[1]; b[1] < 6.0; b[1] += delta)
        {
            for (b[2] = a[2]; b[2] < 5.0; b[2] += delta)
            {
                res1 = Qa(b);
                if (res1 < res0)
                {
                    res0 = res1;
                    exchange(b, a0);
                    //cout << res0 << endl;
                }
            }
        }
        //cout << a0[0] << " , " << a0[1] << " , " << a0[2] << endl;
    }
    exchange(a0, a);
    return;
}

void output(double a[3])
{
    cout << Qa(a) << endl;
    for (int i = 0; i < 3; i++)
    {
        cout << "a" << i << " is : " << a[i] << endl;
    }
}

void exchange(double a[3], double b[3])
{
    for (int i = 0; i < 3; i++)
    {
        b[i] = a[i];
    }

    return;
}