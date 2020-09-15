//Powered By Walter
//四阶Runge-Kutta方法数值解二阶微分方程程序的主函数
//include: "RungeKutta.hpp"
//outfile: "homework7.txt"

#include "RungeKutta.hpp"
#include <fstream>
#include <iostream>
#include <math.h>

#define count 101 //定义循环次数

using namespace std;

double fx(double a);                                      //二阶微分方程的解析计算
void output(double a[101], double b[101], double c[101]); //输出结果到文件

int main()
{
    double x[101];  //自变量
    double y[101];  //解析结果
    double y1[101]; //计算结果
    double z1[101]; //一阶计算结果
    RungeKutta *rungekutta = NULL;

    y1[0] = 4;  //初值
    z1[0] = -2; //一阶初值
    x[0] = 0;   //自变量取值[0,1],步长0.01
    for (int i = 1; i < count; i++)
    {
        x[i] = x[i - 1] + 0.01;
    }

    for (int i = 0; i < count; i++)
    {
        y[i] = fx(x[i]);
    }

    for (int i = 1; i < count; i++)
    {
        rungekutta = new RungeKutta(y1[i - 1], z1[i - 1]);

        y1[i] = rungekutta->getYn1();
        z1[i] = rungekutta->getZn1();

        delete rungekutta;
    }

    output(x, y, y1);

    system("pause");

    return 0;
}

double fx(double a)
{
    double res;
    res = exp(-a) * (4 + 2 * a);
    return res;
}

void output(double a[101], double b[101], double c[101])
{
    ofstream outfile;

    outfile.open("homework7.txt");
    for (int i = 0; i < count; i++)
    {
        outfile << a[i] << " ";
    }
    outfile << endl;
    for (int i = 0; i < count; i++)
    {
        outfile << b[i] << " ";
    }
    outfile << endl;
    for (int i = 0; i < count; i++)
    {
        outfile << c[i] << " ";
    }
    outfile.close();

    cout << "All done." << endl;

    return;
}