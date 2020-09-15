//Powered By Walter
//期末试题第二题：求解二阶常微分方程
//include: "FinalExam2.hpp"
//outfile: "FinalExam2-1.txt","FinalExam2-2.txt"
//借用matlab绘制图像："FinalExam2.m"
#include "FinalExam2.hpp"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;

void output1(double a[251], double b[251], double c[251]); //输出测试结果到文件
void output2(double a[251], double b[251]);                //输出题目结果到文件
double fx(double a);                                       //测试用二阶微分方程的解析计算
double maxdif(double a[251], double b[251]);               //求最大值

int main()
{
    system("chcp 65001");
    double x[251];  //自变量
    double y[251];  //解析结果
    double y1[251]; //计算结果
    double z1[251]; //一阶计算结果
    double maxdiff; //最大误差
    RungeKutta *rungekutta = NULL;

    //测试
    y1[0] = 4;  //初值
    z1[0] = -2; //一阶初值
    x[0] = 0;   //自变量取值[0,5],步长0.02
    for (int i = 1; i < 251; i++)
    {
        x[i] = x[i - 1] + 0.02;
    }
    for (int i = 0; i < 251; i++)
    {
        y[i] = fx(x[i]);
    }
    for (int i = 1; i < 251; i++)
    {
        rungekutta = new RungeKutta(x[i - 1], y1[i - 1], z1[i - 1], 1);

        y1[i] = rungekutta->getYn1();
        z1[i] = rungekutta->getZn1();

        delete rungekutta;
    }
    maxdiff = maxdif(y, y1);
    output1(x, y, y1);
    cout << "测试阶段中的最大误差为：" << maxdiff << endl;

    //计算
    y1[0] = 0.99752738;    //初值
    z1[0] = -0.0073995279; //一阶初值
    for (int i = 1; i < 251; i++)
    {
        rungekutta = new RungeKutta(x[i - 1], y1[i - 1], z1[i - 1], 2);

        y1[i] = rungekutta->getYn1();
        z1[i] = rungekutta->getZn1();

        delete rungekutta;
    }
    output2(x, y1);

    system("pause");

    return 0;
}

void output1(double a[251], double b[251], double c[251])
{
    ofstream outfile;

    outfile.open("FinalExam2-1.txt");
    for (int i = 0; i < 251; i++)
    {
        outfile << a[i] << " ";
    }
    outfile << endl;
    for (int i = 0; i < 251; i++)
    {
        outfile << b[i] << " ";
    }
    outfile << endl;
    for (int i = 0; i < 251; i++)
    {
        outfile << c[i] << " ";
    }
    outfile.close();

    cout << "测试结束，下面进行题目计算。" << endl;

    return;
}

void output2(double a[251], double b[251])
{
    ofstream outfile;

    outfile.open("FinalExam2-2.txt");
    for (int i = 0; i < 251; i++)
    {
        outfile << a[i] << " ";
    }
    outfile << endl;
    for (int i = 0; i < 251; i++)
    {
        outfile << b[i] << " ";
    }
    outfile.close();

    cout << "题目计算结束。" << endl;

    return;
}

double fx(double a)
{
    double res;
    res = exp(-a) * (4 + 2 * a);
    return res;
}

double maxdif(double a[251], double b[251])
{
    double res;
    double dif[251];

    for (int i = 0; i < 251; i++)
    {
        dif[i] = abs(a[i] - b[i]);
    }
    res = dif[0];
    for (int i = 1; i < 251; i++)
    {
        if (dif[i] > res)
        {
            res = dif[i];
        }
    }

    return res;
}