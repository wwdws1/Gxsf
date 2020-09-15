//Powered By Walter
//期末试题第一题：用Lagrange插值的结果进行Simpon积分
//include: "Interpolation.hpp", "Integral.hpp"
//infile: "FinalExam1.txt"
#include "Integral.hpp"
#include "Interpolation.hpp"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;

void input(double a[19], double b[19]);
void test(double a[19], double b[19]);
double fx(double a);
double Fx(double a);

int main()
{
    system("chcp 65001");

    double x0[19], y0[19]; //初始数据点
    double y[6001];        //插值数据点--中间值
    double x1[4], y1[4];   //抽取的4组初始数据点
    double sum1 = 0.0;     //累计积分值
    double sum2, dif1;     //测试用

    Interpolation *interpolation = NULL; //定义指向插值对象的指针
    Integral *integral = NULL;           //定义指向积分对象的指针

    test(x0, y0);
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x1[j] = x0[3 * i + j];
            y1[j] = y0[3 * i + j];
        }

        interpolation = new Interpolation(x1, y1);
        interpolation->getXY(y);

        integral = new Integral(y);
        sum1 += integral->getSn();

        delete interpolation;
        delete integral;
    }

    sum2 = Fx(3.0) - Fx(1.2);
    dif1 = abs(sum2 - sum1);
    //cout << sum1 << endl;
    //cout << sum2 << endl;
    //cout << dif1 << endl;
    cout << "测试部分：" << endl;
    cout << "测试函数：f(x) = x^3 + x^2 + x + 1" << endl;
    cout << "计算结果：" << sum1 << endl;
    cout << "解析结果：" << sum2 << endl;
    cout << "误差为：" << dif1 << endl;

    input(x0, y0);
    sum1 = 0.0;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x1[j] = x0[3 * i + j];
            y1[j] = y0[3 * i + j];
        }

        interpolation = new Interpolation(x1, y1);
        interpolation->getXY(y);

        integral = new Integral(y);
        sum1 += integral->getSn();

        delete interpolation;
        delete integral;
    }

    //cout << sum1 << endl;
    cout << "根据题目数据得出的结果为：" << sum1 << endl;

    system("pause");

    return 0;
}

void input(double a[19], double b[19])
{
    ifstream infile;

    infile.open("FinalExam1.txt");
    for (int i = 0; i < 19; i++)
    {
        infile >> a[i];
    }
    for (int i = 0; i < 19; i++)
    {
        infile >> b[i];
    }

    infile.close();

    return;
}

void test(double a[19], double b[19])
{
    double h = 0.1;

    for (int i = 0; i < 19; i++)
    {
        a[i] = 1.2 + h * i;
        b[i] = fx(a[i]);
    }

    return;
}

double fx(double a)
{
    double res;

    res = pow(a, 3) + pow(a, 2) + a + 1;

    return res;
}

double Fx(double a)
{
    double res;

    res = pow(a, 4) / 4.0 + pow(a, 3) / 3.0 + pow(a, 2) / 2.0 + a;

    return res;
}