//Powered By Walter
//期末试题第一题的Simpon积分部分的具体实现
#include <iostream>

using namespace std;

class Integral
{
public:
    double getSn();                //返回积分值
    Integral(double y[6001]);      //构造函数
    Integral(const Integral &obj); //拷贝构造函数
    ~Integral();                   //析构函数

private:
    double _h = 0.0001; //步长
    double _y[6001];    //f(x)的值
    double _sn;         //积分的结果

    double Tn();                              //计算Tn的值
    double Hn();                              //计算Hn的值
    void setSn();                             //计算Sn的值
    void Copy(double a[], double b[], int n); //数组拷贝
};

Integral::Integral(double y[6001])
{
    Copy(y, _y, 6001);
    setSn();
}

Integral::Integral(const Integral &obj)
{
}

double Integral::Tn()
{
    double res = (_y[0] + _y[6001]) / 2.0; //返回值
    for (int i = 1; i < 3000; i++)
    {
        res += _y[2 * i];
    }
    res *= _h;

    return res;
}

double Integral::Hn()
{
    double res = 0.0; //返回值
    for (int i = 0; i < 3000; i++)
    {
        res += _y[2 * i + 1];
    }

    res *= _h;

    return res;
}

void Integral::setSn()
{
    _sn = (1.0 / 3.0) * (Tn() + 2.0 * Hn());

    return;
}

void Integral::Copy(double a[], double b[], int n)
{
    for (int i = 0; i < n; i++)
    {
        b[i] = a[i];
    }

    return;
}

double Integral::getSn()
{
    return _sn;
}

Integral::~Integral()
{
}