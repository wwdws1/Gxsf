//Powered By Walter
//期末试题第一题的Lagrange插值部分的具体实现
#include <iostream>

using namespace std;

class Interpolation
{
public:
    void getXY(double y[6001]);                //返回插值出来的值
    Interpolation(double x0[4], double y0[4]); //构造函数
    Interpolation(const Interpolation &obj);   //拷贝构造函数
    ~Interpolation();                          //析构函数

private:
    double _h = 0.00005; //步长
    double _x0[4];       //数据点
    double _y0[4];       //数据点
    double _x[6001];     //插值点
    double _y[6001];     //插值点

    void setX();                              //初始化数组x
    void Ln();                                //插值主函数
    void Copy(double a[], double b[], int n); //数组拷贝
};

Interpolation::Interpolation(double x0[4], double y0[4])
{
    Copy(x0, _x0, 4);
    Copy(y0, _y0, 4);
    setX();
    Ln();
}

Interpolation::Interpolation(const Interpolation &obj)
{
}

void Interpolation::setX()
{
    for (int i = 0; i < 6001; i++)
    {
        _x[i] = _x0[0] + _h * i;
    }

    return;
}

void Interpolation::Ln()
{
    double pro1 = 1.0; //分子的乘积
    double pro2 = 1.0; //分母的乘积
    double sum1 = 0.0; //4项的和
    for (int i = 0; i < 6001; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (k == j)
                    continue;
                pro1 *= (_x[i] - _x0[k]);
                pro2 *= (_x0[j] - _x0[k]);
            }
            sum1 += (pro1 / pro2) * _y0[j];
            pro1 = 1.0;
            pro2 = 1.0;
        }
        _y[i] = sum1;
        sum1 = 0.0;
    }

    return;
}

void Interpolation::Copy(double a[], double b[], int n)
{
    for (int i = 0; i < n; i++)
    {
        b[i] = a[i];
    }

    return;
}

void Interpolation::getXY(double y[6001])
{
    Copy(_y, y, 6001);

    return;
}

Interpolation::~Interpolation()
{
}