//Powered By Walter
//四阶Runge-Kutta方法数值解二阶微分方程程序的具体实现

#include <iostream>
#include <math.h>

using namespace std;

static double h = 0.01; //设置步长

class RungeKutta
{
public:
    RungeKutta(double yn, double zn);  //构造函数
    RungeKutta(const RungeKutta &obj); //拷贝构造函数
    double getYn1();                   //返回下一个y值
    double getZn1();                   //返回下一个z值
    ~RungeKutta();                     //析构函数
private:
    double _yn;      //当前的y值
    double _yn1;     //下一个y值
    double _zn;      //当前的z值
    double _zn1;     //下一个z值
    double _k[2][4]; //8个K参数
    void setKs();    //计算各个K的值
    void setYZ1();   //计算下一个y和z的值
};

RungeKutta::RungeKutta(double yn, double zn)
{
    _yn = yn;
    _zn = zn;
    setKs();
    setYZ1();
}

RungeKutta::RungeKutta(const RungeKutta &obj)
{
}

double RungeKutta::getYn1()
{
    return _yn1;
}

double RungeKutta::getZn1()
{
    return _zn1;
}

void RungeKutta::setKs()
{
    _k[0][0] = _zn;
    _k[1][0] = -_yn - (2 * _zn);
    _k[0][1] = _zn + (h / 2) * _k[1][0];
    _k[1][1] = -1 * (_yn + (h / 2) * _k[0][0]) - 2 * (_zn + (h / 2) * _k[1][0]);
    _k[0][2] = _zn + (h / 2) * _k[1][1];
    _k[1][2] = -1 * (_yn + (h / 2) * _k[0][1]) - 2 * (_zn + (h / 2) * _k[1][1]);
    _k[0][3] = _zn + h * _k[1][2];
    _k[1][3] = -1 * (_yn + h * _k[0][2]) - 2 * (_zn + h * _k[1][2]);

    return;
}

void RungeKutta::setYZ1()
{
    _yn1 = _yn + (h / 6) * (_k[0][0] + 2 * _k[0][1] + 2 * _k[0][2] + _k[0][3]);
    _zn1 = _zn + (h / 6) * (_k[1][0] + 2 * _k[1][1] + 2 * _k[1][2] + _k[1][3]);
}

RungeKutta::~RungeKutta()
{
}