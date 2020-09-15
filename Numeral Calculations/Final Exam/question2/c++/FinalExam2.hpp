//Powered By Walter
//期末试题第二题的具体实现
#include <iostream>
#include <math.h>

using namespace std;

class RungeKutta
{
public:
    double getYn1();                                      //返回下一个y值
    double getZn1();                                      //返回下一个z值
    RungeKutta(double xn, double yn, double zn, int num); //构造函数
    RungeKutta(const RungeKutta &obj);                    //拷贝构造函数
    ~RungeKutta();                                        //析构函数
private:
    int _num;         //选择函数用的
    double _xn;       //当前的x值
    double _yn;       //当前的y值
    double _zn;       //当前的z值
    double _yn1;      //下一个y值
    double _zn1;      //下一个z值
    double _k[2][4];  //8个K参数
    double _h = 0.02; //步长

    void setK1s();                            //计算测试函数各个K的值
    void setK2s();                            //计算题目函数各个K的值
    void setYZ1();                            //计算下一个y和z的值
    double fx1(double x, double y, double z); //测试用的函数
    double fx2(double x, double y, double z); //题目中的函数
};

RungeKutta::RungeKutta(double xn, double yn, double zn, int num)
{
    _xn = xn;
    _yn = yn;
    _zn = zn;
    switch (num)
    {
    case 1:
        setK1s();
        break;

    default:
        setK2s();
        break;
    }
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

void RungeKutta::setK1s()
{
    _k[0][0] = _zn;
    _k[1][0] = fx1(_xn, _yn, _zn);
    _k[0][1] = _zn + 0.5 * _h * _k[1][0];
    _k[1][1] = fx1((_xn + 0.5 * _h), (_yn + 0.5 * _h * _k[0][0]), (_zn + 0.5 * _h * _k[1][0]));
    _k[0][2] = _zn + 0.5 * _h * _k[1][1];
    _k[1][2] = fx1((_xn + 0.5 * _h), (_yn + 0.5 * _h * _k[0][1]), (_zn + 0.5 * _h * _k[1][1]));
    _k[0][3] = _zn + _h * _k[1][2];
    _k[1][3] = fx1((_xn + _h), (_yn + _h * _k[0][2]), (_zn + _h * _k[1][2]));

    return;
}

void RungeKutta::setK2s()
{
    _k[0][0] = _zn;
    _k[1][0] = fx2(_xn, _yn, _zn);
    _k[0][1] = _zn + 0.5 * _h * _k[1][0];
    _k[1][1] = fx2((_xn + 0.5 * _h), (_yn + 0.5 * _h * _k[0][0]), (_zn + 0.5 * _h * _k[1][0]));
    _k[0][2] = _zn + 0.5 * _h * _k[1][1];
    _k[1][2] = fx2((_xn + 0.5 * _h), (_yn + 0.5 * _h * _k[0][1]), (_zn + 0.5 * _h * _k[1][1]));
    _k[0][3] = _zn + _h * _k[1][2];
    _k[1][3] = fx2((_xn + _h), (_yn + _h * _k[0][2]), (_zn + _h * _k[1][2]));

    return;
}

void RungeKutta::setYZ1()
{
    _yn1 = _yn + (_h / 6.0) * (_k[0][0] + 2.0 * _k[0][1] + 2.0 * _k[0][2] + _k[0][3]);
    _zn1 = _zn + (_h / 6.0) * (_k[1][0] + 2.0 * _k[1][1] + 2.0 * _k[1][2] + _k[1][3]);
}

double RungeKutta::fx1(double x, double y, double z)
{
    double res;

    x = x;
    res = -y - 2 * z;

    return res;
}

double RungeKutta::fx2(double x, double y, double z)
{
    double res;

    z = z;
    res = ((18.0 * exp(6.0 * x - 12.0)) / pow((1.0 + exp(3.0 * x - 6.0)), 2)) * y - (9.0 * exp(3.0 * x - 6.0)) / pow((1.0 + exp(3.0 * x - 6.0)), 2);

    return res;
}

RungeKutta::~RungeKutta()
{
}