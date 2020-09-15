//Powered by Walter
//Designed for "Homework4.cpp"

#include <iostream>
#include <math.h>

using namespace std;

class Simpson
{
public:
    Simpson(int n);
    Simpson(const Simpson &obj);
    double getSn(void);
    ~Simpson();

private:
    int _n;
    double _a = 0.0;
    double _b = 5.0;
    double _n0;
    double _delta;
    double _sn;

    void setN(int n);
    void setDelta(int n);
    double Fx(double x);
    double Tn(void);
    double Hn(void);
    void Sim(void);
};

Simpson::Simpson(int n)
{
    setN(n);
    setDelta(n);
}

Simpson::Simpson(const Simpson &obj)
{
}

void Simpson::setN(int n)
{
    _n = n;
    _n0 = static_cast<double>(n);
}

void Simpson::setDelta(int n)
{
    _delta = 5.0 / static_cast<double>(n);
}

double Simpson::Fx(double x)
{
    double res;
    res = (6.0 - 10.0 * x + 5.0 * pow(x, 2)) * exp(-1.5 * x);

    return res;
}

double Simpson::Tn(void)
{
    double sum, res, xi;
    sum = 0;

    for (int i = 1; i < _n; i++)
    {
        xi = _a + (i * _delta);
        sum += Fx(xi);
    }
    res = ((_b - _a) / (2.0 * _n0)) * (Fx(_a) + 2 * sum + Fx(_b));

    return res;
}

double Simpson::Hn(void)
{
    double sum, res, xi;
    sum = 0;

    for (int i = 0; i < _n; i++)
    {
        xi = _a + (i * _delta) + (_delta / 2.0);
        sum += Fx(xi);
    }
    res = ((_b - _a) / _n0) * sum;

    return res;
}

void Simpson::Sim(void)
{
    _sn = (1.0 / 3.0) * (Tn() + 2.0 * Hn());
}

double Simpson::getSn(void)
{
    Sim();
    return _sn;
}

Simpson::~Simpson()
{
}