//Powered by Walter
//Include "Simpson.hpp"

#include "Simpson.hpp"
#include <iostream>
#include <math.h>

using namespace std;

int main()
{
    int par;
    double sn1, sn2, del;
    Simpson *simpson = NULL;
    par = 10;
    del = 0.000001;

    simpson = new Simpson(par);
    sn1 = simpson->getSn();
    delete simpson;
    for (int i = 0; i < 100; i++)
    {
        Simpson *simpson = NULL;
        par *= 2;

        simpson = new Simpson(par);
        sn2 = simpson->getSn();
        delete simpson;
        if (pow((sn2 - sn1), 2) < pow(del, 2))
        {
            sn1 = sn2;
            cout << par << endl;
            break;
        }
        else
        {
            sn1 = sn2;
        }
    }
    cout << "The result is : " << sn1 << endl;

    system("pause");

    return 0;
}