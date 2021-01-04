#include <iostream>
#include <math.h>
#include <string>
void significance()
{
    // if s fit , PARA = 1; if likelihood fit , PARA = 2
    int para, n;
    double abs_diff, diff, prob, significance;

    cout << "If result is value of s, please set para = 1, else if result is likelihood, set para = 2" << endl;
    cout << "Please enter the value of para:" << endl;
    cin >> para;

    if (para != 1 && para != 2)
    {
        cerr << "ERROR: YOU ENTER THE WRONG NUMBER." << endl;
        return;
    }

    cout << "Please enter the number of parameters:" << endl;

    cin >> n;

    if (para == 1)
    {
        cout << "Please enter the value of Delta s:" << endl;
    }
    else if (para == 2)
    {
        cout << "Please enter the value of Delta likelihood:" << endl;
    }

    cin >> diff;

    abs_diff = fabs(diff);

    prob = TMath::Prob(abs_diff * para, n);

    significance = RooStats::PValueToSignificance(prob * 0.5);

    cout << "Your sets:" << endl;
    if (para == 1)
    {
        cout << "para = " << para << ", number of parameters = " << n << ", Delta s = " << diff << endl;
    }
    else if (para == 2)
    {
        cout << "para = " << para << ", number of parameters = " << n << ", Delta likelihood = " << diff << endl;
    }
    cout << "Significance WHICH_PARAM : " << significance << endl;
}