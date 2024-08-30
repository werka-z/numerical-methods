#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

float posZero(float value) {
    if (value == 0.0f) {
        return 0.0f;
    }
    return value;
}

int main(){
    
    int n;
    float sum, prod, q;
    float x1, x2, x3;

    cin >> n;
    
    float b, delta;

    for (int i = 0; i < n; i++){
        cin >> prod >> sum;

        // q^2 + (1-(sum/V^3(prod))) + 1
        b = 1 - sum/cbrt(prod);

        // delta = (b-2)*(b+2);
        x2 = cbrt(prod);
        delta = (-1 - sum/cbrt(prod))*(3 - sum/cbrt(prod));

        if (!isnan(delta) && !isinf(delta) && delta >= 0){

            q = -b;

            if (q > 0) {
                q += sqrt(delta);
            }

            else{ //q<0
                q -= sqrt(delta);
            }

            q/=2;

            x1 = x2*q;
            x3 = x2/q;

            if (x1<x3){
                swap(x1, x3);
            }
        }
        else {
            x1 = 0, x2 = 0, x3 = 0;
        }

        cout << scientific << setprecision(10) << posZero(x1) << " " << posZero(x2) << " " << posZero(x3) << endl;
    }
}