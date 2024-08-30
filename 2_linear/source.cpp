#include <cmath>
#include <iostream>
using namespace std;

double findZero(
    double (*f)(double),  // funkcja której zera szukamy w [a, b] 
    double a,             // lewy koniec przedziału
    double b,             // prawy koniec przedziału
    int M,                // maksymalna dozwolona liczba wywołań funkcji f
    double eps,           // spodziewana dokładność zera
    double delta          // wystarczający błąd bezwzględny wyniku
){
// mniej niz polowa krokow - bisekcje (M), czy sa roznych znakow
    if (M <= 0) return 0;
    if (M==1 || M==2) return a;

    double fa = f(a);
    double fb = f(b), c=0, fc = 0;
    bool flag = false;

    for (int i = 0; i < M-2; i++){
        //pierwsza polowa bisekcji, o ile sa roznych znakow, potem sieczne
        if (i < M/2 && fa*fb <= 0){
            c = (a+b)/2;
            fc = f(c);
            if (fc*fa <= 0){
                fb = fc;
                b = c;
            }
            else{
                fa = fc;
                a = c;
            }
        }
        else {//sieczne
            c = a - (fa * (b - a)) / (fb - fa);
            fc = f(c); 

            if (flag){
                b = c;
                fb = fc;
            }

            else{
                a = c;
                fa = fc;
            }

            flag = !flag;
        }

        if (fabs(fc) <= eps || fabs(b-a) <= delta) return c;
        
    }
    return c;
}