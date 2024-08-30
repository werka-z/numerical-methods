#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;

typedef void (*FuncPointer)(const double* x, double* y, double* Df);

void printVector(const double* x, unsigned N){
  for(unsigned i=0;i<N;++i)
    printf("%17.17f ",x[i]);
  printf("\n");
}

int findCurve(FuncPointer f, double* x, unsigned k, double h) {
  double y[2];
  double Df[6];
  double delta[2];
  double x_original[3] = {x[0], x[1], x[2]};
  double d = 0;

  for (int i = 1; i <= k; i++) {
        x[2] = x_original[2] + h*i;
        f(x, y, Df);
        
        do {
            d = Df[0] * Df[4] - Df[1] * Df[3];

            delta[0] = (Df[4] * y[0] - Df[1] * y[1]) / d;
            delta[1] = (Df[0] * y[1] - Df[3] * y[0]) / d;

            x[0] -= delta[0];
            x[1] -= delta[1];

            f(x, y, Df);

            if (fabs(x[0]) >= fabs(x_original[0]) + 500 || fabs(x[1]) >= fabs(x_original[1]) + 500) return i;
            
        } while (fabs(y[0]) > 1e-14 || fabs(y[1]) > 1e-14);

        printVector(x, 3);
  }
  return 0;
}

int findSurface(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2) {
    double y;
    double Df[3];
    double x_original[3] = {x[0], x[1], x[2]};
    double delta;

    for (int i = 0; i < k1; i++) {
        for (int j = 0; j < k2; j++) {
            x[0] = x_original[0]; // reset
            x[1] = x_original[1] + h1 * (i + 1);
            x[2] = x_original[2] + h2 * (j + 1);
            
            int iter = 0;
            do {
                f(x, &y, Df);
                
                if (fabs(y) < 1e-14) break;
                if (fabs(Df[0]) < 1e-14) return i * k1 + j + 1;
                
                delta = -y / Df[0];
                x[0] += delta;
                
                if (iter > 500) return i * k1 + j + 1;
                iter++;
                
            } while (fabs(delta) >= 1e-14);
            
            printVector(x, 3);
        }
        cout << endl;
    }
    return 0;
}

int findFixedPoints(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2) {
  double y[2], Df[8], delta[2];
  double x_original[4] = {x[0], x[1], x[2], x[3]};
  double d = 0;

  for (int i = 0; i < k1; i++) {
    for (int j = 0; j < k2; j++) {
      x[0] = x_original[0]; // reset
      x[1] = x_original[1];
      x[2] = x_original[2] + h1 * (i+1);
      x[3] = x_original[3] + h2 * (j+1);
      
      int iter = 0;
      do {
        f(x, y, Df);

        y[0] -= x[0];
        y[1] -= x[1];

        Df[0] -= 1;
        Df[5] -= 1;
        
        d = Df[0] * Df[5] - Df[1] * Df[4];

        if (fabs(d) < 10e-14) return i * k1 + j + 1;
        
        delta[0] = (Df[5] * y[0] - Df[1] * y[1]) / d;
        delta[1] = (-Df[4] * y[0] + Df[0] * y[1]) / d;

        x[0] -= delta[0];
        x[1] -= delta[1];

      } while (fabs(y[0]) > 1e-14 || fabs(y[1]) > 1e-14); 
    
      printVector(x, 4);
    }
    cout << endl;
  }
  return 0;
}