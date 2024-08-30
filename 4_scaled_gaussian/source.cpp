#include "vectalg.h"
#include <cmath>
using namespace std;


Vector GaussEliminationSolve(Matrix &A, Vector &b){
    int n = A.size();
    Vector scales(n); // skale każdego wiersza

    // 1. scales - wektor do którego wkładamy normy wierszy
    for (int i = 0; i < n; i++) {
        Vector row(n);
        for (int j = 0; j < n; j++) {
            row[j] = A(i,j);
        }
        scales[i] = row.max_norm();
    }

    // 2. główna pętla: w k-tym kroku bierzemy element a(i,k) (i-wiersz, k-kolumna), którego (a(ik))/scales[i] jest największe, 
    // czyli z k-tej kolumny kolumny wybieramy najwiekszy po skalowaniu element
    for (int k = 0; k < n; k++) {
        int maxIdx = k;
        double maxVal = 0;
        for (int i = k; i < n; i++) {
            double val = abs(A(i,k)) / scales[i];
            if (val > maxVal) {
                maxVal = val;
                maxIdx = i;
            }
        }

        //3. swapujemy wiersz k i wybrany X w macierzy, w wektorze też przestawiamy element k i X
        for (int j = 0; j < n; j++) {
            double swp = A(maxIdx, j);
            A(maxIdx, j) = A(k, j);
            A(k, j) = swp;
        }
        
        double swp = b[k];
        b[k] = b[maxIdx];
        b[maxIdx] = swp;


        //4. zerujemy wszystko pod tym X (który aktualnie jest na przekątnej), itd.
        for (int i = k + 1; i < n; i++) {
                double factor = A(i, k) / A(k,k);
                
                for (int j = k; j < n; ++j) {
                    A(i, j) -= factor * A(k,j);
                }
                b[i] -= factor * b[k]; //no i w tym wektorze przydałoby się też
                A(i, k) = 0;
        }
    }

    Vector x(n); // rozwiazanie układu
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A(i,j) * x[j];
        }
        x[i] = sum / A(i,i);
    }
    return x;
}

Vector solveEquations(const Matrix & M, const Vector & l, double eps){
    Matrix A = M;
    Vector b = l;
    
    //Ax = b, 
    int n = M.size();

    Vector x = GaussEliminationSolve(A, b);
    // iteracyjne poprawianie
    Vector r = residual_vector(A, b, x);
    while (r.max_norm() >= eps) {
        A = M;
        Vector e = GaussEliminationSolve(A, r); // Ae = r
        for (int i = 0; i < n; ++i) {
            x[i] += e[i];
        }
        r = residual_vector(A, b, x);
    }
    return x;
}
