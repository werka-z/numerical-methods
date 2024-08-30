#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

unsigned long long factorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

void printTable(vector<vector<double>>& M){
    for (int i=0; i<M.size(); i++){
        for (int j=0; j< M[i].size(); j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

void hermiteCoefficients(vector<vector<double>>& M, vector<double>& coefficients, vector<double>& all_values) {
    int n = M.size();
    //i - columns, j  - rows
    for (int i=2; i < n+1; i++){ // right
    int curr_newidx = 0;
        for (int j=0; j < n - i + 1; j++){ // down
            double val1 = M[j][i-1];
            double val2 = M[j+1][i-1];
            double result;

            int beg = j;
            int end = j+i-1;

            if (M[beg][0] == M[end][0]){ // same x values
                result = all_values[curr_newidx+i-1]/factorial(end-beg); //last value
                M[j][i] = result;
                //printTable(M);
            }
            else{
                result = (val2-val1)/(M[end][0] - M[beg][0]);
                M[j][i] = result;
            }
            if (M[j+1][0] != M[j][0]) curr_newidx = j+1;
        }
    }

    for (int i = 0; i < coefficients.size(); ++i) {
        coefficients[i] = M[0][i+1];
    }
}


long double calcPolynomial(const vector<vector<double>>& M, double point, const vector<double>& coef) {
    long double result = 0;

    for (int i=0; i< coef.size(); i++){

        long double prod = 1;
        for (int j=0; j<i; j++){
            prod *= (point-M[j][0]);
        }

        result += coef[i]*prod;
    }
    return result;
}

int main() {
    int m, N;
    cin >> m >> N;
    
    vector<vector<double>> M(m);
    vector<double> all_values(m);

    for (int i = m+1; i > 1; i--) {
        M[m-i+1] = vector<double>(i);
    }

    for (int i = 0; i < m; ++i) {
        cin >> M[i][0]; //points in the first column
    }
    
    for (int i = 0; i < m; ++i) {
        cin >> all_values[i];
    }    
    
    double value = M[0][0];
    for (int i = 0; i < m; ++i) {
        if (i > 0 && M[i][0]==M[i-1][0]){}
        else{
            value = all_values[i];
        }
        M[i][1] = value; // repeated values in the second column
    }


    vector<double> end_points(N);
    for (int i = 0; i < N; ++i) {
        cin >> end_points[i];
    }

    
    vector<double> coef(M[0].size()-1);
    hermiteCoefficients(M, coef, all_values);

    
    for (int i = 0; i < coef.size(); ++i) {
        cout <<  setprecision(17) << coef[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < N; ++i) {
        cout << setprecision(17) << calcPolynomial(M, end_points[i], coef) << " ";
    }
    cout << endl;

    return 0;
}
