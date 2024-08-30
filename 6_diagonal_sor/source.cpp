#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

vector<double> SORwVectors(vector<vector<double>> bands, int n, vector<double> b, vector<double> x0, int iterations, double SORparam, int M)
{
    for (int k = 0; k < iterations; k++)
    {
        for (int i = 0; i < n; i++)
        {
            auto s = b[i];
            ///pod diagonala
            for (int j=max(0, i-M); j<i; j++){ 
                s -= bands[i-j][j] * x0[j];
            }

            //nad diagonala
            for (int j = max(0, i+1); j < i+1+M && j < n; j++)
            {
                s -= bands[j-i][i] * x0[j];
            }
            
            x0[i] = (1 - SORparam) * x0[i] + SORparam * s / bands[0][i];
        }
    }
    return x0;
}

int main()
{
    int N;    // rozmiar
    int offset;    // l wstęg
    int iterations;    // liczba iteracji
    double SORparam; // parametr
    cin >> N >> offset;

    vector<vector<double>> bands(offset + 1); 

    for (int i = offset; i >= 0; i--)
    {
        bands[i].resize(N-i); //coraz mniejsze wstegi
        for (int j=0; j < bands[i].size(); j++)
        {
            cin >> bands[i][j];
        }
    }

    // prawa strona równania
    vector<double> right(N);
    for (int i = 0; i < N; ++i)
    {
        cin >> right[i];
    }

    // wektor początkowy
    vector<double> x0(N);
    for (int i = 0; i < N; ++i)
    {
        cin >> x0[i];
    }

    cin >> SORparam >> iterations;

    vector<double> solution = SORwVectors(bands, N, right, x0, iterations, SORparam, offset);

    for (int i = 0; i < N; ++i)
    {
        cout << setprecision(17) << scientific << solution[i] << endl;
    }
}
