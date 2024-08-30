#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

template <typename T>
T funkcja(const T &x, const T &y)
{
    T w = sin(x * x - 2 * (y + 1)) / exp(-y * y + cos(x * y));
    return w;
}

template <typename T>
class Jet
{
public:
    T values[6];

    Jet()
    {
        for (int i = 0; i < 6; ++i)
            values[i] = 0;
    }

    Jet(const double x)
    {
        for (int i = 1; i < 6; ++i)
            values[i] = 0;
        values[0] = x;
    }

    Jet(const Jet &other)
    {
        for (int i = 0; i < 6; ++i)
            values[i] = other.values[i];
    }

    Jet &operator=(const Jet &other)
    {
        if (this != &other)
        {
            for (int i = 0; i < 6; ++i)
                values[i] = other.values[i];
        }
        return *this;
    }

    Jet operator-() const
    {
        Jet result = Jet();
        for (int i = 0; i < 6; ++i)
        {
            result.values[i] = -values[i];
        }
        return *this;
    }

    // Jet + Jet
    Jet operator+(const Jet &other) const
    {
        Jet result;
        for (int i = 0; i < 6; ++i)
        {
            result.values[i] = values[i] + other.values[i];
        }
        return result;
    }

    Jet operator-(const Jet &other) const
    {
        Jet result;
        for (int i = 0; i < 6; ++i)
        {
            result.values[i] = values[i] - other.values[i];
        }
        return result;
    }

    Jet operator*(const Jet &other) const
    {
        Jet result;
        result.values[0] = values[0] * other.values[0];

        result.values[1] = values[0] * other.values[1] + values[1] * other.values[0];
        result.values[2] = values[0] * other.values[2] + values[2] * other.values[0];

        result.values[3] = values[3] * other.values[0] + 2 * values[1] * other.values[1] + values[0] * other.values[3];
        result.values[4] = values[4] * other.values[0] + values[1] * other.values[2] + values[2] * other.values[1] + values[0] * other.values[4];
        result.values[5] = values[5] * other.values[0] + 2 * values[2] * other.values[2] + values[0] * other.values[5];
        return result;
    }

    Jet operator/(const Jet &other) const
    {
        Jet result;
        T g = other.values[0];
        T g2 = g * g;

        result.values[0] = values[0] / g;

        result.values[1] = (values[1] * g - values[0] * other.values[1]) / g2;
        result.values[2] = (values[2] * g - values[0] * other.values[2]) / g2;

        // dxx = (dxx*g - 2dx*o.dx + f*o.dxx) / g2
        result.values[3] = (values[3] * g - 2 * values[1] * other.values[1] + values[0] * other.values[3]) / g2;

        // dxy = (dxy*g - dx*o.dy - dy*o.dx + f*dxy ) / g2
        result.values[4] = (values[4] * g - values[1] * other.values[2] - values[2] * other.values[1] + values[0] * other.values[4]) / g2;

        // dyy = (dyy*g - 2dy*o.dy + f*o.dyy) / g2
        result.values[5] = (values[5] * g - 2 * values[2] * other.values[2] + values[0] * other.values[5]) / g2;
        return result;
    }

    Jet sin(const Jet &jet) const
    {
        // sin (u, u') = (sin u, cos u * u')
        Jet result;
        result.values[0] = std::sin(values[0]);
        result.values[1] = std::cos(values[0]) * values[1];
        result.values[2] = std::cos(values[0]) * values[2];
        result.values[3] = -std::sin(values[0]) * values[1] * values[1] + std::cos(values[0]) * values[3];
        result.values[4] = -std::sin(values[0]) * values[1] * values[2] + std::cos(values[0]) * values[4];
        result.values[5] = -std::sin(values[0]) * values[2] * values[2] + std::cos(values[0]) * values[5];
        return result;
    }

    Jet cos(const Jet &jet) const
    {
        // cos (u, u') = (cos u, -sin u * u')
        Jet result;
        result.values[0] = std::cos(values[0]);
        result.values[1] = -std::sin(values[0]) * values[1];
        result.values[2] = -std::sin(values[0]) * values[2];
        result.values[3] = -std::cos(values[0]) * values[1] * values[1] - std::sin(values[0]) * values[3];
        result.values[4] = -std::cos(values[0]) * values[1] * values[2] - std::sin(values[0]) * values[4];
        result.values[5] = -std::cos(values[0]) * values[2] * values[2] - std::sin(values[0]) * values[5];
        return result;
    }

    Jet exp(const Jet &jet) const
    {
        // e^(u, u') = (e^u, e^u * u')
        Jet result;
        T e = std::exp(values[0]);
        result.values[0] = e;
        result.values[1] = e * values[1];
        result.values[2] = e * values[2];
        result.values[3] = e * (values[3] + values[1] * values[1]);
        result.values[4] = e * (values[4] + values[1] * values[2]);
        result.values[5] = e * (values[5] + values[2] * values[2]);
        return result;
    }
};

// Jet + const
template <typename T>
Jet<T> operator+(const T &constant) const
{
    Jet result(*this);
    result.values[0] += constant; // tylko do wartosci funkcji?
    return result;
}

template <typename T>
Jet<T> operator-(const T &constant) const
{
    Jet result(*this);
    result.values[0] -= constant;
    return result;
}

template <typename T>
Jet<T> operator*(const T &constant) const
{
    Jet result = Jet();
    result.values[0] = constant;
    for (int i = 1; i < 6; ++i)
    {
        result.values[i] = 0;
    }
    return result;
}

template <typename T>
Jet<T> operator/(const T &constant) const
{
    Jet result - Jet();
    for (int i = 0; i < 6; ++i)
    {
        result.values[i] = values[i] / constant;
    }
    return result;
}

template <typename T>
Jet<T> operator*(const double scalar, const Jet<T> &jet)
{
    Jet<T> result = Jet();
    result.values[0] = scalar;
    for (int i = 1; i < 6; ++i)
    {
        result.values[i] = 0; // I guess that's how multiplication in FAD works
    }
    return result;
}

template <typename T>
Jet<T> operator+(const int scalar, const Jet<T> &jet)
{
    Jet<T> result(jet);
    result.values[0] += scalar;
    return result;
}

template <typename T>
Jet<T> operator+(const Jet<T> &jet, const int scalar)
{
    return scalar + jet;
}

template <typename T>
Jet<T> operator-(const Jet<T> &jet, const int scalar)
{
    return scalar + -jet;
}

template <typename T>
Jet<T> operator/(const Jet<T> &jet, const int scalar)
{
    return jet / scalar;
}

int main()
{
    int M;
    cin >> M;

    for (int i = 0; i < M; i++)
    {
        double x0, y0;
        cin >> x0 >> y0;

        Jet<double> x, y;
        x.values[0] = x0;
        x.values[1] = 1;
        x.values[2] = 0;
        x.values[3] = 0;
        x.values[4] = 0;
        x.values[5] = 0;

        y.values[0] = y0;
        y.values[1] = 0;
        y.values[2] = 1;
        y.values[3] = 0;
        y.values[4] = 0;
        y.values[5] = 0;

        Jet result = funkcja(x, y);

        cout << "Wartość w punkcie (" << x0 << ", " << y0 << "): " << result.values[0] << endl;
        cout << "Pochodne: dx=" << result.values[1] << ", dy=" << result.values[2];
        cout << ", dxx=" << result.values[3] << ", dxy=" << result.values[4] << ", dyy=" << result.values[5] << endl;
    }
}