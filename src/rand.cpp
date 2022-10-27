#include "rand.hpp"
void generateRandomSU2(matrix &m)
{
    double a, b, c, d;
    double mag;
    do
    {
        a = gen_normal(rng);
        b = gen_normal(rng);
        c = gen_normal(rng);
        d = gen_normal(rng);
        mag = std::sqrt(a * a + b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);

    m << a + b * I, c + d * I, -c + d * I, a - b * I;
    m = m / mag;
}
void generateRandomSU2Rot(matrix &m)
{
    matrix temp;
    double a, b, c, d, mag, mag_rot, rot = rot_size * gen(rng);
    do
    {
        b = gen_normal(rng);
        c = gen_normal(rng);
        d = gen_normal(rng);
        mag = std::sqrt(b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);
    mag_rot = rot / mag;
    b = b * mag_rot;
    c = c * mag_rot;
    d = d * mag_rot;
    a = std::sqrt(1.0 - rot * rot);
    temp
        << a + b * I,
        c + d * I, -c + d * I, a - b * I;

    m = temp * m;
}

void generateRandomTracelessHermitian(matrix &m)
{
    double a, b, c;
    matrix tmp;
    const static double normalization = std::sqrt(6.0); // without this <Tr\phi^2> = 6. You can use this to set the expectation value to whatever we want.
    a = gen_normal(rng);
    b = gen_normal(rng);
    c = gen_normal(rng);
    tmp << a, b - c * I, b + c * I, -a;
    m = tmp / normalization;
}
void stepTracelessHermitian(matrix &m)
{
    double a, b, c;
    a = gen_normal(rng);
    b = gen_normal(rng);
    c = gen_normal(rng);
    matrix temp;
    temp << a, b - c * I, b + c * I, -a;
    temp = temp * step_size_higgs;
    m = m + temp;
}