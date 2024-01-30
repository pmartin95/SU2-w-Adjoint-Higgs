#include "rand.hpp"
void generateRandomSU2(std::unique_ptr<Simulation> &sim, matrix &m)
{
    double a, b, c, d;
    double mag;
    do
    {
        a = sim->gen_normal(sim->rng);
        b = sim->gen_normal(sim->rng);
        c = sim->gen_normal(sim->rng);
        d = sim->gen_normal(sim->rng);
        mag = std::sqrt(a * a + b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);

    m << a + b * Simulation::I, c + d * Simulation::I, -c + d * Simulation::I, a - b * Simulation::I;
    m = m / mag;
}
void generateRandomSU2Rot(std::unique_ptr<Simulation> &sim, matrix &m)
{
    matrix temp;
    double a, b, c, d, mag, mag_rot, rot = sim->rot_size * sim->gen(sim->rng);
    do
    {
        b = sim->gen_normal(sim->rng);
        c = sim->gen_normal(sim->rng);
        d = sim->gen_normal(sim->rng);
        mag = std::sqrt(b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);
    mag_rot = rot / mag;
    b = b * mag_rot;
    c = c * mag_rot;
    d = d * mag_rot;
    a = std::sqrt(1.0 - rot * rot);
    temp
        << a + b * Simulation::I,
        c + d * Simulation::I, -c + d * Simulation::I, a - b * Simulation::I;

    m = temp * m;
}

void generateRandomTracelessHermitian(std::unique_ptr<Simulation> &sim, matrix &m)
{
    double a, b, c;
    matrix tmp;
    const static double normalization = std::sqrt(6.0); // without this <Tr\phi^2> = 6. You can use this to set the expectation value to whatever we want.
    a = sim->gen_normal(sim->rng);
    b = sim->gen_normal(sim->rng);
    c = sim->gen_normal(sim->rng);
    tmp << a, b - c * Simulation::I, b + c * Simulation::I, -a;
    m = tmp / normalization;
}
void stepTracelessHermitian(std::unique_ptr<Simulation> &sim, matrix &m)
{
    double a, b, c;
    a = sim->gen_normal(sim->rng);
    b = sim->gen_normal(sim->rng);
    c = sim->gen_normal(sim->rng);
    matrix temp;
    temp << a, b - c * Simulation::I, b + c * Simulation::I, -a;
    temp = temp * sim->step_size_higgs;
    m = m + temp;
}