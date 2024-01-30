#include "simulation.hpp"
#include <memory>
#pragma once
void generateRandomSU2(matrix &m);
void generateRandomSU2Rot(matrix &m);
void generateRandomTracelessHermitian(matrix &m);
void stepTracelessHermitian(matrix &m);