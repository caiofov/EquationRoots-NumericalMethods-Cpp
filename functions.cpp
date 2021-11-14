#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;
#include "functions.h"

float funcao(float x){
  return -(exp(x)/2) + 2*cos(x);
}

float funcaoDerivada(float x){
  return -(exp(x)/2) - x*sin(x);
}

float funcaoModificada (float x){
  return -(exp(x)/2) + 2*cos(x) - M_PI/4;
}

float funcaoPendulo (float x){
  return pow(x,3)-9*x+3;
}
