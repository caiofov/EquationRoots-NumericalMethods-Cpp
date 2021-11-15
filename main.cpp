#include <iostream>
#include <cmath>
using namespace std;
#include "methods.h"

//https://en.wikipedia.org/wiki/Root-finding_algorithms

float funcao2(float x){
  return x - x*log(x);
}
float funcaoDerivada2(float x){
  return -log(x);
}
float funcaoIteracao2(float x){
  return x*log(x);
}
float funcaoIteracao3(float x){
  return x/log(x);
}
float funcaoIteracao4(float x){
  return 2*x - x*log(x);
}

int main() {
  printf("Raiz: %f", pontoFixo(2, funcao2, funcaoIteracao4, 0.00001));
}