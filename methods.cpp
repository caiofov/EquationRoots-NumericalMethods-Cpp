#include <iostream>
#include <cmath>
using namespace std;
#include "methods.h"

float falsePosition(float a, float b, float precisao, float(*function)(float), int count){
  count++;
  float fa = function(a); //f(a)
  float fb = function(b); //f(b)
  
  if( (abs(b-a) > precisao) && (fa*fb < 0)){
    float pm = ((a*fb)-(b*fa))/(fb - fa); //ponto medio
    float fp = function(pm); //f(ponto medio)
        
    if((fp*fa < 0) && (abs(pm - a) > precisao) && (abs(fp) > precisao)){
      return falsePosition(a, pm, precisao, function, count);
    }
    
    else if((fp * fb < 0) && (abs(b-pm) > precisao) && (abs(fp) > precisao)){
      return falsePosition(pm, b, precisao, function, count);
    }
    
    else{
      return pm;
    }
  }
 
  else{
    return -999999;
  }
}

float bissection(float a, float b, float precisao, float(*function)(float), int count){
  count++;
  float fa = function(a); //f(a)
  float fb = function(b); //f(b)

  if(fa == 0){ //checks if it's an exact root
    return a;
  }
  if(fb == 0){ //checks if it's an exact root
    return b;
  }
  
  if( (abs(b-a) > precisao) && (fa*fb < 0)){
    float pm = (a+b)/2; //ponto medio
    float fp = function(pm); //f(ponto medio)

    if((fp*fa < 0) 
    && (abs(pm - a) > precisao) 
    && (abs(fp) > precisao)){
      return bissection(a, pm, precisao, function, count);
    }
    
    else if((fp * fb < 0) 
    && (abs(b-pm) > precisao) 
    && (abs(fp) > precisao)){
      return bissection(pm, b, precisao, function, count);
    }
    
    else{
      if(fp == 0){ //checks if it's an exact root
      }
      return pm;
    }
  }
 
  else{
    printf("Error \n");
    return -999999;
  }
}

float newtonRaphson(float x0, float precisao, float (*function)(float), float (*functionderivative)(float), int count){
  count++;
  printf("\n%dª Iteração \n", count);
  float fx0 = function(x0);
  float phiX = x0 - (fx0/ functionDerivative(x0));
  
  if(fx0 == 0){ //checks if it's an exact root
    return x0;}; 
  
  if(abs(fx0) > precisao){
    return newtonRaphson(phiX, precisao, function, functionDerivative, count);
  }
  else{
    return x0;
  }

}

float secant(float x0, float x1, float precisao, float (*function)(float), float (*functionderivative)(float), int count){
  count++;
  
  float fx0 = function(x0);
  float fx1 = function(x1);
  
  float phiX2 = ((x0*fx1) - x1*fx0)/(fx1 - fx0);

  if(fx0 == 0){ //checa se é raiz
    return x0;
    }; 
  if(fx1 == 0){ //checa se é raiz
    return x1;
  } 

  if(abs(fx0) > precisao){
    return secant(x1,phiX2,precisao, function, functionDerivative, count);
  }

  else{
    return phiX2;
  }


}

float pontoFixo(float x0, float (*function)(float), float (*functionIteration)(float), float precisao, int count){
  count++;

  float phix0 = functionIteration(x0);
  float fx0 =  function(x0);

  if(fx0 == 0){ //checa se x0 é raiz
    cout << "x0 é raiz" << endl;
    return x0;
  }

  if( abs(fx0) > precisao){
    
    return pontoFixo(phix0, function, functionIteration, precisao, count);
  }

  if(isinf(x0) || isnan(x0)){
    cout << "Tip: try another iteration function" << endl;
  }
  
  return x0;
  
}

float calculoPolinomio(int numCoef, float coeficientes[numCoef], float x){
  int grau = numCoef - 1;
  float z = x;

  float b = coeficientes[grau];
  
  for(int i = grau-1; i > -1; i--){
    b = coeficientes[i] + b*z;
  }

  return b;
}

float calculateDerivative(int numCoef, float coeficientes[numCoef], float x){
  
  int grau = numCoef - 1;
  float z = x;

  float b = coeficientes[grau];
  float c = b;
  
  for(int i = grau-1; i > 0; i--){
    b = coeficientes[i] + b*z;
    c = b + c*z;
  }

  return c;
  
}

float polinomioNewton(int numCoef, float coeficientes[numCoef], float x, float precisao, int maxIter){
  int grau = numCoef - 1;
  float b,c;

  for(int k = 0; k < maxIter + 1; k++){
    
    //inicializa as variaveis B e C
    b = coeficientes[grau];
    c = b;
    
    for(int i = grau-1; i > 0; i--){
      b = coeficientes[i] + b*x;
      c = b + c*x;
    }
    b = coeficientes[0] + b*x;

    //mostrando os valores
    printf("\n%dª Iteração: \nb: %f | c: %f \nx: %f\n", k+1,b,c, x);

    //b = f(x) -> vendo se x já é um bom valor
    if(abs(b) <= precisao){
      return x;
    }

    //atualizando o valor de X pela função de iteração do método de newton
    x = x - b/c;
  }

  cout << "Número máximo de iterações" << endl;
  return -99999;


}


