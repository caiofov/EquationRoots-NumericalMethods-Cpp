#include <iostream>
#include <cmath>
using namespace std;
#include "methods.h"

float falsePosition(float a, float b, float precision, float(*function)(float), int count){
  count++;
  float fa = function(a); //f(a)
  float fb = function(b); //f(b)
  
  if( (abs(b-a) > precision) && (fa*fb < 0)){
    float pm = ((a*fb)-(b*fa))/(fb - fa); //interception point
    float fp = function(pm); //f(interception point)
        
    if((fp*fa < 0) && (abs(pm - a) > precision) && (abs(fp) > precision)){
      return falsePosition(a, pm, precision, function, count);
    }
    
    else if((fp * fb < 0) && (abs(b-pm) > precision) && (abs(fp) > precision)){
      return falsePosition(pm, b, precision, function, count);
    }
    
    else{
      return pm;
    }
  }
 
  else{
    return -999999;
  }
}

float bissection(float a, float b, float precision, float(*function)(float), int count){
  count++;
  float fa = function(a); //f(a)
  float fb = function(b); //f(b)

  if(fa == 0){ //checks if it's an exact root
    return a;
  }
  if(fb == 0){ //checks if it's an exact root
    return b;
  }
  
  if( (abs(b-a) > precision) && (fa*fb < 0)){
    float pm = (a+b)/2; //middle of the interval
    float fp = function(pm); //f(middle)
    if((fp*fa < 0) 
    && (abs(pm - a) > precision) 
    && (abs(fp) > precision)){
      return bissection(a, pm, precision, function, count);
    }
    
    else if((fp * fb < 0) 
    && (abs(b-pm) > precision) 
    && (abs(fp) > precision)){
      return bissection(pm, b, precision, function, count);
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

float newtonRaphson(float x0, float precision, float (*function)(float), float (*functionDerivative)(float), int count){
  count++;
  float fx0 = function(x0);
  float phiX = x0 - (fx0/ functionDerivative(x0));
  
  if(fx0 == 0){ //checks if it's an exact root
    return x0;}; 
  
  if(abs(fx0) > precision){
    return newtonRaphson(phiX, precision, function, functionDerivative, count);
  }
  else{
    return x0;
  }

}

float secant(float x0, float x1, float precision, float (*function)(float), float (*functionDerivative)(float), int count){
  count++;
  
  float fx0 = function(x0);
  float fx1 = function(x1);
  
  float phiX2 = ((x0*fx1) - x1*fx0)/(fx1 - fx0);

  if(fx0 == 0){ //checks if it's an exact root
    return x0;
    }; 
  if(fx1 == 0){ //checks if it's an exact root
    return x1;
  } 

  if(abs(fx0) > precision){
    return secant(x1,phiX2,precision, function, functionDerivative, count);
  }

  else{
    return phiX2;
  }


}

float pontoFixo(float x0, float (*function)(float), float (*functionIteration)(float), float precision, int count){
  count++;

  float phix0 = functionIteration(x0);
  float fx0 =  function(x0);

  if(fx0 == 0){ //checks if it's an exact root
    return x0;
  }

  if( abs(fx0) > precision){
    
    return pontoFixo(phix0, function, functionIteration, precision, count);
  }

  if(isinf(x0) || isnan(x0)){
    printf("Tip: try another iteration function");
  }
  
  return x0;
  
}

float calculatePolynomial(int coefNum, float coeficients[], float x){
  int degree = coefNum - 1;
  float z = x;

  float b = coeficients[degree];
  
  for(int i = degree-1; i > -1; i--){
    b = coeficients[i] + b*z;
  }

  return b;
}

float calculateDerivative(int coefNum, float coeficientes[], float x){
  
  int degree = coefNum - 1;
  float z = x;

  float b = coeficientes[degree];
  float c = b;
  
  for(int i = degree-1; i > 0; i--){
    b = coeficientes[i] + b*z;
    c = b + c*z;
  }

  return c;
  
}

float polynomialNewton(int coefNum, float coeficients[], float x, float precision, int maxIter){
  int degree = coefNum - 1;
  float b,c;

  for(int k = 0; k < maxIter + 1; k++){
    
    //initializing B and C variables
    b = coeficients[degree];
    c = b;
    
    for(int i = degree-1; i > 0; i--){
      b = coeficients[i] + b*x;
      c = b + c*x;
    }
    b = coeficients[0] + b*x;


    //b = f(x) -> checks if x is already a good value for a root
    if(abs(b) <= precision){
      return x;
    }

    //updating X by the iteration function of the original newton method
    x = x - b/c;
  }

  cout << "Max number of iterations reached" << endl;
  return -99999;


}


