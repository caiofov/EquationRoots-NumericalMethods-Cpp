#include <iostream>
#include <cmath>
using namespace std;
#include "methods.h"
#include "functions.h"

float posicaoFalsa(float a, float b, float erro, int count){
  count++;
  float fa = funcao(a); //f(a)
  float fb = funcao(b); //f(b)
  
  cout << "\n" <<count << "ª Iteração \n";

  if( (abs(b-a) > erro) && (fa*fb < 0)){
    float pm = ((a*fb)-(b*fa))/(fb - fa); //ponto medio
    float fp = funcao(pm); //f(ponto medio)

    cout << "\n a = " << a << " | f(a) =" << fa << " | b = " << b << " | f(b) = " << fb << " \n p = " << pm << " | f(p) = " << fp << endl;
    
    if((fp*fa < 0) && (abs(pm - a) > erro) && (abs(fp) > erro)){
      return posicaoFalsa(a, pm, erro, count);
    }
    
    else if((fp * fb < 0) && (abs(b-pm) > erro) && (abs(fp) > erro)){
      return posicaoFalsa(pm, b, erro, count);
    }
    
    else{
      return pm;
    }
  }
 
  else{
    cout << "Erro \n";
    return -999999;
  }
}

float bissecao(float a, float b, float erro, int count){
  count++;
  float fa = funcao(a); //f(a)
  float fb = funcao(b); //f(b)
  cout << "\n" <<count << "ª Iteração \n";

  if(fa == 0){ //checa se é raiz
    cout << "a é raiz" << endl;
    return a;
  }
  if(fb == 0){ //checa se é raiz
    cout << "b é raiz" << endl;
    return b;
  }
  
  if( (abs(b-a) > erro) && (fa*fb < 0)){
    float pm = (a+b)/2; //ponto medio
    float fp = funcao(pm); //f(ponto medio)

    cout << "\n a = " << a << " | f(a) =" << fa << " | b = " << b << " | f(b) = " << fb << " \n p = " << pm << " | f(p) = " << fp << endl;
    
    if((fp*fa < 0) 
    && (abs(pm - a) > erro) 
    && (abs(fp) > erro)){
      return bissecao(a, pm, erro, count);
    }
    
    else if((fp * fb < 0) 
    && (abs(b-pm) > erro) 
    && (abs(fp) > erro)){
      return bissecao(pm, b, erro, count);
    }
    
    else{
      if(fp == 0){ //checa se é raiz
        cout << "Ponto médio é raiz" << endl;
      }
      return pm;
    }
  }
 
  else{
    cout << "Erro \n";
    return -999999;
  }
}



float newtonRaphson(float x0, float erro, int count){
  count++;
  cout << "\n" <<count << "ª Iteração \n";
  float fx0 = funcao(x0);
  float phiX = x0 - (fx0/funcaoDerivada(x0));

  cout << "x0 = "<<x0<<" | F(x0) = "<<fx0<<" | Phi(x0) = " <<phiX << endl;
  
  if(fx0 == 0){
    cout << "x0 é raiz" << endl;
    return x0;}; //checa se é raiz
  
  if(abs(fx0) > erro){
    return newtonRaphson(phiX, erro,count);
  }
  else{
    return x0;
  }

}

float secante(float x0, float x1, float erro, int count){
  count++;
  cout << "\n" <<count << "ª Iteração \n";
  
  float fx0 = funcao(x0);
  float fx1 = funcao(x1);
  
  float phiX2 = ((x0*fx1) - x1*fx0)/(fx1 - fx0);
  cout << "X(k-1) = " << x0 << " | Xk = " << x1 <<"\nF(xk-1) = "<< fx0 << " | F(xk) = " << fx1 << "\nPhi(x2) = " << phiX2 << endl;

  if(fx0 == 0){ //checa se é raiz
    cout << "x0 é raiz" << endl;
    return x0;
    }; 
  if(fx1 == 0){ //checa se é raiz
    cout << "x1 é raiz" << endl;
    return x1;
  } 

  if(abs(fx0) > erro){
    return secante(x1,phiX2,erro,count);
  }

  else{
    return phiX2;
  }


}


