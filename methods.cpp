#include <iostream>
#include <cmath>
using namespace std;
#include "methods.h"

float posicaoFalsa(float a, float b, float precisao, float(*funcao)(float), int count){
  count++;
  float fa = funcao(a); //f(a)
  float fb = funcao(b); //f(b)
  
  printf("\n %dª Iteração \n", count);

  if( (abs(b-a) > precisao) && (fa*fb < 0)){
    float pm = ((a*fb)-(b*fa))/(fb - fa); //ponto medio
    float fp = funcao(pm); //f(ponto medio)
    
    cout << "\n a = " << a << " | f(a) =" << fa << " | b = " << b << " | f(b) = " << fb << " \n p = " << pm << " | f(p) = " << fp << endl;
    
    if((fp*fa < 0) && (abs(pm - a) > precisao) && (abs(fp) > precisao)){
      return posicaoFalsa(a, pm, precisao, funcao, count);
    }
    
    else if((fp * fb < 0) && (abs(b-pm) > precisao) && (abs(fp) > precisao)){
      return posicaoFalsa(pm, b, precisao, funcao, count);
    }
    
    else{
      return pm;
    }
  }
 
  else{
    printf("precisao \n");
    return -999999;
  }
}

float bissecao(float a, float b, float precisao, float(*funcao)(float), int count){
  count++;
  float fa = funcao(a); //f(a)
  float fb = funcao(b); //f(b)
  printf("\n %dª Iteração \n", count);

  if(fa == 0){ //checa se é raiz
    printf("a é raiz\n");
    return a;
  }
  if(fb == 0){ //checa se é raiz
    printf("b é raiz\n");
    return b;
  }
  
  if( (abs(b-a) > precisao) && (fa*fb < 0)){
    float pm = (a+b)/2; //ponto medio
    float fp = funcao(pm); //f(ponto medio)

    cout << "\n a = " << a << " | f(a) =" << fa << " | b = " << b << " | f(b) = " << fb << " \n p = " << pm << " | f(p) = " << fp << endl;
    
    if((fp*fa < 0) 
    && (abs(pm - a) > precisao) 
    && (abs(fp) > precisao)){
      return bissecao(a, pm, precisao, funcao, count);
    }
    
    else if((fp * fb < 0) 
    && (abs(b-pm) > precisao) 
    && (abs(fp) > precisao)){
      return bissecao(pm, b, precisao, funcao, count);
    }
    
    else{
      if(fp == 0){ //checa se é raiz
        printf("Ponto médio é raiz");
      }
      return pm;
    }
  }
 
  else{
    printf("Erro \n");
    return -999999;
  }
}

float newtonRaphson(float x0, float precisao, float (*funcao)(float), float (*funcaoDerivada)(float), int count){
  count++;
  printf("\n%dª Iteração \n", count);
  float fx0 = funcao(x0);
  float phiX = x0 - (fx0/funcaoDerivada(x0));

  cout << "x0 = "<<x0<<" | F(x0) = "<<fx0<<" | Phi(x0) = " <<phiX << endl;
  
  if(fx0 == 0){
    cout << "x0 é raiz" << endl;
    return x0;}; //checa se é raiz
  
  if(abs(fx0) > precisao){
    return newtonRaphson(phiX, precisao, funcao, funcaoDerivada, count);
  }
  else{
    return x0;
  }

}

float secante(float x0, float x1, float precisao, float (*funcao)(float), float (*funcaoDerivada)(float), int count){
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

  if(abs(fx0) > precisao){
    return secante(x1,phiX2,precisao, funcao, funcaoDerivada, count);
  }

  else{
    return phiX2;
  }


}

float pontoFixo(float x0, float (*funcao)(float), float (*funcaoIteracao)(float), float precisao, int count){
  count++;
  printf("\n %dª Iteração \n", count);

  float phix0 = funcaoIteracao(x0);
  float fx0 =  funcao(x0);
  cout << "x0 = " << x0 << " | F(x0) = " << fx0 << "\nPhi(x0) = " << phix0 << endl;
  if(fx0 == 0){ //checa se x0 é raiz
    cout << "x0 é raiz" << endl;
    return x0;
  }

  if( abs(fx0) > precisao){
    
    return pontoFixo(phix0, funcao, funcaoIteracao, precisao, count);
  }

  if(isinf(x0) || isnan(x0)){
    cout << "Tente outra função de iteração" << endl;
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

float calculoDerivada(int numCoef, float coeficientes[numCoef], float x){
  
  int grau = numCoef - 1;
  float z = x;

  float b = coeficientes[grau];
  float c = b;
  
  for(int i = grau-1; i > 0; i--){
    b = coeficientes[i] + b*z;
    c = b + c*z;
  }
  //b = coeficientes[0] + b*z;

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


