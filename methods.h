//questao 1
float posicaoFalsa(float a, float b, float precisao, float (*funcao)(float), int count = 0);

float bissecao(float, float, float, float (*funcao)(float), int count = 0);

//questao 2
float newtonRaphson(float x0, float precisao, float(*funcao)(float), float (*funcaoDerivada)(float), int count = 0);

float secante(float x0, float x1, float precisao, float (*funcao)(float), float (*funcaoDerivada)(float), int count = 0);

//questao 3
float pontoFixo (float, float (*funcao)(float), float (*funcaoIteracao)(float), float, int count = 0);

float calculoPolinomio (int, float[], float);

float calculoDerivada (int, float[], float);

float polinomioNewton(int, float coeficientes[], float x, float, int maxIter = 50);