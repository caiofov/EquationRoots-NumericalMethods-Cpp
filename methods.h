//  ### False Positon method - 
//[a,b] - interval | precision - precision for the method. It's used to stop the iterations | function - the function which we want to find the roots | count - counter for iterations
float falsePosition(float a, float b, float precision, float (*function)(float), int count = 0);


//  ### Bissection method - 
// [a,b] - interval
//precision - precision for the method. It's used to stop the iterations | function - the function which we want to find the roots | count - counter for iterations
float bissection(float a, float b, float precision, float (*function)(float), int count = 0);

// ### Newton-Raphson method
// x0 - the first guess for the root
// __precision__ - precision for the method. It's used to stop the iterations | function - the function which we want to find the roots | functionDerivative - the derivative for the function used | count - counter for iterations
float newtonRaphson(float x0, float precision, float(*function)(float), float (*functionDerivative)(float), int count = 0);

float secant(float x0, float x1, float precision, float (*function)(float), float (*functionDerivative)(float), int count = 0);

//questao 3
float pontoFixo(float, float (*function)(float), float (*functionIteration)(float), float, int count = 0);

float calculoPolinomio (int, float[], float);

float calculateDerivative (int, float[], float);

float polinomioNewton(int, float coeficientes[], float x, float, int maxIter = 50);