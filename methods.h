/*  ### False Positon method - 
- __[a,b]__ - interval
- __ precision__ - precision for the method. It's used to stop the iterations
- __function__ - the function which we want to find the roots
- __count__ - counter for iterations (default 0)*/
float falsePosition(float a, float b, float precision, float (*function)(float), int count = 0);


/*  ### Bissection method - 
 - __[a,b]__ - interval
- __precision__ - precision for the method. It's used to stop the iterations
- __function__ - the function which we want to find the roots
- __count__ - counter for iterations (default 0)*/
float bissection(float a, float b, float precision, float (*function)(float), int count = 0);

/* ### Newton-Raphson method
- __x0__ - the first guess for the root
- __precision__ - precision for the method. It's used to stop the iterations 
- __function__ - the function which we want to find the roots
- __functionDerivative__ - the derivative for the function used count - counter for iterations (default 0)
- __count__ - counter for iterations (default 0) */
float newtonRaphson(float x0, float precision, float(*function)(float), float (*functionDerivative)(float), int count = 0);

float secant(float x0, float x1, float precision, float (*function)(float), float (*functionDerivative)(float), int count = 0);

//questao 3
float pontoFixo(float, float (*function)(float), float (*functionIteration)(float), float, int count = 0);

float calculoPolinomio (int, float[], float);

/*Returns the derivative of a function for a specific value.
- __coefNum__ - Number of coeficients of the given function
- __coeficients__ - An array with all coeficients (including 0's)
- __x__ - the specific value for calculating the function
*/
float calculateDerivative (int coefNum, float coeficients[], float x);

float polinomioNewton(int, float coeficientes[], float x, float, int maxIter = 50);