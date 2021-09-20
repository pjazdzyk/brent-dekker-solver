## BRENT-DECKER SOLVER

### Numerical solver for finding function roots. An improved algorithm proposed by Zhengqiu Zhang with an experimental counterpart points evaluation procedure.

> AUTHOR: <b>PIOTR JAŻDŻYK</b> <br>
> VERSION: 1.1.0 <br>
> DATE: 2021.09.17<br>
> EMAIL: jazdzyk@gmail.com <br>
> LINKEDIN: https://www.linkedin.com/in/pjazdzyk/ <br>
> COMPANY: SYNERSET | https://www.synerset.com | info@synerset.com

### SHORT DESCRIPTION:
This is the single variable equation solver which allows finding roots for any type of function within a defined function argument range.
The solver algorithm uses a combination of Secant and Bisection numerical schemes and inverse quadratic interpolation to reduce calculation time as much
as possible. The solver algorithm is based on a paper published by Zhengqiu Zhang in 2011. The new procedure is simpler and provides
faster convergence than a classical Brent-Decker approach. To mitigate the problematic issue with providing valid counterpart points <a,b>
for which f(a) and f(b) will have opposite signs (meaning, that the root must be in the range between these points) the automatic
point evaluation algorithm was proposed, in which case solution will be found even for invalid points, but only if they are
close enough to the root and have a similar order of magnitude as the expected function root. However, this may not work properly for strong
non-linearity's or if the root is far away from any of the proposed points. The evaluation algorithm is an experimental feature
and will be subject of further development and improvement. Any suggestions on how to make it better are most welcome.

### LEGAL DISCLAIMER:
This code is my intellectual property. Please have respect for this. You can use it freely in any academic or non-commercial use if you 
properly include the source and author in your references or documentation. For commercial use, please contact me first.

### CURRENT FUNCTIONALITY:
1. Iterative solver based on Brent-Decker numerical scheme for solving any type of single variable equations.
2. Improved convergence speed by implementing algorithm modifications proposed by Zhengqiu Zhang.
3. Added automatic counterpart evaluation algorithm to make solver less susceptible to invalid initial guess.

### IMPROVEMENT REQUEST:
1. Code review and comments, suggestions how to make this better.
2. Any convergence speedup proposals.
3. Modifications, improvements or new ideas for counterpart point evaluation procedure.

### USER GUIDE:
STEP 1: Brent-Decker Solver instance creation:<br>
```BrentSolver solver = new BrentSolver();```

STEP 2: Provide equation as double function. Equation must be in form ie: a*x+b = 0:<br>
* I. Linear function: 2x + 10 = 0 <br>
```DoubleFunction<Double> func1 = x -> 2*x + 10;```<br>
* II. Quadratic function: 2x^2 +5x -3 = 0 <br>
```DoubleFunction<Double> func2 = x -> 2*x*x + 5*x -3;```      
* III. Complex function with nested argument under Logarithm <br>
```DoubleFunction<Double> func3 = x -> 93.3519196629417 - (-237300 * Math.log(0.001638 * x) / (1000 * Math.log(0.001638 * x) - 17269));```  

STEP 3: Provide counterpart points a and b, as a scope, between a solution is expected to be found. 
Function value for that points must result in opposite signs. If this fails, solver is equipped with the weak evaluation procedure, 
which will attempt to adjust your points to meet the criteria. If it does not succeed an exception is thrown.
Let's try the first function. Default counterpart points are +50 / -50, lets try with these:
```
var result1 = solver.calcForFunction(func1);
System.out.println("result1= " + result1);   //Outputs -5.0 
```
Solver outputs -5.0 what is the expected solution.<br>
Now the more difficult part, the quadratic function. This function has two roots: -3.0, 0.5. Lets he how solver will behave:
```
var result2a = solver.calcForFunction(func2);
System.out.println("result2a= " + result2a); //Outputs -3.0
```
It is important to understand what happened above. The solver outputs the first root it finds and only one.
If the found root is not the one you have expected - you may need to adjust the scope in which desired solution can be found. 
Let's try to adjust counterpart points to get the other root: */
```
solver.setCounterpartPoints(-1,2);              // We can expect that second root is between 2 and 1
var result2b = solver.calcResult();             // Desired function is already set, so let's simply invoke calculation method.
System.out.println("result2b= " + result2b);    // Outputs 0.5, the second root
```
Last example is to test if solver will manage to find root for a complex coupled argument, for an example inside a logarithm:
```
var result3a = solver.calcForFunction(func3);
System.out.println("result3a= " + result3a);
```
If you run the code above, you will get an exception with information that NAN value was detected, what means that solution is not converged.
Expected root for this function is a value of 80 000. The counterpart we last set are -1 and 2, which are several orders of magnitude
smaller than any expected root. For such differences, automatic evaluation procedure will not manage to overcome it. Therefore,
we need to change solution range scope to more appropriate threshold, for an example: 20 000 - 200 000. Instead of invoking
a separate counterpart setter, you can use overloaded version of calcForFunction method:
```
var result3b = solver.calcForFunction(func3,20000,200000);
System.out.println("result3b= " + result3b);    // Outputs: 79999.99999999991, with some small omittable numerical error
```
Even if the provided counterpart points are invalid as shown above - the automatic evaluation algorithm will adjust them properly using linear extrapolation and
few cycles of tries and error approach. This will work as long provided points are closed to the root and are of the similar order of magnitude.
Is an extremely useful feature extending the usability of this solver. It may seem a bit of defensive approach, but my experience says that in many cases
it will not be possible to automatically provide correct points each time.
To see each step calculation you can turn on the diagnostic output. */
```
solver.setShowDiagnostics(true);
var result3c = solver.calcForFunction(func3,100000,200000);     // for these points f_a and f_b will not have an opposite signs
System.out.println("\nresult3c=" + result3c); //Outputs: 79999.99999999991
```
Despite tha fact that solution was not in the scope of provided counterpart points, solver still managed to converge.
You can observe the diagnostic output and see that evaluation procedure algorithm was launched at the beginning and determined
that second counterpart point needs to be changed to 74939.747 for solver to work. 

As it was demonstrated, this solver can handle almost any equation in a very few iterations steps and manages to
provide a solution even if counterpart points are invalid (but still they need to be of proper order of magnitude).
Evaluation procedure is not bulletproof, I am looking forward for any improvement proposals from the community. 

### REFERENCE SOURCE:
* [1] BRENT-DECKER ITERATIVE SOLVER - MODIFIED ALGORITHM PROPOSED BY Zhengqiu Zhang / International Journal of Experimental Algorithms (IJEA), Volume (2) : Issue (1) : 2011
