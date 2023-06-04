## BRENT-DEKKER SOLVER

---

### Numerical solver for finding function roots based on the improved algorithm proposed by Zhengqiu Zhang with new experimental counterpart points evaluation procedure.

> AUTHOR: PIOTR JAŻDŻYK <br>
> LINKEDIN: https://www.linkedin.com/in/pjazdzyk <br>

[![Maven Central](https://maven-badges.herokuapp.com/maven-central/com.synerset/brent-dekker-solver/badge.svg)](https://maven-badges.herokuapp.com/maven-central/com.synerset/brent-dekker-solver) &nbsp;
[![Vulnerabilities](https://sonarcloud.io/api/project_badges/measure?project=pjazdzyk_unitility&metric=vulnerabilities)](https://sonarcloud.io/summary/new_code?id=pjazdzyk_brent-dekker-solver) &nbsp;
[![Security Rating](https://sonarcloud.io/api/project_badges/measure?project=pjazdzyk_unitility&metric=security_rating)](https://sonarcloud.io/summary/new_code?id=pjazdzyk_brent-dekker-solver) &nbsp;
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=pjazdzyk_unitility&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=pjazdzyk_brent-dekker-solver) &nbsp;

### SHORT DESCRIPTION
This is a Java implementation of a math algorithm for finding roots of continuous and single-variable functions, based on a modified version of the Brent-Dekker algorithm.
This solver uses a combination of Secant and Bisection numerical schemes and inverse quadratic interpolation to reduce calculation time as much as possible.
The algorithm code is based on the paper published by Zhengqiu Zhang in 2011. The new procedure has lower complexity and provides faster convergence than a classical
Brent-Dekker algorithm.
To mitigate the problematic issue by providing valid counterpart points <a,b> for which f(a) and f(b) will have opposite signs (meaning, that the root must be in the range between these points) the
automatic point evaluation algorithm was proposed, in which case solution will be found even for invalid points, but only if they
are close enough to the root and have a similar order of magnitude as the expected function root. However, this may not work
correctly for strong non-linearity or if the root is far away from any of the proposed points. The evaluation algorithm is an experimental
feature and will be subject to further development and improvement. Any suggestions on how to make it better are most welcome.

### TECH
Core: <br>
![image](https://img.shields.io/badge/17-Java-orange?style=for-the-badge) &nbsp;
![image](https://img.shields.io/badge/apache_maven-C71A36?style=for-the-badge&logo=apachemaven&logoColor=white) &nbsp;
![image](https://img.shields.io/badge/apache_LOG4J-C71A36?style=for-the-badge&logo=apachemaven&logoColor=white) &nbsp;

Testing:<br>
![image](https://img.shields.io/badge/Junit5-25A162?style=for-the-badge&logo=junit5&logoColor=white) &nbsp;

CI/CD:<br>
![image](https://img.shields.io/badge/GitHub_Actions-2088FF?style=for-the-badge&logo=github-actions&logoColor=white) &nbsp;
![image](https://img.shields.io/badge/Sonar%20cloud-F3702A?style=for-the-badge&logo=sonarcloud&logoColor=white) &nbsp;

### INSTALLATION

Just copy Maven dependency provided below to your pom.xml file, and you are ready to go. For other package managers, check maven central repository: 
[Brent-Dekker-Solver](https://search.maven.org/artifact/com.synerset/brent-dekker-solver/1.1.3/jar).

```xml
<dependency>
  <groupId>com.synerset</groupId>
  <artifactId>brent-dekker-solver</artifactId>
  <version>1.1.3</version>
</dependency>
```

### LEGAL DISCLAIMER

This code is my intellectual property. Please have respect for this. You can use it freely in any academic or
non-commercial use if you
properly include the source and author in your references or documentation. For commercial use, please contact me first.

### CURRENT FUNCTIONALITY

1. Iterative solver based on Brent-Dekker numerical scheme for solving any type of single variable equations.
2. Improved convergence speed by implementing algorithm modifications proposed by Zhengqiu Zhang in his paper.
3. Automatic counterpart evaluation algorithm to make solver less susceptible to invalid initial guess.

### COLLABORATION

1. Code review and comments, suggestions how to make this better.
2. Convergence speedup and complexity reduction ideas.
3. Counterpart points evaluation procedure improvement.

### USER GUIDE WITH EXAMPLES

STEP 1: Brent-Dekker solver instance creation:<br>
```java
BrentSolver solver = new BrentSolver();
```

STEP 2: Provide equation as double function. Equation must be in form ie: a * x + b = 0:<br>

* I. Linear function: 2x + 10 = 0 <br>

```java
DoubleFunction<Double> linear = x -> 2 * x + 10;
```

* II. Quadratic function: 2x^2 +5x -3 = 0
```java
DoubleFunction<Double> quadratic = x -> 2 * x * x + 5 * x -3;
```

* III. Complex function with nested argument under Logarithm <br>
```java
DoubleFunction<Double> logNested = x -> 93.3519196629417 - (-237300 * Math.log(0.001638 * x) / (1000 * Math.log(0.001638 * x) - 17269));
```

STEP 3: Provide counterpart points a and b, as a scope, between a solution is expected to be found.
Function value for that points must result in opposite signs. If this fails, solver is equipped with the weak evaluation
procedure,
which will attempt to adjust your points to meet the criteria. If it does not succeed an exception is thrown.
Let's try the first function. Default counterpart points are +50 / -50, lets try with these:

```java
var resultLinear = solver.calcForFunction(linear);
System.out.println("Linear function root = " + resultLinear);    
    //Outputs -5.0 
```

solver outputs -5.0 what is the expected solution.<br>
Now the more difficult part, the quadratic function. This function has two roots: -3.0, 0.5. Lets he how solver will
behave:

```java
var resultQuadratic1stRoot = solver.calcForFunction(quadratic);
System.out.println("Quadratic function 1st root = " + resultQuadratic1stRoot);  
    //Outputs -3.0
```

It is important to understand what happened above. The solver outputs the first root it finds and only one.
If the found root is not the one you have expected - you may need to adjust the scope in which the desired solution can
be found.
Let's try to adjust counterpart points to get the other root:

```java
solver.setCounterpartPoints(-1,2);                                              
    // We can expect that second root is between 2 and 1
var resultQuadratic2ndRoot = solver.findRoot();                                 
    // Desired function is already set, so can simply invoke calculation method.
System.out.println("Quadratic function 2nd root = " + resultQuadratic2ndRoot);  
    // Outputs 0.5, the second root
```

Last example is to test if solver will manage to find root for a complex coupled argument, for an example inside a
logarithm:

```java
var resultLogNested = solver.findRoot();
System.out.println("Nested log function = " + resultLogNested);
```

If you run the code above, you will get an exception with information that NAN value was detected, which means that
solution is not converged.
The expected root for this function is a value of 80 000. The counterpart we last set are -1 and 2, which are several
orders of magnitude
smaller than any expected root. For such differences, automatic evaluation procedures will not manage to overcome it.
Therefore,
we need to change the solution range scope to a more appropriate threshold, for example, 20 000 - 200 000. Instead of
invoking
a separate counterpart setter, you can use an overloaded version of calcForFunction method:

```java
var resultLogNested = solver.calcForFunction(logNested,20000,200000);
System.out.println("Nested log function root = " + resultLogNested);   
    // Outputs: 79999.99999999991, with some small numerical error
```

Even if the provided counterpart points are invalid as shown above - the automatic evaluation algorithm will adjust them
properly using linear extrapolation and
few cycles of tries and error approach. This will work as long provided points are closed to the root and are of the
similar order of magnitude.
Is an extremely useful feature extending the usability of this solver. It may seem a bit of defensive approach, but my
experience says that in many cases
it will not be possible to automatically provide correct points each time.
To see each step calculation you can turn on the diagnostic output.

```java
solver.setShowDiagnostics(true);                                
resultLogNested = solver.calcForFunction(logNested,100000,200000);        
    // for these points f_a and f_b will not have an opposite signs                        
System.out.println("\nNested log function root = "  + resultLogNested);   
    // despite invalid points, AB evaluation algorithms works fine
```

Despite the fact that the solution was not in the scope of provided counterpart points, the solver still managed to
converge.
You can observe the diagnostic output and see that the evaluation procedure algorithm was launched at the beginning and
determined
that second counterpart point needs to be changed to 74939.747 for the solver to work.

For single result, you can quickly use static method ofFunction() and get result immediately.

```java
var functionRoot = BrentSolver.ofFunction(x-> 2*x+10, -10,10);
System.out.println(functionRoot);    
    //Outputs: -5.0
```

As it was demonstrated, this solver can handle almost any equation in a very few iterations steps and manages to
provide a solution even if counterpart points are invalid (but still they need to be of the proper order of magnitude).
The evaluation procedure is not bulletproof, I am looking forward to any improvement proposals from the community.

### REFERENCE SOURCE

* [1] BRENT-DEKKER ITERATIVE SOLVER - MODIFIED ALGORITHM PROPOSED BY Zhengqiu Zhang / International Journal of
  Experimental Algorithms (IJEA), Volume (2) : Issue (1) : 2011

## Acknowledgments
* [Shields.io](https://img.shields.io)
* [Badges 4 README.md](https://github.com/alexandresanlim/Badges4-README.md-Profile)
