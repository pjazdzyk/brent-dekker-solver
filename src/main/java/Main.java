import Solver.BrentSolver;

import java.util.function.DoubleFunction;

public class Main {

    public static void main(String[] args) {

        //SIMPLE USER GUIDE
        /*
        * This is an equation solver which allows finding roots for any type of function within a defined argument range. The solver uses
        * a combination of Secant/Bisection numerical schemes and inverse quadratic interpolation to reduce calculation time as much
        * as possible. The solver algorithm was based on a paper published by Zhengqiu Zhang in 2011, new procedure is simpler and provides
        * faster convergence than a classical Brent-Decker approach. To mitigate the problematic issue with providing valid counterpart points a,b
        * for which f(a) and f(b) will have opposite signs (meaning, that the root must be in the range between these points) - the automatic
        * point evaluation algorithm was proposed, in which case solution will be found for even for invalid points only if they
        * close enough and have a similar order of magnitude as the expected function root. However, this may not work for strongly non-linear
        * functions or if the root is far away from any of the proposed points. The evaluation algorithm is an experimental feature
        * and will be improved in the future. Any suggestions on how to make it better are most welcome.
        *
        * */

        /* Step 1: Brent-Decker Solver instance creation */

        BrentSolver solver = new BrentSolver();

        /* Step 2: Provide equation as double function. Equation must be in form ie: a*x+b = 0. */

        DoubleFunction<Double> func1 = x -> 2*x + 10;             //linear function: 2x + 10 = 0
        DoubleFunction<Double> func2 = x -> 2*x*x + 5*x -3;       //quadratic function: 2x^2 +5x -3 = 0
        DoubleFunction<Double> func3 = x -> 93.3519196629417      //complex function with nested argument under Logarithm
                         - (-237300 * Math.log(0.001638 * x)
                 / (1000 * Math.log(0.001638 * x) - 17269));

        /* Step 3: Provide counterpart points a and b, as a scope, between a solution is expected to be found.
         * Function value for that points must result in opposite signs. If this fails, solver is equipped with
         * a weak evaluation procedure, which will attempt to adjust your points to meet the criteria. If it does not
         * succeed an exception is thrown. */

        // Let's try the first function. Default counterpart points are +50 / -50, lets try with these:
        var result1 = solver.calcForFunction(func1);
        System.out.println("result1= " + result1);   //Outputs -5.0

        // Now the more difficult part. This function has two roots: -3.0, 0.5. Lets he how solver will behave:
        var result2a = solver.calcForFunction(func2);
        System.out.println("result2a= " + result2a); //Outputs -3.0

        /* It is important to understand what happened above. The solver outputs the first root it finds and only one.
        * If the found root is not the one you have expected - you need to adjust the scope in which desired
        * solution can be found. Let's try to adjust counterpart points to get the other root: */

        solver.setCounterpartPoints(-1,2);    // Second root 0.5 is between 2 and 1
        var result2b = solver.calcResult();        // Desired function is already set, so can simply invoke calculation method.
        System.out.println("result2b= " + result2b);      // Outputs 0.5

        /* Last example is to present how test if complex nested argument function can be handled properly.*/
        // var result3a = solver.calcForFunction(func3);
        // System.out.println("result3a= " + result3a);

        /* If you uncomment above, you will get an exception that NAN value was detected, what means that solution is not converged.
        * Expected root for this function is 80000. The counterpart we last set are -1 and 2, which are several orders of magnitude
        * smaller than any expected root. For such differences, automatic evaluation procedure will not manage to overcome it. Therefore,
        * we need to change solution search scope to more appropriate threshold, for an example: 20000 - 200000. Instead of invoking
        * a separate counterpart setter, we can use overloaded version of calcForFunction method: */

        var result3b = solver.calcForFunction(func3,20000,200000);
        System.out.println("result3b= " + result3b); // Outputs: 79999.99999999991, with some small insignificant numerical error

        /* in case the order of magnitude for counterpart points are correct, even if you are not sure if provided points will result
        * in opposite sign of the function value - the automatic evaluation algorithm will adjust it using linear extrapolation and
        * few cycles of tries and error. Is an extremely useful feature extending the usability of this solver. To see each step
        * of calculation you can turn on the diagnostic output. */

        solver.setShowDiagnostics(true);
        var result3c = solver.calcForFunction(func3,100000,200000); // for these points f_a and f_b will not have an opposite signs
        System.out.println("\nresult3c=" + result3c); //Outputs: 79999.99999999991

        /* Despite tha fact that solution was not in the scope of provided counterpart points, solver still managed to converge
        * You can see in the diagnostic output, that Evaluation procedure algorithm was launched at first, and determined that
        * second counterpart point needs to be changed to 74939.747 for solver to work. */

        /* As it was demonstrated, proposed solver can handle almost any equation in a very few iterations steps and manages to
        provide a solution even if counterpart points are invalid (but still they need to be of proper order of magnitude).
        Evaluation procedure is not perfect, and I am looking forward for any improvement proposals. */

        /*
        * Author: MScEng Piotr Jażdżyk
        * Mail: info@synerset.com
        * Web: www.synerset.com
        * LinkedIn: https://www.linkedin.com/in/pjazdzyk
        */

    }

}
