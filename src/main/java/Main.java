import solver.BrentSolver;

import java.util.function.DoubleFunction;

public class Main {

    public static void main(String[] args) {

        BrentSolver solver = new BrentSolver();

        DoubleFunction<Double> func1 = x -> 2*x + 10;
        DoubleFunction<Double> func2 = x -> 2*x*x + 5*x -3;
        DoubleFunction<Double> func3 = x -> 93.3519196629417 - (-237300 * Math.log(0.001638 * x) / (1000 * Math.log(0.001638 * x) - 17269));

        var result1 = solver.calcForFunction(func1);
        System.out.println("result1= " + result1);                 // Outputs -5.0

        var result2a = solver.calcForFunction(func2);
        System.out.println("result2a= " + result2a);               // Outputs -3.0

        solver.setCounterpartPoints(-1,2);           // Second root 0.5 is between 2 and 1
        var result2b = solver.calcResult();                // Desired function is already set, so can simply invoke calculation method.
        System.out.println("result2b= " + result2b);               // Outputs 0.5

        // var result3a = solver.calcForFunction(func3);
        // System.out.println("result3a= " + result3a);

        /* If you uncomment above, you will get an exception that NAN value was detected, what means that solution is not converged.
        * Expected root for this function is 80000. The counterpart we last set are -1 and 2, which are several orders of magnitude
        * smaller than any expected root. For such differences, automatic evaluation procedure will not manage to overcome it. Therefore,
        * we need to change solution search scope to more appropriate threshold, for an example: 20000 - 200000. Instead of invoking
        * a separate counterpart setter, we can use overloaded version of calcForFunction method: */

        var result3b = solver.calcForFunction(func3,20000,200000);
        System.out.println("result3b= " + result3b); // Outputs: 79999.99999999991, with some small insignificant numerical error

        solver.setShowDiagnostics(true);
        var result3c = solver.calcForFunction(func3,100000,200000);         // for these points f_a and f_b will not have an opposite signs
        System.out.println("\nresult3c=" + result3c); //Outputs: 79999.99999999991          // despite invalid points, solver managed to converge

        // Quick instance example
        var result = BrentSolver.ofFunction(x-> 2*x+10, -10,10);
        System.out.println(result);

    }

}
