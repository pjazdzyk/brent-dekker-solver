package com.github.pjazdzyk;

import com.github.pjazdzyk.solver.BrentSolver;

import java.util.function.DoubleFunction;

public class Examples {

    public static void main(String[] args) {
        runUserGuideExamples();
    }

    public static void runUserGuideExamples() {
        BrentSolver solver = new BrentSolver();

        // Tested functions in this presentation
        DoubleFunction<Double> linear = x -> 2 * x + 10;
        DoubleFunction<Double> quadratic = x -> 2 * x * x + 5 * x - 3;
        DoubleFunction<Double> logNested = x -> 93.3519196629417 - (-237300 * Math.log(0.001638 * x) / (1000 * Math.log(0.001638 * x) - 17269));

        var resultLinear = solver.calcForFunction(linear);
        System.out.println("Linear function root = " + resultLinear);                      // Outputs -5.0

        var resultQuadratic1stRoot = solver.calcForFunction(quadratic);
        System.out.println("Quadratic function 1st root = " + resultQuadratic1stRoot);     // Outputs -3.0

        solver.setCounterpartPoints(-1, 2);                                   // Second root 0.5 is between 2 and 1
        var resultQuadratic2ndRoot = solver.findRoot();                            // Desired function is already set, so can simply invoke calculation method.
        System.out.println("Quadratic function 2nd root = " + resultQuadratic2ndRoot);     // Outputs 0.5

        // Counterpart point significance:

        // var resultLogNested = solver.findRoot();
        // System.out.println("Nested log function = " + resultLogNested);

        /* If you uncomment above, you will get an exception that NAN value was detected, what means that solution is not converged.
         * Expected root for this function is 80000. The counterpart we last set are -1 and 2, which are several orders of magnitude
         * smaller than any expected root. For such differences, automatic evaluation procedure will not manage to overcome it. Therefore,
         * we need to change solution search scope to more appropriate threshold, for an example: 20000 - 200000. Instead of invoking
         * a separate counterpart setter, we can use overloaded version of calcForFunction method: */

        var resultLogNested = solver.calcForFunction(logNested, 20000, 200000);
        System.out.println("Nested log function root = " + resultLogNested);                                    // Outputs: 79999.99999999991, with some small numerical error

        solver.setShowDiagnostics(true);
        resultLogNested = solver.calcForFunction(logNested, 100000, 200000);                              // for these points f_a and f_b will not have an opposite signs
        System.out.println("\nNested log function root = " + resultLogNested); //Outputs: 79999.99999999991    // despite invalid points, AB evaluation algorithms works fine

        // Quick instance example
        var functionRoot = BrentSolver.ofFunction(x -> 2 * x + 10, -10, 10);
        System.out.println(functionRoot);   //Outputs: -5.0

    }

}

