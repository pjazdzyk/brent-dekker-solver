package com.synerset.brentsolver;

import java.util.function.DoubleUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * BRENT-DEKKER ITERATIVE SOLVER - MODIFIED ALGORITHM PROPOSED BY Zhengqiu Zhang / International Journal of Experimental<br>
 * Algorithms (IJEA), Volume (2) : Issue (1) : 2011
 * <br>
 * Single variable equation solver for finding roots for any type of function within defined argument range.
 * The algorithm uses a combination of following numerical schemes: Secant method, Bisection method and inverse quadratic interpolation to reduce calculation steps as much as possible.
 * Implemented improvements are based on a scientific paper published by Zhengqiu Zhang in 2011.
 * User function must be rearranged to the form of [expression = 0].
 * Variable name convention is intended to be as much similar as possible to the actual equations from math textbooks or in this case scientific papers
 * which was the reference source to construct this algorithm. Variables with zero (0) suffix stands for initial guess values. All variables starting with f_x are an output value from the user function
 * in point "x". Variables "a" and "b" stands for already checked and adjusted values which are put into the root-finding algorithm. Values of "c" stores previous iteration result.
 * Value of "s" stores the result of an applied interpolation algorithm. Root is found when function value in for the "b" point outputs 0 or is smaller than specified solver accuracy.
 * <br>
 * <p>AUTHOR: Piotr Jażdżyk, MScEng</p>
 * CONTACT:
 * <a href="https://www.linkedin.com/in/pjazdzyk">LinkedIn</a> |
 * <a href="http://synerset.com/">www.synerset.com</a><br>
 */

public class BrentSolver {

    // IDENTITY DATA
    private final String id;

    // POINTS AND THEIR FUNCTION VALUES
    private double a0 = -50;                                    // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    private double b0 = 50;                                     // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    private double a, f_a;                                      // first point and its function value f(a), must have opposite sign than b
    private double b, f_b;                                      // second point and its function value f(b), absolute value of b should be lower than a
    private double c, f_c;                                      // previous iteration point and f(c), initially = a
    private double s, f_s;                                      // new approximate value based on secant method or bisection method and f_s
    private int counter;
    private DoubleUnaryOperator userFunction;                   // function to compute, to be provided by user

    // SOLVER SOLUTION CONTROL
    private boolean runFlag = true;
    private double iterationsLimit = 100;
    private double accuracy = 0.00001;

    // AB CONDITION EVAL ALGORITHM CONTROL
    private int evalCycles = 2;
    private int evalX2Divider = 2;
    private int evalXDivider = 2;

    // DIAGNOSTICS OUTPUT CONTROL
    private boolean showDiagnostics = false;

    private static final Logger LOGGER = Logger.getLogger(BrentSolver.class.getName());

    /**
     * Initializes solver instance with function output set as 0, with default name.
     */
    public BrentSolver() {
        this("DefaultSolver");
    }

    /**
     * Initializes solver instance with function output set as 0.
     */
    public BrentSolver(String id) {
        this.id = id;
        this.userFunction = val -> 0.0;
    }

    /**
     * Initializes solver instance with function output set as 0 and sets custom p2 and p3 coefficients.
     * P2 and P3 coefficients are used to tune the counterpart point evaluation algorithm. Default values of p2=2, p3=2
     * works properly for a typical temperature range (-100,100). For pressures variable range, it is recommended to use p2=2, p3=0.
     * For other cases, if the use of point evaluation procedure is expected - p2 and p3 points have to be determined empirically.
     *
     * @param evalX2Divider AB evaluation procedure second point coefficient
     * @param evalXDivider  AB evaluation procedure third point coefficient
     */
    public BrentSolver(String id, int evalX2Divider, int evalXDivider) {
        this(id);
        setEvalX2Divider(evalX2Divider);
        setEvalXDivider(evalXDivider);
    }

    /**
     * Initializes solver instance with function provided by user.
     *
     * @param userFunction tested function (use lambda expression or method reference)
     */
    public BrentSolver(String id, DoubleUnaryOperator userFunction, int evalX2Divider, int evalXDivider) {
        this(id, evalX2Divider, evalXDivider);
        this.userFunction = userFunction;
    }

    /**
     * Returns root value based on Brent-Dekker numerical scheme. Calculation procedure will use
     * expressions provided in userFunction.
     *
     * @return actual root value
     */
    public final double findRoot() {
        if (showDiagnostics) {
            LOGGER.log(Level.INFO, "CALCULATION IN PROGRESS: {}", id);
        }
        // To check and set value b as being closer to the root
        checkSetAndSwapABPoints(a0, b0);
        // In case provided by user point "a" or "b" is actually a root
        if (Math.abs(f_b) < accuracy) {
            if (showDiagnostics) {
                LOGGER.log(Level.INFO, "CALCULATION COMPLETE. Initial value is a root.");
            }
            return b;
        }
        // If solver were stopped
        if (!runFlag) {
            if (showDiagnostics) {
                LOGGER.log(Level.WARNING, "CALCULATION STOPPED: {}", id);
            }
            return b;
        }
        // Checking if Brent AB condition is not met to launch automatic AB points evaluation procedure
        if (initialABConditionIsNotMet()) {
            evaluateValidCondition();
        }
        // In case evaluation procedure will output b as a root
        if (Math.abs(f_b) < accuracy) {
            return b;
        }
        // If at this stage proper A-B condition is not achievable - an exception is thrown.
        if (initialABConditionIsNotMet()) {
            String errorMsg = id + ": EVALUATION PROCEDURE FAILED: f(a) i f(b) must have an opposite signs. Current values:"
                              + String.format(" a = %.3f, b = %.3f,  f(a)= %.3f, f(b)=%.3f", a, b, f_a, f_b);
            LOGGER.log(Level.SEVERE, errorMsg);
            throw new BrentSolverException(errorMsg);
        }
        printSolverResultsDiagnostics(id + ": INITIAL VALUES: ", "");

        /*--------BEGINNING OF ITERATIVE LOOP--------*/
        while (runFlag) {
            counter++;

            // New additional condition proposed by Zhengqiu Zhang
            c = (a + b) / 2;
            f_c = userFunction.applyAsDouble(c);

            //Determining better interpolation: inverse quadratic or else - secant method
            if ((f_a != f_c) && (f_b != f_c))
                s = inverseQuadraticInterpolation(a, b, c, f_a, f_b, f_c);
            else
                s = secantMethod(a, b, f_a, f_b);

            // Difference between a and b, to be checked against the target (initialized by default value)
            double currentDifference = Math.abs(b - a);
            if (c > s) {
                var tempC = c;
                var tempS = s;
                c = tempS;
                s = tempC;
            }
            f_c = userFunction.applyAsDouble(c);
            f_s = userFunction.applyAsDouble(s);
            if (f_c * f_s < 0) {
                a = s;
                b = c;
            } else {
                if (f_s * f_b < 0) {
                    a = c;
                } else {
                    b = s;
                }
            }
            f_a = userFunction.applyAsDouble(a);
            f_b = userFunction.applyAsDouble(b);

            // Calculating current difference after this iteration cycle
            printSolverResultsDiagnostics(id + ": ITERATION: " + counter + " ", "Diff= " + currentDifference);
            if (currentDifference < accuracy) {
                runFlag = false;
            } else if (counter > iterationsLimit) {
                runFlag = false;
            } else if (f_b == 0) {
                runFlag = false;
            }

            // Exception will be thrown if NaN or Infinite values are detected
            checkForInfiniteOrNaN(f_a, f_b, f_c, f_s);

            /*-----------END OF ITERATIVE LOOP-----------*/
        }

        // b is always closer to the root
        if (showDiagnostics) {
            LOGGER.log(Level.INFO, "CALCULATIONS FINISHED: {}", id);
        }
        return b;
    }

    private boolean initialABConditionIsNotMet() {
        return (f_a * f_b) >= 0;
    }

    private void checkSetAndSwapABPoints(double pointA, double pointB) {
        f_a = userFunction.applyAsDouble(pointA);
        f_b = userFunction.applyAsDouble(pointB);
        if (Math.abs(f_a) < Math.abs(f_b)) {
            a = pointB;
            b = pointA;
            f_a = userFunction.applyAsDouble(a);
            f_b = userFunction.applyAsDouble(b);
        } else {
            a = pointA;
            b = pointB;
        }
    }

    private void checkForInfiniteOrNaN(double... values) {
        for (double num : values) {
            if (Double.isInfinite(num))
                throw new BrentSolverException(id + ": Solution error. Infinite number detected.");
            if (Double.isNaN(num))
                throw new BrentSolverException(id + ": Solution error. NaN value detected.");
        }
    }

    private void evaluateValidCondition() {
        /*This method attempts to evaluate valid Brent Solver counterpart point condition.
        PROCEDURE EXPLANATION
        Using linear extrapolation to determine opposite sign of the "b" value.
        1. Linear extrapolations needs a pair of points: P1(x1,f_x1), P2(x2,f_x2) and f_y to determine x.
        2. First point is already provided (b, f_b) and it is already evaluated to be closer to the root
        3. Second point P2 is created from b/evalX2Divider and its function value.
        4. f_x is an x value we expect to have an opposite sign to -b, and to make it closer to the root - it is divided by evalXDivider
        5. In some cases initial values of evalX2Divider and evalXDivider may be adjusted by user. For an example if your initial guess is very close to the root
        you are looking for small evalXDivider values, like 2 or even 1. */

        printEvaluationDiagnostics(id + " EVALUATION PROCEDURE \nINITIAL:");
        double x, f_x, x1, f_x1, x2, f_x2, f_xExact;
        if (evalXDivider > evalCycles)
            evalCycles = evalXDivider;
        for (int i = 0; i <= evalCycles; i++) {
            if (evalXDivider - i == 0.0)
                continue;
            //Point 1, using b as a starting point
            x1 = b;
            f_x1 = f_b;
            //Point 2, creating second point from the b
            x2 = b / evalX2Divider;
            f_x2 = userFunction.applyAsDouble(x2);
            //Point 3, searching for a negative value of -f_b. Further division is to get result as close to the root as possible.
            f_x = -f_b / (evalXDivider - i);
            x = linearInterpolationFromValue(x1, f_x1, x2, f_x2, f_x);
            // When x is determined - to check if it really gives negative value (it may not occur for strong non-linearity)
            f_xExact = userFunction.applyAsDouble(x);
            checkSetAndSwapABPoints(b, x);
            printEvaluationDiagnostics("STEP " + i + ":");
            if (f_xExact * f_x1 < 0)
                break;
        }
    }

    public final void stopSolver() {
        this.runFlag = false;
    }

    /**
     * Sets custom initial counterpart points for Brent-Dekker method condition.
     *
     * @param pointA - first initial counterpart point
     * @param pointB - second initial counterpart point
     */
    public final void setCounterpartPoints(double pointA, double pointB) {
        this.a0 = pointA;
        this.b0 = pointB;
        resetSolverRunFlags();
    }

    /**
     * Resets flags and iteration counter
     */
    public final void resetSolverRunFlags() {
        this.runFlag = true;
        this.counter = 0;
    }

    /**
     * Resets counter part points to default values (+50,-50);
     */
    public final void resetCounterPartPoints() {
        a0 = -50;
        b0 = 50;
    }

    /**
     * Resets evaluation procedure coefficients to default values.
     */
    public final void resetEvaluationCoefficients() {
        evalCycles = 5;
        evalX2Divider = 2;
        evalXDivider = 10;
    }

    /**
     * Sets function to be solved.
     *
     * @param func tested function (use lambda expression or method reference)
     */
    public void setFunction(DoubleUnaryOperator func) {
        this.userFunction = func;
        resetSolverRunFlags();
    }

    /**
     * Returns a function root (calculation result) based on provided user-function.
     *
     * @param func tested function (use lambda expression or method reference)
     * @return function root (Double)
     */
    public double calcForFunction(DoubleUnaryOperator func) {
        setFunction(func);
        return findRoot();
    }

    /**
     * Returns a function root (calculation result) based on provided user-function and counterpart points.
     *
     * @param func tested function (use lambda expression or method reference)
     * @param a0   first point
     * @param b0   second point
     * @return function root (Double)
     */
    public double calcForFunction(DoubleUnaryOperator func, double a0, double b0) {
        setCounterpartPoints(a0, b0);
        return calcForFunction(func);
    }

    /**
     * Sets diagnostic output mode. (true = show output, false = no output).
     *
     * @param showDiagnostics true = on, false = off
     */
    public void setShowDiagnostics(boolean showDiagnostics) {
        this.showDiagnostics = showDiagnostics;
    }

    // INTERPOLATION ALGORITHMS

    /**
     * Linear interpolation/extrapolation. Returns a function value f_x for provided argument of x, based on provided two points P1(x1,f_x1), P2(x2,f_x2).
     * Can be used for interpolation or extrapolation for linear functions.
     *
     * @param x1   first point x value
     * @param f_x1 function value in x1
     * @param x2   second point x value
     * @param f_x2 function value in x2
     * @param x    external point x value
     * @return f_x value of function for argument x.
     */
    public static double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x) {
        return f_x1 + ((x - x1) / (x2 - x1)) * (f_x2 - f_x1);
    }

    /**
     * Linear interpolation/extrapolation. Returns a function argument x for provided argument value f_x, based on provided two points P1(x1,f_x1), P2(x2,f_x2).
     * Can be used for interpolation or extrapolation for linear functions.
     *
     * @param x1   first point x value
     * @param f_x1 function value in x1
     * @param x2   second point x value
     * @param f_x2 function value in x2
     * @param f_x  function value for external point x
     * @return x external point x value
     */
    public static double linearInterpolationFromValue(double x1, double f_x1, double x2, double f_x2, double f_x) {
        return ((x1 - x2) / (f_x1 - f_x2)) * (f_x - f_x1 + x1 * (f_x1 - f_x2) / (x1 - x2));
    }

    /**
     * Inverse quadratic interpolation. Returns a function argument x for f_x = 0.0, based on three provided points P2(x2,f_x2), P1(x1,f_x1), PN(xn,f_xn)
     * It attempts to fit y-based parabola to intersect with axis X as potential function root. It is faster than secant method but more sensitive to initial guesses.
     * It should be used for points very close to the root.
     *
     * @param x2   second point x value
     * @param x1   first point x value
     * @param xn   any other point x value
     * @param f_x2 function value in x2
     * @param f_x1 function value in x1
     * @param f_xn function value in xn
     * @return estimated function root x as argument for f_x = 0.
     */
    public static double inverseQuadraticInterpolation(double x2, double x1, double xn, double f_x2, double f_x1, double f_xn) {
        return x2 * f_x1 * f_xn / ((f_x2 - f_x1) * (f_x2 - f_xn))
               + x1 * f_x2 * f_xn / ((f_x1 - f_x2) * (f_x1 - f_xn))
               + xn * f_x2 * f_x1 / ((f_xn - f_x2) * (f_xn - f_x1));
    }

    /**
     * Secant method. Returns a function value for X-axis secant intersection between points P1(x1,f_x1) and P2(x2,f_x2).
     *
     * @param x2   second point x value
     * @param x1   first point x value
     * @param f_x2 function value in x2
     * @param f_x1 function value in x1
     * @return intersection point for f_x=0
     */
    public static double secantMethod(double x2, double x1, double f_x2, double f_x1) {
        return x1 - f_x1 * (x1 - x2) / (f_x1 - f_x2);
    }

    // DIAGNOSTIC OUTPUT
    private void printEvaluationDiagnostics(String titleMsg) {
        if (showDiagnostics) {
            String formattedMsg = String.format("%s\t EVAL VALUES: a = %.3f, b = %.3f, f(a)= %.3f, f(b)=%.3f", titleMsg, a, b, f_a, f_b);
            LOGGER.log(Level.INFO, formattedMsg);
        }
    }

    private void printSolverResultsDiagnostics(String titleMsg, String endMsg) {
        if (showDiagnostics) {
            String formattedMsg = String.format("%s\t s= %.5f, a= %.5f, f(a)= %.5f, b= %.5f, f(b)= %.5f, c= %.5f, f(c)= %.5f \t %s", titleMsg, s, a, f_a, b, f_b, c, f_c, endMsg);
            LOGGER.log(Level.INFO, formattedMsg);
        }
    }

    // QUICK INSTANCE

    /**
     * Method for obtaining quick and single result for a provided function and expected result range.
     *
     * @param func   provided eqn = 0 as an lambda expression: value -> f(value)
     * @param rangeA first point of the expected result range
     * @param rangeB second point of the expected result range
     * @return calculated root
     */
    public static double ofFunction(DoubleUnaryOperator func, double rangeA, double rangeB) {
        BrentSolver solver = new BrentSolver();
        return solver.calcForFunction(func, rangeA, rangeB);
    }

    // GETTERS & SETTERS

    /**
     * Returns Brent-Dekker accuracy.
     *
     * @return accuracy level
     */
    public double getAccuracy() {
        return accuracy;
    }

    /**
     * Sets accuracy if other than default is required.
     *
     * @param accuracy accuracy
     */
    public final void setAccuracy(double accuracy) {
        this.accuracy = Math.abs(accuracy);
    }

    /**
     * Returns current iteration limit value.
     *
     * @return max iteration limit (int)
     */
    public double getIterationsLimit() {
        return iterationsLimit;
    }

    /**
     * Sets maximum iteration limit
     *
     * @param iterationsLimit maximum iteration limit
     */
    public void setIterationsLimit(double iterationsLimit) {
        this.iterationsLimit = iterationsLimit;
    }

    /**
     * Returns current maximum AB point evaluation procedure cycles
     *
     * @return current maximum evaluation cycles
     */
    public int getEvalCycles() {
        return evalCycles;
    }

    /**
     * Sets maximum AB point evaluation procedure cycles
     *
     * @param evalCycles maximum evaluation cycles
     */
    public void setEvalCycles(int evalCycles) {
        this.evalCycles = evalCycles;
    }

    /**
     * Returns Second point division coefficient (AB evaluation procedure)
     *
     * @return second point div coefficient
     */
    public int getEvalX2Divider() {
        return evalX2Divider;
    }

    /**
     * Sets second point division coefficient (AB evaluation procedure)
     *
     * @param evalX2Divider second point div coefficient
     */
    public void setEvalX2Divider(int evalX2Divider) {
        this.evalX2Divider = evalX2Divider;
    }

    /**
     * Returns second point division coefficient (AB evaluation procedure)
     *
     * @return third point div coefficient
     */
    public int getEvalXDivider() {
        return evalXDivider;
    }

    /**
     * Sets third point division coefficient (AB evaluation procedure calibration)
     *
     * @param evalXDivider third point div coefficient
     */
    public void setEvalXDivider(int evalXDivider) {
        this.evalXDivider = evalXDivider;
    }

    /**
     * Returns number of iterations.
     *
     * @return iteration count
     */
    public int getCounter() {
        return counter;
    }

}