package com.synerset.brentsolver;

import java.util.function.DoubleUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;

import static com.synerset.brentsolver.BrentSolverValidators.requireNonInfiniteAndNonNANResults;
import static com.synerset.brentsolver.BrentSolverValidators.requireNonNull;
import static com.synerset.brentsolver.InterpolationAlgorithms.inverseQuadraticInterpolation;
import static com.synerset.brentsolver.InterpolationAlgorithms.linearInterpolationFromValue;
import static com.synerset.brentsolver.InterpolationAlgorithms.secantMethod;

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

    private final String name;

    private double a0;                         // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    private double b0;                         // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    private double a, f_a;                     // first point and its function value f(a), must have opposite sign than b
    private double b, f_b;                     // second point and its function value f(b), absolute value of b should be lower than a
    private double c, f_c;                     // previous iteration point and f(c), initially = a
    private double s, f_s;                     // new approximate value based on secant method or bisection method and f_s
    private int counter;

    private DoubleUnaryOperator userFunction;  // function to compute, to be provided by user

    private boolean runFlag;
    private double iterationsLimit;
    private double accuracy;

    private int evalCycles;
    private int evalDividerX2;
    private int evalDividerX2Value;

    private boolean showDiagnostics;
    private boolean showSummary;

    private static final Logger LOGGER = Logger.getLogger(BrentSolver.class.getName());

    /**
     * Initializes solver instance with function output set as 0, with default name.
     */
    public BrentSolver() {
        this(BrentSolverDefaults.DEF_NAME, x -> 0, BrentSolverDefaults.DEF_A0, BrentSolverDefaults.DEF_B0);
    }

    /**
     * Initializes solver instance with function output set as 0, with specified name
     */
    public BrentSolver(String name) {
        this(name, x -> 0, BrentSolverDefaults.DEF_A0, BrentSolverDefaults.DEF_B0);
    }

    /**
     * Initializes solver instance with a name and user provided function and specified counterpart points. Counterpart points
     * should enclose a range within only one root exists.
     *
     * @param name              solver name
     * @param functionToCompute user function representing an equation to be computed
     * @param a0                lower bound counterpart point, also referred as "a"
     * @param b0                upper bound counterpart point, also referred as "b"
     */
    public BrentSolver(String name, DoubleUnaryOperator functionToCompute, double a0, double b0) {
        requireNonNull("functionToCompute", functionToCompute);
        this.name = name == null ? BrentSolverDefaults.DEF_NAME : name;
        this.userFunction = functionToCompute;
        this.a0 = a0;
        this.b0 = b0;
        this.runFlag = true;
        this.iterationsLimit = BrentSolverDefaults.DEF_ITERATIONS;
        this.accuracy = BrentSolverDefaults.DEF_ACCURACY;
        this.evalCycles = BrentSolverDefaults.DEF_EVAL_CYCLES;
        this.evalDividerX2Value = BrentSolverDefaults.DEF_EVAL_X2_COEF;
        this.evalDividerX2 = BrentSolverDefaults.DEF_EVAL_X2_VALUE_COEF;
    }

    // Calculation methods

    /**
     * Root finding algorithm. Returns root value based on Brent-Dekker numerical scheme. Calculation procedure will use
     * expressions provided in userFunction in each iteration step.
     *
     * @return actual root value meeting accuracy criteria
     */
    public double findRoot() {
        long startTime = System.nanoTime();

        // Initializing Solver input variables, checking counterpart points conditions
        initializeAndCheckConditions();

        /*--------BEGINNING OF ITERATIVE LOOP--------*/
        log("Starting calculations....");
        logCurrentSolutionStatus("INITIAL VALUES: ", "");
        while (runFlag) {
            counter++;

            // New additional condition proposed by Zhengqiu Zhang
            c = (a + b) / 2;
            f_c = userFunction.applyAsDouble(c);

            // Determining interpolation type for this iteration: inverse quadratic or secant method
            if ((f_a != f_c) && (f_b != f_c)) {
                s = inverseQuadraticInterpolation(a, b, c, f_a, f_b, f_c);
            } else {
                s = secantMethod(a, b, f_a, f_b);
            }

            // Checking current results and swap if needed
            if (c > s) {
                double tempC = c;
                double tempS = s;
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
            double currentDifference = Math.abs(b - a);


            logCurrentSolutionStatus(String.format("ITERATION: %-5s", counter), "Diff = " + currentDifference);

            // Checking conditions if they meet solution criteria
            if (currentDifference < accuracy) {
                log("Solver stopped, calculated solution meets accuracy requirement of {0}.", String.valueOf(accuracy));
                runFlag = false;
            } else if (counter > iterationsLimit) {
                runFlag = false;
                log("Solver stopped, reached iterations limit of {0}. Solution is not converged.", iterationsLimit);
            } else if (f_b == 0) {
                runFlag = false;
            }

            // Exception will be thrown if NaN or Infinite values are detected
            requireNonInfiniteAndNonNANResults(name, f_a, f_b, f_c, f_s);

            /*-----------END OF ITERATIVE LOOP-----------*/
        }

        long endTime = System.nanoTime();
        long durationNano = endTime - startTime;
        double durationMillis = durationNano / 1_000_000.0;
        double durationSeconds = durationNano / 1_000_000_000.0;

        // b is always closer to the root
        logSummary("CALCULATIONS COMPLETE: Root found: [{0}] in {1} iterations. Completed in: {2} millis or {3} seconds. Target accuracy: {4}.",
                String.valueOf(b), counter, durationMillis, durationSeconds, String.valueOf(accuracy));

        return b;

    }

    /**
     * Returns a function root (calculation result) based on provided user-function.
     *
     * @param func tested function (use lambda expression or method reference)
     * @return function root (Double)
     */
    public double findRoot(DoubleUnaryOperator func) {
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
    public double findRoot(DoubleUnaryOperator func, double a0, double b0) {
        setCounterpartPoints(a0, b0);
        return findRoot(func);
    }

    // Solver control

    /**
     * Resets flags and iteration counter
     */
    public void resetSolverRunFlags() {
        this.runFlag = true;
        this.counter = 0;
    }

    /**
     * Resets counter part points to default values (+50,-50);
     */
    public void resetCounterPartPoints() {
        a0 = BrentSolverDefaults.DEF_A0;
        b0 = BrentSolverDefaults.DEF_B0;
    }

    /**
     * Resets evaluation procedure coefficients to default values.
     */
    public void resetEvaluationCoefficients() {
        evalCycles = BrentSolverDefaults.DEF_EVAL_CYCLES;
        evalDividerX2 = BrentSolverDefaults.DEF_EVAL_X2_COEF;
        evalDividerX2Value = BrentSolverDefaults.DEF_EVAL_X2_VALUE_COEF;
    }

    // Getters, setters

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
     * Sets custom initial counterpart points for Brent-Dekker method condition.
     *
     * @param pointA - first initial counterpart point
     * @param pointB - second initial counterpart point
     */
    public void setCounterpartPoints(double pointA, double pointB) {
        this.a0 = pointA;
        this.b0 = pointB;
        resetSolverRunFlags();
    }

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
    public void setAccuracy(double accuracy) {
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
    public int getEvalDividerX2() {
        return evalDividerX2;
    }

    /**
     * Sets second point division coefficient (AB evaluation procedure)
     *
     * @param evalDividerX2 second point div coefficient
     */
    public void setEvalDividerX2(int evalDividerX2) {
        this.evalDividerX2 = evalDividerX2;
    }

    /**
     * Returns second point division coefficient (AB evaluation procedure)
     *
     * @return third point div coefficient
     */
    public int getEvalDividerX2Value() {
        return evalDividerX2Value;
    }

    /**
     * Sets third point division coefficient (AB evaluation procedure calibration)
     *
     * @param evalDividerX2Value third point div coefficient
     */
    public void setEvalDividerX2Value(int evalDividerX2Value) {
        this.evalDividerX2Value = evalDividerX2Value;
    }

    /**
     * Returns number of iterations.
     *
     * @return iteration count
     */
    public int getCounter() {
        return counter;
    }

    // Initial conditions initializer

    private void initializeAndCheckConditions() {
        log("Starting condition evaluation procedure.");

        // To check and set value b as being closer to the root
        f_a = userFunction.applyAsDouble(a0);
        f_b = userFunction.applyAsDouble(b0);
        checkSetAndSwapABPoints(a0, b0);

        // Checking if Brent AB condition is not met to launch automatic AB points evaluation procedure
        if (initialABConditionIsNotMet()) {
            runCounterPartConditionEvaluator();
        }

        // In case provided by user point "a" or "b" is actually a root
        if (accuracyRequirementIsMet()) {
            log("CALCULATION COMPLETE. Initial value is a root: {}.", b);
            runFlag = false;
            return;
        }

        // If at this stage proper A-B condition is not achievable - an exception is thrown.
        if (initialABConditionIsNotMet()) {
            String errorMsg = String.format("%s: EVALUATION PROCEDURE FAILED: f(a) i f(b) must have an opposite signs. " +
                                            "Current values: a = %.3f, b = %.3f, f(a)= %.3f, f(b)=%.3f", name, a, b, f_a, f_b);
            log(Level.SEVERE, errorMsg);
            throw new BrentSolverException(errorMsg);
        }

        log("Condition evaluation successful.");
    }

    // Debugging solution with diagnostic output

    /**
     * Sets diagnostic output mode. (true = show output, false = no output).
     *
     * @param showDebugLogs true = on, false = off
     */
    public void showDebugLogs(boolean showDebugLogs) {
        this.showDiagnostics = showDebugLogs;
    }

    /**
     * Sets solver calculations summary. (true = show output, false = no output).
     *
     * @param showSummaryLogs true = on, false = off
     */
    public void showSummaryLogs(boolean showSummaryLogs) {
        this.showSummary = showSummaryLogs;
    }

    // Solution state requirements

    private boolean initialABConditionIsNotMet() {
        return (f_a * f_b) >= 0;
    }

    private boolean accuracyRequirementIsMet() {
        return Math.abs(f_b) <= accuracy;
    }

    private void checkSetAndSwapABPoints(double pointA, double pointB) {
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

    // Counterpart points evaluation procedure

    /**
     * This method attempts to evaluate valid Brent Solver counterpart point condition. <p>
     * PROCEDURE EXPLANATION <br>
     * Using linear extrapolation to determine opposite sign of the "b" value. <br>
     * 1. Linear extrapolations needs a pair of points: P1(x1,f_x1), P2(x2,f_x2) and f_y to determine x. <br>
     * 2. First point is already provided (b, f_b) and it is already evaluated to be closer to the root <br>
     * 3. Second point P2 is created from b/evalX2Divider and its function value. <br>
     * 4. f_x is an x value we expect to have an opposite sign to -b, and to make it closer to the root - it is divided by evalXDivider <br>
     * 5. In some cases initial values of evalX2Divider and evalXDivider may be adjusted by user. For an example if your initial guess is very close to the root
     * you are looking for small evalXDivider values, like 2 or even 1.
     */
    private void runCounterPartConditionEvaluator() {
        log("EVALUATION PROCEDURE:");
        logCurrentSolutionStatus("INITIAL");

        double x, f_x, x1, f_x1, x2, f_x2, f_xExact;

        if (evalDividerX2Value > evalCycles) {
            evalCycles = evalDividerX2Value;
        }

        for (int i = 0; i <= evalCycles; i++) {

            if (evalDividerX2Value - i == 0.0) {
                continue;
            }

            // Point 1, using b as a starting point
            x1 = b;
            f_x1 = f_b;

            // Point 2, creating second point from the b
            x2 = b / evalDividerX2;
            f_x2 = userFunction.applyAsDouble(x2);

            // Point 3, searching for a negative value of -f_b. Further division is to get result as close to the root as possible.
            f_x = -f_b / (evalDividerX2Value - i);
            x = linearInterpolationFromValue(x1, f_x1, x2, f_x2, f_x);

            // When x is determined - to check if it really gives negative value (it may not occur for strong non-linearity)
            f_xExact = userFunction.applyAsDouble(x);
            f_a = userFunction.applyAsDouble(b);
            f_b = userFunction.applyAsDouble(x);
            checkSetAndSwapABPoints(b, x);

            logCurrentSolutionStatus("STEP " + i + ":");

            if(accuracyRequirementIsMet()){
                return;
            }

            if (f_xExact * f_x1 < 0)
                break;
        }
    }

    // Loggers and diagnostic outputs

    private void log(String msg, Object... msgParams) {
        log(Level.INFO, msg, msgParams);
    }

    private void log(Level level, String msg, Object... msgParams) {
        if (showDiagnostics) {
            LOGGER.log(level, String.format("[%s] - %s", name, msg), msgParams);
        }
    }

    private void logSummary(String msg, Object... msgParams) {
        if (showSummary || showDiagnostics) {
            LOGGER.log(Level.INFO, String.format("[%s] - %s", name, msg), msgParams);
        }
    }

    private void logCurrentSolutionStatus(String titleMsg) {
        String formattedMsg = String.format("%s EVAL VALUES: a = %.3f, b = %.3f, f(a)= %.3f, f(b)=%.3f",
                titleMsg, a, b, f_a, f_b);

        log(formattedMsg);
    }

    private void logCurrentSolutionStatus(String titleMsg, String endMsg) {
        String formattedMsg = String.format("%s s= %.5f, a= %.5f, f(a)= %.5f, b= %.5f, f(b)= %.5f, c= %.5f, f(c)= %.5f %s",
                titleMsg, s, a, f_a, b, f_b, c, f_c, endMsg);

        log(formattedMsg);
    }

    // Static factory methods

    /**
     * Method for obtaining quick and single result for a provided function and expected result range.
     *
     * @param func provided eqn = 0 as an lambda expression: value -> f(value)
     * @param a0   first point of the expected result range
     * @param b0   second point of the expected result range
     * @return calculated root
     */
    public static double findRootOf(DoubleUnaryOperator func, double a0, double b0) {
        BrentSolver solver = BrentSolver.of(func, a0, b0);
        return solver.findRoot();
    }

    public static BrentSolver of() {
        return new BrentSolver();
    }

    public static BrentSolver of(String name) {
        return new BrentSolver(name, x -> 0, BrentSolverDefaults.DEF_A0, BrentSolverDefaults.DEF_B0);
    }

    public static BrentSolver of(DoubleUnaryOperator functionToCompute, double a0, double b0) {
        return new BrentSolver(BrentSolverDefaults.DEF_NAME, functionToCompute, a0, b0);
    }

}