package com.synerset.brentsolver;

public class InterpolationAlgorithms {

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

}
