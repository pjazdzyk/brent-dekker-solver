package com.synerset.brentsolver;

class BrentSolverDefaults {

    private BrentSolverDefaults() {
        throw new IllegalStateException("Utility class");
    }

    public static final String DEF_NAME = "DefaultBrentSolver"; // default name
    public static final double DEF_A0 = -50;                    // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    public static final double DEF_B0 = 50;                     // initial guess or arbitrarily assumed value to get opposite (negative and positive) result from tested equation
    public static final double DEF_ITERATIONS = 100;            // default limit of iterations
    public static final double DEF_ACCURACY = 1E-11;            // expected accuracy level
    public static final int DEF_EVAL_CYCLES = 5;                // number of evaluation cycles for counterpart evaluation procedure
    public static final int DEF_EVAL_X2_COEF = 2;               // divider coefficient used in counterpart evaluation
    public static final int DEF_EVAL_X2_VALUE_COEF = 2;         // divider coefficient used in counterpart evaluation

}