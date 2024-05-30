package com.synerset.brentsolver;

class BrentSolverValidators {

    private BrentSolverValidators() {
        throw new IllegalStateException("Utility class");
    }

    public static void requireNonInfiniteAndNonNANResults(String name, double... values) {
        for (double num : values) {
            if (Double.isInfinite(num))
                throw new BrentSolverException(name + ": Solution error. Infinite number detected.");
            if (Double.isNaN(num))
                throw new BrentSolverException(name + ": Solution error. NaN value detected.");
        }
    }

    public static void requireNonNull(String variableName, Object object) {
        if (object == null) {
            throw new BrentSolverException("Argument [" + variableName + "] cannot be null.");
        }
    }

}
