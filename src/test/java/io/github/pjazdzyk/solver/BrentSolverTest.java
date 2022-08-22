package io.github.pjazdzyk.solver;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.function.DoubleFunction;
import java.util.stream.Stream;

public class BrentSolverTest {

    @Test
    @DisplayName("Returns root for linear single variable function")
    public void calcForFunction_singleVariableFunction_returnsRoot() {

        // Arrange
        BrentSolver solver = new BrentSolver("Test-SOLVER");
        DoubleFunction<Double> func = p -> (p + 10) / 20;
        var expected = -10;

        // Act
        var result = solver.calcForFunction(func);

        // Assert
        Assertions.assertEquals(expected, result, solver.getAccuracy());

    }

    @ParameterizedTest
    @MethodSource("polyTestInlineData")
    @DisplayName("Returns root for polynomial nested single variable function")
    public void calcForFunction_polynomialFunction_returnRoot(double a, double b) {

        // Arrange
        BrentSolver solver = new BrentSolver("TEST-SOLVER");
        DoubleFunction<Double> func = p -> 93.3519196629417 - (-237300 * Math.log(0.001638 * p) / (1000 * Math.log(0.001638 * p) - 17269));
        var expected = 80000;

        //Act
        var actual = solver.calcForFunction(func, a, b);

        //Assert
        Assertions.assertEquals(actual, expected, solver.getAccuracy());

    }

    public static Stream<Arguments> polyTestInlineData() {
        return Stream.of(
                Arguments.of(50000, 120000),
                Arguments.of(80000, 200000),
                Arguments.of(80000, 80000),
                Arguments.of(20000, 80000),
                Arguments.of(10000, 20000)
        );
    }

}
