package com.synerset.brentsolver;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.function.DoubleUnaryOperator;
import java.util.stream.Stream;

class BrentSolverTest {

    private BrentSolver solver;

    @BeforeEach
    void setUp() {
        solver = BrentSolver.of("Test-SOLVER");
        solver.showDebugLogs(true);
    }

    @Test
    @DisplayName("should return root for simple linear function")
    void findRoot_givenSingleVariableFunction_returnsRoot() {
        // Arrange
        DoubleUnaryOperator func = p -> (p + 10) / 20;
        var expectedRoot = -10;

        // Act
        var actualRoot = solver.findRoot(func);

        // Assert
        Assertions.assertEquals(expectedRoot, actualRoot, 1E-10);
    }

    @Test
    @DisplayName("should return one of two roots within specified solution boundary")
    void findRoot_givenQuadraticFunction_returnRoot() {
        // Arrange
        DoubleUnaryOperator quadraticFunction = x -> 2 * x * x + 5 * x - 3;
        var expectedFirstRoot = -3;
        var expectedSecondRoot = 0.5;

        // Act
        var actualFirstRoot = solver.findRoot(quadraticFunction);
        solver.setCounterpartPoints(-1, 2);
        var actualSecondRoot = solver.findRoot(quadraticFunction);

        // Assert
        Assertions.assertEquals(expectedFirstRoot, actualFirstRoot, 1E-10);
        Assertions.assertEquals(expectedSecondRoot, actualSecondRoot, 1E-10);
    }

    @ParameterizedTest
    @MethodSource("polyTestInlineData")
    @DisplayName("should return root for nested log function for series of counterpart points which brakes brent-decker counterpart points condition")
    void findRoot_givenPolynomialFunction_returnRoot(double pointA, double pointB) {
        // Arrange
        DoubleUnaryOperator func = p -> 93.3519196629417 - (-237300 * Math.log(0.001638 * p) / (1000 * Math.log(0.001638 * p) - 17269));
        var expectedRoot = 80000;

        //Act
        var actualRoot = solver.findRoot(func, pointA, pointB);

        //Assert
        Assertions.assertEquals(actualRoot, expectedRoot, 1E-9);
    }

    static Stream<Arguments> polyTestInlineData() {
        return Stream.of(
                Arguments.of(50000, 120000),
                Arguments.of(80000, 200000),
                Arguments.of(80000, 80000),
                Arguments.of(20000, 80000),
                Arguments.of(10000, 20000)
        );
    }

    @Test
    @DisplayName("should throw an exception if point evaluation procedure fails to determine valid counterpart points")
    void findRoot_givenAcosFunction_throwsSolverResultException() {
        // Arrange
        solver.setCounterpartPoints(10, 5);
        DoubleUnaryOperator func = x -> Math.acos(x / 2);

        // Assert
        Assertions.assertThrows(BrentSolverException.class, () -> solver.findRoot(func));
    }

}
