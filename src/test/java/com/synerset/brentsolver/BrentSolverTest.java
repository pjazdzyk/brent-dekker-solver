package com.synerset.brentsolver;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.function.DoubleFunction;
import java.util.stream.Stream;

class BrentSolverTest {

    private double solver_accuracy;
    private BrentSolver solver;

    @BeforeEach
    void setUp() {
        solver = new BrentSolver("Test-SOLVER");
        solver_accuracy = solver.getAccuracy();
        solver.setShowDiagnostics(true);
    }

    @Test
    @DisplayName("should return root for simple linear function")
    void findRoot_givenSingleVariableFunction_returnsRoot() {
        // Arrange
        DoubleFunction<Double> func = p -> (p + 10) / 20;
        var expectedRoot = -10;

        // Act
        var actualRoot = solver.calcForFunction(func);

        // Assert
        Assertions.assertEquals(expectedRoot, actualRoot, solver_accuracy);
    }

    @Test
    @DisplayName("should return one of two roots within specified solution boundary")
    void findRoot_givenQuadraticFunction_returnRoot() {
        // Arrange
        DoubleFunction<Double> quadraticFunction = x -> 2 * x * x + 5 * x - 3;
        var expectedFirstRoot = -3;
        var expectedSecondRoot = 0.5;

        // Act
        var actualFirstRoot = solver.calcForFunction(quadraticFunction);
        solver.setCounterpartPoints(-1, 2);
        var actualSecondRoot = solver.calcForFunction(quadraticFunction);

        // Assert
        Assertions.assertEquals(expectedFirstRoot, actualFirstRoot, solver_accuracy);
        Assertions.assertEquals(expectedSecondRoot, actualSecondRoot, solver_accuracy);
    }

    @ParameterizedTest
    @MethodSource("polyTestInlineData")
    @DisplayName("should return root for nested log function for series of counterpart points which brakes brent-decker counterpart points condition")
    void findRoot_givenPolynomialFunction_returnRoot(double pointA, double pointB) {
        // Arrange
        DoubleFunction<Double> func = p -> 93.3519196629417 - (-237300 * Math.log(0.001638 * p) / (1000 * Math.log(0.001638 * p) - 17269));
        var expectedRoot = 80000;

        //Act
        var actualRoot = solver.calcForFunction(func, pointA, pointB);

        //Assert
        Assertions.assertEquals(actualRoot, expectedRoot, solver_accuracy);
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
        DoubleFunction<Double> func = x -> Math.acos(x / 2);

        // Assert
        Assertions.assertThrows(BrentSolverException.class, () -> solver.calcForFunction(func));
    }

}
