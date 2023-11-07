# Interior Point Algorithm

The interior point algorithm is a powerful numerical optimization method used to solve linear programming problems. It aims to find the optimal solution of a linear programming problem by iteratively approaching the optimal point while staying within the feasible region.

## Algorithm Steps

1. **Initialize Decision Variables**: Start with an initial feasible solution (trial solution) within the feasible region. Typically, this is done by choosing non-negative values for decision variables.

2. **Calculate Diagonal Matrix (D)**: Form a diagonal matrix 'D' using the current values of decision variables. Each element of 'D' is the value of the corresponding decision variable.

3. **Form Matrices (AA and cc)**:
   - Calculate 'AA' by multiplying the coefficient matrix 'A' of constraint functions with the diagonal matrix 'D'.
   - Calculate 'cc' by multiplying 'D' with the coefficients of the objective function 'C'.

4. **Construct Identity Matrix (I)**: Create an identity matrix 'I' with appropriate dimensions (the number of decision variables).

5. **Calculate Matrix 'F' and Inverse 'FI'**:
   - Compute matrix 'F' by multiplying 'AA' with the transpose of 'AA'.
   - Calculate the inverse of 'F' and store it in 'FI' (FI = F^(-1)).

6. **Calculate Matrix 'H'**:
   - Compute matrix 'H' by multiplying the transpose of 'AA' with 'FI'.

7. **Construct Matrix 'P'**:
   - Form matrix 'P' by subtracting 'H' multiplied by 'AA' from the identity matrix 'I' (P = I - H * AA).

8. **Calculate Vector 'cp'**:
   - Compute vector 'cp' by multiplying matrix 'P' with vector 'cc' (cp = P * cc).

9. **Calculate Minimum Value (nu)**:
   - Determine the minimum value in vector 'cp' and assign it to 'nu'.

10. **Update Decision Variables (y)**:
    - Calculate vector 'y' using the formula: y = (1 + (alpha / nu) * cp).

11. **Update Trial Solution (x)**:
    - Update the current values of decision variables ('x') using the formula: x = D * y.

12. **Convergence Check**:
    - Check if the difference between the new decision variables and the previous ones is less than a predefined approximation accuracy ('epsilon'). If the difference is smaller than 'epsilon', stop the algorithm.

13. **Iteration**:
    - Repeat steps 3 to 11 iteratively until convergence is achieved or a predefined maximum number of iterations is reached.

## Output
The algorithm provides the optimal solution (decision variables 'x*') that satisfies all constraints and the maximum (or minimum) value of the objective function, depending on the nature of the optimization problem.
