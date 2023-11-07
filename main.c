#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to calculate the dot product of two vectors
void dotProduct(double *result, double *vector1, double *vector2, int size)
{
  *result = 0;
  for (int i = 0; i < size; i++)
  {
    *result += vector1[i] * vector2[i];
  }
}

int main()
{
  int n, m;

  printf("Enter the number of decision variables (n): ");
  scanf("%d", &n);

  printf("Enter the number of constraints (m): ");
  scanf("%d", &m);

  // Allocate memory for the coefficients of the objective function (C)
  double *C = (double *)malloc(n * sizeof(double));

  printf("Enter the coefficients of the objective function (C):\n");
  for (int i = 0; i < n; i++)
  {
    scanf("%lf", &C[i]);
  }

  // Allocate memory for the coefficients of constraint function (A)
  double **A = (double **)malloc(m * sizeof(double *));
  for (int i = 0; i < m; i++)
  {
    A[i] = (double *)malloc(n * sizeof(double));
  }

  printf("Enter the coefficients of the constraint function (A):\n");
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      scanf("%lf", &A[i][j]);
    }
  }

  // Allocate memory for the right-hand side numbers (b)
  double *b = (double *)malloc(m * sizeof(double));

  printf("Enter the right-hand side numbers (b):\n");
  for (int i = 0; i < m; i++)
  {
    scanf("%lf", &b[i]);
  } 

  double epsilon;

  printf("Enter the approximation accuracy (epsilon): ");
  scanf("%lf", &epsilon);

  double alpha = 0.5; // Set alpha to 0.5 by default

  printf("Choose alpha (0.5 or 0.9): ");
  scanf("%lf", &alpha);

  if (alpha != 0.5 && alpha != 0.9)
  {
    printf("Invalid choice for alpha. Please choose 0.5 or 0.9.\n");
    return 1;
  }

  int maxIterations = 10000; // Maximum number of iterations

  // Allocate memory for decision variables (x)
  double *x = (double *)malloc(n * sizeof(double));
  
  // Free dynamically allocated memory
  free(C);
  for (int i = 0; i < m; i++)
  {
    free(A[i]);
  }
  free(A);
  free(b);
  free(x);

  return 0;
}
