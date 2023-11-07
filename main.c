#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
  int size;
  double *data;
} Vector;

typedef struct
{
  int rows;
  int cols;
  double **data;
} Matrix;

Vector createVector(int size)
{
  Vector v;
  v.size = size;
  v.data = (double *)malloc(size * sizeof(double));
  return v;
}

Matrix createMatrix(int rows, int cols)
{
  Matrix m;
  m.rows = rows;
  m.cols = cols;
  m.data = (double **)malloc(rows * sizeof(double *));
  for (int i = 0; i < rows; i++)
  {
    m.data[i] = (double *)malloc(cols * sizeof(double));
  }
  return m;
}

Matrix identityMatrix(int rows)
{
  Matrix m;
  m.rows = rows;
  m.cols = rows;
  m.data = (double **)malloc(rows * sizeof(double *));
  for (int i = 0; i < rows; i++)
  {
    m.data[i] = (double *)malloc(rows * sizeof(double));
  }
  for (int i = 0; i < rows; i++)
  {
    m.data[i][i] = 1;
  }
  return m;
}

void freeVector(Vector v)
{
  free(v.data);
}

void freeMatrix(Matrix *m)
{
  for (int i = 0; i < m->rows; i++)
  {
    free(m->data[i]);
  }
  free(m->data);
}

Vector scalarMultiplyVector(double scalar, Vector v)
{
  Vector result = createVector(v.size);
  for (int i = 0; i < v.size; i++)
  {
    result.data[i] = scalar * v.data[i];
  }
  return result;
}

Vector subtractVectors(Vector a, Vector b)
{
  Vector result = createVector(a.size);
  for (int i = 0; i < a.size; i++)
  {
    result.data[i] = a.data[i] - b.data[i];
  }
  return result;
}

double vectorNorm(Vector v)
{
  double sum = 0;
  for (int i = 0; i < v.size; i++)
  {
    sum += v.data[i] * v.data[i];
  }
  return sqrt(sum);
}

Matrix subtractMatrices(Matrix *m1, Matrix *m2)
{
  Matrix result = createMatrix(m1->rows, m1->cols);
  for (int i = 0; i < m1->rows; i++)
  {
    result.data[i] = m2->data[i] - m2->data[i];
  }
  return result;
}

// Function to create a diagonal matrix from a vector
Matrix diagonalMatrix(Vector v)
{
  Matrix m = createMatrix(v.size, v.size);
  for (int i = 0; i < v.size; i++)
  {
    for (int j = 0; j < v.size; j++)
    {
      m.data[i][j] = (i == j) ? v.data[i] : 0.0;
    }
  }
  return m;
}

Matrix multiplyMatrices(Matrix *m1, Matrix *m2)
{
  Matrix result = createMatrix(m1->cols, m2->rows);
  if (m1->cols != m2->rows)
  {
    printf("Error: incompatible matrix sizes for multiplication.\n");
    return;
  }

  for (int i = 0; i < m1->rows; i++)
  {
    for (int j = 0; j < m2->cols; j++)
    {
      result.data[i][j] = 0.0;
      for (int k = 0; k < m1->cols; k++)
      {
        result.data[i][j] += m1->data[i][k] * m2->data[k][j];
      }
    }
  }
  return result;
}

Matrix transposeMatrix(Matrix *m)
{
  Matrix transposed = createMatrix(m->cols, m->rows);

  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      transposed.data[j][i] = m->data[i][j];
    }
  }

  return transposed;
}

Matrix invertMatrix(Matrix *m)
{
  if (m->rows != m->cols)
  {
    printf("Error: only square matrices can be inverted.\n");
    exit(-1);
  }

  // Create the augmented matrix with the identity matrix
  Matrix augmented = createMatrix(m->rows, m->cols * 2);
  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      augmented.data[i][j] = m->data[i][j];
      augmented.data[i][j + m->cols] = (i == j) ? 1.0 : 0.0;
    }
  }

  // Perform Gauss-Jordan elimination
  for (int i = 0; i < m->rows; i++)
  {
    // TODO: Add pivoting and check for zero pivot
    double pivot = augmented.data[i][i];
    if (fabs(pivot) < 1e-10)
    {
      printf("Error: matrix is singular and cannot be inverted.\n");
      freeMatrix(&augmented);
      exit(-1);
    }

    for (int j = 0; j < 2 * m->cols; j++)
    {
      augmented.data[i][j] /= pivot;
    }

    for (int k = 0; k < m->rows; k++)
    {
      if (k != i)
      {
        double factor = augmented.data[k][i];
        for (int j = 0; j < 2 * m->cols; j++)
        {
          augmented.data[k][j] -= factor * augmented.data[i][j];
        }
      }
    }
  }
  Matrix result = createMatrix(m->rows, m->cols);
  // Extract the inverse matrix from the augmented matrix
  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      result.data[i][j] = augmented.data[i][j + m->cols];
    }
  }

  freeMatrix(&augmented);

  return result;
}

Matrix scalarDivideMatrix(Matrix *m, double scalar)
{
  Matrix result = createMatrix(m->rows, m->cols);
  if (scalar == 0)
  {
    printf("Division by zero error.\n");
    return;
  }

  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      result.data[i][j] = m->data[i][j] / scalar;
    }
  }
  return result;
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
  Vector x = createVector(n);    // Example initialization
  Matrix A = createMatrix(m, n); // Example initialization
  Matrix C = createMatrix(n, 1); // Example initialization

  Matrix D = diagonalMatrix(x); // A function to create a diagonal matrix from x
  Matrix ATA, ATAI, H, P, cp;
  double nu, alpha = 0.5;
  int iteration = 0;

  while (iteration < maxIterations)
  {
    Matrix AA = multiplyMatrices(&A, &D); // AA = A * D
    Matrix cc = multiplyMatrices(&D, &C); // cc = D * C
    Matrix I = identityMatrix(n);         // I = Identity matrix
    Matrix AAT = transposeMatrix(&AA);

    ATA = multiplyMatrices(&AA, &AAT); // ATA = AA * AA'
    ATAI = invertMatrix(&ATA);         // ATAI = (ATA)^(-1)
    H = multiplyMatrices(&AAT, &ATAI); // H = AA' * ATAI

    Matrix HAA = multiplyMatrices(&H, &AA);
    P = subtractMatrices(&I, &HAA); // P = I - H * AA
    cp = multiplyMatrices(&P, &cc); // cp = P * cc

    nu = absoluteMin(&cp);                                                          // Function to find the absolute minimum value in cp
    Matrix y = scalarDivideMatrix(addMatrices(identityMatrix(n), &cp), alpha / nu); // y = 1 + (alpha / nu) * cp

    Matrix yy = multiplyMatrices(&D, &y); // yy = D * y

    // Check for convergence
    if (matrixNorm(subtractMatrices(&yy, &x)) < epsilon)
    {
      break;
    }

    // Update x for the next iteration
    copyMatrix(&x, &yy);

    iteration++;

    // Free temporary matrices
    freeMatrix(&AA);
    freeMatrix(&cc);
    freeMatrix(&ATA);
    freeMatrix(&ATAI);
    freeMatrix(&H);
    freeMatrix(&P);
    freeMatrix(&cp);
    freeMatrix(&yy);
  }

  printMatrix(&x); // Print the final result

  // Free all matrices
  freeVector(x);
  freeMatrix(&A);
  freeMatrix(&C);
  freeMatrix(&D);

  // Free dynamically allocated memory
  free(C);
  for (int i = 0; i < m; i++)
  {
    free(A[i]);
  }
  return 0;
}
