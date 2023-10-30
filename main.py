'''
z = x1 + 2x2
x1 + x2 <= 8

A = [1, 1]
C = [1, 2]
B = [8]

'''

def form_d(D, x):
	i = 0
	while i < len(D):
		D[i][i] = x[i][0]
		i += 1

def mult(A, B, R):
	for i in range(len(A)):
		for j in range(len(B[0])):
			s = 0
			for k in range(len(B)):
				s += A[i][k] * B[k][j]
			R[i][j] = s

def sub(A, B, R):
	for i in range(len(A)):
		for j in range(len(A[0])):
			R[i][j] = A[i][j] - B[i][j]

def add(A, B, R):
	for i in range(len(A)):
		for j in range(len(A[0])):
			R[i][j] = A[i][j] + B[i][j]

def scale(A, k, R):
	for i in range(len(A)):
		for j in range(len(A[0])):
			R[i][j] = A[i][j] * k

def transpose(A, R):
	for i in range(len(A)):
		for j in range(len(A[0])):
			R[j][i] = A[i][j]

def init(A):
	for i in range(len(A)):
		A[i][i] = 1

def invert(A, R):
	for i in range(len(R)):
		R[i][i] = 1

	for i in range(len(A)):
		div_line(R, i, A[i][i])
		div_line(A, i, A[i][i])
		for j in range(i + 1, len(A)):
			sub_line(R, i, i, j, A[j][i])
			sub_line(A, i, i, j, A[j][i])

	for i in range(len(A) - 1, 0, -1):
		for j in range(i - 1, -1, -1):
			sub_line(R, i, i, j, A[j][i])
			sub_line(A, i, i, j, A[j][i])

def sub_line(A, i, j, k, c):
	for h in range(len(A)):
		A[k][h] -= A[i][h] * c

def div_line(A, i, c):
	for j in range(len(A)):
		A[i][j] /= c

def get_max_negative(A):
	m = 0
	for i in range(len(A)):
		if A[i][0] < m:
			m = A[i][0]
	return m

def print_matrix(A):
	for i in range(len(A)):
		for j in range(len(A[0])):
			print(A[i][j], end='')
			print(' ', end='')
		print('\n', end='')
	print_line()

def print_line():
	print('------------')

def nullify(A):
	for i in range(len(A)):
		for j in range(len(A[0])):
			if abs(A[i][j]) < 0.0000000000001:
				A[i][j] = 0

A = [[1, 1, 1, 0],
	 [1, 2, 0, 1]]
C = [[1], [2], [0], [0]]
B = [8, 5]
X = [[2, 2], [2, 2], [B[0] - 2 - 2, 0], [0, B[1] - 2 - 2]]

HEAP = {}

for i in range(100000):
	D = [[0 for _ in range(3)] for _ in range(3)]
	form_d(D, X)
	# print_matrix(D)
	A_ = [[0 for _ in range(3)] for _ in range(1)]
	mult(A, D, A_)
	# print_matrix(A_)
	C_ = [[0 for _ in range(1)] for _ in range(3)]
	mult(D, C, C_)
	# print_matrix(C_)
	A_T = [[0 for _ in range(1)] for _ in range(3)]
	transpose(A_, A_T)
	A_A_T = [[0 for _ in range(1)] for _ in range(1)]
	mult(A_, A_T, A_A_T)
	# print_matrix(A_A_T)
	A_A_TI = [[0 for _ in range(1)] for _ in range(1)]
	invert(A_A_T, A_A_TI)
	A_TA_A_TI = [[0 for _ in range(1)] for _ in range(3)]
	mult(A_T, A_A_TI, A_TA_A_TI)
	# print_matrix(A_TA_A_TI)
	A_TA_A_TIA_ = [[0 for _ in range(3)] for _ in range(3)]
	mult(A_TA_A_TI, A_, A_TA_A_TIA_)
	# print_matrix(A_TA_A_TIA_)
	I = [[0 for _ in range(3)] for _ in range(3)]
	init(I)
	P = [[0 for _ in range(3)] for _ in range(3)]
	sub(I, A_TA_A_TIA_, P)
	# print_matrix(P)
	C_P = [[0 for _ in range(1)] for _ in range(3)]
	mult(P, C_, C_P)
	# print_matrix(C_P)
	X_ = [[1] for _ in range(3)]
	# print_matrix(X_)
	alpha = 0.5
	v = get_max_negative(C_P)
	v = abs(v)
	if v == 0:
		break
	X__ = [[0] for _ in range(3)]
	add(X__, C_P, X__)
	scale(X__, alpha / v, X__)
	add(X_, X__, X_)
	X = [[0] for _ in range(3)]
	mult(D, X_, X)
	nullify(X)
	print_matrix(X)

print_matrix(X)