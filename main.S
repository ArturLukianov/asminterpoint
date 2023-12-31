section .data
  enterVarNum   db "Enter number of variables: ", 0
  enterConNum   db "Enter number of constraints: ", 0
  enterConCoeff db "Enter the constraints' coefficients:", 0Ah, 0
  enterVarCoeff db "Enter the variables' coefficients:", 0Ah, 0
  enterConB     db "Enter the constraints' B:", 0Ah, 0

section .text
	global _start

; Externs
;  print.S
extern print
extern print_matrix
;  input.S
extern read_int
;  matrix.S
extern C_matrix
extern create_matrix

; Macros
%macro print 1
  push %1
  call print
  add  rsp, 8
%endmacro

_start:
  print enterVarNum
  call  read_int
  push  rax
  print enterConNum
  call  read_int
  pop   rbx

  add  rax, rbx
  mov  rbx, 1
  call create_matrix
  mov  QWORD [C_matrix], rax
  call print_matrix

exit:
	mov rax, 60
	mov rdi, 0
	syscall ; exit(0)
