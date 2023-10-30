section .data
  global C_matrix
  C_matrix resq 1

section .text
  global print_matrix

; Externs

; Macros
%macro prologue 0
  push rbp
  mov  rbp, rsp
%endmacro

%macro epilogue 0
  mov  rsp, rbp
  pop  rbp
  ret
%endmacro