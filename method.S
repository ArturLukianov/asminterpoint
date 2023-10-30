section .data
  global C_table
  C_table resq 1

section .text
  global print_table

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