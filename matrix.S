section .text
  global create_matrix

; Externs
;  memory.S
extern malloc

%macro prologue 0
  push rbp
  mov  rbp, rsp
%endmacro

%macro epilogue 0
  mov  rsp, rbp
  pop  rbp
  ret
%endmacro


; matrix* create_matrix(uint64_t width: rax, uint64_t height:rbx):rax
; Initializes a matrix of size width x height on heap
; Table structure:
;   uint64 width
;   uint64 height
;   float  data[width][height] -> {row_0}{row_1}...
; Uses malloc
; SIZE EQUAL TO ZERO NOT SUPPORTED
create_matrix:
  prologue

  ; size = 8 + 8 + width * height
  push rax
  push rbx
  mul  rbx
  add  rax, 0x10
  call malloc
  pop rbx
  pop rcx
  mov QWORD [rax], rcx
  mov QWORD [rax+8], rbx
  add rax, 0x10
  imul rbx, rcx
  xor rcx, rcx
  create_matrix_loop:
    mov QWORD [rax+rcx*8], 0
    inc rcx
    cmp rcx, rbx
    jl  create_matrix_loop
  sub rax, 0x10

  epilogue