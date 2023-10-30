section .text
  global create_table

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


; table* create_table(uint64_t width: rax, uint64_t height:rbx):rax
; Initializes a table of size width x height on heap
; Table structure:
;   uint64 width
;   uint64 height
;   float  data[width][height] -> {row_0}{row_1}...
; Uses malloc
; SIZE EQUAL TO ZERO NOT SUPPORTED
create_table:
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
  create_table_loop:
    mov QWORD [rax+rcx*8], 0
    inc rcx
    cmp rcx, rbx
    jl  create_table_loop
  sub rax, 0x10

  epilogue