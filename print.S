section .data
	NaN db "NaN", 0
	Inf db "INF", 0

section .text
	global print
	global print_float
	global print_int
	global putchar
	global print_table

extern strlen
extern malloc
extern free
extern memcpy
extern memcpy_reverse

%macro prologue 0
	push rbp
	mov  rbp, rsp
%endmacro

%macro epilogue 0
	mov  rsp, rbp
	pop  rbp
%endmacro

print:
	call strlen
	mov rdi, 1         ; stdout
	mov rsi, [rsp+8]   ; buffer
	mov rdx, rax       ; size
	mov rax, 1         ; write
	syscall
	ret


;; Prints float
print_float:
	push rbp
	mov rbp, rsp

	push 256
	call malloc
	mov  rbx, rax

	mov BYTE [rbp-0x01], 0

	mov rdx, [rbp+0x10]
	mov r8, 0x8000000000000000
	mov rcx, rdx
	and rcx, r8
	cmp rcx, r8
	jne  print_float_positive
	mov r8, 0x7fffffffffffffff
	and rdx, r8
	mov BYTE [rbp-0x01], 1
	mov BYTE [rbx], '-'
	inc rbx
	print_float_positive:

	mov r8, 0x7fffffffffffffff
	cmp rdx, r8
	je  print_float_nan

	mov r8, 0x7ff0000000000000
	and rcx, r8
	cmp rcx, r8
	je print_float_inf

	push rdx
	push 256
	call malloc
	mov r9, rax
	pop rdx
	pop r11
	push r11
	fld1
	fld QWORD [rsp]
	fprem

	; Integer part
	fld    QWORD [rsp]
	push   rax
	fisttp QWORD [rbp-0x11]
	mov    rax, [rbp-0x11]
	mov    rcx, 10
	mov    QWORD [rbp-0x09], 0

	print_float_int_loop_begin:
	xor    rdx, rdx
	div    rcx
	lea    r8, [rdx+0x30]
	mov    BYTE [r9], r8b
	inc    r9
	inc    QWORD [rbp-0x09]
	cmp    rax, 0
	jne print_float_int_loop_begin
	
	sub  r9, QWORD [rbp-0x09]
	push QWORD [rbp-0x09]
	push r9
	push rbx
	call memcpy_reverse
	pop rbx
	
	; Decimal point
	mov r9, QWORD [rbp-0x09]
	mov BYTE [rbx+r9], '.'
	inc r9

	add r9, rbx
	; dec r9

	mov QWORD [rbp-0x09], 0

	; Fractional part
	print_float_frac_loop_begin:
	mov  rax, 10
	push  rax
	fild QWORD [rsp]
	fmulp
	fst  QWORD [rsp]
	fld  QWORD [rsp]
	fisttp QWORD [rsp]
	pop rdx
	lea r8, [rdx+0x30]
	mov BYTE [r9], r8b
	cmp QWORD [rbp-0x09], 6
	jge print_float_frac_loop_end
	ftst
	inc r9
	inc QWORD [rbp-0x09]
	fld1
	fxch ST0, ST1
	fprem
	fxch
	fstp QWORD [rsp]
	jnz print_float_frac_loop_begin
	print_float_frac_loop_end:

	; while ([r9] == 0x30): r9--
	print_float_remove_zeros_begin:
	cmp BYTE [r9], 0x30
	jne print_float_remove_zeros_end
	dec r9
	jmp print_float_remove_zeros_begin
	print_float_remove_zeros_end:
	inc r9
	mov BYTE [r9], 0

	jmp print_float_print
	
	print_float_inf:
	push 4
	push Inf
	push rbx
	call memcpy
	pop  rbx
	jmp print_float_print

	print_float_nan:
	push 4
	push NaN
	push rbx
	call memcpy
	pop  rbx
	jmp print_float_print

	print_float_print:
	cmp BYTE [rbp-0x01], 1
	jne print_float_print_positive
	dec rbx
	print_float_print_positive:
	push rbx
	call print
	call free

	jmp print_float_end

	print_float_end:

	push rax
	fstp QWORD [rsp]
	fstp QWORD [rsp]

	mov rsp, rbp
	pop rbp
	ret

; void print_int(int64: rax)
print_int:
	prologue
	mov  rbx, rax
	push 256
	call malloc

	test rbx, rbx
	jns   print_int_positive
	neg  rbx
	push rax
	mov  rax, '-'
	call putchar
	pop  rax
	print_int_positive:

	xor  rsi, rsi
	mov  r8, rax
	mov  rax, rbx
	mov  rcx, 10

	print_int_loop:
		xor rdx, rdx
		div rcx
		add rdx, 0x30

		mov BYTE [r8+rsi], dl
		inc rsi

		cmp rax, 0
		je print_int_end
		jmp print_int_loop

	print_int_end:

	push r8
	push rsi
	push 256
	call malloc
	pop rsi
	pop rsi
	pop r8

	push rsi
	push r8
	push rax
	call memcpy_reverse
	call print
	call free
	pop rax
	call free

	epilogue

putchar:
	push rax
	mov rdi, 1         ; stdout
	mov rsi, rsp-1     ; buffer
	mov rdx, 1         ; size
	mov rax, 1         ; write
	syscall
	pop rax
	ret


; void print_table(table*:rax)
print_table:
  prologue

  mov  rax, QWORD [rax]
  call print_int

  epilogue