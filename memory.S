section .text
	global memcpy
	global memcpy_reverse

	global malloc
	global free

;; memcpy_reverse(void* dst, void* src, uint64_t count)
;; > dst   (rdx)
;;   src   (rcx)
;;   count (r8)
memcpy_reverse:
	push rbp
	mov  rbp, rsp

	mov  r8, [rbp+0x20]
	mov  rcx, [rbp+0x18]
	mov  rdx, [rbp+0x10]
	add  rcx, r8
	dec  rcx
	xor  rax, rax
	memcpy_reverse_loop_begin:
	cmp  rax, r8
	jge  memcpy_reverse_loop_end
	mov  bl, [rcx]
	mov  [rdx], bl
	dec  rcx
	inc  rdx
	inc  rax
	jmp  memcpy_reverse_loop_begin
	memcpy_reverse_loop_end:

	mov rsp, rbp
	pop rbp

	ret

;; memcpy(void* dst, void* src, uint64_t count)
;; > dst   (rdx)
;;   src   (rcx)
;;   count (r8)
memcpy:
	push rbp
	mov  rbp, rsp

	mov  r8, [rbp+0x20]
	mov  rcx, [rbp+0x18]
	mov  rdx, [rbp+0x10]
	xor  rax, rax
	memcpy_loop_begin:
	cmp  rax, r8
	jge  memcpy_loop_end
	mov  bl, [rcx]
	mov  [rdx], bl
	inc  rcx
	inc  rdx
	inc  rax
	jmp  memcpy_loop_begin
	memcpy_loop_end:

	mov rsp, rbp
	pop rbp

	ret

;; Chunk structure:
;; +------------+
;; | CHUNK SIZE | QWORD
;; +------------+
;; | CHUNK DATA | LENGTH: SIZE - 8
;;      ....
;; 

;; void* malloc(uint64_t count: rax): rax
malloc:
	push rbp
	mov rbp, rsp

	push rax

	; Call mmap
	mov rdi, 0x00
	add rax, 8
	mov rsi, rax ; SIZE
	mov rdx, 0x03 ; READ|WRITE
	mov r10, 0x22 ; MAP_ANONYMOUS|MAP_PRIVATE
	mov r8,  -1   ; fd
	mov r9,  0x00 ; offset
	mov rax, 0x09	; mmap
	syscall
	pop rdx
	mov QWORD [rax], rdx
	add rax, 0x08

	mov rsp, rbp
	pop rbp

	ret

;; free(void* ptr)
free:
	push rbp
	mov rbp, rsp

	mov rdi, QWORD [rbp+0x10] ; ptr
	sub rdi, 8

	mov rsi, QWORD [rdi] ; size
	
	mov rax, 0x0B	; munmap
	syscall

	mov rsp, rbp
	pop rbp

	ret