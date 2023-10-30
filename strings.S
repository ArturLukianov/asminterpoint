section .text
	global strlen

strlen:
	add rsp, 8
	mov rax, [rsp+8]
strlen_beginloop:
	cmp byte [rax], 0
	je  strlen_endloop
	inc rax
	jmp strlen_beginloop
strlen_endloop:
	sub rax, [rsp+8]
	sub rsp, 8
	ret

