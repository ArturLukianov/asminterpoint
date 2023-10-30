section .text
    global read_float
    global read_int

extern strlen
extern malloc
extern free
extern memcpy
extern memcpy_reverse

getchar:
    push rbp
    mov  rbp, rsp

    push QWORD 0
    
    mov  rax, 0
    mov  rdi, 0
    mov  rsi, rsp  ; 
    mov  rdx, 1    ; size
    syscall

    pop  rax

    mov  rsp, rbp
    pop  rbp
    ret

;; Read float from stdin to RAX
read_float:
    push rbp
    mov rbp, rsp

    fldz
    sub rsp, 0x8

    call getchar
    mov  r12, rax

    xor r14b, r14b

    cmp r12b, '-'
    jne integer_part_input_skip
    inc r14b
    jmp integer_part_input

integer_part_input: ; Input part before dot (integer part)
    call getchar
    mov  r12, rax
    integer_part_input_skip:

    ; check if it's line feed
    cmp r12b, 0Ah
    je read_line_into_float_end

    ; check if it's space
    cmp r12b, 20h
    je read_line_into_float_end

    ; check if it's dot
    cmp r12b, 2eh
    je char_dot

    ; Process digit
    push 10
    fild QWORD [rsp]
    fmulp 

    ; Get digit from ASCII
    sub r12b, '0'

    push r12
    fild QWORD [rsp]
    faddp

    inc rsi
    jmp integer_part_input

char_dot: ; If we found dot, switch to inputing decimal part
    xor r9, r9
    inc r9

    jmp fractional_part_input

fractional_part_input: ; Input the part after the dot
    call getchar
    mov  r12, rax

    ; check if it's line feed
    cmp r12b, 0Ah
    je fractional_part_input_end

    ; check if it's space
    cmp r12b, 20h
    je fractional_part_input_end

    sub r12b, '0'
    push QWORD 10
    fild QWORD [rsp]
    fmulp

    push r12
    fild QWORD [rsp]
    faddp

    imul r9, QWORD 10

    jmp fractional_part_input

fractional_part_input_end:
    push r9
    fild QWORD [rsp]
    fdivp

read_line_into_float_end:
    fstp QWORD [rsp]
    pop rax
    cmp r14b, 1
    jne read_float_return

    mov r9, 0x8000000000000000
    or rax, r9

    read_float_return:
    mov rsp, rbp
    pop rbp
    ret

; Reads an unsigned 64-bit int from STDIN to RAX
read_int:
    push rbp
    mov  rbp, rsp
    
    xor  rbx, rbx
    xor  rax, rax

    read_int_loop:
        call getchar
        cmp  al, 0Ah
        je   read_int_end
        cmp  al, 20h
        je   read_int_end
        
        sub  al, 0x30

        imul rbx, 10
        add  rbx, rax
        jmp  read_int_loop

    read_int_end:
    mov  rax, rbx
    mov  rsp, rbp
    pop  rbp
    ret