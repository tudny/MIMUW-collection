
global so_emul

section .text

; R8b is the code counter
so_emul:

        push    r12
        push    r13
        push    r14
        push    r15

        or      rdi, rdi        ; uint16_t* code
        or      rsi, rsi        ; uint8_t*  data
        or      rdx, rdx        ; size_t    steps
        or      rcx, rcx        ; size_t    core

        mov     [rel code], rdi
        mov     [rel data], rsi
        mov     [rel steps], rdx
        mov     [rel core], rcx

        xchg    rcx, rdx

;        mov     BYTE [rel A], 3
;        mov     BYTE [rel D], 5
;        mov     BYTE [rel X], 7
;        mov     BYTE [rel Y], 9
;        mov     BYTE [rel PC], 11
;        ; unused
;        mov     BYTE [rel C], 1
;        mov     BYTE [rel Z], 1

        ; Set R8 to 0 (PC register)

        test    rcx, rcx
        jz      .end_of_code

        xor     r8, r8
        mov     r8b, [rel PC]

.code_loop:

        xor     rax, rax ; TODO remove
        mov     ax, [rdi + r8 * 2]
        call    handle_instruction

        inc     r8b
        loop    .code_loop

.end_of_code:
        mov     [rel PC], r8b
        mov     rax, [rel A]

        pop     r15
        pop     r14
        pop     r13
        pop     r12

        ret


; Modifies R9
handle_instruction:

        ; TODO ax
        xor     r9, r9
        mov     r9w, ax
        shr     r9w, 14

        cmp     r9w, 0x00
        je      handle_group_0
        cmp     r9w, 0x01
        je      handle_group_1
        cmp     r9w, 0x02
        je      handle_group_2
        cmp     r9w, 0x03
        je      handle_group_3

; ============================================================================
;                   GET      ARG        1
; ============================================================================
; Code in R10, result in R12
; R14 <- temp
get_arg_1:
        xor     r12, r12
        xor     r14, r14

        ; TODO można krócej 0,1,2,3
        cmp     r10, 0
        jne     .not_0
        mov     r12b, [rel A]
        ret

.not_0:
        cmp     r10, 1
        jne     .not_1
        mov     r12b, [rel D]
        ret

.not_1:
        cmp     r10, 2
        jne     .not_2
        mov     r12b, [rel X]
        ret

.not_2:
        cmp     r10, 3
        jne     .not_3
        mov     r12b, [rel Y]
        ret

.not_3:
        cmp     r10, 4
        jne     .not_4
        mov     r14b, [rel X]
        mov     r12b, [rsi + r14]
        ret

.not_4:
        cmp     r10, 5
        jne     .not_5
        mov     r14b, [rel Y]
        mov     r12b, [rsi + r14]
        ret

.not_5:
        cmp     r10, 6
        jne     .not_6
        mov     r14b, [rel X]
        add     r14b, [rel D]
        mov     r12b, [rsi + r14]
        ret

.not_6:
        ; must be 7 then
        mov     r14b, [rel Y]
        add     r14b, [rel D]
        mov     r12b, [rsi + r14]
        ret


; ============================================================================
;                   GET      ARG        2
; ============================================================================
; Code in R11, result in R13
; R14 <- temp
get_arg_2:
        xor     r13, r13
        xor     r14, r14

        ; TODO można krócej 0,1,2,3
        cmp     r11, 0
        jne     .not_0
        mov     r13b, [rel A]
        ret

.not_0:
        cmp     r11, 1
        jne     .not_1
        mov     r13b, [rel D]
        ret

.not_1:
        cmp     r11, 2
        jne     .not_2
        mov     r13b, [rel X]
        ret

.not_2:
        cmp     r11, 3
        jne     .not_3
        mov     r13b, [rel Y]
        ret

.not_3:
        cmp     r11, 4
        jne     .not_4
        mov     r14b, [rel X]
        mov     r13b, [rsi + r14]
        ret

.not_4:
        cmp     r11, 5
        jne     .not_5
        mov     r14b, [rel Y]
        mov     r13b, [rsi + r14]
        ret

.not_5:
        cmp     r11, 6
        jne     .not_6
        mov     r14b, [rel X]
        add     r14b, [rel D]
        mov     r13b, [rsi + r14]
        ret

.not_6:
        ; must be 7 then
        mov     r14b, [rel Y]
        add     r14b, [rel D]
        mov     r13b, [rsi + r14]
        ret
; ============================================================================

; ============================================================================
;                   SET      ARG        1
; ============================================================================
; Code in R10, value in R12
; R14 <- temp
set_arg_1:
        xor     r14, r14
        ; TODO można krócej 0,1,2,3
        cmp     r10, 0
        jne     .not_0
        mov     [rel A], r12b
        ret

.not_0:
        cmp     r10, 1
        jne     .not_1
        mov     [rel D], r12b
        ret

.not_1:
        cmp     r10, 2
        jne     .not_2
        mov     [rel X], r12b
        ret

.not_2:
        cmp     r10, 3
        jne     .not_3
        mov     [rel Y], r12b
        ret

.not_3:
        cmp     r10, 4
        jne     .not_4
        mov     r14b, [rel X]
        mov     [rsi + r14], r12b
        ret

.not_4:
        cmp     r10, 5
        jne     .not_5
        mov     r14b, [rel Y]
        mov     [rsi + r14], r12b
        ret

.not_5:
        cmp     r10, 6
        jne     .not_6
        mov     r14b, [rel X]
        add     r14b, [rel D]
        mov     [rsi + r14], r12b
        ret

.not_6:
        ; must be 7 then
        mov     r14b, [rel Y]
        add     r14b, [rel D]
        mov     [rsi + r14], r12b
        ret



; TODO SETcc
update_z:
        jnz     .zero_not_set
        mov     BYTE [rel Z], 1
        ret
.zero_not_set:
        mov     BYTE [rel Z], 0
        ret


; TODO SETcc
update_c:
        jnc     .carry_not_set
        mov     BYTE [rel C], 1
        ret
.carry_not_set:
        mov     BYTE [rel C], 0
        ret


set_c:
        cmp     BYTE [rel C], 1
        jne     .c_not_set
        stc
        ret
.c_not_set:
        clc
        ret


set_z:
        cmp     BYTE [rel Z], 1
        ret


; ============================================================================
;                                    HANDLERS
; ============================================================================

; Load case in r10 result in r13
get_addr:
        xor     r14, r14

        cmp     r10, 0
        jne     .not_0
        lea     r13, [rel A]
        ret

.not_0:
        cmp     r10, 1
        jne     .not_1
        lea     r13, [rel D]
        ret

.not_1:
        cmp     r10, 2
        jne     .not_2
        lea     r13, [rel X]
        ret

.not_2:
        cmp     r10, 3
        jne     .not_3
        lea     r13, [rel Y]
        ret

.not_3:
        cmp     r10, 4
        jne     .not_4
        mov     r14b, [rel X]
        lea     r13, [rsi + r14]
        ret

.not_4:
        cmp     r10, 5
        jne     .not_5
        mov     r14b, [rel Y]
        lea     r13, [rsi + r14]
        ret

.not_5:
        cmp     r10, 6
        jne     .not_6
        mov     r14b, [rel X]
        add     r14b, [rel D]
        lea     r13, [rsi + r14]
        ret

.not_6:
        ; must be 7 then
        mov     r14b, [rel Y]
        add     r14b, [rel D]
        lea     r13, [rsi + r14]
        ret


xchg_handler:

        call    get_addr                ; r10 -> r13
        mov     r12, r13                ;     -> r12
        mov     r10, r11
        call    get_addr                ; r11 -> r13

        ; swap r12, r13 via r14 and r15
        mov     r14b, [r12]
        mov     r15b, [r13]
        mov     [r12], r15b
        mov     [r13], r14b

        ret

; Modifies R9, R10, R11
handle_group_0:
        ; TODO group 00

        mov     r9w, ax
        and     r9w, 0x00ff             ; instruction code
        and     ax, 0xff00              ; args code

        movzx   r10, ax                 ; get arg1
        and     r10, 0x0700
        shr     r10, 8
        movzx   r11, ax                 ; get arg2
        and     r11, 0x3800
        shr     r11, 11

        ; XCHG has its own handler
        cmp     r9w, 0x08
        jnz     .not_xchg
        call    xchg_handler
.not_xchg:

        call    get_arg_1               ; loads arg1 (r10) into r12b
        call    get_arg_2               ; loads arg2 (r11) into r13b

        ; ======================== MOV ========================
        ; Handle MOV
        cmp     r9w, 0x00
        jnz     .not_move
        ; the meat
        mov     r12b, r13b
        jmp     .end
.not_move:

        ; ======================== OR ========================
        ; Handle OR
        cmp     r9w, 0x02
        jnz     .not_or
        ; the meat
        or      r12b, r13b
        call    update_z
        jmp     .end
.not_or:

        ; ======================== ADD ========================
        ; Handle ADD
        cmp     r9w, 0x04
        jnz     .not_and
        ; the meat
        add     r12b, r13b
        call    update_z
        jmp     .end
.not_and:

        ; ======================== SUB ========================
        ; Handle SUB
        cmp     r9w, 0x05
        jnz     .not_sub
        ; the meat
        sub     r12b, r13b
        call    update_z
        jmp     .end
.not_sub:

        ; ======================== ADC ========================
        ; Handle ADC
        cmp     r9w, 0x06
        jnz     .not_adc
        ; the meat
        call    set_c
        adc     r12b, r13b
        call    update_z
        call    update_c
        jmp     .end
.not_adc:

        ; ======================== SBB ========================
        ; Handle SBB
        cmp     r9w, 0x07
        jnz     .not_sbb
        ; the meat
        call    set_c
        sbb     r12b, r13b
        call    update_z
        call    update_c
        jmp     .end
.not_sbb:

.end:
        call    set_arg_1                ; loads arg1 (r12b) into register
        ret
; ============================================================================

handle_group_1:
        ; TODO group 01

        mov     r9w, ax
        and     r9w, 0x3800             ; instruction code
        shr     r9w, 11
        and     ax, 0x07ff              ; args code

        movzx   r10, ax                 ; get arg1
        and     r10, 0x0700
        shr     r10, 8

        movzx   r11, ax                 ; get imm8 (r11b)
        and     r11, 0x00ff

        call    get_arg_1               ; loads arg1 (r10) into r12b

        ; ======================== MOVI ========================
        ; Handle MOVI
        cmp     r9w, 0x00
        jnz     .not_movei
        ; the meat
        mov     r12b, r11b
        jmp     .end
.not_movei:

        ; ======================== XORI ========================
        ; Handle XORI
        cmp     r9w, 0x3
        jnz     .not_xori
        ; the meat
        xor     r12b, r11b
        call    update_z
        jmp     .end
.not_xori:

        ; ======================== ADDI ========================
        ; Handle ADDI
        cmp     r9w, 0x4
        jnz     .not_addi
        ; the meat
        add     r12b, r11b
        call    update_z
        jmp     .end
.not_addi:

        ; ======================== CMPI ========================
        ; Handle CMPI
        cmp     r9w, 0x5
        jnz     .not_cmpi
        ; the meat
        cmp     r12b, r11b
        call    update_z
        call    update_c
        jmp     .end
.not_cmpi:

        ; ======================== RCR ========================
        ; Handle RCR
        cmp     r9w, 0x6
        jnz     .not_rcr
        ; the meat
        call    set_c
        rcr     r12b, 1
        call    update_c
        jmp     .end
.not_rcr:

.end:

        call    set_arg_1               ; loads arg1 (r12b) into register
        ret

handle_group_2:

        and     ax, 0x0100              ; get the value (0 or 1)
        mov     [rel C], ah             ; set the C flag

        ret

handle_group_3:
        ; TODO group 11

        cmp     ax, 0xffff              ; BRK special case
        je      so_emul.end_of_code

        mov     r9w, ax
        and     ax, 0x00ff              ; load imm8 into ax - can be skipped using al
        and     r9w, 0x3f00
        shr     r9w, 8                  ; r9w is the instruction code

        ; ======================== JMP ========================
        ; Handle JMP
        cmp     r9w, 0x00
        jnz     .not_jmp
        ; the meat
        add     r8b, al
        jmp     .end
.not_jmp:

        ; ======================== JNC ========================
        ; Handle JNC
        cmp     r9w, 0x02
        jnz     .not_jnc
        ; the meat
        call    set_c
        jc      .end
        add     r8b, al
        jmp     .end
.not_jnc:

        ; ======================== JC ========================
        ; Handle JC
        cmp     r9w, 0x03
        jnz     .not_jc
        ; the meat
        call    set_c
        jnc     .end
        add     r8b, al
        jmp     .end
.not_jc:

        ; ======================== JNZ ========================
        ; Handle JNZ
        cmp     r9w, 0x04
        jnz     .not_jnz
        ; the meat
        call    set_z
        jz      .end
        add     r8b, al
        jmp     .end
.not_jnz:

        ; ======================== JZ ========================
        ; Handle JZ
        cmp     r9w, 0x05
        jnz     .not_jz
        ; the meat
        call    set_z
        jnz     .end
        add     r8b, al
        jmp     .end
.not_jz:

.end:
        ret


section .bss

A:      resb 1
D:      resb 1
X:      resb 1
Y:      resb 1
PC:     resb 1
unused: resb 1
C:      resb 1
Z:      resb 1

code:   resq 1
data:   resq 1
steps:  resq 1
core:   resq 1
