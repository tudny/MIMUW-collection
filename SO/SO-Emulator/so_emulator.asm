; so emulator solution

; invariants:
; 1* r8b is the steps counter (0, ..., 255) PC
; 2* r9 is the pointer to REGISTER for this core
; 3* rdi is the pointer to the original code (uint16_t*)
; 4* rsi is the pointer to the data (uint8_t*)
; 5* rcx is the steps count (decreasing) (size_t)
; 6* rax is mostly used as temporary register
; 7* each instruction is being held in rdx (after value of cores is used to calculate memory)

global so_emul

section .text

so_emul:
        push    r12
        push    r13
        push    r14

        ; Initiate
        xchg    rcx, rdx                    ; invariant 5 and 6
        lea     r9, [rel REGISTERS]         ; invariant 2
        lea     r9, [r9 + rdx * 8]          ; += core * 8
        mov     r10, 8                      ; find PC value for this core
        call    get_addr_case
        xor     r8, r8                      ; we need to clear first 64-8 bits
        mov     r8b, [r10]                  ; initiate PC register
        ;

        test    rcx, rcx                    ; zero steps to be done (almost-infinite loop possible)
        jz      .end_of_so_emul

.main_loop:

        xor     rdx, rdx
        mov     dx, [rdi + r8 * 2]
        inc     r8b                         ; we use first 8-bit, so it will overflow properly
        jmp     handle_instruction
.back_to_loop:

        loop    .main_loop

.end_of_so_emul:

        ; Save the result into RAX
        mov     r10, 0                      ; get A
        call    get_addr_case               ; r10 -> A register address
        mov     [r10 + 4], r8b              ; load current instruction counter (r8b) to PC (&A + 4)
        mov     rax, [r10]                  ; load all registers into RAX (64-bit)
        ;

        pop     r14
        pop     r13
        pop     r12

        ret


; Changes R10, R11
; R10 <- load case (0, 1, 2, 3, 4, 5, 6, 7)
; We extend all task cases for more useful cases
; 0  <- A
; 1  <- D
; 2  <- X
; 3  <- Y
; 4  <- [X]
; 5  <- [Y]
; 6  <- [X + D]
; 7  <- [Y + D]
; -----
; 8  <- PC
; 9  <- C
; 10 <- Z
;
; R10 -> result (address)
get_addr_case:
        cmp     r10, 4                      ; A, D, X, Y
        jnb     .not_first_registers        ; diff then /\
        ; If I remove this empty line next jmp is never done (WHY???)
        jmp     .ret_reg_addr               ; REGISTERS + core*8 + offset

.not_first_registers:
        cmp     r10, 8                      ; PC, C, Z
        jnae    .not_second_registers       ; memory access via [*]
        je      .pc_reg
        sub     r10, 3  ; 9->6, 10->7       ; C, Z
        jmp     .ret_reg_addr

.pc_reg:
        mov     r10, 4                      ; PC
        ; jmp     .ret_reg_addr, but no need to

.ret_reg_addr:
        lea     r10, [r9 + r10]             ; depends of value of core
        ret

.not_second_registers:                      ; [X], [Y], [X + D], [Y + D]
        cmp     r10, 5
        mov     r11, 0                      ; We will store the value of D if [* + D], 0 otherwise
        jna      .delta_zero
        push    r10                         ; [X + D], [Y + D] - we save r10
        mov     r10, 1                      ; call to get the pointer to D
        call    get_addr_case               ; r10 -> pointer to D
        mov     r11b, [r10]                 ; value of D
        pop     r10                         ; restore value or r10
        sub     r10, 2                      ; [X + D](6 -> 4 [X]), [Y + D](7 -> 5 [Y])

.delta_zero:
        ; [X](4 -> 2 X), [Y](5 -> 3 Y)
        sub     r10, 2
        call    get_addr_case               ; R10 -> pointer to X or Y
        movzx   r10, BYTE [r10]             ; we get the value of X or Y
        add     r10b, r11b                  ; we add D (or 0 if not relevant)
        lea     r10, [rsi + r10]            ; &data[r10 = (X|Y + D)]
        ret


load_z:
        cmp     BYTE [r9 + 7], 1
        ret

load_c:
        mov     rax, 0
        cmp     al, BYTE [r9 + 6]
        ret

save_z:
        jnz     .zero_not_set
        mov     BYTE [r9 + 7], 1
        ret
.zero_not_set:
        mov     BYTE [r9 + 7], 0
        ret

save_c:
        mov     BYTE [r9 + 6], 0
        rcl     BYTE [r9 + 6], 1
        ret



; The instruction code is in DX (16-bit)
handle_instruction:
        cmp     dx, 0xffff                  ; handle BRK
        je      so_emul.end_of_so_emul

        movzx   r12, dx
        and     r12, 0x3800
        shr     r12, 11                     ; First part (first 3-bit)

        movzx   r13, dx
        and     r13, 0x0700
        shr     r13, 8                      ; Second part (first 3-bit)

        movzx   r14, dx
        and     r14, 0x00ff                 ; Third part (first 8-bit)

        shr     dx, 14                      ; Group (first 2-bit)

        je      handle_group_0

        cmp     dx, 0x01
        je      handle_group_1

        cmp     dx, 0x02
        je      handle_group_2

        cmp     dx, 0x03
        je      handle_group_3

        jmp     so_emul.back_to_loop


handle_group_0:
        mov     r10, r13                    ; Get arg1
        call    get_addr_case
        mov     r13, r10

        mov     r10, r12                    ; Get arg2
        call    get_addr_case

        mov     dl, BYTE [r10]                  ; Move to rdx to make operations

        cmp     r14, 0x00                   ; Check MOV
        jnz     .not_mov
        ; ============== MOV ==============
        mov     [r13], dl
        jmp     .end
.not_mov:

        cmp     r14, 0x02                   ; Check OR
        jnz     .not_or
        ; ============== OR ==============
        or      [r13], dl
        call    save_z
        jmp     .end
.not_or:

        cmp     r14, 0x04                   ; Check ADD
        jnz     .not_add
        ; ============== ADD ==============
        add     [r13], dl
        call    save_z
        jmp     .end
.not_add:

        cmp     r14, 0x05                   ; Check SUB
        jnz     .not_sub
        ; ============== SUB ==============
        sub     [r13], dl
        call    save_z
        jmp     .end
.not_sub:

        cmp     r14, 0x06                   ; Check ADC
        jnz     .not_adc
        ; ============== ADC ==============
        call    load_c
        adc     [r13], dl
        call    save_z
        call    save_c
        jmp     .end
.not_adc:

        cmp     r14, 0x07                   ; Check SBB
        jnz     .not_sbb
        ; ============== SBB ==============
        call    load_c
        sbb     [r13], dl
        call    save_z
        call    save_c
        jmp     .end
.not_sbb:

        cmp     r14, 0x08                   ; Check XCHG
        jnz     .not_xchg
        ; ============== XCHG ==============

        mov     al, [r10]
        lock xchg    [r13], al
        mov     [r10], al

        ; jmp     .end, but no need
.not_xchg:

.end:
        jmp     so_emul.back_to_loop


handle_group_1:
        mov     r10, r13                    ; Get arg1
        call    get_addr_case

        cmp     r12, 0x0                    ; Check MOVI
        jnz     .not_movi
        ; ============== MOVI ==============
        mov     [r10], r14b
        jmp     .end
.not_movi:

        cmp     r12, 0x3                    ; Check XORI
        jnz     .not_xori
        ; ============== XORI ==============
        xor     [r10], r14b
        call    save_z
        jmp     .end
.not_xori:

        cmp     r12, 0x4                    ; Check ADDI
        jnz     .not_addi
        ; ============== ADDI ==============
        add     [r10], r14b
        call    save_z
        jmp     .end
.not_addi:

        cmp     r12, 0x5                    ; Check CMPI
        jnz     .not_cmpi
        ; ============== CMPI ==============
        cmp     [r10], r14b
        call    save_z
        call    save_c
        jmp     .end
.not_cmpi:

        cmp     r12, 0x6                    ; Check RCR
        jnz     .not_rcr
        ; ============== RCR ==============
        call    load_c
        rcr     BYTE [r10], 1
        call    save_c
        ; jmp     .end, but no need
.not_rcr:

.end:
        jmp     so_emul.back_to_loop


handle_group_2:
        mov     [r9 + 6], r13b
        jmp     so_emul.back_to_loop


handle_group_3:
        cmp     r13, 0x0                    ; Check JMP
        jne     .not_jmp
        ; ============== JMP ==============
        add     r8b, r14b
        jmp     .end
.not_jmp:

        cmp     r13, 0x2
        jne     .not_jnc
        ; ============== JNC ==============
        call    load_c
        jc      .end
        add     r8b, r14b
        jmp     .end
.not_jnc:

        cmp     r13, 0x3
        jne     .not_jc
        ; ============== JC ==============
        call    load_c
        jnc     .end
        add     r8b, r14b
        jmp     .end
.not_jc:

        cmp     r13, 0x4
        jne     .not_jnz
        ; ============== JNZ ==============
        call    load_z
        jz      .end
        add     r8b, r14b
        jmp     .end
.not_jnz:

        cmp     r13, 0x5
        jne     .not_jz
        ; ============== JZ ==============
        call    load_z
        jnz     .end
        add     r8b, r14b
        ; jmp     .end, but no need
.not_jz:

.end:
        jmp     so_emul.back_to_loop


section .bss

; Registers are ordered like that:
; A, D, X, Y, PC, unused, C, Z
; They hold 64-bit space of each core
; for core M they start at CORES + M and are one after another
REGISTERS: resq CORES
