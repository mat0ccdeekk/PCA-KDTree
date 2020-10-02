%include "sseutils.nasm"

section .data      ; Sezione contenente dati inizializzati

section .bss      ; Sezione contenente dati non inizializzati

section .text      ; Sezione contenente il codice macchina

global trasposta
  
trasposta:
  
    ; ------------------------------------------------------------
    ; Sequenza di ingresso nella funzione
    ; ------------------------------------------------------------
    push    rbp        ; salva il Base Pointer
    mov    rbp, rsp      ; il Base Pointer punta al Record di Attivazione corrente
    pushaq            ; salva i registri generali

    ; ------------------------------------------------------------
    ; I parametri sono passati nei registri
    ; ------------------------------------------------------------
    ; rdi = indirizzo della struct input

    ;RDI, RSI, RDX, RCX, R8 ed R9.

    ; rdi = Indirizzo di A
    ; rsi = Indirizzo di B
    ; rdx = Indirizzo di x
    ; rcx = Indirizzo di y

  ;------------------------------------------------------------
  ; Corpo della funzione
  ;------------------------------------------------------------
  mov rax, 0   ; i=0
fori:
  cmp rax, rdx  ; if(i==x)
  je theend
  
  mov rbx, 0  ; j=0
forj:
  cmp rbx, rcx  ; if(j==y)
  je fineforj
  
  mov rbp, rcx  ; ebp = y
  imul rbp, 4    ; y*4
  imul rbp, rax  ; y*4*i
  
ind_a:
  add rbp, rdi  ; a[y*4*i][]
  
  vmovss xmm0, [rbp+4*rbx]    ; a[y*4*i][j*dim]
  
  mov rbp, rdx  ; ebp = x
  imul rbp, 4    ; x*4
  imul rbp, rbx  ; x*4*j
ind_b:
  add rbp, rsi  ; b[x*4*j]
  
  vmovss [rbp+4*rax], xmm0    ; b[x*4*j][i*dim]
  
  inc rbx
  jmp forj
fineforj:
  inc rax
  jmp fori
  
theend:
  
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

    popaq                    ; ripristina i registri generali
    mov        rsp, rbp        ; ripristina lo Stack Pointer
    pop        rbp                ; ripristina il Base Pointer
    ret                        ; torna alla funzione C chiamante
    