; ---------------------------------------------------------
; KDTREE con istruzioni AVX a 64 bit
; ---------------------------------------------------------
;
;
;
;Sottrazione tra due matrici elemento per elemento
;riceviamo una matrice Media matr tmp 
;Loop Vectorization con p=8, Loop Unroll con UNROLL=4


%include "sseutils64.nasm"

section .data            ; Sezione contenente dati inizializzati
p    equ    8
UNROLL  equ    4
    
section .bss            ; Sezione contenente dati non inizializzati



section .text            ; Sezione contenente il codice macchina
; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;    getmem    <size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;    fremem    <address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro    getmem    2
    mov    eax, %1
    push   eax
    mov    eax, %2
    push   eax
    call   get_block
    add    esp, 8
%endmacro

%macro    fremem    1
    push   %1
    call   free_block
    add    esp, 4
%endmacro


global SottraiMatrici

SottraiMatrici:
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

; rdi = Indirizzo di M
; rsi = Indirizzo di TMP
; rdx = dimensione

; Per comodita sposto i registri per farli combaciare con la versione a 32
;mov rax,rdi
;mov rbx,rsi

mov r9,rdx ;copio il valore intero che mi indica la dimensione
mov rax,0


;------------------------------------------------------------
; Inizio della funzione
;------------------------------------------------------------

;selezioniamo quanti cicli devono essere fatti
; a tal fine facciamo una divisione parte intera per 32
    shr r9,5
    jz provaLoopV

loopUNROLLQ:
    vmovups ymm0,[rdi+rax]
    vmovups ymm1,[rdi+rax+32]
    vmovups ymm2,[rdi+rax+64]
    vmovups ymm3,[rdi+rax+96]

    vmovups ymm4,[rsi+rax]
    vmovups ymm5,[rsi+rax+32]
    vmovups ymm6,[rsi+rax+64]
    vmovups ymm7,[rsi+rax+96]

    vsubps ymm0,ymm4
    vsubps ymm1,ymm5
    vsubps ymm2,ymm6
    vsubps ymm3,ymm7

    vmovups [rdi+rax],ymm0
    vmovups [rdi+rax+32],ymm1
    vmovups [rdi+rax+64],ymm2
    vmovups [rdi+rax+96],ymm3

    add rax,128 ;p*unroll*dimensione 8*4*4
    dec r9
    jnz loopUNROLLQ
    
    
lUResto_avvio:
; per calcolarci i resti utilizziamo l'operazione di shift per selezionare la parte intera
    mov r9,rdx
    shr r9,5
    shl r9,5
    ;troviamo efettivamente quanti cicli dobbiamo fare
    ;sottraiamo la dimensione al totale appena calcolato
    mov r8,r9
    mov r9,rdx
    sub r9,r8
    jz  fine


lUResto:
    vmovss xmm0,[rdi+rax]
    vmovss xmm1,[rsi+rax]

    vsubss xmm0,xmm1

    add rax, qword 4
    dec r9
    jnz lUResto
    jmp fine

provaLoopV:
	
    mov r9,rdx
;divisione parte intera per 8 per capire quanti cicli facciamo
    shr r9,3
    jz  loopVect_R_init

loopVecQ:
    vmovups ymm0,[rdi+rax]
    vsubps  ymm0,[rsi+rax]

    vmovups [rdi+rax],ymm0
    add rax,8*4 ; loop Vectorization*Dim
    dec r9
    jnz loopVecQ
    
loopVect_R_init:
; per calcolarci i resti utilizziamo l'operazione di shift per selezionare la parte intera
    mov rcx,rdx
    shr r9,2
    shl r9,2
    ;troviamo efettivamente quanti cicli dobbiamo fare
    ;sottraiamo la dimensione al totale appena calcolato
    mov r8,r9
    mov r9,rdx
    sub r9,r8
    jz  fine


loopVect_R:
    vmovss xmm0,[rdi+rax]
    vsubss xmm0,[rsi+rax]
    vmovss [rdi+rax],xmm0

    add rax,qword 4
    dec r9
    jnz loopVect_R

fine:

; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

popaq                    ; ripristina i registri generali
mov        rsp, rbp        ; ripristina lo Stack Pointer
pop        rbp                ; ripristina il Base Pointer
ret                        ; torna alla funzione C chiamante
    