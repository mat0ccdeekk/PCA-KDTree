; ---------------------------------------------------------
; KDTREE con istruzioni SSE a 32 bit
; ---------------------------------------------------------
;
;
;
;Sottrazione tra due matrici elemento per elemento
;riceviamo una matrice Media matr tmp 
;Loop Vectorization con p=4, Loop Unroll con UNROLL=4


%include "sseutils.nasm"

section .data            ; Sezione contenente dati inizializzati
;Posizione dei parametri nel record di attivazione della funzione (i primi 8 byte sono occupati dall'indirizzo di ritorno a e da ebp)
    posM    equ 8        ;puntatore a float, occupa 32 bit (=4 byte)
    posTMP  equ 12       ;puntatore a float, occupa 32 bit (=4 byte)
    posDim  equ 16       ;puntatore a float, occupa 32 bit (=4 byte)
    
section .bss            ; Sezione contenente dati non inizializzati
    m resd 1          ;Usato per memorizzare la posizione di m
    tmp resd 1          ;Usato per memorizzare la posizione della matrice tmp 
    dim resd 1          ;Usato per memorizzare la dimensione (prodotto tra righe*colonne)


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


;numero di cicli per vectorization o loop unroll ecx=(d/p*UNROLL)
%macro  nCicli   2  
  mov %1,[ebp+posDim]
  shr %1,%2
%endmacro

;calcola il numero di cicli di resto per la code vectorization
%macro  nResti 3  
  mov %2,[ebp+posDim]  
  mov %1,%2
  shr %2,%3
  shl %2,%3
  sub %1,%2
%endmacro

global SottraiMatrici

SottraiMatrici:

; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------

    push ebp            ; salva il Base Pointer
    mov  ebp, esp       ; il Base Pointer punta al Record di Attivazione corrente
    push ebx            ; salva i registri da preservare
    push esi
    push edi

;------------------------------------------------------------
; Lettura dei parametri dal Record di attivazione
;------------------------------------------------------------

    mov eax, [ebp+posM]        ; Matrice Media 
    mov ebx, [ebp+posTMP]      ; Matrice Temporanea
    mov ecx, [ebp+posDim]      ; intero dimensione Prodotto righe * colonne

    mov [dim],ecx 

;------------------------------------------------------------
; Inizio della funzione
;------------------------------------------------------------

;selezioniamo quanti cicli devono essere fatti
    nCicli ecx,4
    jz provaLoopV

loopUNROLLQ:
    movups xmm0,[eax+edi]
    movups xmm1,[eax+edi+16]
    movups xmm2,[eax+edi+32]
    movups xmm3,[eax+edi+48]

    movups xmm4,[ebx+edi]
    movups xmm5,[ebx+edi+16]
    movups xmm6,[ebx+edi+32]
    movups xmm7,[ebx+edi+48]

    subps xmm0,xmm4
    subps xmm1,xmm5
    subps xmm2,xmm6
    subps xmm3,xmm7

    movups [eax+edi],xmm0
    movups [eax+edi+16],xmm1
    movups [eax+edi+32],xmm2
    movups [eax+edi+48],xmm3

    add edi,4*4*4 ;p*unroll*dimensione
    dec ecx 
    jnz loopUNROLLQ
    
    
lUResto_avvio:    
    nResti ecx,esi,4
    jz  fine 

lUResto:
    movss xmm0,[eax+edi]
    movss xmm1,[ebx+edi]

    subss xmm0,xmm1

    add edi,4
    dec ecx 
    jnz lUResto
    jmp fine

provaLoopV:
    xor edi,edi

    nCicli ecx,2
    jz  loopVect_R_init

loopVecQ:
    movups xmm0,[eax+edi]
    subps  xmm0,[ebx+edi]

    movups [eax+edi],xmm0
    add edi,4*4
    dec ecx
    jnz loopVecQ

loopVect_R_init:
    nResti ecx,esi,2
    jz  fine

loopVect_R:
    movss xmm0,[eax+edi]
    subss xmm0,[ebx+edi]
    movss [eax+edi],xmm0

    add ebx,4
    dec ecx
    jnz loopVect_R

fine:
    ; ------------------------------------------------------------
    ; Sequenza di uscita dalla funzione
    ; ------------------------------------------------------------

    pop  edi            ; ripristina i registri da preservare
    pop  esi
    pop  ebx
    mov  esp, ebp      ; ripristina lo Stack Pointer
    pop  ebp            ; ripristina il Base Pointer
    ret                ; torna alla funzione C chiamante
    