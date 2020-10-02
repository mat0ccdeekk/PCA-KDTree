; ---------------------------------------------------------
; KDTREE con istruzioni SSE a 32 bit
; ---------------------------------------------------------
;
;
;
;Questo file contiene il prodotto scalare di un vettore, riceviamo come parametro un vettore di float e uno scalare, 
;poi uno scalare come risultato


%include "sseutils.nasm"

section .data            ; Sezione contenente dati inizializzati
;Posizione dei parametri nel record di attivazione della funzione (i primi 8 byte sono occupati dall'indirizzo di ritorno a e da ebp)
    posA equ 8        ;puntatore a float, occupa 32 bit (=4 byte)
    lenA equ 12       ;intero a 32 bit
    punRis equ 16     ;puntatore al risultato

section .bss            ; Sezione contenente dati non inizializzati
    a resd 1          ;Usato per memorizzare la posizione di A
    len resd 1          ;Usato per memorizzare la posizione di B
    ris resd 1          ;Usato per memorizzare la posizione di C

    quoK resd 1       ;Usato per memorizzare il valore quoziente della divisione tra dimensione della matrice e UNROLL*p


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

;ci calcoliamo il numero di cicli interi di loop vectorization ecx=(d/p*UNROLL)
%macro  nCicli   2  
  mov %1,[ebp+lenA]
  shr %1,%2
%endmacro

%macro  nResti 3  ;calcola il numero di cicli di resto per la code vectorization
  mov %2,[ebp+lenA]  
  mov %1,%2
  shr %2,%3
  shl %2,%3
  sub %1,%2
%endmacro

global prodottoScalare1

; il prodotto scalare lo eseguiamo prendendo il vettore 4 elementi per volta utilizzando quindi la loop Vectorization
; in più usiamo la loop UNROLL 

prodottoScalare1:

; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------

    push ebp            ; salva il Base Pointer
    mov  ebp, esp       ; il Base Pointer punta al Record di Attivazione corrente
    push ebx            ; salva i registri da preservare
    push esi
    push edi
    push edx

;------------------------------------------------------------
; Lettura dei parametri dal Record di attivazione
;------------------------------------------------------------

    mov eax, [ebp+posA]        ; Vettore A
    mov ebx, [ebp+lenA]        ; Lunghezza di A
    mov ecx, [ebp+punRis]      ; puntatore risultato


    mov [a], eax            ; A
    mov [len], ebx          ; len(A) 
    mov [ris], ecx          ; risultato

    push ebp                ; si mette ebp nuovamente sullo stack

    

;------------------------------------------------------------
; Inizio della funzione
;------------------------------------------------------------

;prima di tutto verifichiamo quanti cicli possiamo fare per intero 

    nCicli ecx,4
    jz unoXVolta
    
;puntiamo al Vettore
    mov ebx, [a]

;in xmm0 memorizziamo la somma delle moltiplicazioni 
    xorps xmm0,xmm0
    mov edi,0

    

forKQuoziente:
; controlliamo se fare ancora un ciclo
    cmp ecx,0
    je unoXVolta
    
    mov eax, edi
    imul eax, 4
    
    movups xmm1, [ebx+eax]
    movups xmm2, [ebx+eax+16]
    movups xmm3, [ebx+eax+32]
    movups xmm4, [ebx+eax+48]
    mulps xmm1,xmm1
    mulps xmm2,xmm2
    mulps xmm3,xmm3
    mulps xmm4,xmm4

    addps xmm1,xmm2
    addps xmm3,xmm4
    addps xmm1,xmm3

    haddps xmm1,xmm1
    haddps xmm1,xmm1

    addps xmm0,xmm1

    add edi,16
    dec ecx 
    jmp forKQuoziente
    
forKResto:
    cmp edx,0
    je fine

    mov eax, edi
    imul eax, 4
    
    movss xmm1, [ebx+eax]
    mulss xmm1,xmm1

    addss xmm0,xmm1

    dec edx 
    inc edi 
    jmp forKResto
    

unoXVolta:
    ;calcoliamo quanti resti 
    nResti ecx,esi,2
    jz fine


calcolo:

    cmp ecx,0
    je fine 

    mov eax,edi 
    imul eax,4

    movss xmm1, [ebx+eax]
    mulss xmm1,xmm1

    addss xmm0,xmm1

    dec ecx 
    inc edi 
    jmp calcolo


fine:
  movss [ris],xmm0
  
  ; ------------------------------------------------------------
  ; Sequenza di uscita dalla funzione
  ; ------------------------------------------------------------
  pop ebp      ; ripiglio ebp
  pop  edi      ; ripristina i registri da preservare
  pop  esi
  pop  ebx
  mov  esp, ebp  ; ripristina lo Stack Pointer
  pop  ebp      ; ripristina il Base Pointer
  ret        ; torna alla funzione C chiamante
    