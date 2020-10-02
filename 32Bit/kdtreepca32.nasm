; ---------------------------------------------------------
; KDTREE con istruzioni SSE a 32 bit
; ---------------------------------------------------------
;
;
;
;Questo file contiene il prodotto matriciale tra la matrice A e la matrice B
;il tutto memorizzato sulla matrice C, altri parametri necessari
;sono l'indice n che identifica la dimensione delle righe di A indice m la dimensione delle colonne di A
;e delle righe di B e l'indice o che rappresenta le colonne di B
;Loop Vectorization con p=4, Loop Unroll con UNROLL=4


%include "sseutils.nasm"

section .data            ; Sezione contenente dati inizializzati
;Posizione dei parametri nel record di attivazione della funzione (i primi 8 byte sono occupati dall'indirizzo di ritorno a e da ebp)
    posA equ 8        ;puntatore a float, occupa 32 bit (=4 byte)
    posB equ 12       ;puntatore a float, occupa 32 bit (=4 byte)
    posC equ 16       ;puntatore a float, occupa 32 bit (=4 byte)
    posN equ 20       ;intero a 32 bit
    posM equ 24       ;intero a 32 bit
    posO equ 28       ;intero a 32 bit

section .bss            ; Sezione contenente dati non inizializzati
    a resd 1          ;Usato per memorizzare la posizione di A
    b resd 1          ;Usato per memorizzare la posizione di B
    c resd 1          ;Usato per memorizzare la posizione di C
    n resd 1          ;Usato per memorizzare la posizione di n
    o resd 1          ;Usato per memorizzare la posizione di o
    m resd 1          ;Usato per memorizzare la posizione di m
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


global moltiplicaMatrici

; il prodotto Matriciale viene eseguito moltiplicando le righe di una matrice per le colonne dell'altra
; poi si sommano i risultati, il primo ciclo lo si effettua sulle colonne della seconda matrice

moltiplicaMatrici:

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

    mov eax, [ebp+posA]        ; Matrice A
    mov ebx, [ebp+posB]        ; Matrice B
    mov ecx, [ebp+posC]        ; Matrice C
    mov edx, [ebp+posN]        ; intero n
    mov edi, [ebp+posO]        ; intero o
    mov esi, [ebp+posM]        ; intero m

    mov [a], eax            ; A
    mov [b], ebx            ; B
    mov [c], ecx            ; C
    mov [n], edx            ; n
    mov [o], edi            ; o
    mov [m], esi            ; m

    push ebp                ; si mette ebp nuovamente sullo stack

;------------------------------------------------------------
; Inizio della funzione
;------------------------------------------------------------
mov ebp, 0        ; j=0
forj:
;se arriviamo all'ultima riga usciamo
    cmp ebp, [o]
    je uscita

;altrimenti inizializziamo gli indici che servono per iterare sulle matrici
    mov edx, 0
    mov esi, 0        ; inizializziamo la i


;analizziamo gli elementi che sono multipli di 16

forIQuoziente:
;verifichiamo che non siamo arrivati alla fine della seconda matrice
    cmp esi, [n]
    je fineForIQuoz


;in xmm0 memorizziamo la somma delle moltiplicazioni di un'intera riga per l'altra martice
    xorps xmm0,xmm0

    mov edi, 0        ; k=0
    mov edx, 0
    mov ecx, 16        ; p*unroll = 16
    mov eax, [m]

;idiv ECX fa la divisione del registro EAX per il registro ECX, memorizzandone la parte intera in EAX e
;in EDX il resto
    idiv ecx        ; [m]/16

;in EAX abbiamo il numero di cicli interi da fare tramite Loop vectorization e unroll
;ora bisogna capire quanti effettivamente sono i cicli da fare per intero-->k= eax*16
    imul eax,16

;si memorizza il valore in quoK
    mov [quoK], eax
    mov ebx, [a]
    mov ecx, [b]

forKQuoziente:
    cmp edi, [quoK]
    je fineForKQuoziente

;qui settiamo correttamente il valore associato alle righe,
;moltiplicando il valore di i per dimensione, poi per m cioè il numero di colonne della prima matrice
    mov edx, esi
    imul edx, 4
    imul edx, [m]

;tramite la add ci spostiamo sulla matrice A
    add edx, ebx    ; A[i*dim*m]
    mov eax, edi
    imul eax, 4

;possiamo leggere 4 volari di A, da notare che non usiamo xmm0 che sarà usato per la somma
;inoltre i valori non sono allineati
    movups xmm1, [edx+eax]
    movups xmm2, [edx+eax+16]
    movups xmm3, [edx+eax+32]
    movups xmm4, [edx+eax+48]


    mov edx, ebp
    imul edx, 4
    imul edx, [m]

;tramite la add ci spostiamo sulla matrice B
    add edx, ecx

;prendiamo ora 3 valori, poiche i registri sono finiti
    movups xmm5, [edx+eax]
    movups xmm6, [edx+eax+16]
    movups xmm7, [edx+eax+32]

;si fa una moltiplicazione per liberare un registro
    mulps xmm1, xmm5

;prendiamo cosi l'ultimo valore che serve
    movups xmm5, [edx+eax+48]

;facciamo le moltiplicazioni restanti
    mulps xmm2, xmm6
    mulps xmm3, xmm7
    mulps xmm4, xmm5

;sommiamo tutti i valori in xmm0
    addps xmm0, xmm4
    addps xmm0, xmm3
    addps xmm0, xmm2
    addps xmm0, xmm1

;aumentiamo k di 16 per passare al blocco successivo
    add edi, 16
    jmp forKQuoziente

fineForKQuoziente:
;qui verifichiamo se possiamo applicare solo la loop Vectorization e non la unroll
;edx lo utilizziamo per muoverci sulla matrice A quindi per il momento lo si azzera
    mov edx, 0
    mov ecx, 4

    mov eax, [m]
    idiv ecx
    imul eax,4

    mov [quoK], eax
    mov ebx, [a]
    mov ecx, [b]

forKQuoSec:
;inizia qui effettivamente il ciclo sugli elementi a blocchi di 4
    cmp edi, [quoK]
    je fineForKQuoSec

    mov edx, esi    ; edx = i
    imul edx, 4
    imul edx, [m]

    add edx, ebx
    mov eax, edi    ; eax = k
    imul eax, 4

;leggiamo un solo gruppo da 4 elementi sulla matrice A
    movups xmm1, [edx+eax]

    mov edx, ebp    ; edx = j
    imul edx, 4
    imul edx, [m]
    add edx, ecx

;leggiamo un solo gruppo di 4 elementi sulla matrice B
    movups xmm5, [edx+eax]

;moltiplichiamo i due gruppi e li sommiamo in xmm0
    mulps xmm5, xmm1
    addps xmm0, xmm5

;incrementiamo k di 4 elementi dato che si sta facendo solo la loop vectorization
    add edi, 4
    jmp forKQuoSec

fineForKQuoSec:

; Da qui poi si inzia ad analizare il resto, fossero ad esempio 39 elementi avremo:
; 2 cicli di forIQuoziente che analizza 16 valori per volta, poi 1 ciclo da 4 (forKQuoSec) e poi il resto che segue
forKResto:
    cmp edi, [m]
    je fineForKResto

    mov edx, esi    ; edx = i
    imul edx, 4
    imul edx, [m]
    add edx, ebx
    mov eax, edi

    imul eax, 4

;prendiamo il singolo valore su A
    movss xmm3, [edx+eax]

    mov edx, ebp    ; edx = j
    imul edx, 4
    imul edx, [m]
    add edx, ecx

;prendiamo il singolo valore su B
    movss xmm7, [edx+eax]
    ;moltiplichiamo tra loro i valori e accumuliamo i valori di somma in xmm0
    mulss xmm3, xmm7
    addss xmm0, xmm3

;incrementiamo per passare all'elemento successivo
    inc edi
    jmp forKResto

fineForKResto:

;si fa puntare ebx alla matrice C
    mov ebx, [c]

;facciamo la doppia riduzione per ottenere 4 volte ripetuto il valore che ci interessa
    haddps xmm0, xmm0
    haddps xmm0, xmm0

;troviamo il punto dove inserire il valore in C
    mov edx, esi    ; edx = i
    imul edx, 4
    imul edx, [o]
    add edx, ebx

;spostiamo il calcolo fatto nella posizione che gli compete
    movss [edx+ebp*4], xmm0

;ci spostiamo alla colonna successiva in B aumentando i
    inc esi
    jmp forIQuoziente


fineForIQuoz:
;passiamo alla riga successiva in A aumentando j
    inc ebp
    jmp forj

uscita:
;si sono analizzate tutte le righe della matrice A quindi possiamo chiudere
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------
    pop    ebp            ; ripiglio ebp
    pop    edi            ; ripristina i registri da preservare
    pop    esi
    pop    ebx
    mov    esp, ebp       ; ripristina lo Stack Pointer
    pop    ebp            ; ripristina il Base Pointer
    ret                   ; torna alla funzione C chiamante