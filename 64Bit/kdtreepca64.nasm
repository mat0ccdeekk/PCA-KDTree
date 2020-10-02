
; ---------------------------------------------------------
; KDTREE con istruzioni AVX a 64 bit
; ---------------------------------------------------------
;
;
;
;Questo file contiene il prodotto matriciale tra la matrice A e la matrice B
;il tutto memorizzato sulla matrice C, altri parametri necessari
;sono l'indice n che identifica la dimensione delle righe di A indice m la dimensione delle colonne di A
;e delle righe di B e l'indice o che rappresenta le colonne di B
;Loop Vectorization con p=4, Loop Unroll con UNROLL=4



section .data      ; Sezione contenente dati inizializzati
p    equ    8
UNROLL  equ    4
;
;align 32
;vec1:    dd    1.0, 2.0, 3.0, 4.0

section .bss      ; Sezione contenente dati non inizializzati
a: resq 1
b: resq 1
c: resq 1
n: resq 1
m: resq 1
o: resq 1
q: resq 1


;alignb 32

section .text      ; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;  getmem  <size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;  fremem  <address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro  getmem  2
mov  rdi, %1
mov  rsi, %2
call  get_block
%endmacro

%macro  fremem  1
mov  rdi, %1
call  free_block
%endmacro


global moltiplicaMatrici


moltiplicaMatrici:
; ------------------------------------------------------------
; Sequenza di ingresso nella funzione
; ------------------------------------------------------------
push    rbp        ; salva il Base Pointer
mov    rbp, rsp      ; il Base Pointer punta al Record di Attivazione corrente
;pushaq            ; salva i registri generali

; ------------------------------------------------------------
; I parametri sono passati nei registri
; ------------------------------------------------------------
; rdi = indirizzo della struct input

;RDI, RSI, RDX, RCX, R8 ed R9.

; rdi = Indirizzo di A
; rsi = Indirizzo di B
; rdx = Indirizzo di C
; rcx = Indirizzo di n
; r8  = indirizzo di m
; r9  = Indirizzo di o

    mov [a], rdi            ; A
    mov [b], rsi            ; B
    mov [c], rdx            ; C
    mov [n], rcx            ; n
    mov [o], r8             ; o
    mov [m], r9             ; m
; ------------------------------------------------------------
; Inizio della funzione
; ------------------------------------------------------------

mov rbp, 0  ; j=0
forj:
;se arriviamo all'ultima riga usciamo
cmp rbp,[o]
je uscita

;altrimenti inizializziamo gli indici che servono per iterare sulle matrici
mov rdx,0 ; edx
mov rsi,0 ;inizializziamo la i = 0(esi)

;analizziamo gli elementi che sono multipli di 32

forIQuoziente:
;verifichiamo che non siamo arrivati alla fine della seconda matrice
cmp rsi,[n]
je fineForIQuoz

;in ymm0 memorizziamo la somma delle moltiplicazioni di un'intera riga per l'altra martice
vxorps ymm0,ymm0

mov rdi,0 ;k = 0
mov rdx,0 ;edx
mov rax,0

;per capire quanti cicli trasliamo in rcx la dimensione e facciamo lo shift destro di 5
mov rcx,0 ;ecx
mov rcx,[m]
shr rcx,5

mov rax,rcx
imul rax,32

;si memorizza il valore in quoK
mov [q], rax
mov rbx, [a]
mov rcx, [b]

forKQuoziente:
cmp rdi,[q]
je fineForKQuoziente


;qui settiamo correttamente il valore associato alle righe,
;moltiplicando il valore di i per dimensione, poi per m cioè il numero di colonne della prima matrice
mov rdx, rsi
imul rdx, 4
imul rdx, [m]

;tramite la add ci spostiamo sulla matrice A
add rdx, rbx
mov rax, rdi
imul rax, 4

;possiamo leggere 32 volari di A, da notare che non usiamo xmm0 che sarà usato per la somma
;inoltre i valori non sono allineati
vmovups ymm1, [rdx+rax]
vmovups ymm2, [rdx+rax+32]
vmovups ymm3, [rdx+rax+64]
vmovups ymm4, [rdx+rax+96]
mov rdx, rbp
imul rdx, 4
imul rdx,[m]

;tramite la add ci spostiamo sulla matrice B
add rdx, rcx

;prendiamo ora 3 valori, poiche i registri sono finiti
vmovups ymm5, [rdx+rax]
vmovups ymm6, [rdx+rax+32]
vmovups ymm7, [rdx+rax+64]

;si fa una moltiplicazione per liberare un registro
vmulps ymm1, ymm5

;prendiamo cosi l'ultimo valore che serve
vmovups ymm5, [rdx+rax+96]

;facciamo le moltiplicazioni restanti
vmulps ymm2, ymm6
vmulps ymm3, ymm7
vmulps ymm4, ymm5

;sommiamo tutti i valori in xmm0
vaddps ymm0, ymm4
vaddps ymm0, ymm3
vaddps ymm0, ymm2
vaddps ymm0, ymm1

;aumentiamo k di 32 per passare al blocco successivo
add rdi, 32
jmp forKQuoziente


fineForKQuoziente:
;qui verifichiamo se possiamo applicare solo la loop Vectorization e non la unroll
;rdx lo utilizziamo per muoverci sulla matrice A quindi per il momento lo si azzera
mov rdx, 0
mov rcx, 8

mov rax, [m]

;divido per 8
shr rax,3
imul rax,8

mov [q], rax
mov rbx, [a]
mov rcx, [b]

forKQuoSec:
;inizia qui effettivamente il ciclo sugli elementi a blocchi di 4
cmp rdi, [q]
je fineForKQuoSec

mov rdx, rsi    ; edx = i
imul rdx, 4
imul rdx, [m]

add rdx, rbx
mov rax, rdi    ; eax = k
imul eax, 4

;leggiamo un solo gruppo da 8 elementi sulla matrice A
vmovups ymm1, [rdx+rax]

mov rdx, rbp    ; rdx = j
imul rdx, 8
imul rdx, [m]
add rdx, rcx

;leggiamo un solo gruppo di 8 elementi sulla matrice B
vmovups ymm5, [rdx+rax]

;moltiplichiamo i due gruppi e li sommiamo in xmm0
vmulps ymm5, ymm1
vaddps ymm0, ymm5

;incrementiamo k di 8 elementi dato che si sta facendo solo la loop vectorization
add edi, 8
jmp forKQuoSec

fineForKQuoSec:
forKResto:
cmp rdi, [m]
je fineForKResto

mov rdx, rsi    ; edx = i
imul rdx, 8
imul rdx, [m]
add rdx, rbx
mov rax, rdi

imul rax, 4
;prendiamo il singolo valore su A
vmovss xmm3, [rdx+rax]

mov rdx, rbp    ; edx = j
imul rdx, 8
imul edx, [m]
add rdx, rcx

;prendiamo il singolo valore su B
vmovss xmm7, [rdx+rax]
;moltiplichiamo tra loro i valori e accumuliamo i valori di somma in xmm0
vmulss xmm3, xmm7
vaddss xmm0, xmm3

;incrementiamo per passare all'elemento successivo
inc edi
jmp forKResto


fineForKResto:

;si fa puntare ebx alla matrice C
mov rbx, [c]


;facciamo la doppia riduzione per ottenere 4 volte ripetuto il valore che ci interessa
vhaddps ymm0, ymm0
vhaddps ymm0, ymm0

;troviamo il punto dove inserire il valore in C
mov rdx, rsi    ; edx = i
imul rdx, 4 ; CONTROLLARE
imul rdx, [o]
add rdx, rbx

;spostiamo il calcolo fatto nella posizione che gli compete
vmovss [rdx+rbp*4], xmm0


;ci spostiamo alla colonna successiva in B aumentando i
inc rsi
jmp forIQuoziente

fineForIQuoz:
;passiamo alla riga successiva in A aumentando j
inc rbp
jmp forj

uscita:
; ------------------------------------------------------------
; Sequenza di uscita dalla funzione
; ------------------------------------------------------------

;popaq            ; ripristina i registri generali
mov    rsp, rbp      ; ripristina lo Stack Pointer
pop    rbp          ; ripristina il Base Pointer
ret              ; torna alla funzione C chiamante