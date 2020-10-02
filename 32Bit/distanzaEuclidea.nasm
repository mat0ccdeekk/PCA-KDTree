;Funzione distance in assembly x86-32+SSE - Località degli accessi, Loop Vectorization con p=4, Loop Unroll con UNROLL=4
;VERSIONE CON ACCESSI NON ALLINEATI ---> PRESTAZIONI MOLTO INFERIORI A asmdistance2 (loop vectorization e basta con accessi allineati)
%include "sseutils.nasm"
section .data
  p  equ  4
  UNROLL  equ  4

section .bss
  
section .text

global distanzaEuclidea  ;rende la funzione visibile all'esterno

  ;Posizione dei parametri nel record di attivazione della funzione (i primi 8 byte sono occupati dall'indirizzo di ritorno a e da ebp)
  P  equ  8  ;puntatore a float, occupa 32 bit (=4 byte)
  Q  equ   12  ;puntatore a float, occupa 32 bit (=4 byte)
  d   equ  16  ;intero a 32 bit
  dist  equ  20  ;puntatore a float del risultato
distanzaEuclidea:
    ; ------------------------------------------------------------
    ; Sequenza di ingresso nella funzione
    ; ------------------------------------------------------------
    push  ebp        ; salva il Base Pointer
    mov  ebp, esp      ; il Base Pointer punta al Record di Attivazione corrente
    push  ebx        ; salva i registri da preservare
    push  esi
    push  edi
  
    ; ------------------------------------------------------------
    ; legge i parametri dal Record di Attivazione corrente
    ; ------------------------------------------------------------
    mov eax,[ebp+P]      ;legge P
    mov edx,[ebp+Q]      ;legge Q
    mov esi,[ebp+d]      ;legge d
    
    xorps xmm1,xmm1      ;azzero la somma    
    xor edi,edi
    mov ebx,0      ;offset
    mov ecx,esi      ;ecx=d
  bp1:  shr ecx,4      ;contatore numero di cicli ecx=(d/p*UNROLL)
    jz  try_LoopVect
    
    
loopUnroll_Q:   
    movups xmm0,[eax+ebx]    ;sposta p[i..i+p-1] in xmm0 partendo con offset 0
    movups xmm2,[edx+ebx]
  u0:  subps  xmm0,xmm2    ;xmm0[i..i+p-1]=p[i..i+p-1]-q[i..i+p-1]
    mulps  xmm0,xmm0    ;xmm0[i..i+p-1]=xmm0[i..i+p-1]*xmm0[i..i+p-1]
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive. In xmm0[0] c'è la somma di tutte e 4 le sottrazioni
    addss  xmm1,xmm0    ;aggiungo la somma delle 4 differenze in s
    
    movups xmm0,[eax+ebx+16]    ;sposta p[i..i+p-1] in xmm0 partendo con offset 0
    movups xmm2,[edx+ebx+16]
  u1:  subps  xmm0,xmm2    ;xmm0[i..i+p-1]=p[i..i+p-1]-q[i..i+p-1]
    mulps  xmm0,xmm0    ;xmm0[i..i+p-1]=xmm0[i..i+p-1]*xmm0[i..i+p-1]
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive. In xmm0[0] c'è la somma di tutte e 4 le sottrazioni
    addss  xmm1,xmm0    ;aggiungo la somma delle 4 differenze in s
    
    movups xmm0,[eax+ebx+32]    ;sposta p[i..i+p-1] in xmm0 partendo con offset 0
    movups xmm2,[edx+ebx+32]
  u2:  subps  xmm0,xmm2    ;xmm0[i..i+p-1]=p[i..i+p-1]-q[i..i+p-1]
    mulps  xmm0,xmm0    ;xmm0[i..i+p-1]=xmm0[i..i+p-1]*xmm0[i..i+p-1]
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive. In xmm0[0] c'è la somma di tutte e 4 le sottrazioni
    addss  xmm1,xmm0    ;aggiungo la somma delle 4 differenze in s
    
    movups xmm0,[eax+ebx+48]    ;sposta p[i..i+p-1] in xmm0 partendo con offset 0
    movups xmm2,[edx+ebx+48]
  u3:  subps  xmm0,xmm2    ;xmm0[i..i+p-1]=p[i..i+p-1]-q[i..i+p-1]
    mulps  xmm0,xmm0    ;xmm0[i..i+p-1]=xmm0[i..i+p-1]*xmm0[i..i+p-1]
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive. In xmm0[0] c'è la somma di tutte e 4 le sottrazioni
    addss  xmm1,xmm0    ;aggiungo la somma delle 4 differenze in s
    
    add ebx,64        ;puntatore base p*UNROLL*dim
    dec ecx          ;decremento il contatore del numero di cicli
    jnz loopUnroll_Q
loopUnroll_R_init:    
    ;DEVO FARE: mov ecx, d-(d/p*UNROLL)*(p*UNROLL)
    mov ecx,esi        ;ecx=d
    shr esi,4        ;esi=d/(p*UNROLL)(intero)
    shl esi,4        ;esi=(d/p*UNROLL)*(p*UNROLL)(intero)
    sub ecx,esi        ;ecx=d-(d/p*UNROLL)*(p*UNROLL)
    jz  end  
loopUnroll_R:  
    movss xmm0,[eax+ebx]
    subss xmm0,[edx+ebx]
    mulss xmm0,xmm0
    addss xmm1,xmm0
    add ebx,dword 4
    dec ecx
    jnz loopUnroll_R
    jmp end

;Nel caso in cui l'UNROLL non sia applicabile, proviamo ad applicare solo la Loop Vectorization

try_LoopVect:  mov ecx,esi      ;ecx=d
    shr ecx,2      ;contatore numero di cicli ecx=(d/p)
    jz  loopVect_R_init
    xorps xmm1,xmm1      ;azzero la somma    
    xor edi,edi
    
loopVect_Q:   
    movups xmm0,[eax+ebx]    ;sposta p[i..i+p-1] in xmm0 partendo con offset 0
    subps  xmm0,[edx+ebx]    ;xmm0[i..i+p-1]=p[i..i+p-1]-q[i..i+p-1]
    mulps  xmm0,xmm0    ;xmm0[i..i+p-1]=xmm0[i..i+p-1]*xmm0[i..i+p-1]
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive
    haddps xmm0,xmm0    ;sommo gli elementi di xmm0 in posizioni consecutive. In xmm0[0] c'è la somma di tutte e 4 le sottrazioni
    addss  xmm1,xmm0    ;aggiungo la somma delle 4 differenze in s
    
    add ebx,16      ;puntatore base p*dim
    dec ecx        ;decremento il contatore del numero di cicli
    jnz loopVect_Q
loopVect_R_init:    
    ;DEVO FARE: mov ecx, d-(d/p)*p
    mov ecx,esi        ;ecx=d
    shr esi,2        ;esi=d/p(intero)
    shl esi,2        ;esi=(d/p)*p(intero)
    sub ecx,esi        ;ecx=d-(d/p)*p

loopVect_R:  jz  end   
    movss xmm0,[eax+ebx]
    subss xmm0,[edx+ebx]
    mulss xmm0,xmm0
    addss xmm1,xmm0
    add ebx,dword 4
    dec ecx
    jnz loopVect_R

end:
    sqrtss xmm1,xmm1
    mov eax,[ebp+dist]
    movss [eax],xmm1

    ; ------------------------------------------------------------
    ; Sequenza di uscita dalla funzione
    ; ------------------------------------------------------------

    pop  edi            ; ripristina i registri da preservare
    pop  esi
    pop  ebx
    mov  esp, ebp          ; ripristina lo Stack Pointer
    pop  ebp            ; ripristina il Base Pointer
    ret              ; torna alla funzione C chiamante