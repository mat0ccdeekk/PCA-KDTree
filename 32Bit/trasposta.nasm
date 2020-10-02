%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati
 	inda equ 8
	indb equ 12
	indx equ 16
	indy equ 20
section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global trasposta
	
trasposta:
	
	; ------------------------------------------------------------
	; Sequenza di ingresso nella funzione
	; ------------------------------------------------------------

	push ebp			; salva il Base Pointer
	mov	ebp, esp		; il Base Pointer punta al Record di Attivazione corrente
	push ebx			; salva i registri da preservare
	push esi
	push edi
	
	;------------------------------------------------------------
	; Lettura dei parametri dal Record di attivazione
	;------------------------------------------------------------
	
	mov	eax, [ebp+inda]		; A
	mov ebx, [ebp+indb]		; B
	mov ecx, [ebp+indx]		; x
	mov edx, [ebp+indy]		; y	
	push ebp				; salvo ebp nello stack

	;------------------------------------------------------------
	; Corpo della funzione
	;------------------------------------------------------------
	mov esi, 0 	; i=0
fori:
	cmp esi, ecx	; if(i==x)
	je fine
	
	mov edi, 0	; j=0
forj:
	cmp edi, edx	; if(j==y)
	je fineforj
	
	mov ebp, edx	; ebp = y
	imul ebp, 4		; y*4
	imul ebp, esi	; y*4*i
	
indexa:
	add ebp, eax	; a[y*4*i][]
	
	movss xmm0, [ebp+4*edi]		; a[y*4*i][j*dim]
	
	mov ebp, ecx	; ebp = x
	imul ebp, 4		; x*4
	imul ebp, edi	; x*4*j
indexb:
	add ebp, ebx	; b[x*4*j]
	
	movss [ebp+4*esi], xmm0		; b[x*4*j][i*dim]
	
	inc edi
	jmp forj
fineforj:
	inc esi
	jmp fori
	
fine:
	
	; ------------------------------------------------------------
	; Sequenza di uscita dalla funzione
	; ------------------------------------------------------------
	pop ebp			; ripiglio ebp
	pop	edi			; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp	; ripristina lo Stack Pointer
	pop	ebp			; ripristina il Base Pointer
	ret				; torna alla funzione C chiamante
