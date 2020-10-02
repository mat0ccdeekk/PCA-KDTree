/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2018/19
* 
* Progetto dell'algoritmo di Product Quantization for Nearest Neighbor Search
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 kdtreepca32.nasm && gcc -O0 -m32 -msse kdtreepca32.o kdtreepca32c.c -o kdtreepca32c && ./kdtreepca32c
* 
* oppure
* 
* ./runkdtreepca32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

struct NodoAlbero{
	int id_vicino;
	float* vetRadice;
	float* vetCoppie;
	struct NodoAlbero *sinistro, *destro;
};//Struttura per costruire il KDTREE
typedef struct NodoAlbero *Albero;

#define	MATRIX		float*
#define	KDTREE		float* // modificare con il tipo di dato utilizzato

typedef struct {
    char* filename; //nome del file, estensione .ds per il data set, estensione .qs per l'eventuale query set
    MATRIX ds; //data set
    MATRIX qs; //query set
    int n; //numero di punti del data set
    int k; //numero di dimensioni del data/query set
    int nq; //numero di punti del query set
    int h; //numero di componenti principali da calcolare 0 se PCA non richiesta
    int kdtree_enabled; //1 per abilitare la costruzione del K-d-Tree, 0 altrimenti
    KDTREE kdtree; //riferimento al K-d-Tree, NULL se costruzione non richiesta
    float r; //raggio di query, -1 se range query non richieste
    int silent; //1 per disabilitare le stampe, 0 altrimenti
    int display; //1 per stampare i risultati, 0 altrimenti
    MATRIX dsMedia; //matrcie del data set centrato rispetto alla media
    MATRIX dsT; //matrice trasposta del data/query set centrato rispetto alla media
    MATRIX U; //matrice U restituita dall'algoritmo PCA
    MATRIX V; //matrice V restituita dall'algoritmo PCA
    MATRIX pcaQS; //query set con PCA
    float* vetScore; //vettore degli score
    float* vetLoad;  //vettore dei load
    //STRUTTURE OUTPUT MODIFICABILI
    int* QA; //risposte alle query in forma di coppie di interi (id_query, id_vicino)
    int nQA; //numero di risposte alle query
} params;



struct Nodo{
	int id_vicino;
	int id_query;
	struct Nodo *succ;
};//Struttura per memorizzare la coppia (id_query, id_vicino)
typedef struct Nodo *Lista;

static Lista L, inizio; //lista per output id query

// PROCEDURE ASSEMBLY
extern void trasposta(float* a, float* b, int x, int y);
extern void moltiplicaMatrici(float* matrice, float* vettore, float* vettColonna, int r, int c, int a);
extern void SottraiMatrici(float* matriceA, float* matriceB, int x);

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/


void* get_block(int size, int elements) { 
    return _mm_malloc(elements*size,16);
}

void free_block(void* p) { 
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(float),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
*
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/

MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status;
    fp = fopen(filename, "rb");
    if (fp == NULL){
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }
    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    MATRIX data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(float), rows*cols, fp);
    fclose(fp);
    *n = rows;
    *k = cols;
    return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
        fwrite(&n, 4, 1, fp);
        fwrite(&k, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, 4, k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += 4*k;
        }
    }
    fclose(fp);
}

//--------------------PCA-----------------//

void centrareDataset(params* input, float* ds, int n, int k){
	int i=0, col=0;
	float * media = alloc_matrix(1, k);
	input->dsMedia = alloc_matrix(n, k);
	if(!media){
		printf("Problema allocazione vettore dei valori medi");
	}
	for(i = 0; i < k*n; i++){
		if(i % k == 0 && i != 0) col = 0;
		media[col] = media[col] + ds[i];
		col++;
	}//somma colonne

	for(i = 0; i < k; i++){
		media[i] = media[i] / n;
	}//mediaColonne
	col=0;
	for(i = 0; i < k*n; i++){
		if(i % k == 0 && i != 0) col = 0;
		input->dsMedia[i] = ds[i] - media[col];
		col++;
	}//nuovoDataset
	//dealloc_matrix(media);
}//centrareDataset


void traspostaDataSet(params* input){
	int j=0, offset=0;
	trasposta(input->dsMedia, input->dsT, input->n, input->k);

}//traspostaDataSet

void estraiVettoreScore(params* input){
	int j = 0, z = 0;
	for(j = 0; j < input->n*input->k; j=j+input->k ){
		input->vetScore[z] = input->dsMedia[j];
		z++;
	}
	//printf("\n");
}//estraiVettoreScore


float prodottoScalare(float* vettore, int lunghezza){
	float  p = 0;
	for( int i = 0; i < lunghezza; i++)
		p = (vettore[i] * vettore[i]) + p;
	return p;
}//prodottoScalare

void calcoloVettoreLoad(params* input){
	float* numeratore = alloc_matrix(input->k, 1);
	float  denominatore = 0;
	moltiplicaMatrici(input->dsT, input->vetScore, numeratore, input->k, input->n, 1);
	//numeratore = prodottoMatriciale(input->dsT, input->vetScore, input->k, input->n);
	denominatore = prodottoScalare(input->vetScore, input->n);
	for(int i = 0; i < input->k; i++){
		input->vetLoad[i] = numeratore[i] / denominatore;
	}
	//dealloc_matrix(numeratore);
}//calcoloVettoreLoad

void normalizzaVettoreLoad(params* input){
	int i;
	float p = 0, radice = 0;
	p = prodottoScalare(input->vetLoad, input->k);
	radice = sqrt(p);
	for(i = 0; i < input->k; i++){
		input->vetLoad[i] = input->vetLoad[i] / radice;
	}
}//normalizzaVettoreLoad

void aggiornaVettoreScore(params* input){
	float* numeratore1 = alloc_matrix(input->n, 1);
	float  denominatore = 0;
	moltiplicaMatrici(input->dsMedia, input->vetLoad, numeratore1, input->n, input->k, 1);
	//numeratore1 = prodottoMatriciale(input->dsMedia, input->vetLoad, input->n, input->k);
	denominatore = prodottoScalare(input->vetLoad, input->k);
	for(int i = 0; i < input->n; i++){
		input->vetScore[i] = numeratore1[i] / denominatore;
	}
	dealloc_matrix(numeratore1);
}//aggiornaVettoreScore

void aggiornaMatriceScore(params*input , int iColonna){
	int i, iRiga = 0;
	for(i = 0; i < input->n; i++){
		input->U[iColonna + iRiga] = input->vetScore[i];
		iRiga = iRiga + input->h;
	}
}//aggiornaMatriceScore

void aggiornaMatriceLoad(params*input , int iColonna){
	int i, iRiga = 0;
	for(i = 0; i < input->k; i++){
		input->V[iColonna + iRiga] = input->vetLoad[i];
		iRiga = iRiga + input->h;
	}
}//aggiornaMatriceLoad

void aggiornaDataSet(params* input, float *matriceTemp){
	int i, j = 0, z = 0, col = 0, riga = 0;
    //prodotto vettore colonna * vettore riga = matrice
	moltiplicaMatrici(input->vetScore, input->vetLoad, matriceTemp, input->n, 1, input->k);
	//sottrazione
	SottraiMatrici(input->dsMedia, matriceTemp, input->n*input->k);
	//dealloc_matrix(matriceTemp);
}//aggiornaDataSet

/*
*	PCA
* 	=====================
*/
void pca(params* input) {
    int iColonna=0;
    float t=0, t1=0;
    float teta = 1e-6;
	float* matriceTemp = alloc_matrix(input->n,input->k);
	input->vetScore = alloc_matrix(input->n,1);
    input->vetLoad = alloc_matrix(input->k,1);
	input->U = alloc_matrix(input->n,input->h);
	input->V = alloc_matrix(input->k,input->h);
	input->dsT = alloc_matrix(input->k, input->n);
	centrareDataset(input, input->ds, input->n, input->k);
	estraiVettoreScore(input);
	//NIPALS
	for(int i = 0; i < input->h; i++ ){
    	traspostaDataSet(input);
		ricalcola:
				   calcoloVettoreLoad(input);
		           normalizzaVettoreLoad(input);
		           t = prodottoScalare(input->vetScore, input->n);
		           aggiornaVettoreScore(input);
		           t1 = prodottoScalare(input->vetScore, input->n);
			           if ((abs(t1 - t)) >= t1*teta ){
		        	   goto ricalcola;
		           }
		           aggiornaMatriceScore(input, iColonna);
		           aggiornaMatriceLoad(input, iColonna);
		           iColonna++;
		           aggiornaDataSet(input, matriceTemp);
	}//NIPALS
    dealloc_matrix(input->vetLoad);
    dealloc_matrix(input->vetScore);
    dealloc_matrix(input->dsT);
	dealloc_matrix(matriceTemp);
}//PCA

//--------------------KDTree-------------------//

int calcolaMediano(float* ds, int n, int k, int colonna, float* vettP[], int vettI[]){
	float* vettC = alloc_matrix(n, 1);
	int med;
	int z = 0;
	//estrai colonna
	for(int j = colonna; j < n*k; j=j+k ){
		vettC[z] = ds[j];
		vettP[z] = &ds[j];
		z++;
	}
	//ordiniamo colonna
	int j;
	float val;
	float* val1;
	int val2;
	for(int i=1;i<n;i++){
		val = vettC[i];
		val1 = vettP[i];
		val2 = vettI[i];
		j=i-1;
		for(j;j>=0 && vettC[j]>val; j--){
			vettC[j+1] = vettC[j];
			vettP[j+1] = vettP[j];
			vettI[j+1] = vettI[j];
		}
		vettC[j+1] = val;
		vettP[j+1] = val1;
		vettI[j+1] = val2;
	}
	if( n % 2 == 0 )
		 med = ( n / 2 );
	else
		 med = (n - 1) / 2;
	//printf("STAMPA MEDIANO \n");
	return med;
}//calcolaMediano

float* costruisciDS1(float* vettP[], int med, int n, int k, int c, int vettI[], int* vettI1){
	//printf("STAMPA DS1 \n");
	float* ds1 = alloc_matrix(n, k);
	int f=0;
	float* riga;
	for(int i = 0; i < med; i++){
		riga = vettP[i] - c;
		vettI1[i] = vettI[i];
		for(int z = 0; z < k; z++){
			ds1[z+f] = *(riga+z);
		}
		f=f+k;
	}
	return ds1;
}//costruisciDS1

float* costruisciDS2(float* vettP[], int med, int n, int k, int c, int vettI[], int* vettI2){
	float* ds2 = alloc_matrix(n, k);
	int f=0;
	int cont = 1;
	float* riga;
	for(int i = 0; i < n; i++){
		vettI2[i] = vettI[med+cont];
		riga = vettP[med+cont] - c;
		for(int j = 0; j < k; j++){
			ds2[j+f] = *(riga);
			riga = riga + 1;
		}
		cont++;
		f = f + k;
	}
	return ds2;
}//costruisciDS2

int stampaAlbero(Albero A, int k){
	if(A == NULL) {return 0;}
	if(A->id_vicino != -1){
		printf("stampa albero \n");
		for(int i = 0; i < k; i++){
			printf(" %f ",A->vetRadice[i]);
		}
		printf(" id vicino: %d ",A->id_vicino);
		printf("\n");
	}
	stampaAlbero(A->sinistro, k);
	stampaAlbero(A->destro, k);
	return 0;
}//stampaAlbero


void memorizzaRadice(float* vettP[], int med, int n, int k, int c, Albero A){
	A->vetRadice = alloc_matrix(1,k);
	float* riga = vettP[med] - c;
	for(int i = 0; i < k; i++){
		A->vetRadice[i] = *riga;
		riga = riga + 1;
	}
}//memorizzaRadice

//regione indotta, massimo e minimo per colonna
void coppieMaxMin(float* ds, int n, int k, Albero A){
	float max = ds[0], min = ds[0];
	int cont=0;
	A->vetCoppie = alloc_matrix(k,2);
	for(int j = 0; j < k; j++){
		max = ds[j];
		min = ds[j];
		for(int i = j; i < n*k; i=i+k){
			if( ds[i] > max) max = ds[i];
			if( ds[i] < min) min = ds[i];
		}

		A->vetCoppie[j+cont] = min;
		cont++;
		A->vetCoppie[j+cont] = max;
	}
}//coppieMaxMin

Albero buildKdtree(float* ds, int n, int k, int livello, float* vettP[], int vettI[]) {
	if (n == 0) {
		return NULL;
	}
	Albero A = malloc(sizeof(struct NodoAlbero));
	int med;
	float *ds1;
	float *ds2;
	int n1=0, n2=0;
	int c;
	c = livello % k;
	med = calcolaMediano(ds, n, k, c, vettP, vettI);
	A->id_vicino = vettI[med];
	coppieMaxMin(ds, n, k, A);
	if ( n % 2 == 0)
		{n1 = n / 2; n2 = (n / 2) - 1; }
	else
		{ n1 = (n - 1) / 2; n2 = (n - 1) / 2;}
	int *vettI1 = (int*)malloc(sizeof(int)*n1);
	int *vettI2 = (int*)malloc(sizeof(int)*n2);
	ds1 = costruisciDS1(vettP, med, n1, k, c, vettI, vettI1);
	ds2 = costruisciDS2(vettP, med, n2, k, c, vettI, vettI2);
	livello++;
	memorizzaRadice(vettP, med, n, k, c, A);
	A->sinistro = buildKdtree(ds1, n1, k, livello, vettP, vettI1);
	A->destro = buildKdtree(ds2, n2, k, livello, vettP, vettI2);
	return A;
}//kdtree

Albero kdtree(params* input){
    float* vettP[input->n];
  	int* vettI = (int*)malloc(sizeof(int)*input->n);
  	int livello = 0;
    Albero A1 = malloc(sizeof(struct NodoAlbero));
  	for(int y = 0; y < input->n; y++)
  		vettI[y] = y;
  	if(input->h > 0){
  		dealloc_matrix(input->ds); //DA PROVARE
  		A1 = buildKdtree(input->U, input->n, input->h, livello, vettP, vettI);
  		dealloc_matrix(input->U);
  	}else{
  		A1 = buildKdtree(input->ds, input->n, input->k, livello, vettP, vettI);
  		dealloc_matrix(input->ds);
  	}
  	return A1;
}//kdtree

//--------------------RangeQuery-----------------//

float distanzaEuclidea1(float* p, float* query, int k){
	float d = 0;
	float somma = 0;
	float x = 0;
	for(int i = 0; i < k; i++){
		x = p[i] - query[i];
		somma = x*x + somma;
	}
	d = sqrt(somma);
	return d;
}//distanzaEuclidea


float calcoloDistanza(float* query, float* coppie, int k){
	int cont = 0;
	float d;
	float* p = alloc_matrix(k, 1);
	for( int i = 0; i < k; i++){
		if (query[i] <= coppie[cont])
			p[i] = coppie[cont];
		else if (query[i] >= coppie[cont+1])
			p[i] = coppie[cont+1];
		else p[i] = query[i];
		cont = cont + 2;
	}
	//distanzaEuclidea(p, query, k, d);
	d = distanzaEuclidea1(p, query, k);
	dealloc_matrix(p);
	return d;
}//calcoloDistanza

int stampaLista(Lista L, int nQA){
	if(L == NULL || nQA == 0) return 0;
	if(L->id_vicino != -1){	printf("%d , %d \n",L->id_query,L->id_vicino);}
	stampaLista(L->succ, nQA-1);
	return 0;
}

void add(Albero A, int p){
    L = malloc(sizeof(struct Nodo));
    L->id_query = p;
    L->id_vicino = A->id_vicino;
    L->succ = inizio;
    inizio = L;
}

int calcolaQuery(Albero A, float* query, float r, int k, int p, params* input){
	if( calcoloDistanza(query, A->vetCoppie, k) > r ) {
		return 0;
	}

	float d=0;
	float* punto = alloc_matrix(k,1);
	for(int j = 0; j < k; j++)
		punto[j] = A->vetRadice[j];

	//distanzaEuclidea(p, query, k, d);
	d=distanzaEuclidea1(punto, query, k);
	if( d <= r){
		input->nQA++;
		add(A, p);
	}
	dealloc_matrix(punto);
	if( A->sinistro != NULL ){
		calcolaQuery(A->sinistro, query, r, k, p, input);
	}
	if( A->destro != NULL ){
		calcolaQuery(A->destro, query, r, k, p, input);
	}
	return 0;
}

Lista* range_query(Albero A, float* qs, float r, int k, int nq, params* input){
	float* query = alloc_matrix(k,1);
	Lista* B2 = (Lista*)malloc(sizeof(Lista)*nq);
	int dim=0;
	for(int j = 0; j < nq; j++){
		for(int i = 0; i < k; i++)
			query[i] = qs[i+dim];
		dim = dim+k;
		inizio = NULL;
	    L = malloc(sizeof(struct Nodo));
		L->id_query = j;
	    L->id_vicino = -1;
	    L->succ = NULL;
		calcolaQuery(A, query, r, k, j, input);
		B2[j] = L;
	}
	dealloc_matrix(query);
	return B2;
}//rangeQuery

void prodottoMatrici(params* input){
	input->pcaQS = alloc_matrix(input->nq, input->h);
	int rigaqs = 0, rigav = 0, rigaqsm = 0;
	float mul=0;
	for(int i = 0; i < input->nq; i++){
		for(int j = 0; j < input->h; j++){
			for(int z = 0; z < input->k; z++){
				mul = (input->dsMedia[rigaqsm+z] * input->V[rigav+j]) + mul;
				rigav = rigav + input->h;
			}
			input->pcaQS[rigaqs+j] = mul;
			rigav = 0;
			mul = 0;
		}
		rigaqsm = rigaqsm + input->k;
		rigaqs = rigaqs + input->h;
	}
}

int salvaLista(Lista L, FILE *f, int nQA){
	if(L == NULL || nQA == 0){ return 0; }
	if(f && L->id_vicino != -1){
		fprintf(f,"%d %d \n",L->id_query, L->id_vicino);
	}
	salvaLista(L->succ, f, nQA-1);
}//salvaLista

int save_dataL(Lista L, FILE *f, int nQA){
	if(L == NULL || nQA == 0){ return 0; }
	if(f && L->id_vicino != -1){
        fwrite(&L->id_query, sizeof(int), 1, f);
        fwrite(&L->id_vicino, sizeof(int), 1, f);
	}
	save_dataL(L->succ, f, nQA-1);
}//salvaLista

int main(int argc, char** argv){

    char fname[256];
    int i, j, k;
    clock_t t;
    char* dsname;
    float time;
    //
    // Imposta i valori di default dei parametri
    //
    params* input = malloc(sizeof(params));

    input->filename = NULL;
    input->h = 0;
    input->kdtree = NULL;
    input->r = -1;
    input->silent = 0;
    input->display = 0;
    input->QA = NULL;
    input->nQA = 0;

    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //
    
    if (argc <= 1 && !input->silent) {
        printf("Usage: %s <data_name> [-pca <h>] [-kdtree [-rq <r>]]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\t-d: display query results\n");
        printf("\t-s: \n");
        printf("\t-pca <h>: h-component PCA enabled\n");
        printf("\t-kdtree: kdtree building enabled\n");
        printf("\t-rq <r>: range query search with radius r enabled\n");
        printf("\n");
        //exit(0);
    }
    
    //
    // Legge i valori dei parametri da riga comandi
    //

    int par = 1;
    while (par < argc) {
        if (par == 1) {
            input->filename = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-pca") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing h value!\n");
                exit(1);
            }
            input->h = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-kdtree") == 0) {
            input->kdtree_enabled = 1;
            par++;
            if (par < argc && strcmp(argv[par],"-rq") == 0) {
                par++;
                if (par >= argc) {
                    printf("Missing radius value!\n");
                    exit(1);
                }
                input->r = atof(argv[par]);

                if(input->r < 0){
                    printf("Range query radius must be non-negative!\n");
                    exit(1);
                }
                par++;
            }
        } else{
            printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
            par++;
        }
    }
    
    //
    // Legge i dati e verifica la correttezza dei parametri
    //


    if(input->filename == NULL || strlen(input->filename) == 0){
        printf("Missing input file name!\n");
        exit(1);
    }
    sprintf(fname, "%s.ds", input->filename);
    dsname = basename(strdup(input->filename));
    input->ds = load_data(fname, &input->n, &input->k);

	if(input->h < 0){
        printf("Invalid value of PCA parameter h!\n");
        exit(1);
    }
    if(input->h > input->k){
        printf("Value of PCA parameter h exceeds data set dimensions!\n");
        exit(1);
    }
    if(input->r >= 0){
        sprintf(fname, "%s.qs", input->filename);
        input->qs = load_data(fname, &input->nq, &k);
        if(input->k != k){
            printf("Data set dimensions and query set dimensions are not compatible!\n");
            exit(1);
        }
    }
    
    //
    // Visualizza il valore dei parametri
    //
    
    if(!input->silent){
        printf("Input file name: '%s'\n", input->filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [k]: %d\n", input->k);
        if(input->h > 0){
            printf("PCA search enabled\n");
            printf("Number of principal components [h]: %i\n",input->h);
        }else{
            printf("PCA search disabled\n");
        }
        if(input->kdtree_enabled){
            printf("Kdtree building enabled\n");
            if(input->r >= 0){
                printf("Range query search enabled\n");
                printf("Range query search radius [r]: %f\n",input->r);
            }else{
                printf("Range query search disabled\n");
            }
        }else{
            printf("Kdtree building disabled\n");
        }
    }
    //
    // Calcolo PCA
    //
    if(input->h > 0){
        t = clock();
        pca(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
        sprintf(fname, "%s.U", dsname);
        save_data(fname, input->U, input->n, input->h);
        sprintf(fname, "%s.V", dsname);
        save_data(fname, input->V, input->k, input->h);
    }else
        time = -1;
    if (!input->silent)
    	printf("\nPCA time = %.3f sec\n", time);
    else
    	printf("%.3f\n", time);

    //
    // Costruzione K-d-Tree
    //
    Albero A1 = malloc(sizeof(struct NodoAlbero));
    if(input->kdtree_enabled){
        t = clock();
        A1 = kdtree(input);
        t = clock() - t;
        time = ((float)t)/CLOCKS_PER_SEC;
    }else
        time = -1;
    if (!input->silent)
        printf("\nIndexing time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    //stampaAlbero(A1, input->h);
    //
    // Range query search
    //
	Lista* B3 = (Lista*)malloc(sizeof(Lista)*input->nq);

	if(input->r >= 0){
		t = clock();
		if(input->h > 0){
			centrareDataset(input, input->qs, input->nq, input->k);
			prodottoMatrici(input);
			dealloc_matrix(input->qs);
			B3 = range_query(A1, input->pcaQS, input->r, input->h, input->nq, input);
			dealloc_matrix(input->pcaQS);
		}else{
			B3 = range_query(A1, input->qs, input->r, input->k, input->nq, input);
			dealloc_matrix(input->qs);
		}
		//printf("stampa albero \n");
		t = clock() - t;
		time = ((float)t)/CLOCKS_PER_SEC;
	}else
		time = -1;
	if (!input->silent)
		printf("\nQuerying time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

    //
    // Salva il risultato delle query
    // da modificare se si modifica il formato delle matrici di output
    //moltiplicaMatrici
    
    if(input->r >= 0){
        if(!input->silent && input->display) {
            //NB: il codice non assume che QA sia ordinata per query, in caso lo sia ottimizzare il codice
            printf("\nQuery Answer:\n");
            for(i = 0; i < input->nq; i++){
             //   printf("query %d: [ ", i);
                for(j = 0; j < input->nQA; j++)
                    if(input->QA[j*2] == i)
                        printf("%d ", input->QA[j*2+1]);
                printf("]\n");
            }
        }

         sprintf(fname, "%s.qa", dsname);
        //save_data(fname, input->QA, input->nQA, 2);
/*
        for(int i = 0; i < input->nq; i++){
    		printf("query %d \n",i);
    		stampaLista(B3[i], input->nQA);
    		printf("\n");
    	}
*/
        int c=2;
        FILE *f;
        f = fopen(fname,"wb");
	    fwrite(&c, sizeof(int), 1, f);
		fwrite(&input->nQA, sizeof(int), 1, f);
    	for(int i = 0; i < input->nq; i++){
           	save_dataL(B3[i], f, input->nQA);
    	}
    	fclose(f);

        FILE *f1;
        f1 = fopen("lista32bit.txt","w");
    	for(int i = 0; i < input->nq; i++){
    		salvaLista(B3[i], f1, input->nQA);
    	}
    	fclose(f1);

    }
    if (!input->silent)
        printf("\nDone.\n");
    return 0;
}

