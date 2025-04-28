#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h> //struct timeval t1, t2; gettimeofday(&t1, NULL); gettimeofday(&t2, NULL); 

typedef int indice;
typedef float proba;

int M, C, L; // Nb elements non nuls, colonnes, lignes
proba alpha = 0.85; // Proba du surfeur aléatoire
proba invAlpha = 1 - 0.85;

struct elem {
    indice i, j; // Les indexes
    proba val; // La valeur
};

struct elem *P; // Matrice stochastique
int *f; // Vecteur ligne de taille N tel que f[i] = 1 si la ligne i de P ne contient que des zéros et sinon f[i] = 0
proba *w;
proba *x;
proba *y;
proba *z;

void initvec(proba *z, float val) {
    indice i; for (i = 0; i<C; i++) { z[i] = val;}
}

void initvecStoc(proba *z, float val) {
    indice i; for (i = 0; i<C; i++) { z[i] = val/(float) C;}
}

void calculeF(int *f, struct elem *P) {
    if (f == NULL) { exit(50); }
    for (indice i = 0; i<C; i++) {
        f[i] = 1;
    }
    int ligne;
    for (indice k = 0; k < M; k++) {
        ligne = P[k].i;
        if (f[ligne]){
            f[ligne] = 0;
        }
    }
}

void affichevec(proba *z) {
    indice i; for (i = 0; i<C; i++) {printf("%f ", z[i]);}
    printf("\n");
}

void affichevecInt(int *f) {
    printf("Affichage de vecteur de Int\n");
    indice i; for (i = 0; i<C; i++) {printf("%d ", f[i]);}
    printf("\n");
}

void affichemat(struct elem *P) {
    indice k; for (k = 0; k < M; k++) {{printf("%d %d %f\n", P[k].i, P[k].j, P[k].val);}
    }
}

struct elem* lectureMatrice(char *nom_fic) {
    indice i, j, k;
    struct elem *P;
    double v;
    FILE *F;
    F = fopen(nom_fic, "r");
    fscanf(F, "%d", &(L));
    fscanf(F, "%d", &(C));
    fscanf(F, "%d", &(M));

    P = malloc(M * sizeof(struct elem));
    if (P == NULL) { exit(10); };

    for (k = 0; k < M; k++) {
        struct elem e;
        fscanf(F, "%d", &(i));
        fscanf(F, "%d", &(j));
        fscanf(F, "%le", &(v));
        e.i = i - 1;
        e.j = j - 1;
        e.val = v;
        P[k] = e;
    }
    fclose(F);
    return P;
}

proba* diffVect(proba *x, proba *y) {
    indice k;
    proba *z = malloc(L * sizeof(float));
    for (k = 0; k < L; k++) {
        z[k] = x[k] - y[k];
    }
    return z;
}

float normeVect(proba *x, proba *y) {
    indice k;
    proba* delta = diffVect(x, y);
    float norme = 0;
    for (k = 0; k < L; k++) {
        norme += fabs(delta[k]);
    }
    free(delta);
    return norme;
}

// Addition entre deux vecteurs de proba (ajoute dans y)
void addVect(proba *x, proba *y) {
    indice i; for (i = 0; i < C; i++) {y[i] += x[i];}
}

// Multiplie un vecteur de proba par un float
proba* multFloatProba(float a, proba *x) {
    indice i;
    proba *res = malloc(C * sizeof(proba)); if (res == NULL) {exit(40);}
    for (i = 0; i < C; i++) {
        res[i] = a * x[i];
    }
    return res;
}

// Multiplie un vecteur de proba (horizontal) avec un vecteur de int (vertical)
float multVectPVectI(proba *x, int *f) {
    indice i;
    float res = 0.;
    for (i = 0; i<C; i++) {
        res += x[i] * f[i];
    }
    return res;
}

// Multiplie à gauche une matrice par un vecteur
void multVecMat(proba *x, struct elem *P, proba *y) {
    indice i, j, k;
    for (k = 0; k < M; k++) {
        i = P[k].i;
        j = P[k].j;
        float val = P[k].val;
        y[j] += x[i] * val;
    }
}

void metZero(proba *y) {
    indice i;
    for (i = 0; i < L; i++) {
        y[i] = 0.0;
    }
}

void recopie(proba *x, proba *y) {
    indice i;
    for (i = 0; i < L; i++) {
        x[i] = y[i];
    }
}

void iterer(proba *x, proba *y, proba *w) {
    proba *aX = multFloatProba(alpha, x);
    multVecMat(aX, P, y);

    proba gauche = invAlpha / (float) C;
    proba droite = (alpha / (float) C) * multVectPVectI(x, f);

    initvec(w, gauche + droite);
    addVect(w, y);
}

void puissances(struct elem *P, proba *x) {
    
    float delta = 1.0;
    float epsilon = 10e-6;
    int cnt = 0;
    f = malloc(C * sizeof(int));
    w = malloc(L * sizeof(proba));
    y = malloc(L * sizeof(proba));
   
    calculeF(f, P);
    // x : PI(n)
    // y : PI(n+1)
    while (delta > epsilon) {
        cnt ++;
        metZero(y);
        metZero(w);
        iterer(x, y, w);
        printf("pi %d :\n", cnt);
        affichevec(y);
        delta = normeVect(x, y);
        recopie(x, y);
    }

    free(f);
    free(y);
}

int main(int argc, char **argv) {
    
    char nom[] = "matriceCreuse.txt";
    P = lectureMatrice(nom);
    x = malloc(L * sizeof(proba));

    if (x == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    initvecStoc(x, 1); // x0 = (1/N)e

    affichevec(x);

    puissances(P, x);

    free(x);
    free(P);
    return 0;
}