#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h> //struct timeval t1, t2; gettimeofday(&t1, NULL); gettimeofday(&t2, NULL); 

typedef int indice;
typedef float proba;

    int M; // Nb elements non nuls
int C; // Nb colonnes
proba alpha = 0.85; // Proba du surfeur aléatoire
proba invAlpha = 1 - 0.85;

struct elem {
    indice i, j; // Les indexes
    proba val; // La valeur
};

struct elem *P; // Matrice stochastique
struct timeval t1, t2;

int *f; // Vecteur ligne de taille N tel que f[i] = 1 si la ligne i de P ne contient que des zéros et sinon f[i] = 0
proba *w;
proba *x;
proba *y;
proba *z;

void initvec(proba *z, float val) {
    indice i; for (i = 0; i<C; i++) { z[i] = val;}
}

void initvecStoc(proba *z, float val) {
    initvec(z, (double) C);
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
    indice i; for (i = 0; i<C; i++) {printf("%d ", f[i]);}
    printf("\n");
}

void affichemat(struct elem *P) {
    indice k; for (k = 0; k < M; k++) {{printf("%d %d %f\n", P[k].i, P[k].j, P[k].val);}
    }
}

struct elem* lectureMatrice(char *nom_fic) {
    printf("Lecture Matrice\n");
    indice i, j;
    int nb;
    struct elem *P;
    double v;
    FILE *F;
    F = fopen(nom_fic, "r");
    if (F == NULL) { exit(100); }
    printf("Lecture header\n");
    fscanf(F, "%d", &(C));
    fscanf(F, "%d", &(M));
    printf("Lecture header done: %d, %d\n", C, M);
    P = malloc(M * sizeof(struct elem));
    if (P == NULL) { exit(10); };

    indice courant = 0;
    for (int k = 0; k < C; k++) {
        fscanf(F, "%d", &(i));
        fscanf(F,"%d",  &(nb)); 
        printf("%d\n", i);
        printf("nombre attendu : %d\n", nb);
        for (int z = 0; z < nb; z++) {
            struct elem e;

            fscanf(F, "%d", &(j));
            fscanf(F, "%le", &(v));
            printf("elem: %d, %f\n", j, v);
            e.i = i - 1;
            e.j = j - 1;
            e.val = v;
            P[courant] = e;
            courant++;
            if (courant > M) {
                exit(202);
            }
            printf("\n");
        }
    }
    printf("Nombre d'elements attendu: %d, nombre réel: %d\n", M, courant);
    fclose(F);
    return P;
}

proba* diffVect(proba *x, proba *y) {
    indice k;
    proba *z = malloc(C * sizeof(float));

    for (k = 0; k < C; k++) {
        z[k] = x[k] - y[k];
    }
    return z;
}

float normeVect(proba *x, proba *y) {
    indice k;
    proba* delta = diffVect(x, y);
    float norme = 0;
    for (k = 0; k < C; k++) {
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
    for (i = 0; i < C; i++) {
        y[i] = 0.0;
    }
}

void recopie(proba *x, proba *y) {
    indice i;
    for (i = 0; i < C; i++) {
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

proba* puissances(char nom[], proba *x) {
    P = lectureMatrice(nom);
    if (!x) {

        x = malloc(C * sizeof(proba));
        initvec(x, 1/(double) C);  // x0 = (1/N)e

    } else {  // adapter en fonction du nombre de nouveaux elements
        
        z = malloc(C * sizeof(proba));

        for (indice i = 0; i < C; i++) {
            if (x[i]) {
                z[i] = x[i];
            } else {
                z[i] = 0.0;
            }
        }
        printf("Vecteur fourni: \n");
        x = z;
    }
    affichevec(x);

    float delta = 1.0;
    float epsilon = 10e-6;
    int cnt = 0;

    f = malloc(C * sizeof(int));
    w = malloc(C * sizeof(proba));
    y = malloc(C * sizeof(proba));
   
    if (f == NULL || w == NULL || y == NULL) {
        printf("Erreur d'allocation");
        exit(40);
    }

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
    return x;
}

int main(int argc, char *argv[]) {

    //FILE *gnuplot = popen("gnuplot -persistent", "w");
    //if (gnuplot) {
    //    
    //    fprintf(gnuplot, "set title 'Temps de convergence en fonction de alpha'\n");
    //    fprintf(gnuplot, "set xlabel 'alpha'\n");
    //    fprintf(gnuplot, "set ylabel 'Temps (microsecondes)'\n");
    //    fprintf(gnuplot, "plot '-' with linespoints title 'PageRank'\n");
//
    //    double pas = 10.;
    //    for (int i = 0; i < pas; i++) {
    //        alpha = 0 + (i * (1 / pas));
    //        invAlpha = 1 - alpha;
    //        long moyenne_temps = 0;
//
    //        double moyenne = 100.;
    //        for (int index = 0; index < moyenne; index++) {
    //            gettimeofday(&t1, NULL);
    //            puissances("Matrices/StanfordBerkeley.txt", 0);
    //            gettimeofday(&t2, NULL);
    //            moyenne_temps += t2.tv_usec - t1.tv_usec;
    //        }
    //        if (moyenne_temps > 0) {
    //            fprintf(gnuplot, "%f %f\n", alpha, moyenne_temps / moyenne);
    //        }
    //    }
    //}

//    printf("\nitération 2\n");
//    x = puissances("matriceCreuseV2.txt", x);
//
//    printf("\nitération 3\n");
//    x = puissances("matriceCreuseV3.txt", x);

    puissances("Matrices/StanfordBerkeley.txt", 0);
    free(x);
    free(P);
    return 0;
}