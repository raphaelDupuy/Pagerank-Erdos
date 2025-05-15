#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include "calcul.c"
#include "affichage.c"

struct resultat_puissances {  // Sortie de la fonction Puissances
    struct vecteur *V;        // Vecteur de proba stationnaire
    int n;                    // Nombre d'étapes avant la convergence
};

struct timeval t1, t2;

int *f; // Vecteur ligne de taille N tel que f[i] = 1 si la ligne i de P ne contient que des zéros et sinon f[i] = 0

void initvec(struct vecteur *Vec, float val, int taille) {
    indice i; for (i = 0; i < taille; i++) { Vec->v[i] = val; }
}

void calculeF(int *f, struct matrice *Mat) {
    for (indice i = 0; i<Mat->C; i++) {
        f[i] = 1;
    }
    int ligne;
    for (indice k = 0; k < Mat->M; k++) {
        ligne = Mat->P[k].i;
        if (f[ligne]){
            f[ligne] = 0;
        }
    }
}

struct matrice *lectureMatrice(const char *nom_fic) {
    indice i, j, nb;
    int C, M;
    double v;
    FILE *F = fopen(nom_fic, "r");
    if (!F) {
        perror("ouverture fichier");
        exit(100);
    }

    if (fscanf(F, "%d %d", &C, &M) != 2) {
        fprintf(stderr, "Format de header invalide\n");
        exit(101);
    }
    printf("Lecture header done: C = %d, M = %d\n", C, M);

    struct elem *P = malloc(M * sizeof *P);
    if (!P) {
        perror("malloc P");
        exit(102);
    }

    //Boucle de lecture jusqu'à avoir M éléments ou EOF
    int courant = 0;
    while (courant < M && fscanf(F, "%d %d", &i, &nb) == 2) {
        if (i < 1 || i > C) {
            fprintf(stderr, "Index de ligne hors bornes: %d\n", i);
            exit(103);
        }

        // Lecture des 'nb' paires (j, v)
        for (int z = 0; z < nb; z++) {
            if (fscanf(F, "%d %lf", &j, &v) != 2) {
                fprintf(stderr, "Erreur de format à la ligne %d, voisin %d\n", i, z);
                exit(104);
            }
            if (courant >= M) {
                fprintf(stderr, "Trop d'éléments lus (>%d)\n", M);
                exit(105);
            }
            P[courant].i   = i - 1;
            P[courant].j   = j - 1;
            P[courant].val = (proba)v;
            courant++;
        }
    }
    fclose(F);

    // Vérification
    if (courant != M) {
        fprintf(stderr,
          "Nombre d'éléments lus (%d) != M attendu (%d)\n", courant, M);
        exit(106);
    }
    printf("Lecture matrice terminée, %d éléments lus\n", courant);

    struct matrice *Mat = malloc(sizeof *Mat);
    if (!Mat) {
        perror("malloc Mat");
        exit(107);
    }
    Mat->C = C;
    Mat->M = M;
    Mat->P = P;
    return Mat;
}

void recopie(struct vecteur *X, struct vecteur *Y) {
    indice i;
    int C = X->C;
    if (C != Y->C) { exit(320); }
    for (i = 0; i < C; i++) {
        X->v[i] = Y->v[i];
    }
}

// Mat0 -> Matrice à partir de laquelle on génère une nouvelle
// n    -> Proportion de noeuds à rajouter (%)
// p    -> Proba d'avoir un arc entre un nouveau noeud et un ancien
struct matrice *genereErdosStochastique(const struct matrice *Mat0, double n, proba p){
    indice C0 = Mat0->C;
    indice C  = C0 + (C0 * (n/100));
    int M0 = Mat0->M;

    // on recopie M0 + n * C arcs au plus
    int maxM = M0 + n * C;
    struct elem *tmp = malloc(maxM * sizeof *tmp);
    if (!tmp) {
        perror("malloc tmp");
        exit(EXIT_FAILURE);
    }

    // Copie de la matrice initiale
    memcpy(tmp, Mat0->P, M0 * sizeof *tmp);
    int k = M0;

    // Buffer temporaire pour stocker les voisins d’un nouveau sommet
    indice *voisins = malloc(C * sizeof *voisins);
    if (!voisins) {
        perror("malloc voisins");
        exit(EXIT_FAILURE);
    }

    for (indice di = 0; di < n; di++) {
        indice i = C0 + di;
        int deg = 0;

        // Collecte des voisins selon probabilité p
        for (indice j = 0; j < C0; j++) {
            if (j == i) continue;
            if ((double)rand() / RAND_MAX < p) {
                voisins[deg++] = j;
            }
        }

        if (deg == 0) {
            for (indice j = 0; j < C0; j++) {
                if (j == i) continue;
                voisins[deg++] = j;
            }
        }

        // Ecriture immédiate de ces deg arcs avec valeur 1/deg
        proba v = 1.0f / deg;
        for (int t = 0; t < deg; t++) {
            tmp[k].i   = i;
            tmp[k].j   = voisins[t];
            tmp[k].val = v;
            k++;
        }
    }

    free(voisins);

    // Construction de la matrice finale de taille k
    struct matrice *Mat = malloc(sizeof *Mat);
    if (!Mat) {
        perror("malloc Mat");
        exit(EXIT_FAILURE);
    }
    Mat->C = C;
    Mat->M = k;
    Mat->P = malloc(k * sizeof *Mat->P);
    if (!Mat->P) {
        perror("malloc Mat->P");
        exit(EXIT_FAILURE);
    }
    memcpy(Mat->P, tmp, k * sizeof *tmp);
    free(tmp);

    printf("généré Erdos-Rényi stochastique : C=%d, M=%d\n", Mat->C, Mat->M);
    return Mat;
}

void iterer(struct vecteur *X, struct vecteur *Y, struct vecteur *W, struct matrice *Mat, float alpha) {
    struct vecteur *aX = multFloatProba(alpha, X);

    multVecMat(aX, Mat, Y);

    proba gauche = (1 - alpha) / (float) Mat->C;
    proba droite = (alpha / (float) Mat->C) * multVectPVectI(X, f);

    initvec(W, gauche + droite, Mat->C);
    addVect(W, Y);
}

struct resultat_puissances *puissances(struct matrice *Mat, struct vecteur *Pi, float alpha) {
    int C = Mat->C;
    if (!Pi) {
        Pi = alloueVecteur(C);
        initvec(Pi, 1/(double) C, C);  // x0 = (1/N)e
    } else {  // adapter en fonction du nombre de nouveaux elements
        if (C < Pi->C) { printf("%d > %d\n", Pi->C, C); exit(64); }
        struct vecteur *Z = alloueVecteur(C);

        for (indice i = 0; i < C; i++) {
            if (Pi->v[i]) {
                Z->v[i] = Pi->v[i];
            } else {
                Z->v[i] = 0.0;
            }
        }

        free(Pi->v);
        free(Pi);
        Pi = Z;
    }

    float delta = 1.0;
    float epsilon = 10e-6;
    int cnt = 0;

    f = malloc(Mat->C * sizeof(int));
    struct vecteur *y = alloueVecteur(Mat->C);
    struct vecteur *w = alloueVecteur(Mat->C);
    if (f == NULL || y == NULL || w == NULL) {
        printf("Erreur d'allocation");
        exit(40);
    }

    calculeF(f, Mat);
    // x : PI(n)
    // y : PI(n+1)
    while (delta > epsilon) {
        cnt ++;
        metZero(y);
        metZero(w);
        iterer(Pi, y, w, Mat, alpha);
        delta = normeVect(Pi, y);
        recopie(Pi, y);
    }

    free(f);
    free(y);
    struct resultat_puissances *res = malloc(sizeof *res);
    if (!res) {perror("malloc res"); exit(27);}
    res->V = Pi;
    res->n = cnt;
    return res;
}

struct vecteur **plot(struct matrice *Mat, struct vecteur **Pi_tab, char nom[]) {

    FILE *gnuplot = initialiseGNU(nom);
    int pas = 20;
    struct vecteur **vecteurs = malloc(pas * sizeof *vecteurs);
    if (vecteurs == NULL) { exit(501); }
    if (gnuplot) {
        printf("Starting\n");
        for (int i = 0; i < pas - 1; i++) {
            float alpha = 0 + (i * (1 / (double) pas));
            long long moyenne_iter = 0;

            struct vecteur *init = (Pi_tab != NULL) ? Pi_tab[i] : NULL;
            struct resultat_puissances *res = puissances(Mat, init, alpha);
            moyenne_iter += res->n;
            vecteurs[i] = res->V;

            fprintf(gnuplot, "%f %lld\n", alpha, moyenne_iter);
            printf("%f, %lld\n", alpha, moyenne_iter);
        }
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
    }
    return vecteurs;
}

int main(int argc, char *argv[]) {

    char nom[] = "Matrices/StanfordBerkeley.txt";

    struct matrice *Mat = lectureMatrice(nom);

    struct vecteur **vecteurs = plot(Mat, 0, "plot_depart");
    struct matrice *erdos = genereErdosStochastique(Mat, 2, 1);
    plot(erdos, vecteurs, "tst");   
    //for (int nb = 1; nb < 1000; nb += 50) {
    //    printf("%d\n", nb);
    //    struct matrice *erdos = genereErdosStochastique(Mat, nb, 0.2);
//
    //    char nom_png_N[64];
    //    snprintf(nom_png_N, sizeof nom_png_N, "plot_erdos_N_%03d", nb);
//
    //    char nom_png_Pi[64];
    //    snprintf(nom_png_Pi, sizeof nom_png_Pi, "plot_erdos_Pi_%03d", nb);
    //    plot(erdos, 0, nom_png_N);
    //    plot(erdos, vecteurs, nom_png_Pi);
    //}
    return 0;
}