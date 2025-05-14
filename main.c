#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/time.h> //struct timeval t1, t2; gettimeofday(&t1, NULL); gettimeofday(&t2, NULL); 

typedef int indice;
typedef float proba;

struct elem {
    indice i, j; // Les indexes
    proba val;   // La valeur
};

struct matrice {    // Matrice carrée
    struct elem *P; // Ensemble des (i, j, val)
    indice C;       // Nombre de colonnes (et lignes, matrice carrée)
    int M;          // Nombre d'elements non nuls (de *P)
};

struct vecteur {
    proba *v;
    indice C;
};

struct resultat_puissances {  // Sortie de la fonction Puissances
    struct vecteur *V;        // Vecteur de proba stationnaire
    int n;                    // Nombre d'étapes avant la convergence
};

struct timeval t1, t2;

int *f; // Vecteur ligne de taille N tel que f[i] = 1 si la ligne i de P ne contient que des zéros et sinon f[i] = 0

struct vecteur *alloueVecteur(int taille) {
    struct vecteur *Z = malloc(sizeof(struct vecteur)); if(!Z) { exit(420); }
    Z->C = taille;; if (!Z->C) { exit(422); }
    Z->v = malloc(taille * sizeof(proba)); if (!Z->v) { exit(421); }
    return Z;
}

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

void affichevec(struct vecteur *Vec) {
    printf("Affichage vecteur\n");
    if (!Vec || !Vec->v) { exit(6); }
    indice i; for (i = 0; i < Vec->C; i++) {printf("%f ", Vec->v[i]);}
    printf("\n");
}

void affichemat(struct matrice *Mat) {
    indice k; for (k = 0; k < Mat->M; k++) {{printf("%d %d %f\n", Mat->P[k].i, Mat->P[k].j, Mat->P[k].val);}
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

proba *diffVect(struct vecteur *X, struct vecteur *Y) {
    indice k;
    int C = X->C;
    if (C != Y->C) { exit(320); }
    proba *Z = malloc(C * sizeof(proba));
    for (k = 0; k < C; k++) {
        Z[k] = X->v[k] - Y->v[k];
    }
    return Z;
}

float normeVect(struct vecteur *X, struct vecteur *Y) {
    indice k;
    int C = X->C;
    if (C != Y->C) { exit(320); }

    proba* delta = diffVect(X, Y);
    float norme = 0;
    for (k = 0; k < C; k++) {
        norme += fabs(delta[k]);
    }
    free(delta);
    return norme;
}

// Addition entre deux vecteurs de proba (ajoute dans y)
void addVect(struct vecteur *X, struct vecteur *Y) {
    int C = X->C;
    if (C != Y->C) { exit(320); }
    indice i; for (i = 0; i < C; i++) {Y->v[i] += X->v[i];}
}

// Multiplie un vecteur de proba par un float
struct vecteur *multFloatProba(float a, struct vecteur *Vec) {
    int taille = (int) Vec->C;
    struct vecteur *res = alloueVecteur(taille);
    for (indice i = 0; i < taille; i++) {
        res->v[i] = a * Vec->v[i];
    }
    return res;
}

// Multiplie un vecteur de proba (horizontal) avec un vecteur de int (vertical)
float multVectPVectI(struct vecteur *Vec, int *f) {
    indice i;
    float res = 0.;
    for (i = 0; i < Vec->C; i++) {
        res += Vec->v[i] * f[i];
    }
    return res;
}

// Multiplie à gauche une matrice par un vecteur
void multVecMat(struct vecteur *X, struct matrice *Mat, struct vecteur *Y) {
    indice i, j, k;
    if ((X->C != Y->C) ||(Y->C != Mat->C)) { exit(321); }
    for (k = 0; k < Mat->M; k++) {
        i = Mat->P[k].i;
        j = Mat->P[k].j;
        float val = Mat->P[k].val;
        Y->v[j] += X->v[i] * val;
    }
}

void metZero(struct vecteur *Vec) {
    indice i;
    for (i = 0; i < Vec->C; i++) {
        Vec->v[i] = 0.0;
    }
}

void recopie(struct vecteur *X, struct vecteur *Y) {
    indice i;
    int C = X->C;
    if (C != Y->C) { exit(320); }
    for (i = 0; i < C; i++) {
        X->v[i] = Y->v[i];
    }
}


struct matrice *genereErdosStochastique(const struct matrice *Mat0, indice n, proba p) {
    
    indice C0 = Mat0->C;
    indice C  = C0 + n;             // nouvelle taille
    int M0    = Mat0->M;
    
    // Compte pour chaque nouveau sommet i = C0 .. C-1 le nombre d'arcs
    int *deg_out = calloc(n, sizeof *deg_out);
    if (!deg_out) {
        perror("calloc deg_out");
        exit(EXIT_FAILURE);
    }
    for (int di = 0; di < n; di++) {
        indice i = C0 + di;
        int    cnt = 0;
        for (indice j = 0; j < C; j++) {
            if (j == i) continue;
            if ((double)rand() / RAND_MAX < p) {
                cnt++;
            }
        }
        deg_out[di] = cnt;
    }

    // Calcule pour chaque di de la taille effective et l'offset dans P
    //    taille[di] = deg_out[di] > 0 ? deg_out[di] : C
    int *taille   = malloc(n * sizeof *taille);
    int *offset   = malloc(n * sizeof *offset);
    if (!taille || !offset) {
        perror("malloc taille/offset");
        exit(EXIT_FAILURE);
    }

    for (int di = 0; di < n; di++) {
        taille[di] = (deg_out[di] > 0 ? deg_out[di] : C);
    }
    offset[0] = M0;
    for (int di = 1; di < n; di++) {
        offset[di] = offset[di - 1] + taille[di - 1];
    }
    int M_new = offset[n - 1] + taille[n - 1];  // nombre total d'éléments

    struct elem *P = malloc(M_new * sizeof *P);
    if (!P) {
        perror("malloc P");
        exit(EXIT_FAILURE);
    }
    struct matrice *Mat = malloc(sizeof *Mat);
    if (!Mat) {
        perror("malloc Mat");
        exit(EXIT_FAILURE);
    }
    Mat->C = C;
    Mat->M = M_new;
    Mat->P = P;

    printf("Copie de la matrice initiale\n");
    memcpy(P, Mat0->P, M0 * sizeof *P);

    printf("Remplissage des offset\n");
    // Remplir les blocs pour chaque nouveau sommet
    for (int di = 0; di < n; di++) {
        indice i     = C0 + di;
        int    start = offset[di];
        int    written = 0;
        if (deg_out[di] > 0) {
            // val stochastique = 1/deg_out[di]
            proba v = 1.0f / deg_out[di];
            for (indice j = 0; j < C; j++) {
                if (j == i) continue;
                if ((double)rand() / RAND_MAX < p) {
                    P[start + written].i   = i;
                    P[start + written].j   = j;
                    P[start + written].val = v;
                    written++;
                }
            }
        }
    }

    free(deg_out);
    free(taille);
    free(offset);

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
    struct resultat_puissances *res = malloc(sizeof(int) + sizeof(struct vecteur));
    res->V = Pi;
    res->n = cnt;
    return res;
}

FILE *initialiseGNU(char nom[]) {

    if (!nom) { nom = "plot"; }
    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot) {

        fprintf(gnuplot, "reset\n");    
        fprintf(gnuplot, "set terminal pngcairo size 1000,700 font 'Helvetica,12'\n");
        fprintf(gnuplot, "set output '%s.png'\n", nom);

        fprintf(gnuplot, "set lmargin 12\n");
        fprintf(gnuplot, "set rmargin 4\n");
        fprintf(gnuplot, "set bmargin 5\n");
        fprintf(gnuplot, "set tmargin 3\n");

        fprintf(gnuplot, "set xlabel 'alpha'\n");
        fprintf(gnuplot, "set ylabel 'Nombre etapes avant la convergence'\n");
        fprintf(gnuplot, "set title 'Nombre etapes avant la convergence en fonction de alpha'\n");
        fprintf(gnuplot, "set grid\n");

        fprintf(gnuplot, "plot '-' with linespoints title 'Convergence'\n");

        return gnuplot;
    } else { return 0; }
}

struct vecteur **plot(struct matrice *Mat, struct vecteur **Pi_tab, char nom[]) {

    FILE *gnuplot = initialiseGNU(nom);
    double pas = 20.;
    struct vecteur **vecteurs = malloc(pas * sizeof(struct vecteur *));
    if (vecteurs == NULL) { exit(501); }

    if (gnuplot) {
        printf("Starting\n");
        for (int i = 0; i < pas - 1; i++) {
            float alpha = 0 + (i * (1 / pas));
            long long moyenne_iter = 0;

            struct vecteur *init = (Pi_tab != NULL) ? Pi_tab[i] : NULL;
            struct resultat_puissances *res = puissances(Mat, init, alpha);
            moyenne_iter += res->n;
            vecteurs[i] = res->V;

            fprintf(gnuplot, "%f %lld\n", alpha, moyenne_iter);
            printf("%f, %lld\n", alpha, moyenne_iter);
        }
    printf("Done\n");
    }
    return vecteurs;
}

int main(int argc, char *argv[]) {

    char nom[] = "Matrices/StanfordBerkeley.txt";

    struct matrice *Mat = lectureMatrice(nom);

    struct vecteur **vecteurs = plot(Mat, 0, "plot_depart");
    struct matrice *erdos = genereErdosStochastique(Mat, 600, 0.8);
    //plot(erdos, 0, "plot_erdos_N");
    //affichemat(erdos);
    plot(erdos, vecteurs, "plot_erdos_pi");
    return 0;
}