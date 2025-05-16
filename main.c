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

struct resultat_calculs {
    struct vecteur **vecteurs;
    int *Nb_iterations; 
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

    // lecture jusqu'à avoir M éléments ou EOF
    int courant = 0;
    while (courant < M && fscanf(F, "%d %d", &i, &nb) == 2) {
        if (i < 1 || i > C) {
            fprintf(stderr, "Index de ligne hors bornes: %d\n", i);
            exit(103);
        }

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

    // Verif
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
// n -> Proportion de noeuds à rajouter (%)
// p -> Proba d'avoir un arc entre un nouveau noeud et un ancien
struct matrice *genereErdosStochastique(const struct matrice *Mat0, float new_nodes, proba p) {
    int C0 = Mat0->C;
    int C  = C0 + new_nodes;
    int M0 = Mat0->M;

    // capacité initiale
    int cap = M0;
    struct elem *tmp = malloc(cap * sizeof *tmp);
    if (!tmp) { perror("malloc tmp"); exit(EXIT_FAILURE); }
    // copie de la matrice de départ
    memcpy(tmp, Mat0->P, M0 * sizeof *tmp);
    int len = M0;

    // buffer temporaire pour la collecte de voisins
    int *vois = malloc(C * sizeof *vois);
    if (!vois) { perror("malloc vois"); exit(EXIT_FAILURE); }

    // génération Erdős–Rényi pour les nouveaux sommets
    for (int di = 0; di < new_nodes; di++) {
        int i   = C0 + di, deg = 0;
        for (int j = 0; j < C0; j++) {
            if ((double)rand()/RAND_MAX < p) {
                vois[deg++] = j;
            }
        }
        if (deg == 0) {  // éviter division par 0
            for (int j = 0; j < C0; j++) vois[deg++] = j;
        }
        proba w = 1.0f / deg;
        for (int t = 0; t < deg; t++) {
            if (len >= cap) {
                cap *= 2;
                struct elem *x = realloc(tmp, cap * sizeof *tmp);
                if (!x) { perror("realloc tmp"); free(tmp); exit(EXIT_FAILURE); }
                tmp = x;
            }
            tmp[len].i   = i;
            tmp[len].j   = vois[t];
            tmp[len].val = w;
            len++;
        }
    }
    free(vois);

    struct matrice *Mat = malloc(sizeof *Mat);
    if (!Mat) { perror("malloc Mat"); free(tmp); exit(EXIT_FAILURE); }
    Mat->C = C;
    Mat->M = len;
    Mat->P = tmp;

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
    free(w);
    struct resultat_puissances *res = malloc(sizeof *res);
    if (!res) {perror("malloc res"); exit(27);}
    res->V = Pi;
    res->n = cnt;
    return res;
}

struct resultat_calculs *calcule_prob_stat(struct matrice *Mat, struct vecteur **Pi_tab, int pas) {

    printf("C = %d, M = %d\n", Mat->C, Mat->M);
    int *Nb_iterations = malloc(pas * sizeof(int));
    if (Nb_iterations == NULL) { exit(500); }

    struct vecteur **vecteurs = malloc(pas * sizeof(*vecteurs));
    if (vecteurs == NULL) { exit(501); }

    for (int i = 0; i < pas; i++) {
        float alpha = 0 + (i * (1 / (double) pas));
        struct vecteur *init = (Pi_tab != NULL) ? Pi_tab[i] : NULL;
        struct resultat_puissances *res = puissances(Mat, init, alpha);
        Nb_iterations[i] = res->n;
        vecteurs[i] = res->V;
        printf("%f, %d\n", alpha, Nb_iterations[i]);
    }
    
    struct resultat_calculs *res = malloc(sizeof(struct resultat_calculs));
    res->Nb_iterations = Nb_iterations;
    res->vecteurs = vecteurs;
    return res;
}

int main(int argc, char *argv[]) {

    srand(time(NULL));
    char nom[] = "Matrices/StanfordBerkeley.txt";

    struct matrice *Mat = lectureMatrice(nom);
    int pas = 10;
    struct resultat_calculs *res = calcule_prob_stat(Mat, NULL, pas);
    struct vecteur **vecteurs = res->vecteurs;
    plot(res->Nb_iterations, pas, "Depart");

    struct matrice *erdos = genereErdosStochastique(Mat, 11, 0.5);
    
    struct resultat_calculs *resErdosPi = calcule_prob_stat(erdos, vecteurs, pas);
    struct resultat_calculs *resErdosN = calcule_prob_stat(erdos, NULL, pas);

    plot_diff(resErdosN->Nb_iterations, resErdosPi->Nb_iterations, pas, "Différence");

    return 0;
}