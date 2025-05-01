#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h> //struct timeval t1, t2; gettimeofday(&t1, NULL); gettimeofday(&t2, NULL); 

typedef int indice;
typedef float proba;

proba alpha = 0.85; // Proba du surfeur aléatoire
proba invAlpha = 1 - 0.85;

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

void affichemat(struct matrice Mat) {
    indice k; for (k = 0; k < Mat.M; k++) {{printf("%d %d %f\n", Mat.P[k].i, Mat.P[k].j, Mat.P[k].val);}
    }
}

struct matrice *lectureMatrice(char *nom_fic) {
    printf("Lecture Matrice\n");
    indice i, j;
    int nb, C, M;
    struct matrice *Mat = malloc(sizeof(struct matrice));
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
    printf("Lecture matrice\n");
    for (int k = 0; k < C; k++) {
        fscanf(F, "%d", &(i));
        fscanf(F,"%d",  &(nb)); 
        for (int z = 0; z < nb; z++) {
            struct elem e;

            fscanf(F, "%d", &(j));
            fscanf(F, "%le", &(v));
            e.i = i - 1;
            e.j = j - 1;
            e.val = v;
            P[courant] = e;
            courant++;
            if (courant > M) {
                exit(202);
            }
        }
    }
    printf("Nombre d'elements attendu: %d, nombre réel: %d\n", M, courant);
    fclose(F);

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

// struct elem *genereErdos(indice lignes, proba p) {
//     struct elem *P = malloc(1);
//     for () {
// 
//     }
// }

void iterer(struct vecteur *X, struct vecteur *Y, struct vecteur *W, struct matrice *Mat) {
    struct vecteur *aX = multFloatProba(alpha, X);

    multVecMat(aX, Mat, Y);

    proba gauche = invAlpha / (float) Mat->C;
    proba droite = (alpha / (float) Mat->C) * multVectPVectI(X, f);

    initvec(W, gauche + droite, Mat->C);
    addVect(W, Y);
}

struct vecteur *puissances(struct matrice *Mat, struct vecteur *Pi) {
    int C = Mat->C;
    printf("%d\n", C);
    if (!Pi) {

        Pi = alloueVecteur(C);
        initvec(Pi, 1/(double) C, C);  // x0 = (1/N)e

    } else {  // adapter en fonction du nombre de nouveaux elements
        if (C != Pi->C) { exit(320); }
        struct vecteur *Z = alloueVecteur(C);

        for (indice i = 0; i < C; i++) {
            if (Pi->v[i]) {
                Z->v[i] = Pi->v[i];
            } else {
                Z->v[i] = 0.0;
            }
        }
        printf("Vecteur fourni: \n");

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
        iterer(Pi, y, w, Mat);
        delta = normeVect(Pi, y);
        recopie(Pi, y);
    }

    free(f);
    free(y);
    return Pi;
}

FILE *initialiseGNU() {

    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot) {

        fprintf(gnuplot, "reset\n");    
        fprintf(gnuplot, "set terminal pngcairo size 1000,700 font 'Helvetica,12'\n");
        fprintf(gnuplot, "set output 'plot.png'\n");

        fprintf(gnuplot, "set lmargin 12\n");
        fprintf(gnuplot, "set rmargin 4\n");
        fprintf(gnuplot, "set bmargin 5\n");
        fprintf(gnuplot, "set tmargin 3\n");

        fprintf(gnuplot, "set xlabel 'alpha'\n");
        fprintf(gnuplot, "set ylabel 'Temps (secondes)'\n");
        fprintf(gnuplot, "set title 'Temps de convergence en fonction de alpha'\n");
        fprintf(gnuplot, "set grid\n");

        fprintf(gnuplot, "plot '-' with linespoints title 'Convergence'\n");

        return gnuplot;
    } else { return 0; }
}

int main(int argc, char *argv[]) {

    FILE *gnuplot = initialiseGNU();
    if (gnuplot) {
        
        int cnt = 0;
        double pas = 30.;
        char nom[] = "Matrices/StanfordBerkeley.txt";
        struct matrice *Mat = lectureMatrice(nom);

        for (int i = 0; i < pas - 1; i++) {

            alpha = 0 + (i * (1 / pas));
            invAlpha = 1 - alpha;
            long long moyenne_temps = 0;            
            long moyenne = 1.;

            for (int index = 0; index < moyenne; index++) {

                cnt ++;
                gettimeofday(&t1, NULL);
                puissances(Mat, 0);
                gettimeofday(&t2, NULL);
                moyenne_temps += (t2.tv_sec - t1.tv_sec) * 1000000L + (t2.tv_usec - t1.tv_usec);

            }

            fprintf(gnuplot, "%f %f\n", alpha, (moyenne_temps / moyenne) / 1e6);
            printf("%f, %f\n", alpha, (moyenne_temps/moyenne) / 1e6);     
        }
    }

    fprintf(gnuplot, "e\n");
    fprintf(gnuplot, "unset output\n");

    //    printf("\nitération 2\n");
    //    x = puissances("matriceCreuseV2.txt", x);
    //
    //    printf("\nitération 3\n");
    //    x = puissances("matriceCreuseV3.txt", x);

    // gettimeofday(&t1, NULL);
    // puissances("Matrices/StanfordBerkeley.txt", 0);
    // gettimeofday(&t2, NULL);
    // long temps_micro = (t2.tv_sec - t1.tv_sec) * 1000000L + (t2.tv_usec - t1.tv_usec);
    // printf("Converge en %ld micro sec\n", temps_micro);
    return 0;
}