#include<main.c>

// Multiplie un vecteur de proba par un float
struct vecteur *multFloatProba(float a, struct vecteur *Vec) {
    int taille = (int) Vec->C;
    struct vecteur *res = alloueVecteur(taille);
    for (indice i = 0; i < taille; i++) {
        res->v[i] = a * Vec->v[i];
    }
    return res;
}

// Addition entre deux vecteurs de proba (ajoute dans y)
void addVect(struct vecteur *X, struct vecteur *Y) {
    int C = X->C;
    if (C != Y->C) { exit(320); }
    indice i; for (i = 0; i < C; i++) {Y->v[i] += X->v[i];}
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

// Multiplie un vecteur de proba (horizontal) avec un vecteur de int (vertical)
float multVectPVectI(struct vecteur *Vec, int *f) {
    indice i;
    float res = 0.;
    for (i = 0; i < Vec->C; i++) {
        res += Vec->v[i] * f[i];
    }
    return res;
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

void metZero(struct vecteur *Vec) {
    indice i;
    for (i = 0; i < Vec->C; i++) {
        Vec->v[i] = 0.0;
    }
}