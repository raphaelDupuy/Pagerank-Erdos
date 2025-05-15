#ifndef CALCUL_H
#define CALCUL_H

#include<stdio.h>

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

struct vecteur *multFloatProba(float a, struct vecteur *Vec);
void            addVect(struct vecteur *X, struct vecteur *Y);
void            multVecMat(struct vecteur *X, struct matrice *Mat, struct vecteur *Y);
float           multVectPVectI(struct vecteur *Vec, int *f);
float           normeVect(struct vecteur *X, struct vecteur *Y);
proba          *diffVect(struct vecteur *X, struct vecteur *Y);
void            metZero(struct vecteur *Vec);
struct vecteur *alloueVecteur(int taille);

#endif