#include<stdio.h>
#include"calcul.h"

void affichemat(struct matrice *Mat) {
    indice k; for (k = 0; k < Mat->M; k++) {{printf("%d %d %f\n", Mat->P[k].i, Mat->P[k].j, Mat->P[k].val);}
    }
}

void affichevec(struct vecteur *Vec) {
    printf("Affichage vecteur\n");
    if (!Vec || !Vec->v) { exit(6); }
    indice i; for (i = 0; i < Vec->C; i++) {printf("%f ", Vec->v[i]);}
    printf("\n");
}

FILE *initialiseGNU(const char *ylabel, const char *xlabel, const char *nom) {

    if (nom == NULL) { nom = "plot"; }
    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot) {

        fprintf(gnuplot, "reset\n");
        fprintf(gnuplot, "set terminal pngcairo size 1000,700 font 'Helvetica,12'\n");

        fprintf(gnuplot, "set output '%s.png'\n", nom);
        fprintf(gnuplot, "set lmargin 12\n");
        fprintf(gnuplot, "set rmargin 4\n");
        fprintf(gnuplot, "set bmargin 5\n");
        fprintf(gnuplot, "set tmargin 3\n");

        fprintf(gnuplot, "set xlabel '%s'\n", xlabel);
        fprintf(gnuplot, "set ylabel '%s'\n", ylabel);
        fprintf(gnuplot, "set title 'Nombre etapes avant la convergence en fonction de alpha'\n");
        fprintf(gnuplot, "set grid\n");

        return gnuplot;
    } else { return 0; }
}

void plot(int *moyennes, int pas, char *nom) {
    FILE *gnuplot = initialiseGNU("Nombre etapes avant convergence", "Alpha", nom);
    printf("Plotting\n");
    fprintf(gnuplot, "plot '-' with linespoints title 'Convergence'\n");
    for (int i = 0; i <= pas; i++) {
        
        float alpha = 0 + (i * (1 / (double) pas));
        
        fprintf(gnuplot, "%f %d\n", alpha, moyennes[i]);
    }

    fprintf(gnuplot, "e\n");
    fflush(gnuplot);
}

void plot_diff(int *etapesN, int* etapesPi, int pas, char *nom) {
    FILE *gnuplot = initialiseGNU("Nombre etapes avant convergence", "Alpha", nom);
    printf("Plotting diff\n");
    fprintf(gnuplot, "plot '-' with linespoints title 'Initialisation 1/N', ");
    fprintf(gnuplot, "'-' with linespoints title 'Initialisation Pi[i-1]'\n");
    for (int i = 0; i < pas - 1; i++) {
        float alpha = 0 + (i * (1 / (double) pas));
        fprintf(gnuplot, "%f %d\n", alpha, etapesN[i]);
    }
    fprintf(gnuplot, "e\n");
    fflush(gnuplot);

    for (int i = 0; i < pas - 1; i++) {
        float alpha = 0 + (i * (1 / (double) pas));
        fprintf(gnuplot, "%f %d\n", alpha, etapesPi[i]);
    }
    fprintf(gnuplot, "e\n");
    pclose(gnuplot);
}
