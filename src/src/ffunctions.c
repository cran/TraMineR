#include<R.h>
#include<Rmath.h>

void cLEVEN(int *seq1, int *seq2, double *param, double *scost, double *result) {

    int k,i,j,n,m,nbe,idel;
    //Le vecteur param contient la taille des deux séquences (identique), le nombre d'états, les indels
    double cost,a,b,c,min;
    n = (int)param[0];
    m = (int)param[1];
    nbe = (int)param[3];
    idel = (int)param[2];
    //Rprintf("indel = %f\n", (float)idel);
    //Rprintf("nbetats = %i\n", nbe);
    // Rprintf("n = %i, m = %i\n", n, m);
    n++;
    m++;

    double** leven = (double**) R_alloc(m, sizeof(double*));
    for (i = 0; i < m; i++)
        leven[i] = (double*) R_alloc(n, sizeof(double));

    //Rprintf("Array with n = %i, m = %i\n", n, m);
    leven[0][0] = (double)0;
    //Rprintf("seq1[0] = %i", seq1[0]);
    //Rprintf("seq1[1] = %i", seq1[1]);
    //Rprintf("leven[0][0] = %f\n", leven[0][0]);

    for (k = 1; k < n; k++) {
        leven[0][k] = leven[0][k-1] + idel;
        // Rprintf("%f", leven[0][k]);
    }

    for (k=1;k<m;k++) {
        leven[k][0] = leven[k-1][0] + idel;
        // Rprintf("%f\n", leven[k][0]);
    }
    //Rprintf("Matrix initialized");




    for (i=1;i<m;i++) {
        for (j=1;j<n;j++) {
            if (seq1[j-1] == seq2[i-1]) {
                //Rprintf("Seq1[%i] et Seq2[%i] sont �gaux, %i == %i \n", j-1, i-1, seq1[j-1], seq2[i-1]);
                cost = 0;
            } else {
                int cell = (int)seq1[j-1]+(seq2[i-1]*nbe);

                cost = scost[cell];
//Rprintf("cell number %i, value = %i \n", cell, cost);
                //cost = 2;
            }


            a=leven[i][j-1] + idel;
            b=leven[i-1][j]+idel;
            c=leven[i-1][j-1]+cost;
            min=a;
            if (b<min)
                min=b;
            if (c<min)
                min=c;
            // d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
            //Rprintf("min : %f\n", min);
            leven[i][j]=min;
        }
    }

    for (k = 1; k < n; k++) {
        leven[0][k] = leven[0][k-1] + idel;
        // Rprintf("%f", leven[0][k]);
    }


    // Rprintf("ok : %f", leven[n-1][m-1]);
    *result = leven[m-1][n-1];

}


void cLCS(int *iseq, int *jseq , double *length, int *result) {

    int np,mp,max,a,b,i,j,n,m,temp;
    n = (int)length[1];
    m = (int)length[0];
    np = n + 1;
    mp = m + 1;
    //  Rprintf("m = %i, n = %i \n", m, n);
    //    int** L = (int**) malloc(np * sizeof(int*));
    //    for (k = 0; k< mp;k++) L[k] = (int*)malloc(mp * sizeof(int));
    int L [mp][np];
    //Rprintf("Creation du tableau ok");
    for (i = 0; i<mp; i++) {
        for (j= 0; j<np; j++) {
            L[i][j] = (int)0;
            //	 Rprintf("i = %i, j = %i\n", i, j);
        }
    }
    //    Rprintf("tab : %i \n", L[0][1]);


    for (i = 1; i <= m;i++) {
        temp = iseq[i-1];
        //  Rprintf("i = %i, temp = %i\n", i, temp);
        for (j = 1; j<= n;j++) {
            if (temp == jseq[j-1]) {
                L[i][j] = 1 + L[i-1][j-1];
                //	   Rprintf("temp = %i jseq = %i, L[%i][%i] = %i\n", temp, jseq[j-1], i, j, L[i][j]);
            } else {
                //   Rprintf("�gal 0 pour i %d", i);
                a = L[i-1][j];
                b = L[i][j-1];
                max = a;
                if (a < b) max = b;
                //   Rprintf("max = %i \n", max);
                L[i][j] = max;
            }
        }
    }

    *result = L[m][n];
}



void cLCP(int *iseq, int *jseq , double *length, int *result) {

	int minlength=imin2(length[0], length[1]);
    int i=0;
 

    while (i<minlength && iseq[i]==jseq[i]) {
		i++;
    }

    *result = i;
}

