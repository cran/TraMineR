#include<R.h>
#include<Rmath.h>

// TraMineR 2, Pierre-Alexandre Fonta (2017): not deleted because it used by seqLLCS() in the JSS article
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
                //   Rprintf("égal 0 pour i %d", i);
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


// TraMineR 2, Pierre-Alexandre Fonta (2017): not deleted because it used by seqLLCP() in the JSS article
void cLCP(int *iseq, int *jseq , double *length, int *result) {

	int minlength=imin2(length[0], length[1]);
    int i=0;


    while (i<minlength && iseq[i]==jseq[i]) {
		i++;
    }

    *result = i;
}

