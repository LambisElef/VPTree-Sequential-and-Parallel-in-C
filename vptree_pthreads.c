#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "vptree.h"

typedef struct _dThreadsArg {
    double *X;
    double *vpd;
    struct _vptree *vpn;
    int first;
    int last;
    int d;
} dThreadsArg;

typedef struct _sThreadsArg {
    double *X;
    double *vpd;
    int *idx;
    struct _vptree *vpn;
    int n;
    int d;
    int isInner;
} sThreadsArg;

vptree * myBuildvp(double *X, int n, int d, int *idx);
double quickselectMedian(double *vpd, int vpdSize);
double quickselect(double *vpd, int vpdSize, int k);
void *calcDistance(void *arg);
void *calcSubTree(void *arg);

int dThreadNum = 12;

vptree * buildvp(double *X, int n, int d) {
    // Initialization of tree node.
    vptree *root = (vptree *)malloc(sizeof(vptree));

    // Initialization and setup of idx array.
    int *idx = (int *)malloc(n * sizeof(int));
    for (int i=0; i<n; i++)
        idx[i] = i;

    // Start building the tree.
    root = myBuildvp(X, n, d, idx);

    // Cleanup.
    free(idx);

    // Returns tree node.
    return root;
}

vptree * myBuildvp(double *X, int n, int d, int *idx) {
    // Initialization of tree node.
    vptree *vpn = (vptree *)malloc(sizeof(vptree));

    // If node is empty.
    if (n == 0)
        return NULL;

    // Initialization of vantage point and setting it as the first element of X.
    vpn->vp = (double *)malloc(d * sizeof(double));
    for (int i=0; i<d; i++)
        vpn->vp[i] = X[i];

    // Setup of vantage point id in the original X Array.
    vpn->idx = idx[0];

    // If node is a leaf.
    if (n == 1) {
        vpn->inner = NULL;
        vpn->outer = NULL;
        vpn->md = 0;
        return vpn;
    }

    // Calculation of distances between vantage point and points below.
    double *vpd = (double *)malloc(n * sizeof(double));

    // Parallel code if trigger is met.
    if (n > 256) {
        // Sets the number of distances each thread will calculate.
        int dThreadElementNum = (n - n%dThreadNum) / dThreadNum;

        // Initialization of threads.
        pthread_t *dThreads = (pthread_t *)malloc(dThreadNum * sizeof(*dThreads));
        dThreadsArg *dThreadArg = (dThreadsArg *)malloc(dThreadNum * sizeof(dThreadsArg));

        // Initialization of thread arguments and creation of threads.
        for (int i=0; i<dThreadNum; i++) {
            dThreadArg[i].X = X;
            dThreadArg[i].vpn = vpn;
            dThreadArg[i].vpd = vpd;
            dThreadArg[i].first = i * dThreadElementNum;
            dThreadArg[i].last = (i+1) * dThreadElementNum;
            if ((i+2)*dThreadElementNum >= n)
                dThreadArg[i].last = (n-1);
            dThreadArg[i].d = d;
            pthread_create(&dThreads[i], NULL, calcDistance, (void *)(&dThreadArg[i]));
        }

        // Wait for all threads to finish.
        for (int i=0; i<dThreadNum; i++) {
            pthread_join(dThreads[i], NULL);
        }
    }

    // Sequential code.
    else {
        for (int i=0; i<(n-1); i++) {
        vpd[i] = 0;
        for (int j=0; j<d; j++)
            vpd[i] += pow(vpn->vp[j] - X[(i+1)*d+j],2);
        vpd[i] = sqrt(vpd[i]);
        }
    }

    // Calculation of median distance between vantage point and points.
    vpn->md = quickselectMedian(vpd, (n-1));

    // Parallel code if trigger is met.
    if (n > 256) {
        // Initialization of inner and outer sub tree threads.
        pthread_t innerThread, outerThread;

        // Inner Thread Arguments and Creation.
        sThreadsArg sInnerThreadArg;
        sInnerThreadArg.X = X;
        sInnerThreadArg.vpd = vpd;
        sInnerThreadArg.idx = idx;
        sInnerThreadArg.vpn = vpn;
        sInnerThreadArg.n = n;
        sInnerThreadArg.d = d;
        sInnerThreadArg.isInner = 1;
        pthread_create (&innerThread, NULL, calcSubTree, (void *)(&sInnerThreadArg));

        // Outer Thread Arguments and Creation.
        sThreadsArg sOuterThreadArg;
        sOuterThreadArg.X = X;
        sOuterThreadArg.vpd = vpd;
        sOuterThreadArg.idx = idx;
        sOuterThreadArg.vpn = vpn;
        sOuterThreadArg.n = n;
        sOuterThreadArg.d = d;
        sOuterThreadArg.isInner = 0;
        pthread_create (&outerThread, NULL, calcSubTree, (void *)(&sOuterThreadArg));

        // Wait for the two threads to finish.
        pthread_join(innerThread, NULL);
        pthread_join(outerThread, NULL);
    }

    // Sequential code.
    else {
        // Dynamic initialization of inner array, outer array, inner IDs and outer IDs.
        int inn = -1;
        int out = -1;

        double *I = (double *)malloc(d * sizeof(double));
        double *O = (double *)malloc(d * sizeof(double));

        int *Iidx = (int *)malloc(sizeof(int));
        int *Oidx = (int *)malloc(sizeof(int));

        // Dynamic split of points to inner and outer array.
        for(int i=0; i<(n-1); i++) {
            if (vpd[i] <= vpn->md) {
                inn++;
                I = (double *)realloc(I, (inn+1)*d * sizeof(double));
                for (int j=0; j<d; j++)
                    I[inn*d+j] = X[(i+1)*d+j];

                // Setup of point id in the original X Array.
                Iidx = (int *)realloc(Iidx, (inn+1) * sizeof(int));
                Iidx[inn] = idx[i+1];
            }
            else {
                out++;
                O = (double *)realloc(O, (out+1)*d * sizeof(double));
                for (int j=0; j<d; j++)
                    O[out*d+j] = X[(i+1)*d+j];

                // Setup of point id in the original X Array.
                Oidx = (int *)realloc(Oidx, (out+1) * sizeof(int));
                Oidx[out] = idx[i+1];
            }
        }

        vpn->inner = myBuildvp(I, inn+1, d, Iidx);
        vpn->outer = myBuildvp(O, out+1, d, Oidx);

        // Cleanup.
        free(I);
        free(O);

        free(Iidx);
        free(Oidx);
    }

    // Returns tree node.
    return vpn;
}

double quickselectMedian(double *vpd, int vpdSize) {
    if (vpdSize % 2 == 1)
        return quickselect(vpd, vpdSize, vpdSize/2);
    else
        return 0.5 * (quickselect(vpd, vpdSize, vpdSize/2 - 1) + quickselect(vpd, vpdSize, vpdSize/2));
}

double quickselect(double *vpd, int vpdSize, int k) {
    if (vpdSize == 1)
        return vpd[0];

    double pivot = vpd[0];

    double *lows = (double *)malloc(sizeof(double));
    double *highs = (double *)malloc(sizeof(double));
    double *pivots = (double *)malloc(sizeof(double));

    int l = -1;
    int h = -1;
    int p = -1;

    for (int i=0; i<vpdSize; i++) {
        if (vpd[i] < pivot) {
            l++;
            lows = (double *)realloc(lows, (l+1)*sizeof(double));
            lows[l] = vpd[i];
        }
        else if (vpd[i] > pivot) {
            h++;
            highs = (double *)realloc(highs, (h+1)*sizeof(double));
            highs[h] = vpd[i];
        }
        else {
            p++;
            pivots = (double *)realloc(pivots, (p+1)*sizeof(double));
            pivots[p] = vpd[i];
        }
    }

    if (k < (l+1))
        return quickselect(lows, (l+1), k);
    else if (k < (l+1) + (p+1))
        return pivots[0];
    else
        return quickselect(highs, (h+1), k - (l+1) - (p+1));
}

void *calcDistance(void *arg) {

    dThreadsArg *dThreadArg = (dThreadsArg *)arg;
    for (int i=dThreadArg->first; i<dThreadArg->last; i++) {
        dThreadArg->vpd[i] = 0;
        for (int j=0; j<dThreadArg->d; j++)
            dThreadArg->vpd[i] += pow(dThreadArg->vpn->vp[j] - dThreadArg->X[(i+1)*dThreadArg->d+j],2);
        dThreadArg->vpd[i] = sqrt(dThreadArg->vpd[i]);
    }
    return NULL;
}

void *calcSubTree(void *arg) {
    sThreadsArg *sThreadArg = (sThreadsArg *)arg;
    int counter = -1;
    double *A = (double *)malloc(sThreadArg->d * sizeof(double *));
    int *Aidx = (int *)malloc(sizeof(int));

    if (sThreadArg->isInner == 1) {
        for (int i=0; i<((sThreadArg->n)-1); i++) {
            if (sThreadArg->vpd[i] <= sThreadArg->vpn->md) {
                counter++;
                A = (double *)realloc(A, (counter+1)*sThreadArg->d * sizeof(double *));
                for (int j=0; j<sThreadArg->d; j++)
                    A[counter*sThreadArg->d+j] = sThreadArg->X[(i+1)*sThreadArg->d+j];

                // Setup of point id in the original X Array.
                Aidx = (int *)realloc(Aidx, (counter+1) * sizeof(int));
                Aidx[counter] = sThreadArg->idx[i+1];
            }
        }

        sThreadArg->vpn->inner = myBuildvp(A, counter+1, sThreadArg->d, Aidx);

        // Cleanup.
        free(A);
        free(Aidx);
    }
    else {
        for (int i=0; i<((sThreadArg->n)-1); i++) {
            if (sThreadArg->vpd[i] > sThreadArg->vpn->md) {
                counter++;
                A = (double *)realloc(A, (counter+1)*sThreadArg->d * sizeof(double *));
                for (int j=0; j<sThreadArg->d; j++)
                    A[counter*sThreadArg->d+j] = sThreadArg->X[(i+1)*sThreadArg->d+j];

                // Setup of point id in the original X Array.
                Aidx = (int *)realloc(Aidx, (counter+1) * sizeof(int));
                Aidx[counter] = sThreadArg->idx[i+1];
            }
        }

        sThreadArg->vpn->outer = myBuildvp(A, counter+1, sThreadArg->d, Aidx);

        // Cleanup.
        free(A);
        free(Aidx);
    }
    return NULL;
}

int main() {
    int n, d;
    /*
    // (int argc, char* argv[])

    if (argc != 3) {
        printf("Usage: %s n\n  where n is no. of points and d. is no. of dimensions.\n", argv[0]);
        return(1);
    }

    n = atoi(argv[1]);
    d = atoi(argv[2]);
    */
    n = 1000;
    d = 1000;

    srand(time(NULL));

    // Initialization of points in a d-Dimensional space.
    double *X = (double *)malloc(n*d * sizeof(double));

    for (int i=0; i<n; i++)
        for (int j=0; j<d; j++)
            X[i*d+j] = (double)rand();

    // Builds the tree from the X points array.
    vptree *root = buildvp(X, n, d);

    // Cleanup.
    free(X);

    return 0;
}
