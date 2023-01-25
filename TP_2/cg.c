#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "time.h"
#include "omp.h"


double rand01(void)
{
    return (double) rand() / ((double) RAND_MAX);
}
int main(int argc, char **argv) {
    int        i, j, k, n, max_iter;
    double     start, finish, elapsed;
    double    tol;
    double beta, alpha;
    double rr, rr0, rr1, pq;
   
    /* Paramètres */
    n = atoi(argv[1]);  
    tol = atof(argv[2]);
    max_iter = atoi(argv[3]); 
     
    //
    //--------  Données matrice et second-membre ----------------
    //--------     Ne pas paralleliser  -------------------------
    double* A = malloc(n * n * sizeof(double));
    double* B = malloc(n * n * sizeof(double));
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            A[i*n+j] = rand01();
        }
    }
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            B[i*n+j] = A[i*n+j] + A[j*n+i];
        }
    }
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            A[i*n+j] = B[i*n+j];
        }
        A[i*n+i] = n;
    }
    free(B);
    double* b = malloc(n * sizeof(double));  
    for (i = 0; i < n; i++){
        b[i] = 0.0;
        for (j = 0; j < n; j++) b[i] += A[i*n+j];
    }
    //-------------------------------------------------------------
    //--------------------------------------------------------------
    
    
    
    /* Prepare variables for the main loop */
    k = 0;
    
    double* x = malloc(n * sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* p = malloc(n * sizeof(double));
    double* q = malloc(n * sizeof(double));
    
    
    /* Initialisation */
    rr = 0.0;
    for (i = 0; i < n; i++) {
        x[i] = 0.0;
        r[i] = b[i];
        p[i] = r[i];
        rr += r[i]*r[i]; 
    }
    rr0 = rr;
    printf("iter=%d  |g|=%15.8e \n",k,sqrt(rr));

    start = omp_get_wtime();
    
    omp_set_num_threads(4);
    while ((k < max_iter) && ( rr > tol*tol*rr0 )) {
        k++;
        
        /* q=A*p  */
        #pragma omp parallel private(i,j) shared(q, A, p){
            #pragma omp for
            for (i = 0; i < n; i++){
                q[i] = 0.0;
                for (j = 0; j < n; j++){
                    q[i] += A[i*n+j]*p[j];
                }
            }   
        }
        
    
        /* p.q */
        #pragma omp parallel private(i) shared(q, p, pq){
            pq = 0.0;
            #pragma omp for reduction(+:pq)
            for (i = 0; i < n; i++){
                pq += p[i]*q[i];
            }
        }
    
        /* alpha */
        alpha = rr/pq;
    
        /* Updates */
        #pragma omp parallel private(i, alpha) shared(x, p, r, q)
        for (i = 0; i < n; i++){
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * q[i];
        }
    
        /* beta */
        rr1=0.0;
        #pragma omp parallel reduction(+:rr1) private(i) shared(r)
        for (i = 0; i < n; i++){
            rr1 += r[i]*r[i];
        }

        beta = rr1/rr;
    
        /* New p */
        #pragma omp parallel private(i,j, beta) shared(r, p)
        for (i = 0; i < n; i++){
            p[i] = r[i] + beta*p[i];
        }
    
        rr=rr1;
        printf("Iter=%d  err=%15.8f \n",k,sqrt(rr/rr0));
    }

    finish = omp_get_wtime();
    elapsed = finish - start;
    
    
    printf("========= Solver Completed ========= \n");
    printf("The code to be timed took %lf seconds\n", elapsed);
    printf("Number of iterations: %d \n", k);
    
    int verif = 1;
    double epsilon = 0.0001;
    for (i=0; i<n; i++){
    if (abs(1 - x[i]) > epsilon){
        verif = 0;
    }
    }

    if (verif == 1){
        printf("Le résultat est correct\n");
    }
    else{
        printf("Le résultat est érroné\n");
    }

    /* Free all  allocated memory   */
    free(A);
    free(b);

    free(x);
    free(q);
    free(p);
    free(r);
  
    return 0;
} /* main */
