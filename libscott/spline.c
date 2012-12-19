#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "spline.h"
#include "splineimpl.h"


Spline * CreateSpline(double * x,
                      double * y,
                      size_t n){
    Spline * s = malloc(sizeof(*s));
    s->n = n;
    s->x = malloc(n*sizeof(*(s->x)));
    s->y = malloc(n*sizeof(*(s->y)));
    memcpy(s->x,x,n*sizeof(*x));
    memcpy(s->y,y,n*sizeof(*y));
int i = 0;
    for (i=1;i<n;++i){
        if (x[i]<=x[i-1]){
            fprintf(stdout,"Spline inversion %e\t%e (size %zu)\n",
                    x[i-1],x[i],n);
        }
    }

    s->impl = CreateSplineImpl(s);
    return s;
}

void FreeSpline(Spline * s){
    if (s){
       FreeSplineImpl(s->impl);
       free(s->x);
       free(s->y);
    }
    free(s);
}

Spline * CopySpline(Spline * s){
    return CreateSpline(s->x,s->y,s->n);
}

// Returns interpolation at x
double EvalSpline(Spline * s,
        double x){
    if (!s) return 0;
    return EvalSplineImpl(s,x);
}
double EvalSplineDeriv(Spline * s,
        double x){
    if (!s) return 0;
    return EvalSplineDerivImpl(s,x);
}

// Call after changing the vectors
void RearmSpline(Spline * s){
    FreeSplineImpl(s->impl);
    s->impl = CreateSplineImpl(s);
}

Spline * ReadSpline(FILE * stream){
    if (!stream) return NULL;

    size_t n = 256;
    double * x = malloc(n*sizeof(*x));
    double * y = malloc(n*sizeof(*y));
    size_t used = 0;

    double tmpx,tmpy;
    while(fscanf(stream,"%le\t%le\n",&tmpx,&tmpy)==2){
        if (used == n){
            n*=2;
            x = realloc(x,n*sizeof(*x));
            y = realloc(y,n*sizeof(*y));
        }
        x[used] = tmpx;
        y[used] = tmpy;
        ++used;
    }

    Spline * s = CreateSpline(x,y,used);
    free(x);
    free(y);
    return s;
}
