/*
 * sortutils.c
 * Version: 1.1
 *
 * Copyright (C) 2012 Frank Everdij, Delft University of Technology
 * F.P.X.Everdij@tudelft.nl
 *
 * See the COPYRIGHT file for copyright details.
 *
 */
#include "sortutils.h"

int countsortsparse (int *a, int *b, void *c, const int type, const int n)
{
    int i,count,gap,temp;
    int *e,*ec,emax,emin,esize; /*(minval(a):maxval(a))*/
    int *ar,*br;
    double *dc=NULL,*dcr=NULL;
    int *ic=NULL,*icr=NULL;

    ar=ialloc(n,0);
    memcpy(ar,a,n*sizeof(int));
    br=ialloc(n,0);
    memcpy(br,b,n*sizeof(int));

    /* determine maximum and minimum column indices */
    emin=a[0];
    emax=emin;
    for (i=0;i<n;i++) {
        if (a[i]<emin) emin=a[i];
        if (a[i]>emax) emax=a[i];
    }
    esize=emax-emin+1;
    e=ialloc(esize+1,1);
    for (i=0;i<n;i++) e[a[i]-emin]++;
    
    count=0;
    for (i=0; i<=esize;i++) {
        temp=e[i];
    	e[i]=count;
	    count+=temp;
	}

    ec=ialloc(esize+1,0);
    memcpy(ec,e,(esize+1)*sizeof(int));
    
    if (type==DOUBLE) {
        dc=(double *)c;
        dcr=dalloc(n,0);
        memcpy(dcr,dc,n*sizeof(double));
        for (i=0;i<n;i++) {
            temp=e[ar[i]-emin];
        	a[temp]=ar[i];
	        b[temp]=br[i];
  	        dc[temp]=dcr[i];
	        e[ar[i]-emin]=temp+1;
        }
        free(dcr);
        for (i=0;i<esize;i++) {
            gap=ec[i+1]-ec[i];
            dcr=dc+ec[i];
            if (gap<MAXGAP) insortsparse(b+ec[i], (void *)dcr, type, gap);
                else shellsortsparse(b+ec[i], (void *)dcr, type, gap);
        }
    } else {
        ic=(int *)c;
        icr=ialloc(n,0);
        memcpy(icr,ic,n*sizeof(int));
        for (i=0;i<n;i++) {
            temp=e[ar[i]-emin];
        	a[temp]=ar[i];
	        b[temp]=br[i];
  	        ic[temp]=icr[i];
	        e[ar[i]-emin]=temp+1;
        }
        free(icr);
        for (i=0;i<esize;i++) {
            gap=ec[i+1]-ec[i];
            icr=ic+ec[i];
            if (gap<MAXGAP) insortsparse(b+ec[i], (void *)icr, type, gap);
                else shellsortsparse(b+ec[i], (void *)icr, type, gap);
        }
    }

    free(ar);
    free(br);
    free(e);
    free(ec);
   
    return EXIT_SUCCESS;
}

/*
! 
!   Method:	Optimized insertion-sort (ala Jon Bentley)
!
*/
void insortsparse (int *a, void *c, const int type, const int n)
{

    int i,j;
    int temp;
    double *dc,dtemp;
    int *ic,itemp;

	if (n < 2) return;
    if (type==DOUBLE) {
        dc=(double*)c;
        for (i=1; i<n; i++) {
            j = i;
            temp = a[j];
            dtemp = dc[j];
            while ((j>0)&&(a[j-1] > temp)) {
                a[j] = a[j-1];
                dc[j] = dc[j-1];
                j--;
            }
            a[j] = temp;
            dc[j] = dtemp;
        }
    } else {
        ic=(int *)c;
        for (i=1; i<n; i++) {
            j = i;
            temp = a[j];
            itemp = ic[j];
            while ((j>0)&&(a[j-1] > temp)) {
                a[j] = a[j-1];
                ic[j] = ic[j-1];
                j--;
            }
            a[j] = temp;
            ic[j] = itemp;
        }
    }

    return;
}

/*
 * from http://www.softpanorama.org/Algorithms/Sorting/shellsort.shtml
 */
void shellsortsparse(int *a, void *c, const int type, const int n)
{
    int i, j, k, h, v, *ic, iv;
    double *dc, dv;
    int cols[] = {1391376, 463792, 198768, 86961, 33936, 13776, 4592,
                    1968, 861, 336, 112, 48, 21, 7, 3, 1};

    if (type==DOUBLE) {
        dc=(double *)c;
        for (k=0; k<16; k++) {
            h=cols[k];
            for (i=h; i<n; i++) {
                v=a[i];
                dv=dc[i];
                j=i;
                while (j>=h && a[j-h]>v) {
                    a[j]=a[j-h];
                    dc[j]=dc[j-h];
                    j=j-h;
                }
                a[j]=v;
                dc[j]=dv;
            }
        }
    } else {
        ic=(int *)c;
        for (k=0; k<16; k++) {
            h=cols[k];
            for (i=h; i<n; i++) {
                v=a[i];
                iv=ic[i];
                j=i;
                while (j>=h && a[j-h]>v) {
                    a[j]=a[j-h];
                    ic[j]=ic[j-h];
                    j=j-h;
                }
                a[j]=v;
                ic[j]=iv;
            }
        }
    }

    return;
}

/* (C) Copr. 1986-92 Numerical Recipes Software #. */
void heapsortsparse(int *a, void *c, const int type, const int n)
{
	unsigned long i,ir,j,l;
    int *ra,rra;
	double *rb,rrb;
    int *rc,rrc;

	if (n < 2) return;

    ra=a-1;

    if (type==DOUBLE) {
        rb=(double *)c;
        rb--;
    	l=(n >> 1)+1;
    	ir=n;
    	for (;;) {
    		if (l > 1) {
    			rra=ra[--l];
                rrb=rb[l];
    		} else {
    			rra=ra[ir];
                rrb=rb[ir];
    			ra[ir]=ra[1];
                rb[ir]=rb[1];
    			if (--ir == 1) {
    				ra[1]=rra;
                    rb[1]=rrb;
    				break;
    			}
    		}
    		i=l;
    		j=l+l;
    		while (j <= ir) {
    			if (j < ir && ra[j] < ra[j+1]) j++;
    			if (rra < ra[j]) {
    				ra[i]=ra[j];
                    rb[i]=rb[j];
    				i=j;
    				j <<= 1;
    			} else j=ir+1;
    		}
    		ra[i]=rra;
            rb[i]=rrb;
    	}
    } else {
        rc=(int *)c;
        rc--;
    	l=(n >> 1)+1;
    	ir=n;
    	for (;;) {
    		if (l > 1) {
    			rra=ra[--l];
                rrc=rc[l];
    		} else {
    			rra=ra[ir];
                rrc=rc[ir];
    			ra[ir]=ra[1];
                rc[ir]=rc[1];
    			if (--ir == 1) {
    				ra[1]=rra;
                    rc[1]=rrc;
    				break;
    			}
    		}
    		i=l;
    		j=l+l;
    		while (j <= ir) {
    			if (j < ir && ra[j] < ra[j+1]) j++;
    			if (rra < ra[j]) {
    				ra[i]=ra[j];
                    rc[i]=rc[j];
    				i=j;
    				j <<= 1;
    			} else j=ir+1;
    		}
    		ra[i]=rra;
            rc[i]=rrc;
    	}
    }

    return;
}

int countsortint (int *a, const int n)
{
    int i,count,temp,err;
    int *e,emax,emin; /*(minval(a):maxval(a))*/
    int *ar;

    ar=ialloc(n,0);
    memcpy(ar,a,n*sizeof(int));
    
    emin=a[0];
    emax=emin;
    for (i=0;i<n;i++) {
        if (a[i]<emin) emin=a[i];
        if (a[i]>emax) emax=a[i];
    }
    
    e=ialloc(emax-emin+1,1);
    for (i=0;i<n;i++) e[a[i]-emin]++;
    
    count=0;
    for (i=0; i<=(emax-emin);i++) {
        temp=e[i];
    	e[i]=count;
	    count+=temp;
	}
	
    for (i=0;i<n;i++) {
        temp=e[ar[i]-emin];
    	a[temp]=ar[i];
	    e[ar[i]-emin]=temp+1;
    }
    
    free(ar);
    free(e);
    
    err=insortint(a, n);
   
    return EXIT_SUCCESS;
}

/*
! 
!   Method:	Optimized insertion-sort (ala Jon Bentley)
!
*/
int insortint (int *a, const int n)
{

    int i,j;
    int temp;

	if (n < 2) return EXIT_SUCCESS;
    for (i=1; i<n; i++) {
        j = i;
        temp = a[j];
        while ((j>0)&&(a[j-1] > temp)) {
            a[j] = a[j-1];
            j--;
        }
        a[j] = temp;
    }    
    return EXIT_SUCCESS;
}

/* (C) Copr. 1986-92 Numerical Recipes Software #. */
int * heapsortindex(void *b, const int type, const int n)
{
	unsigned long i,ir,j,l;
    int *a,*ra,rra;
	double *rb,rrb;
    int *rc,rrc;

    a=ialloc(n,1);
    for (i=0;i<n;i++) a[i]=i;
	if (n < 2) return a;

    ra=a-1;

    if (type==DOUBLE) {
        rb=(double *)b;
        rb--;
    	l=(n >> 1)+1;
    	ir=n;
    	for (;;) {
    		if (l > 1) {
    			rra=ra[--l];
                rrb=rb[l];
    		} else {
    			rra=ra[ir];
                rrb=rb[ir];
    			ra[ir]=ra[1];
                rb[ir]=rb[1];
    			if (--ir == 1) {
    				ra[1]=rra;
                    rb[1]=rrb;
    				break;
    			}
    		}
    		i=l;
    		j=l+l;
    		while (j <= ir) {
    			if (j < ir && rb[j] < rb[j+1]) j++;
    			if (rrb < rb[j]) {
    				ra[i]=ra[j];
                    rb[i]=rb[j];
    				i=j;
    				j <<= 1;
    			} else j=ir+1;
    		}
    		ra[i]=rra;
            rb[i]=rrb;
    	}
    } else {
        rc=(int *)b;
        rc--;
    	l=(n >> 1)+1;
    	ir=n;
    	for (;;) {
    		if (l > 1) {
    			rra=ra[--l];
                rrc=rc[l];
    		} else {
    			rra=ra[ir];
                rrc=rc[ir];
    			ra[ir]=ra[1];
                rc[ir]=rc[1];
    			if (--ir == 1) {
    				ra[1]=rra;
                    rc[1]=rrc;
    				break;
    			}
    		}
    		i=l;
    		j=l+l;
    		while (j <= ir) {
    			if (j < ir && rc[j] < rc[j+1]) j++;
    			if (rrc < rc[j]) {
    				ra[i]=ra[j];
                    rc[i]=rc[j];
    				i=j;
    				j <<= 1;
    			} else j=ir+1;
    		}
    		ra[i]=rra;
            rc[i]=rrc;
    	}
    }

    return a;
}
