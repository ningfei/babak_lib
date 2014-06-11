#define _hpsort

void hpsort(unsigned long n, float *ra)
{    
        unsigned long i,ir,j,l;
        float rra;
     
        if (n < 2) return;
        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l-1]; 
                } else {
                        rra=ra[ir-1];
                        ra[ir-1]=ra[0];
                        if (--ir == 1) {
                                ra[0]=rra;
                                break;
                        }
                }
                i=l;
                j=l+l;
                while (j <= ir) {
                        if (j < ir && ra[j-1] < ra[j]) j++;
                        if (rra < ra[j-1]) {
                                ra[i-1]=ra[j-1];
                                i=j;
                                j <<= 1;
                        } else j=ir+1;
                }
                ra[i-1]=rra;
        }
}    

void hpsort(unsigned long n, float *ra, int *indx)
{    
        unsigned long i,ir,j,l;
        float rra;
		int rindx;
     
        if (n < 2) return;
        l=(n >> 1)+1;
        ir=n;
        for (;;) {
                if (l > 1) {
                        rra=ra[--l-1]; 
   						rindx=indx[l-1]; 
                } else {
                        rra=ra[ir-1];
                        rindx=indx[ir-1];
                        ra[ir-1]=ra[0];
                        indx[ir-1]=indx[0];
                        if (--ir == 1) {
                                ra[0]=rra;
                                indx[0]=rindx;
                                break;
                        }
                }
                i=l;
                j=l+l;
                while (j <= ir) {
                        if (j < ir && ra[j-1] < ra[j]) j++;
                        if (rra < ra[j-1]) {
                                ra[i-1]=ra[j-1];
                                indx[i-1]=indx[j-1];
                                i=j;
                                j <<= 1;
                        } else j=ir+1;
                }
                ra[i-1]=rra;
                indx[i-1]=rindx;
        }
}    
