//*********************************************************************
// "line2D_norm_constrained" finds the best fitting line of the points given 
// by x0 an x1 arrays that passes through the point (x0a,x1a).
//
// Inputs:
//    x0:      intensitis of the baseline image
//    x1:      intensitis of the follow-up image
//    w:       weights of each point on the histogram which will be updated in each iteration 
//    w_const: weights of each point on the histogram which is constant
//    n:       number of the voxels (size of arrays x0 and x1)
//    u:       the unit vector of the best fitting line
//    d:       the distance of the best fitting line from origin
//    x0a,x1a: coordinates of the constraining point that the line passes from 
void line2D_norm_constrained(double *x0, double *x1, double *w, double *w_const, int n, double *u, double &d, double x0a, double x1a)
{

   double *x0c;  // x0 centered wrt x0a
   double *x1c;  // x1 centered wrt x1a
   double a11, a22, a12;
   double L;     // the smaller eigenvalue of A

   // allocate memory for x0c and x1c
   x0c = (double *)calloc(n, sizeof(double));
   x1c = (double *)calloc(n, sizeof(double));

   // center x0 and x1 
   for(int i=0; i<n; i++)
   {
      x0c[i] = x0[i] - x0a;
      x1c[i] = x1[i] - x1a;
   }

   a11=a12=a22=0.0;
   for(int i=0; i<n; i++)
   {
      a11 += w[i]*w_const[i]*x0c[i]*x0c[i];
      a12 += w[i]*w_const[i]*x0c[i]*x1c[i];
      a22 += w[i]*w_const[i]*x1c[i]*x1c[i];
   }

   // to make the numbers more manageable 
   a11 /= n;
   a12 /= n;
   a22 /= n;

   {
      double a,b,c,D;
      double L1,L2;

      a=1;
      b=-(a11+a22);
      c=a11*a22 - a12*a12;
      D = b*b - 4*a*c;
      L1 = (-b + sqrt(D))/(2*a);
      L2 = (-b - sqrt(D))/(2*a);

      L=L2;
   }

   {
      double u1, u2;

      u1 = sqrt( (a22 - L)/(a11+a22-2*L) );
      u2 = sqrt( 1.0 - u1*u1);

      if( (int)(1000*L*u1)!= (int)(1000*a11*u1+a12*u2) || (int)(1000*L*u2)!=(int)(1000*a12*u1+a22*u2))
      {
         u2 = -u2;
      }

      d = u1*x0a + u2*x1a;

      u[0]=u1;
      u[1]=u2;
   }

   // update weights in preparation for the next iteration
   {
      double r;
      for(int i=0; i<n; i++)
      {
         r = u[0]*x0c[i] + u[1]*x1c[i];
         r = fabs(r);

         // amounts to M-esitmator with Geman-McClure residuals
         w[i] = 1.0/((1.0+r*r)*(1.0+r*r));

      }
   }

   free(x0c);
   free(x1c);
}
