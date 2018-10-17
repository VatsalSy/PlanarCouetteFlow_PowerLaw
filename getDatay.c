/* Title: Janus Drops
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
char filename[80];
scalar D2plot[], shear[];

event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  FILE * fp = ferr;
  int LEVEL = 8;
  foreach(){
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
  }
  boundary ({shear});
  for (double y = 0.; y < 1.0; y += 1./pow(2.,LEVEL)){
    fprintf(ferr, "%g %g %g %g\n", y, interpolate(u.x, L0/2., y), interpolate(shear, L0/2., y),
            interpolate(D2plot, L0/2.0, y));
  }
  fflush (fp);
  fclose (fp);
}

int main(int a, char const *arguments[])
{
  // Boundary condition: periodic right - left
  periodic (right);
  u.t[top] = neumann(0);
  u.n[top] = neumann(0);
  uf.n[top] = neumann(0);
  // bottom is moving wall
  u.n[bottom] = dirichlet(0);
  uf.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  p[top] = neumann(0);
  pf[top] = neumann(0);
  p[bottom] = neumann(0);
  pf[bottom] = neumann(0);
  sprintf (filename, "%s", arguments[1]);
  run();
}
