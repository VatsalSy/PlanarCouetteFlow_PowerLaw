/* Title: Plannar Coeutte flow
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
char filename[80];
double tau_y,mu_0,mumax;
double n;
int imax = 1e5;
#define tconv 1e-3
#define dtmax tconv/10
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  init_grid (1<<6);
  L0 = 1.0;
  origin (0.0, 0.0);
  DT = dtmax;
  stokes = true;
  TOLERANCE = 1e-5;

  /* Values of yeld stress, viscosity and coefficient.
     Newtonian: $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 1
     Power law $\mu_0 = 1.0$; $\tau_y = 0.$ and n = 0.5
     Herschel-Bulkley $\mu_0 = 1.0$; $\tau_y = 0.25$ and n = 0.5
     Bingham $\mu_0 = 1.0; $\tau_y = 0.25$ and n = 1
  */

  mu_0 = 1.0;
  tau_y= 0.0;
  n = 1.0;
  if (a >= 3){
    mu_0 = atof(arguments[2]);
  }
  if (a >= 4){
    tau_y = atof(arguments[3]);
  }
  if (a >= 5){
    n = atof(arguments[4]);
  }
  // The regularisation value of viscosity
  mumax=1000;

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

  run();
}

/* D2plot save the second principal invariant of the shear strain rate tensor
   for post-processing purposes.
*/
scalar un[], D2plot[];
face vector muv[];

// initialization event
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  /**
    presure gradient `mdpdx`
   $$-\frac{dp}{dx} = 1 $$
  */
  const face vector mdpdx[] = {1.0,0.0};
  a = mdpdx;
  /**
   Initialy at rest
  */
  foreach() {
    u.x[] = 0;
    u.y[] = 0;
  }
  foreach(){
    un[] = u.x[];
  }
  dump (file = "start");
}

/**
We look for a stationary solution. */
event logfile (t += tconv; i <= imax) {
  double du = change (u.x, un);
  fprintf(ferr, "err = %g, t = %g\n", du, t);
  if (i > 0 && du < 1e-6){
    dump (file = filename);
    return 1; /* stop */
  }
  if (i>imax-10){
    dump (file = filename);
  }
}

/**
## Implementation of generalized Newtonian viscosity
*/

event properties(i++) {
  scalar muTemp[];
  /*
   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{N-1} + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$
   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = min\left(\mu_{eq}, \mu_{max}\right) $$
  */
    foreach() {
      double D2 = 0.;
      foreach_dimension() {
          double du_x = u.x[1,0] - u.x[-1,0];
          double du_ydv_x = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
          D2 += sq(du_x) + sq(du_ydv_x);
      }
      if (D2 > 0.0) {
        // fprintf(ferr, "%g %g %g\n",t,y, D2);
        D2 = sqrt(D2)/(2.*Delta);
        D2plot[] = D2;
        double temp = tau_y/(sqrt(2.0)*D2) + mu_0*exp((n-1.0)*log(D2/sqrt(2.0)));
        muTemp[] = (temp < mumax ? temp : mumax);
      } else {
        if (tau_y > 0.0 || n < 1.0){
          muTemp[] = mumax;
        } else {
          muTemp[] = (n == 1.0 ? mu_0 : 0.0);
        }
      }
    }
    boundary ({muTemp});
    foreach_face() {
        muv.x[] = (muTemp[] + muTemp[-1,0])/2.;
    }
    boundary ((scalar *){muv});
}
