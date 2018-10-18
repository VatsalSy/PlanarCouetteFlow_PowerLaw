# PlannarCoeutteFlow_PowerLaw
The code implements viscosity for non-Newtonian fluids in Basilisk C (http://basilisk.fr/). This uses the codes and algorithm discussed at http://basilisk.fr/sandbox/M1EMN/Exemples/bingham_simple.c and extends it for Generalized Newtonian fluid. The test case is for planar Couette flow. For further documentation contact vatsalsanjay@gmail.com

In this repository, the second invariant of Velocity vector gradient is calculated at the cell center (to get scalar muTemp[] at the cell center). The coefficient of viscosity is calculated using (muTemp[] + muTemp[-1,0])/2.;
The relevant code snippet is:

```
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
```
