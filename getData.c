/* Title: Janus Drops
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"


char filename[80];
int nx,ny;
double xmin, ymin, xmax, ymax;
scalar * list = NULL;
scalar D2plot[];

event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  list = list_add (list, D2plot);
  // char nameOut[80];
  // sprintf (nameOut, "test.txt");
  // FILE * fp = fopen (nameOut, "w");
  FILE * fp = ferr;
  nx++;
  ny++;
  double Deltax = 0.999999*(xmax-xmin)/(nx - 1);
  double Deltay = 0.999999*(ymax-ymin)/(ny - 1);
  int len = list_len(list);
  double ** field = (double **) matrix_new (nx, ny, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  fprintf (fp, "# 1:x 2:y");
  int i = 3;
  for (scalar s in list){
    fprintf (fp, " %d:%s", i++, s.name);
  }
  fputc ('\n', fp);
  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
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
  xmin = atof(arguments[2]);
  xmax = atof(arguments[3]);
  ymin = atof(arguments[4]);
  ymax = atof(arguments[5]);
  nx = atoi(arguments[6]);
  ny = atoi(arguments[7]);
  list = list_add (list, u.x);
  list = list_add (list, u.y);
  run();
}
