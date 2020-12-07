#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int xcells = 60,ycells=60,ic=(xcells/2),jc=(ycells/2),steps=1000;
  double Ez[xcells][ycells], Hx[xcells][ycells],Hy[xcells][ycells],pulse=0,fh=12.0,a=10.0,t0=40.0,t=0.0;
  FILE *f;
  char file[15]="dat.dat";


  for(int i=0; i<xcells;i++){
    for(int j=0; j<ycells; j++){
      Ez[i][j]= 0.0;
      Hx[i][j]= 0.0;
      Hy[i][j]= 0.0;
    }
  }


  for(int n=0; n < steps; n++){
  	t +=1;

    pulse = a*exp(-0.5 * (pow((t0 - t) / fh, 2.0)));

    sprintf(file,"Data/out%d.dat",n+1);
    f = fopen(file,"w");


    for(int j=1; j<ycells-1;j++){
      for(int i=1; i <xcells-1; i++){
        Ez[i][j] += 0.5*(Hy[i][j] - Hy[i-1][j] - Hx[i][j]+Hx[i][j-1]);
      }
    }

    Ez[ic][jc] += pulse;

    for(int j=0; j< ycells-1; j++){
      for(int i=0; i< xcells; i++){
        Hx[i][j] += 0.5*(Ez[i][j]-Ez[i][j+1]);
      }
    }

    for(int j=0; j < ycells;j++){
      for(int i=0; i< xcells-1; i++){
        Hy[i][j] += 0.5*(Ez[i+1][j]-Ez[i][j]);
      }
    }
    for(int j=0; j<ycells;j++){
      for(int i=0; i < xcells; i++){
        fprintf(f, "%lf ",Ez[i][j]);
      }
      fprintf(f, "\n");
    }
    fflush(f);
    fclose(f);
  }
  f = fopen("gnuE.txt","w");
  fprintf(f, "set pm3d\nset cbrange[-1:1] \nset isosample 11\nset sample 11\nset term gif large animate delay 3 optimize\nset output 'out.gif'\nset zrange [1:-1]\nset print '-'\ndo for [i=1:%d]{\nsplot 'Data/out'.i.'.dat' matrix w l\n print i\n}\nexit",steps);
  fflush(f);
  fclose(f);
  system("gnuplot gnuE.txt");
  return 0;
}
