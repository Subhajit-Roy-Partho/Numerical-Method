#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int cells = 200,kc=(cells/2),steps=500;
  double E[cells], B[cells],pulse=0,fh=12.0,a=1.0,t0=40.0,t=0.0,Eold[2],Bold[2],factor;
  FILE *f;
  f = fopen("Outdata.dat","w");


  for(int i=0; i<cells;i++){
    E[i]= 0.0;
    B[i]= 0.0;
  }


  for(int n=0; n < steps; n++){
    pulse = a*exp(-0.5 * (pow((t0 - t) / fh, 2.0)));
    t+=1;
    Eold[0] = E[1];
    Eold[1] = E[cells-2];
    for(int i=1; i< cells-1; i++){
      E[i]+= 0.5*(B[i-1]-B[i]);
    }

    E[kc] += pulse;
    // Mur Conditions modified (Inelastic boundary)
    E[0] = Eold[0] - (1.0/3.0)*(E[1]-E[0]);
    E[cells-1] = Eold[1] - (1.0/3.0)*(E[cells-2]-E[cells-1]);

    for(int i=0; i<cells-1;i++){
      B[i] += 0.5*(E[i]-E[i+1]);
    }

    for(int i=0;i<cells; i++){
      fprintf(f, "%i\t%f\t%f\n",i+1,E[i],B[i]);
    }
  }
  fclose(f);
  f = fopen("gnuE.txt","w");
  fprintf(f, "set term gif animate delay 3 optimize\nset output 'E.gif'\nset xrange[1:%i]\nset yrange[%lf:%lf]\n",cells,-(a+0.2),(a+0.2));
  fprintf(f, "do for [i=0:%i]{\nplot 'Outdata.dat' every ::(((i-1)*200)+1)::(i*200) u 1:2 w l\n}\nexit",steps-1);
  fclose(f);
  system("gnuplot gnuE.txt");
  f = fopen("gnuM.txt","w");
  fprintf(f, "set term gif animate delay 3 optimize\nset output 'M.gif'\nset xrange[1:%i]\nset yrange[%lf:%lf]\n",cells,-(a+0.2),(a+0.2));
  fprintf(f, "do for [i=0:%i]{\nplot 'Outdata.dat' every ::(((i-1)*200)+1)::(i*200) u 1:3 w l\n}\nexit",steps-1);
  fclose(f);
  system("gnuplot gnuM.txt");
  system("rm -rf gnuE.txt gnuM.txt");
  fflush(f);
  return 0;
}
