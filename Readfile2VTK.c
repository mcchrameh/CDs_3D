#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 128 
#define NY 128
#define NZ 128

#define IMG(cdsp, i, j ,k) (((cdsp) + (NZ*NY*(i))+NZ*j)[k])



int main(int argc, char *argv[])
{
   FILE *f;

   size_t nread;

    if ((f = fopen ("datafile_mpi", "rb")) == NULL)
    {
      fprintf (stderr, "cannot open input file \"%s\".\n", "datafile_mpi");
      return 1;
    }
  printf ("NX=%d, NY=%d, NZ=%d\n", NX, NY, NZ);
  double *cdsp = (double *) malloc (sizeof (double) * NX * NY * NZ);

   if ((nread = fread (cdsp, sizeof (double), NX * NY * NZ, f)) != NX * NY * NZ) //reads data from the file f into the array cdsp
/*  int i1,j1,k1;
   for(i1=0;i1<NX;i1++)
      {
       for(j1=0;j1<NY;j1++)
          {
           for(k1=0;k1<NZ;k1++)
              {
                fscanf(f, "%lf", &cdsp[(i1*NY+j1)*NZ + k1]);
              }
           }
       }*/
    fprintf (stderr, "Warning: only %d/%d doubles read.\n", (int)nread, NX * NY * NZ);
  fclose (f);
  
  FILE *fp = fopen("result.vtk", "w");
  int origins[3]={0,0,0};
  int spacings[3]={1,1,1};
  int float_variable=1;
  /* Write vtk Datafile header */
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "VTK\nASCII\nDATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS                             %d %d %d\n",NX, NY,NZ);
  fprintf(fp, "ORIGIN                                 %d %d %d\n", origins[0],origins[1],origins[2]);
  fprintf(fp, "SPACING                                %d %d %d\n",spacings[0],spacings[1],spacings[2]);
  fprintf(fp, "POINT_DATA                             %d\n", NX*NY*NZ);
  fprintf(fp, "SCALARS  scalars  float                %d\n", float_variable);
  fprintf(fp, "LOOKUP_TABLE default\n");
 
  int i,j,k;
  for ( i = 0; i < NX; i++)
     {
        for( j = 0; j < NY; j++)
           {
          for( k = 0; k < NZ; k++)
           {
            fprintf (fp, "%f\n ", IMG (cdsp, i, j, k));
           }

                    
           }


     }

  fclose(fp);
 free(cdsp);
return 0;


}

