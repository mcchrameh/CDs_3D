
#include "CahnH3D.h"
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <assert.h>
#include "mpi.h"
#include "omp.h"
#define N_Threads 4
#define Random_min  -0.05 //-0.25 // for mim random number generator
#define Random_max  0.05 //0.25
#define PHI(i,j,k)         PHI[i][j][k]
#define PHI_old(i,j,k)     PHI_old[i][j][k]
#define gamma(i,j,k)       gamma[i][j][k]
#define Laplacian2(i,j,k)   Laplacian2[i][j][k]
#define gamma_p(i,j,k)     gamma_p[i][j][k]
#define PHI_old_p(i,j,k)   PHI_old_p[i][j][k]
#define Laplacian2_p(i,j,k) Laplacian2_p[i][j][k]
#define PHI_p(i,j,k)       PHI_p[i][j][k]
#define Nzdir 100
#define Nydir 100
#define MAX_LINE_LENGTH  80 //used for reading the variables line by line in the input file
#define Q(i,j,k)  Q[((k) + (Nzdir) * ((j) + (Nydir) * (i)))] //use to write to file
//#define Q1(i,j,k)  Q1[((k) + (Nzdir) * ((j) + (Nydir) * (i)))] 
using namespace std;



CahnHill3D::CahnHill3D()
      { 
        nx =0;
        ny=0;
         
              
        }
CahnHill3D::CahnHill3D(int Nx, int Ny, int Nz)
        {
                                      
      	  nx=Nx;ny=Ny;nz=Nz;delta_x=1.0;delta_y=1.0; delta_t =0.001;M =1.0;b=M;k=M, u=0.5;
          C1=6.0/80.0; C2=3.0/80.0; C3=1.0/80.0; //C1=1.0/6, C2=1.0/12
          dxx = delta_x;
          B=0.005; D=0.5;A=1.3;v=1.5;tau=0.36;F=0.45;r=0.5;
          x =new double[nx];
          y= new double[ny];
          z=new double[nx];
                                                   
          PHI=new double**[nx];
          PHI_old = new double**[nx];
          gamma = new double**[nx];
          Laplacian2 = new double**[nx];
          Q =new double[nx*ny*nz];
      for(int i=0;i<nx; i++)
             { 
	             PHI[i]=new double*[ny];
	             PHI_old[i]=new double*[ny];
               gamma[i]=new double*[ny]; 
               Laplacian2[i] = new double*[ny];
            for(int j=0;j<ny;j++)
		           {
		            PHI[i][j]=    new double[nz];                                                                                                                                      
		            PHI_old[i][j]=new double[nz]; 
                gamma[i][j]=  new double[nz];                                                                                                                                    
		            Laplacian2[i][j] = new double[nz];

               } 
             }
                                                                                                                                                                                                                                                                                                         
           for(int i=0;i<nx;i++)
              { 
               for(int j=0;j<ny;j++)
                  {   
                   for(int k=0;k<nz;k++)
                     { 
                       PHI[i][j][k]=0.0;
                       PHI_old[i][j][k]=0.0;
                       gamma[i][j][k]=0.0;
                       Laplacian2[i][j][k]=0.0;
                     }
                   }
               }
          for(int i=0;i<nx;i++)
             {
              x[i]=i;
             }

          for(int j=0;j<ny;j++)
             {
               y[j]=j;

              }       
          for(int k=0;k<nz;k++)
             {
               z[k]=k;
             }                                                                                                                                                                                                                                                                                                 
          down1x = new int [nx]; memset(down1x, 0, nx*sizeof(int));
          down1x_2 = new int [nx+2]; memset(down1x_2, 0, (nx+2)*sizeof(int));
          down2x = new int [nx]; memset(down2x, 0, nx*sizeof(int));
          up1x   = new int [nx]; memset(up1x, 0, nx*sizeof(int));
          up1x_2   = new int [nx+2]; memset(up1x_2, 0, (nx+2)*sizeof(int));

          memset(Q, 0, nx*ny*nz*sizeof(double));

          up2x   = new int [nx]; memset(up2x, 0, nx*sizeof(int));
          cout<<"Constructor called"<<endl;
   }

CahnHill3D::~CahnHill3D()
   {
      for(int i=0;i<nx;i++)
        {
         for(int j=0;j<ny;j++)
          {
            delete[] PHI[i][j];
            delete[] PHI_old[i][j];
            delete[] gamma[i][j];
            delete[] Laplacian2[i][j];
          }
            delete[] PHI[i];
            delete[] PHI_old[i];
            delete[] gamma[i];
            delete[] Laplacian2[i];
       }
       delete[] PHI;
       delete[] PHI_old;
       delete[] gamma;
       delete[] Laplacian2;
       delete[] Q;
       delete[] x;
       delete[] y;
       delete[] z;
       delete[] down1x;
       delete[] down1x_2;
       delete[] up1x;
       delete[] up1x_2;

    std::cout<<"Destructor called"<<std::endl;
}



void CahnHill3D:: setIndex(int start, int end)
    {
        for(int s=start; s<end; s++)
           {
                up1x[s]=s+1;
                up2x[s]=s+2;
                down1x[s]=s-1;
                down2x[s]=s-2;
           }



            down1x[start]=end-1;
            down2x[start]=end-2;
            down2x[start+1]=end-1;
            up1x[end-1]=start;
            up2x[end-1]=1;
            up2x[end-2]=0;




    }
            
void CahnHill3D::setIndex2(int start,int end)
   {

      for(int s=start; s<end; s++)
           {    
                up1x_2[s]=s+1;
                
                down1x_2[s]=s-1;
                
           }


            
            down1x_2[start]=end-1;
            //down2x[start]=nx-2;
           // down2x[start+1]=nx-1;
            up1x_2[end-1]=start;
          //  up2x[end-1]=1;
          //  up2x[end-2]=0;



   } 

void CahnHill3D::initialCondition()
{

      for(int i=0; i<nx;i++)
         {
          for(int j=0;j<ny;j++)
            {
             for(int k=0;k<nz;k++)
               {
                   double range =Random_max-Random_min;
                   double div =RAND_MAX / range;
                   PHI_old[i][j][k]=Random_min + (rand()/div);
                     
                }
              }

            }


}
 

double CahnHill3D::g(double phi)
{
    double q=0.0;
    q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    // q = A*phi- (A/3.0)*pow(phi,3); // map from paper
    return q;
}    

void CahnHill3D::setLaplacianBase()
 {
    double AP=0.0, BP=0.0,  ATP=0.0, CP=0.0 ;//AP->nearest neighbor, BP->next nearest neighbor, CP->next next nearest neighbor;
    setIndex(0,nx);
    for(int i=0;i<nx;i++)
         {
          for(int j=0;j<ny;j++)
            {
             for(int k=0;k<nz;k++)
              {
                AP = C1*(PHI_old(up1x[i],j,k) + PHI_old(down1x[i],j,k) 
                       + PHI_old(i,up1x[j],k) + PHI_old(i,down1x[j],k) 
                       + PHI_old(i,j,up1x[k]) + PHI_old(i,j,down1x[k]) );
                      /*if((i==0)&&(j==0)&&(k==0) )
                          {
                             printf("PHI_old(%d,%d,%d)=%lf\n",up1x[i],j,k,PHI_old(up1x[i],j,k));
                             printf("PHI_old(%d,%d,%d)=%lf\n",down1x[i],j,k,PHI_old(down1x[i],j,k));
                             printf("PHI_old(%d,%d,%d)=%lf\n",i,up1x[j],k,PHI_old(i,up1x[j],k));
                             printf("PHI_old(%d,%d,%d)=%lf\n",i,down1x[j],k, PHI_old(i,down1x[j],k));
                             printf("PHI_old(%d,%d,%d)=%lf\n",i,j,up1x[k],  PHI_old(i,j,up1x[k]));
                             printf("PHI_old(%d,%d,%d)=%lf\n",i,j,down1x[k],  PHI_old(i,j,down1x[k]) );
                          }*/
                BP =C2*(PHI_old(down1x[i],up1x[j],k) +   PHI_old(down1x[i],down1x[j],k) 
                      + PHI_old(up1x[i],up1x[j],k)   +   PHI_old(up1x[i],down1x[j],k) 
                      + PHI_old(i,down1x[j],up1x[k]) +   PHI_old(i,down1x[j],down1x[k]) 
                      + PHI_old(i,up1x[j],up1x[k])   +   PHI_old(i,up1x[j],down1x[k]) 
                      + PHI_old(down1x[i],j,up1x[k]) +   PHI_old(down1x[i],j,down1x[k]) 
                      + PHI_old(up1x[i],j,up1x[k])   +   PHI_old(up1x[i],j,down1x[k])) ;
 
                CP=C3*( PHI_old(down1x[i],down1x[j],down1x[k]) + PHI_old(down1x[i],up1x[j],down1x[k])
                      + PHI_old(down1x[i],down1x[j],up1x[k]) + PHI_old(down1x[i],up1x[j],up1x[k])
                      + PHI_old(up1x[i],down1x[j],down1x[k]) + PHI_old(up1x[i],up1x[j], down1x[k])
                      + PHI_old(up1x[i],down1x[j],up1x[k])   + PHI_old(up1x[i],up1x[j],up1x[k]));
                 ATP = AP + BP + CP;
                 gamma(i,j,k)=g(PHI_old(i,j,k))-PHI_old(i,j,k) + D*(ATP-PHI_old(i,j,k));
              //   printf("gamma(%d,%d,%d)=%lf\n",i,j,k,gamma(i,j,k));

                  
                 }
              }
            }


  }




void CahnHill3D::SetSecondLaplacian2()
{   
  double AP=0.0, BP=0.0, CP=0.0 ;
  setIndex(0,nx);
  for(int i=0;i<nx;i++)
         {
          for(int j=0;j<ny;j++)
            {
             for(int k=0;k<nz;k++)
              {
                AP = C1*(gamma(up1x[i],j,k) + gamma(down1x[i],j,k) 
                      +  gamma(i,up1x[j],k) + gamma(i,down1x[j],k)
                      +  gamma(i,j,up1x[k]) + gamma(i,j,down1x[k]));
 
                BP = C2*(gamma(down1x[i],up1x[j],k) + gamma(down1x[i],down1x[j],k)
                      +  gamma(up1x[i],up1x[j],k)   + gamma(up1x[i],down1x[j],k)
                      +  gamma(i,down1x[j],up1x[k]) + gamma(i,down1x[j],down1x[k])
                      +  gamma(i,up1x[j],up1x[k])   + gamma(i,up1x[j],down1x[k])
                      +  gamma(down1x[i],j,up1x[k]) + gamma(down1x[i],j,down1x[k])
                      +  gamma(up1x[i],j,up1x[k])   + gamma(up1x[i],j,down1x[k]));
          
                CP = C3*(gamma(down1x[i],down1x[j],down1x[k]) + gamma(down1x[i],up1x[j],down1x[k])
                      +  gamma(down1x[i],down1x[j],up1x[k])   + gamma(down1x[i],up1x[j],up1x[k])
                      +  gamma(up1x[i],down1x[j],down1x[k])   + gamma(up1x[i],up1x[j],down1x[k])
                      +  gamma(up1x[i],down1x[j],up1x[k])     + gamma(up1x[i],up1x[j],up1x[k]));

                Laplacian2(i,j,k)=AP + BP + CP;
              }
            }
          }
 
}






void CahnHill3D::FiniteDifferenceScheme()
{
  cout<<"in finite difference scheme"<<endl;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
           for(int k=0;k<nz;k++)
             {

               PHI(i,j,k)= PHI_old(i,j,k)-(Laplacian2(i,j,k)-gamma(i,j,k) + B*PHI_old(i,j,k));
             }
        }
    }


} 

void CahnHill3D::UpdateSolution()
{

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
         for(int k=0;k<nz;k++)
            {
              PHI_old[i][j][k]=PHI[i][j][k];
            }
        }

    }


}

void CahnHill3D::result(int count)
{

        FILE *file = fopen(make_output_filename(count).c_str(), "w");
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
          for(int k=0;k<nz;k++)
            {
              //Q(i,j,k)=PHI[i][j][k];
              fprintf(file,"%lf  \t",  PHI(i,j,k));
             }
             // fprintf(file,"\n");
         }
//      fprintf(file,"\n");
     }
          fclose(file);
     save_vtk(PHI, nx, ny, nz);

           
            
 }


void CahnHill3D:: save_vtk(double ***Q1, int Nx, int Ny, int Nz)
 {
    
  FILE *fp = fopen("result.vtk", "w");
  int origins[3]={0,0,0};
  int spacings[3]={1,1,1};
  int float_variable=1;
  /* Write vtk Datafile header */
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "VTK\nASCII\nDATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS                             %d %d %d\n",Nx, Ny,Nz);
  fprintf(fp, "ORIGIN                                 %d %d %d\n", origins[0],origins[1],origins[2]);
  fprintf(fp, "SPACING                                %d %d %d\n",spacings[0],spacings[1],spacings[2]);
  fprintf(fp, "POINT_DATA                             %d\n", Nx*Ny*Nz);
  fprintf(fp, "SCALARS  scalars  float                %d\n", float_variable);
  fprintf(fp, "LOOKUP_TABLE default\n");


 
  for(int i = 0; i < Nx; i++)
    { 
      for(int j= 0; j < Ny; j++)
       {
        for(int k=0;k<Nz;k++)
          {
             fprintf(fp, " %lf", Q1[i][j][k]);
          
          }
  
        }     
    
    }
  fclose(fp);     



 }
         
void CahnHill3D::Solver2()
{

   initialCondition();
    int count =0;
    double t=0.0;
    while(t<102)
    {
        std::cout<<"in while loop solver"<<std::endl;
        setLaplacianBase();
        SetSecondLaplacian2();
        FiniteDifferenceScheme();
        UpdateSolution();
        if (count==100|| count==0)
            result( count);

        t+=1.0;
        count++;
        cout<<"time="<<t<<"\t";

   }



    cout<<""<<endl;



}
void CahnHill3D::ReadFilein(std::ifstream& infile)
  {

    
          //  double *u1;
           //   u1= new double[Nx*Ny*Nz];
              //std::memset(u1, 0, Nx*Ny*Nz*sizeof(double));
            /* u1 = new double**[Nx];
             for(int i=0;i<Nx; i++)
             {
               u1[i]=new double*[Ny];
               for(int j=0;j<Ny;j++)
                  {
                    u1[i][j]  =  new double[Nz];  
                  }
             }
                */ 
          
      
              for(int i=0;i<nx;i++)
                 {
                  for(int j=0;j<ny;j++)
                     {
                      for(int k=0;k<nz;k++)
                        {
                         // infile>>u1[(i*Ny + j)*Nz +k];
                          infile>>PHI_old(i,j,k);
                        }
                      }

                    }
            
    
    
  }
                   

Parallel_CahnHill3D::Parallel_CahnHill3D()//: CahnHill3D(N1x, N1y, N1z)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    C1=6.0/80.0; C2=3.0/80.0; C3=1.0/80.0;
    tag =1;
    double rationals[10];
    int integers[6];
    if(rank==0)
    {
     Read_input_parameters(integers, rationals);
    }
    MPI_Bcast(integers,6,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(rationals,10,MPI_DOUBLE,0,MPI_COMM_WORLD);
     
    Nx = integers[0]; Ny=integers[1]; Nz=integers[2];
    Procx=integers[3]; Procy=integers[4]; Procz =integers[5];
    M = rationals[0]; u=rationals[1]; B = rationals[2]; 
    D = rationals[3]; A=rationals[4]; v = rationals[5];
    tau= rationals[6]; F=rationals[7]; r=rationals[8];
    Max_time_iteration=rationals[9];
      
    int coords[3], dims[3], periods[3];
    MPI_Comm new_comm =CreatCartesianTopology();
    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
    MPI_Cart_get (new_comm,3,dims,periods, coords );
    assert(Procx*Procy*Procz==size);
    printf("dimensions of topology=(%d,%d,%d)\n",dims[0],dims[1],dims[2]);
    if(Nx%Procx!=0)
     {
       if(coords[0]<(dims[0]-1))
         {
           //nlocalx = (int)Nx/std::cbrt(size);
           nlocalx = (int)Nx/Procx;
       //    printf("rank=%d, nlocal=%d\n",rank,nlocalx);
         }
       else
         {

         //  nlocalx = (int)Nx/std::cbrt(size);
           nlocalx = (int)Nx/Procx + (int)Nx%Procx;
      //     startrow_x=coords[0]*nlocalx;
        //   endrow_x=Nx-1;
         //  nlocalx  = endrow_x-startrow_x + 1;
          printf("upper rank=%d, nlocal=%d, coords=(%d,%d,%d)\n",rank,nlocalx, coords[0],coords[1],coords[2]);
         }
     }
   else
    {
       nlocalx = (int)Nx/Procx;
    }

 //divide along y

   if(Ny%Procy!=0)
     {
       if(coords[1]<(dims[1]-1))
         {
     
           nlocaly = (int)Ny/Procy;

     //      printf("rank=%d, nlocaly=%d\n",rank,nlocaly);
         }
       else
         {

           nlocaly = (int)Ny/Procy + (int)Ny%Procy;
        }
    }
   else
    {
     
       nlocaly = (int)Ny/Procy;
    }



  //divide along z
  

    if(Nz%Procz!=0)             
       {                        
          if(coords[2]<(dims[2]-1))        
            {                    
             
             nlocalz = (int)Nz/Procz;
            }
                              
          else                   
           {
             nlocalz = (int)Nz/Procz + (int)Nz%Procz;                    
     
           }
     }            
   else               
    {                 
      nlocalz = (int)Nz/Procz;
    }

 printf("rank=%d, nlocal=(%d,%d,%d)\n",rank,nlocalx,nlocaly,nlocalz);

 /*-----------------------------------------------------------------
 * Create local matrices on each process
 ----------------------------------------------------------------------*/
  matrix_right_oxz= new double [nlocalx*nlocalz];
  matrix_left_oxz = new double [nlocalx*nlocalz];
  matrix_outer_oyx= new double [nlocalx*nlocaly];
  matrix_inner_oyx= new double [nlocalx*nlocaly];
  matrix_north_oyz= new double [nlocaly*nlocalz];
  matrix_south_oyz= new double [nlocaly*nlocalz];
  array_left_oxz  = new double [nlocalz];
  array_right_oxz = new double [nlocalz];
  array_north_oyz = new double [nlocaly];
  array_south_oyz = new double [nlocaly];
  array_inner_oxy = new double [nlocaly];
  array_outer_oxy = new double [nlocaly];
  array_x_right   = new double [nlocalx];
  array_x_left    = new double [nlocalx];
   down1x = new int [Ny]; memset(down1x, 0, Ny*sizeof(int));
   up1x   = new int [Ny]; memset(up1x, 0, Ny*sizeof(int));
 
//  value           = new double [4];
//  value1          = new double [4];
 
  PHI_p       = new double**[nlocalx + 2];
  PHI_old_p   = new double**[nlocalx + 2];
  gamma_p     = new double**[nlocalx + 2]; 
  Laplacian2_p= new double**[nlocalx + 2];
  PHI_local_result  = new double [nlocalx*nlocaly*nlocalz]; 
  PHI_global_result = new double [Nx*Ny*Nz]; 
          for(int i=0;i<nlocalx+2; i++)
             {
               PHI_p[i]       = new double*[nlocaly+2];
               PHI_old_p[i]   = new double*[nlocaly+2];
               gamma_p[i]     = new double*[nlocaly+2];
               Laplacian2_p[i]= new double*[nlocaly+2];
               for(int j=0;j<nlocaly+2;j++)
                  {
                   PHI_p[i][j]  =  new double[nlocalz+2];                                                                                                    
                   PHI_old_p[i][j]=  new double[nlocalz+2];
                   gamma_p[i][j]=  new double[nlocalz+2];
                   Laplacian2_p[i][j] = new double[nlocalz+2];
                 }
             }


           for(int i=0;i<nlocalx+2;i++)
              {
               for(int j=0;j<nlocaly+2;j++)
                  {
                   for(int k=0;k<nlocalz+2;k++)
                     {
                       PHI_p[i][j][k]=0.0;
                       PHI_old_p[i][j][k]=0.0;
                       gamma_p[i][j][k]=0.0;
                       Laplacian2_p[i][j][k]=0.0;
                     }
                   }
               }


 



  FreeCommunicator(new_comm);

}


Parallel_CahnHill3D::~Parallel_CahnHill3D()
   {
      for(int i=0;i<nlocalx+2;i++)
        {
         for(int j=0;j<nlocaly+2;j++)
          {
            delete[] PHI_p[i][j];
            delete[] PHI_old_p[i][j];
            delete[] gamma_p[i][j];
            delete[] Laplacian2_p[i][j];
          }
            delete[] PHI_p[i];
            delete[] PHI_old_p[i];
            delete[] gamma_p[i];
            delete[] Laplacian2_p[i];
       }
       delete[] PHI_p;
       delete[] PHI_old_p;
       delete[] gamma_p;
       delete[] Laplacian2_p;
       delete[] matrix_right_oxz;     
       delete[] matrix_left_oxz;
       delete[] matrix_south_oyz;
       delete[] matrix_north_oyz;
       delete[] matrix_inner_oyx;
       delete[] matrix_outer_oyx;
       delete[] array_left_oxz;
       delete[] array_right_oxz;
       delete[] array_north_oyz;
       delete[] array_south_oyz;
       delete[] array_inner_oxy;
       delete[] array_outer_oxy;
       delete[] array_x_left;
       delete[] array_x_right;
 //      delete[] value;
 //      delete[] value1; 
       delete[] PHI_local_result;
       delete[] PHI_global_result;
       delete[] down1x;
       delete[] up1x;
//        FreeCommunicator(new_comm);

    std::cout<<"Destructor called"<<std::endl;
}
double Parallel_CahnHill3D::g(double phi)
{
    double q=0.0;
    q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    // q = A*phi- (A/3.0)*pow(phi,3); // map from paper
    return q;
}    


void Parallel_CahnHill3D::Initialize_parallel(MPI_Comm new_comm)
{

      MPI_Comm_rank(new_comm, &rank);
      MPI_Comm_size(new_comm, &size);

       srand(time(NULL) + rank);
       printf("Nlocals=(%d,%d,%d)\n",nlocalx,nlocaly,nlocalz);
       std::cout<<"I am rank="<<rank<<std::endl;
       for(int i=1; i<=nlocalx ;i++)
          {
           for(int j=1;j<=nlocaly;j++)
              {
               for(int k=1;k<=nlocalz;k++)
                 {
                    double range =Random_max-Random_min;
                    double div =RAND_MAX / range;
                    PHI_old_p[i][j][k]=Random_min + (rand()/div);
                   // if(rank==2)
                     // {
                      // printf("Rank=%d, Phi_old[%d][%d][%d]=%lf\t",rank, i,j,k,PHI_old_p[i][j][k]);
                     // }
                  }
                }
            }

 std::cout<<" intialize data structure"<<std::endl;

} 

void Parallel_CahnHill3D::Read_input_parameters(int *integers, double *rationals)
  {
    FILE* file;
    char Data[MAX_LINE_LENGTH];
    if((file=fopen("ParameterFile.dat","r"))==NULL)
    {
      printf("Error opening ParameterFile.dat\n");
      return;
    }
    
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[0]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[1]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[2]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[3]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[4]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%d\n",&integers[5]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[0]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[1]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[2]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[3]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[4]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[5]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[6]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[7]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[8]);
    fgets(Data,MAX_LINE_LENGTH,file);
    fscanf(file,"%lf\n",&rationals[9]);
            
    //fgets(Data,MAX_LINE_LENGTH,file);
    //fscanf(file,"%le\n",&conf[0]);
   // fgets(Data,MAX_LINE_LENGTH,file);
    //fscanf(file,"%le",&conf[1]);

  fclose(file);
    
    
} 

MPI_Comm Parallel_CahnHill3D:: CreatCartesianTopology()
    {
     MPI_Comm   new_comm;
     if(Procx*Procy*Procz!=size)
       {
         printf("Number of processors is not factorizable. That is Px*Py*Pz=%d, is not equal to the number of processors=%d\n ",Procx*Procy*Procz,size);
         exit(0);
      
       }
       
      int  ROW=0,COL=1, DEPTH=2;

       int periods[3];
        periods[0]=1;periods[1]=1, periods[2]=1;
       int coords[3];
      int dims[3];
      dims[ROW] = Procx;//std::cbrt(size);
      dims[COL] = Procy;//std::cbrt(size); 
      dims[DEPTH]=Procz;//std::cbrt(size);
     MPI_Comm_size( MPI_COMM_WORLD, &size );
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     MPI_Cart_create( MPI_COMM_WORLD, 3,dims, periods, 0, &new_comm );
     MPI_Comm_rank(new_comm,&my3drank);
     MPI_Cart_coords(new_comm,my3drank,3,coords);

     MPI_Comm_size( new_comm, &size );
     //printf("I am rank=%d, coords=(%d,%d,%d)\n", my3drank,coords[0],coords[1],coords[2]);
    return new_comm;
   }

void Parallel_CahnHill3D::FreeCommunicator(MPI_Comm new_comm)
     {
       MPI_Comm_free(&new_comm);

     }


void Parallel_CahnHill3D::FiniteDifferenceScheme_p(MPI_Comm new_comm)
    {
      // std:: cout<<"in finite difference scheme parallel"<<std::endl;
      MPI_Comm_rank(new_comm, &rank);
      MPI_Comm_size(new_comm, &size);

        for(int i=1;i<=nlocalx;i++)
          {
           for(int j=1;j<=nlocaly;j++)
             {
               for(int k=1;k<=nlocalz;k++)
                 {

                   PHI_p(i,j,k)= PHI_old_p(i,j,k)-(Laplacian2_p(i,j,k)-gamma_p(i,j,k) + B*PHI_old_p(i,j,k));
                   PHI_local_result[((i-1)*nlocaly+(j-1))*nlocalz + k-1]=PHI_p(i,j,k);
                    // PHI(i,j,k)= PHI_old(i,j,k)-(Laplacian2(i,j,k)-gamma(i,j,k) + B*PHI_old(i,j,k));
                 }
             }
          }

       

     }
                                                                                                                                                                                                                                                                                                                                                                        
void Parallel_CahnHill3D:: UpdateSolution_p(MPI_Comm new_comm)
   {

       
      MPI_Comm_rank(new_comm, &rank);
      MPI_Comm_size(new_comm, &size);

      for(int i=1;i<=nlocalx;i++)
        {
         for(int j=1;j<=nlocaly;j++)
           {
            for(int k=1;k<=nlocalz;k++)
              {
                PHI_old_p[i][j][k]=PHI_p[i][j][k];
               }
        }

    }



   }



 void Parallel_CahnHill3D::parallel_solver()
    {
      // Initialize_parallel();
       MPI_Comm new_comm;
       new_comm = CreatCartesianTopology();
      // Initialize_parallel(new_comm);
      // ReadFile_MPI(new_comm);
       int count =0;
       double t =0.0;
       printf("In parallel solver with Maximum iteration=%lf\n",Max_time_iteration);
        while(t<Max_time_iteration)
         {
           ExchangeData_Second( new_comm, PHI_old_p);
           ComputeLaplacianBase_p(new_comm);
           ExchangeData_Second(new_comm,gamma_p);
           setSecond_laplacian(new_comm);
    //       ExchangeData(new_comm,Laplacian2_p);
           FiniteDifferenceScheme_p(new_comm);
           UpdateSolution_p(new_comm);
            if(t==999)
               {
                  WriteToFile(new_comm);
                  WriteToFile_MPI(new_comm);
               }
             t+=1.0;
             printf("time=%lf\n",t);

              count++;
           }
         
       //   WriteToFile_MPI(new_comm);
         FreeCommunicator(new_comm);

     }
 
void Parallel_CahnHill3D::ExchangeData_second(MPI_Comm new_comm, double ***array)
{

       MPI_Status status;
       int coords[3];
       MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);
       int dims[3], periods[3];

       MPI_Cart_get (new_comm,3,dims,periods, coords );

       //-------Exchange along Z direction----------------/
//MPI_Wait(&request,&status);
      MPI_Cart_get (new_comm,3,dims,periods, coords );
     if(dims[2]>1)
       { 
 
       for(int i=1;i<=nlocalx;i++)
         {
          for(int j=1;j<=nlocaly;j++)
             {
               matrix_inner_oyx[(i-1)*nlocaly +(j-1)]=array[i][j][nlocalz];
               matrix_outer_oyx[(i-1)*nlocaly +(j-1)]=array[i][j][1];
             }
          }
      MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);
  //   if((coords[2]>1)&&(coords[2]<dims[2]-1))
  //    {
       MPI_Cart_shift( new_comm, 2, 1, &outer_nbr, &inner_nbr );
     
     
       MPI_Sendrecv_replace(matrix_inner_oyx, nlocalx*nlocaly,  MPI_DOUBLE, inner_nbr, tag,
                             outer_nbr,tag, new_comm,&status);
       MPI_Sendrecv_replace(matrix_outer_oyx, nlocalx*nlocaly,  MPI_DOUBLE, outer_nbr, tag,
                            inner_nbr,tag, new_comm,&status);

    
       //fill ghost points
       for(int i=1;i<=nlocalx;i++)
        {
            for(int j=1;j<=nlocaly;j++)
            {
                array[i][j][0]=matrix_inner_oyx[(i-1)*nlocaly +(j-1)];
                array[i][j][nlocalz+1]=matrix_outer_oyx[(i-1)*nlocaly +(j-1)];
               // if(my3drank==2)
              //  printf("check transfer alomg z: matrix[%d]=%lf\n",(i-1)*nlocaly +(j-1),matrix_inner_oyx[(i-1)*nlocaly +(j-1)]);              
                  
            }
        }
      }
//---------End Exchange along Z----------/

//------Exchange along Y -------------------/

     MPI_Cart_get (new_comm,3,dims,periods, coords );
  if(dims[1]>1)
   {

     
       for(int i=1;i<=nlocalx;i++) 
          {
          
            for(int k=0;k<=nlocalz+1;k++) 
               {
                 matrix_right_oxz[(i-1)*(nlocalz+1) +(k)]=array[i][nlocaly][k];
               }
          }

        for(int i=1;i<=nlocalx;i++)
          {

            for(int k=0;k<=nlocalz+1;k++)
               {
                 matrix_left_oxz[(i-1)*(nlocalz+1) +(k)]=array[i][1][k];
               }
          }
  
     //---------ME send last values to right neighbor and recieve last values from left neighbor along the y axis---------
        MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);
       
        MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr ); //move along the y axis
       
       
       MPI_Sendrecv_replace(matrix_right_oxz, (nlocalz+1)*nlocalx,  MPI_DOUBLE, right_nbr, tag,
                            left_nbr,tag, new_comm,&status);
       MPI_Sendrecv_replace(matrix_left_oxz, (nlocalz+1)*nlocalx,  MPI_DOUBLE, left_nbr, tag,
                            right_nbr,tag, new_comm,&status);
       
    
       
     //---------ME send first values to left neighbor and recieve first values from right neighbor along the y axis---------
         
     
       
    //----Copy recieved data into ghost points ---------/ 
      for(int i=1;i<=nlocalx;i++)
        {
         for(int k=0;k<=nlocalz+1;k++)
           {
             array[i][nlocaly+1][k]=matrix_left_oxz[(i-1)*(nlocalz+1) +(k)];
            array[i][0][k]=matrix_right_oxz[(i-1)*(nlocalz+1) +(k)];
          }
         }

 }

   //------------Exchange along X direction-------------------------------------------------------------/
 if(dims[0]>1)
 { 
  for(int j=0;j<=nlocaly+1;j++)
          {
            
            for(int k=0;k<=nlocalz+1;k++)
               { 
                 matrix_north_oyz[(j)*(nlocalz+1) +(k)]=array[1][j][k];
               }
          }
   
  for(int j=0;j<=nlocaly+1;j++)
          {
 
            for(int k=0;k<=nlocalz+1;k++)
               {
                 matrix_south_oyz[(j)*(nlocalz+1) +(k)]=array[nlocalx][j][k];
               //    if(my3drank==size-1)
                //   printf("array[%d][%d][%d]=%lf\n",nlocalx,j,k,array[nlocalx][j][k]);
               }
          }

     MPI_Comm_rank(new_comm,&my3drank);
     MPI_Cart_coords(new_comm,my3drank,3,coords);
 //  if((coords[0]>1)&&(coords[0]<dims[0]-1))
 //  {
     MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr );

 /*---------ME send first values to north  neighbor and recieve firt values from south neighbor along the x axis---------*/

  
   
   
   MPI_Sendrecv_replace(matrix_north_oyz, (nlocaly+2)*(nlocalz+2),  MPI_DOUBLE, up_nbr, tag,
                       down_nbr,tag, new_comm,&status);
     
   MPI_Sendrecv_replace(matrix_south_oyz, (nlocaly+2)*(nlocalz+2),  MPI_DOUBLE, down_nbr, tag,
                        up_nbr,tag, new_comm,&status);
  

 /*------Copy values into ghost cells---------------*/
  
  for(int j=0;j<=nlocaly+1;j++)                                
        {                                                       
         for(int k=0;k<=nlocalz+1;k++)                            
           {                                                    
            array[nlocalx+1][j][k]=matrix_north_oyz[(j)*(nlocalz+1) +(k)];
            array[0][j][k]=matrix_south_oyz[(j)*(nlocalz+1) +(k)];
                                                                
           }                                                    
         } 


    
  }  
  MPI_Barrier(new_comm);
  //----------------END along X----------------------------------------------------------------




}



void Parallel_CahnHill3D::ExchangeData(MPI_Comm new_comm, double ***array)
{
       MPI_Status status;
     
       int coords[3];

       MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);

     int dims[3], periods[3];

     MPI_Cart_get (new_comm,3,dims,periods, coords );
  if(dims[1]>1)
   {

     
       for(int i=1;i<=nlocalx;i++) 
          {
          
            for(int k=1;k<=nlocalz;k++) 
               {
                 matrix_right_oxz[(i-1)*nlocalz +(k-1)]=array[i][nlocaly][k];
               }
          }

        for(int i=1;i<=nlocalx;i++)
          {

            for(int k=1;k<=nlocalz;k++)
               {
                 matrix_left_oxz[(i-1)*nlocalz +(k-1)]=array[i][1][k];
               }
          }
  
     //---------ME send last values to right neighbor and recieve last values from left neighbor along the y axis---------
        MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);
       
        MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr ); //move along the y axis
       
       
       MPI_Sendrecv_replace(matrix_right_oxz, nlocalz*nlocalx,  MPI_DOUBLE, right_nbr, tag,
                            left_nbr,tag, new_comm,&status);
       MPI_Sendrecv_replace(matrix_left_oxz, nlocalz*nlocalx,  MPI_DOUBLE, left_nbr, tag,
                            right_nbr,tag, new_comm,&status);
       
    
       
     //---------ME send first values to left neighbor and recieve first values from right neighbor along the y axis---------
         
     
       
    //----Copy recieved data into ghost points ---------/ 
      for(int i=1;i<=nlocalx;i++)
        {
         for(int k=1;k<=nlocalz;k++)
           {
             array[i][nlocaly+1][k]=matrix_left_oxz[(i-1)*nlocalz +(k-1)];
            array[i][0][k]=matrix_right_oxz[(i-1)*nlocalz +(k-1)];
          }
         }

 }
    
   /*-------Exchange along Z direction----------------*/
//MPI_Wait(&request,&status);
   MPI_Cart_get (new_comm,3,dims,periods, coords );
  if(dims[2]>1)
   { 
 
    

     
     
         for(int i=1;i<=nlocalx;i++)
         {
             
             for(int j=1;j<=nlocaly;j++)
             {
                 matrix_inner_oyx[(i-1)*nlocaly +(j-1)]=array[i][j][nlocalz];
                 matrix_outer_oyx[(i-1)*nlocaly +(j-1)]=array[i][j][1];
               
             }
          
         }

       MPI_Comm_rank(new_comm,&my3drank);
       MPI_Cart_coords(new_comm,my3drank,3,coords);
  //   if((coords[2]>1)&&(coords[2]<dims[2]-1))
  //    {
       MPI_Cart_shift( new_comm, 2, 1, &outer_nbr, &inner_nbr );
     
     
        MPI_Sendrecv_replace(matrix_inner_oyx, nlocalx*nlocaly,  MPI_DOUBLE, inner_nbr, tag,
                             outer_nbr,tag, new_comm,&status);
       MPI_Sendrecv_replace(matrix_outer_oyx, nlocalx*nlocaly,  MPI_DOUBLE, outer_nbr, tag,
                            inner_nbr,tag, new_comm,&status);

    
       //fill ghost points
       for(int i=1;i<=nlocalx;i++)
        {
            for(int j=1;j<=nlocaly;j++)
            {
                array[i][j][0]=matrix_inner_oyx[(i-1)*nlocaly +(j-1)];
                array[i][j][nlocalz+1]=matrix_outer_oyx[(i-1)*nlocaly +(j-1)];
               // if(my3drank==2)
              //  printf("check transfer alomg z: matrix[%d]=%lf\n",(i-1)*nlocaly +(j-1),matrix_inner_oyx[(i-1)*nlocaly +(j-1)]);              
                  
            }
        }



 }
 //---------------END OF Z---------------------------------------------
 /*---------ME send first values to outer neighbor and recieve firt values from inner neighbor along the z axis---------*/
    
    /*------------Exchange along X direction-------------------------------------------------------------*/
 if(dims[0]>1)
 { 
  for(int j=1;j<=nlocaly;j++)
          {
            
            for(int k=1;k<=nlocalz;k++)
               { 
                 matrix_north_oyz[(j-1)*nlocalz +(k-1)]=array[1][j][k];
               }
          }
   
  for(int j=1;j<=nlocaly;j++)
          {
 
            for(int k=1;k<=nlocalz;k++)
               {
                 matrix_south_oyz[(j-1)*nlocalz +(k-1)]=array[nlocalx][j][k];
               //    if(my3drank==size-1)
                //   printf("array[%d][%d][%d]=%lf\n",nlocalx,j,k,array[nlocalx][j][k]);
               }
          }

     MPI_Comm_rank(new_comm,&my3drank);
     MPI_Cart_coords(new_comm,my3drank,3,coords);
 //  if((coords[0]>1)&&(coords[0]<dims[0]-1))
 //  {
     MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr );

 /*---------ME send first values to north  neighbor and recieve firt values from south neighbor along the x axis---------*/

  
   
   
   MPI_Sendrecv_replace(matrix_north_oyz, nlocaly*nlocalz,  MPI_DOUBLE, up_nbr, tag,
                       down_nbr,tag, new_comm,&status);
     
   MPI_Sendrecv_replace(matrix_south_oyz, nlocaly*nlocalz,  MPI_DOUBLE, down_nbr, tag,
                        up_nbr,tag, new_comm,&status);
  

 /*------Copy values into ghost cells---------------*/
  
  for(int j=1;j<=nlocaly;j++)                                
        {                                                       
         for(int k=1;k<=nlocalz;k++)                            
           {                                                    
            array[nlocalx+1][j][k]=matrix_north_oyz[(j-1)*nlocalz +(k-1)];
            array[0][j][k]=matrix_south_oyz[(j-1)*nlocalz +(k-1)];
                                                                
           }                                                    
         } 


    
  }  
  MPI_Barrier(new_comm);
  //----------------END along X--------------------------------------------------------------------------------------------------------------     
        

  /////Exchange diagonally //////////////
    
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // one dimension array transfer
 //------------------------------------------------------------------------------------------------------------------------------------- 
  /*-----------perpendicular to  the y direction, copy last col     -----------*/
  //-----------------for lower points------------------------------------------------

 if(dims[1]>1)
 {



     for(int k=1;k<=nlocalz;k++)
        {
          array_right_oxz[k-1]= array[0][nlocaly][k];
          array_left_oxz[k-1] = array[0][1][k];
 
        }

    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);

    MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );

    MPI_Sendrecv_replace(array_right_oxz, nlocalz,  MPI_DOUBLE, right_nbr, tag,
                          left_nbr,tag, new_comm,&status);
     
    MPI_Sendrecv_replace(array_left_oxz, nlocalz,  MPI_DOUBLE, left_nbr, tag,
                           right_nbr,tag, new_comm,&status);


   //fill ghost points
   
   for(int k=1;k<=nlocalz;k++)
     {
       array[0][0][k] = array_right_oxz[k-1];
       array[0][nlocaly+1][k]=array_left_oxz[k-1];
     }
MPI_Barrier(new_comm);
  

  //--------for upper points-------------------------------
  std::memset(array_right_oxz,0,nlocalz*sizeof(double));
  std::memset(array_left_oxz,0,nlocalz*sizeof(double));


     for(int k=1;k<=nlocalz;k++)
        {
          array_right_oxz[k-1]= array[nlocalx+1][nlocaly][k];
          array_left_oxz[k-1] = array[nlocalx+1][1][k];
        }
    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
    
    MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
  
    MPI_Sendrecv_replace(array_right_oxz, nlocalz,  MPI_DOUBLE, right_nbr, tag,
                          left_nbr,tag, new_comm,&status);
     
     MPI_Sendrecv_replace(array_left_oxz, nlocalz,  MPI_DOUBLE, left_nbr, tag,
                          right_nbr,tag, new_comm,&status);


       
    //fill ghost points

  for(int k=1;k<=nlocalz;k++)
     {
       array[nlocalx+1][0][k] = array_right_oxz[k-1];
       array[nlocalx+1][nlocaly+1][k]=array_left_oxz[k-1];
     }
 

}
MPI_Barrier(new_comm);
 
/*-----------------------END of along Y-----------------------------------------------*/
//-------------ALONG Z direction-------------------------------------------------------
    /*----Lower points-----*/
 if(dims[2]>1)
 {
    for(int j=1;j<=nlocaly;j++)
        {
          array_inner_oxy[j-1]= array[0][j][nlocalz];
          array_outer_oxy[j-1]= array[0][j][1];

        }

    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
 
    MPI_Cart_shift( new_comm, 2, 1, &outer_nbr, &inner_nbr ); //along z
    
     MPI_Sendrecv_replace(array_inner_oxy, nlocaly,  MPI_DOUBLE, inner_nbr, tag,
                          outer_nbr,tag, new_comm,&status);
     
     MPI_Sendrecv_replace(array_outer_oxy, nlocaly,  MPI_DOUBLE, outer_nbr, tag,
                          inner_nbr,tag, new_comm,&status);
     
 
   //fill ghost points
   

    for(int j=1;j<=nlocaly;j++)
     {
       array[0][j][nlocalz+1] = array_outer_oxy[j-1];
       array[0][j][0]=array_inner_oxy[j-1];
     }
 
 //--------------- for UPPER points---------------------------------------------------------
   std::memset(array_outer_oxy,0,nlocaly*sizeof(double));
   std::memset(array_inner_oxy,0,nlocaly*sizeof(double));

  for(int j=1;j<=nlocaly;j++)
        {
          array_inner_oxy[j-1]= array[nlocalx+1][j][nlocalz];
          array_outer_oxy[j-1]= array[nlocalx+1][j][1];

        }

 
    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
     
 
      MPI_Cart_shift( new_comm, 2, 1, &outer_nbr, &inner_nbr );
     
     MPI_Sendrecv_replace(array_inner_oxy, nlocaly,  MPI_DOUBLE, inner_nbr, tag,
                          outer_nbr,tag, new_comm,&status);
     
     MPI_Sendrecv_replace(array_outer_oxy, nlocaly,  MPI_DOUBLE, outer_nbr, tag,
                          inner_nbr,tag, new_comm,&status);
     

     
     
     //fill ghost points
    
   
    for(int j=1;j<=nlocaly;j++)
     {
       array[nlocalx+1][j][nlocalz+1] = array_outer_oxy[j-1];
       array[nlocalx+1][j][0]=array_inner_oxy[j-1];
     }

 }
MPI_Barrier(new_comm);

//--------------------END ALONG Z-----------------------------------------
 //---For rows at the edges, each processor  will send four and recieve four.---------
if(dims[2]>1)
 {
   std::memset(array_right_oxz,0,nlocalz*sizeof(double));
   std::memset(array_left_oxz,0,nlocalz*sizeof(double));
   //Along y direction
   MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
   for(int i=1;i<=nlocalx;i++)
     {
       array_x_left[i-1]= array[i][nlocaly][0];
       array_x_right[i-1]=array[i][nlocaly][nlocalz+1];
  //     if(my3drank==2)
  //     printf("Left Array:left[%d]=%lf\n",i-1,array_x_left[i-1]);
     }
     
    //if((coords[1]>1)&&(coords[1]<dims[1]-1))
   // {
      MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr ); //along y
     
     MPI_Sendrecv_replace(array_x_left, nlocalx,  MPI_DOUBLE, right_nbr, tag,
                          left_nbr,tag, new_comm,&status);
     
     MPI_Sendrecv_replace(array_x_right, nlocalx,  MPI_DOUBLE, right_nbr, tag,
                          left_nbr,tag, new_comm,&status);
     
    //fill ghost points
   for(int i=1;i<=nlocalx;i++)
      {
        array[i][0][0]=array_x_left[i-1];
        array[i][0][nlocalz+1]=array_x_right[i-1];
      }

//-------------------------------------------------------------
  std::memset(array_x_left,0,nlocalx*sizeof(double));
  std::memset(array_x_right,0,nlocalx*sizeof(double));
  
  for(int i=1;i<=nlocalx;i++)
     {
       array_x_left[i-1]=array[i][1][0];
       array_x_right[i-1]=array[i][1][nlocalz+1];
     }


    MPI_Comm_rank(new_comm,&my3drank);
    MPI_Cart_coords(new_comm,my3drank,3,coords);
//    if((coords[1]>1)&&(coords[1]<dims[1]-1))
//    {
      MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
     
     MPI_Sendrecv_replace(array_x_left, nlocalx,  MPI_DOUBLE, left_nbr, tag,
                          right_nbr,tag, new_comm,&status);
     
     MPI_Sendrecv_replace(array_x_right, nlocalx,  MPI_DOUBLE, left_nbr, tag,
                          right_nbr,tag, new_comm,&status);

  
  //fill ghost points
  
  for(int i=1;i<=nlocalx;i++)
     {

      array[i][nlocaly+1][0] =array_x_left[i-1];
      array[i][nlocaly+1][nlocalz+1]=array_x_right[i-1]; 

     }
 }

MPI_Barrier(new_comm);
    
 /* ------ Implementation of point to point transfer ie Next next nearest neighbors--------*/
if((dims[1]>1)&&(dims[2]>1))
{
 //Exchange bottom edge points
  

  // copy points for exchange
  if(Procz==1)
  {
      
      value1_p1=array[0][1][1];
      value1_p2=array[0][nlocaly][1];
      value1_p3=array[0][1][nlocalz];
      value1_p4=array[0][nlocaly][nlocalz];
      
  }
 else
 {
    value1_p1=array[0][1][0];
    value1_p2=array[0][nlocaly][0];
    value1_p3=array[0][1][nlocalz+1];
    value1_p4=array[0][nlocaly][nlocalz+1];
 }
  MPI_Comm_rank(new_comm,&my3drank);
  MPI_Cart_coords(new_comm,my3drank,3,coords);
 
     MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
  
    
    MPI_Sendrecv_replace(&value1_p1, 1,  MPI_DOUBLE, left_nbr, tag,
                         right_nbr,tag, new_comm,&status);
    
    
    MPI_Sendrecv_replace(&value1_p2, 1,  MPI_DOUBLE, right_nbr, tag,
                         left_nbr,tag, new_comm,&status);

    
    MPI_Sendrecv_replace(&value1_p3, 1,  MPI_DOUBLE, left_nbr, tag,
                         right_nbr,tag, new_comm,&status);
   
    
    MPI_Sendrecv_replace(&value1_p4, 1,  MPI_DOUBLE, right_nbr, tag,
                         left_nbr,tag, new_comm,&status);


     // fill ghost points
       if(Procz==1)
          {
              array[0][0][1]=value1_p2;
              array[0][nlocaly+1][1]=value1_p1;
              array[0][0][nlocalz]=value1_p4;
              array[0][nlocaly+1][nlocalz]=value1_p3;
              
      
          }
       else
         {
           array[0][0][0]=value1_p2;
           array[0][nlocaly+1][0]=value1_p1;
           array[0][0][nlocalz+1]=value1_p4;
           array[0][nlocaly+1][nlocalz+1]=value1_p3;
        }
  

   MPI_Barrier(new_comm);

  
  /*-------------UPPER POINTS-----------------------------------*/
 
  //copy the points to be sent
  if(Procz==1)
  {

      value_p1=array[nlocalx+1][1][1];
      value_p2=array[nlocalx+1][nlocaly][1];
      value_p3=array[nlocalx+1][1][nlocalz];
      value_p4=array[nlocalx+1][nlocaly][nlocalz];

  }
else
{
  value_p1=array[nlocalx+1][1][0];
  value_p2=array[nlocalx+1][nlocaly][0];
  value_p3=array[nlocalx+1][1][nlocalz+1];
  value_p4=array[nlocalx+1][nlocaly][nlocalz+1];
}
 
 
  MPI_Comm_rank(new_comm,&my3drank);
  MPI_Cart_coords(new_comm,my3drank,3,coords);

     MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
  
    
    MPI_Sendrecv_replace(&value_p1, 1,  MPI_DOUBLE, left_nbr, tag,
                         right_nbr,tag, new_comm,&status);
    
    MPI_Sendrecv_replace(&value_p2, 1,  MPI_DOUBLE, right_nbr, tag,
                         left_nbr,tag, new_comm,&status);
    
    MPI_Sendrecv_replace(&value_p3, 1,  MPI_DOUBLE, left_nbr, tag,
                         right_nbr,tag, new_comm,&status);
                         
    MPI_Sendrecv_replace(&value_p4, 1,  MPI_DOUBLE, right_nbr, tag,
                         left_nbr,tag, new_comm,&status);

    
     

     // fill ghost points
       if(Procz==1)
       {
          array[nlocalx+1][0][1]=value_p2;
          array[nlocalx+1][0][nlocalz]=value_p4;
          array[nlocalx+1][nlocaly+1][0]=value_p1; 
          array[nlocalx+1][nlocaly+1][nlocalz]=value_p3;

           
           
       }
       else
       {
  //       if(my3drank==0)
    //      {
//            printf("CHECK POINT ENTERED, (array(%d,%d,0)=%lf,value_p[0]=%lf)\n",nlocalx+1,nlocaly+1,array[nlocalx+1][nlocaly+1][0],value_p1);
      //    }
          array[nlocalx+1][0][0]=value_p2;
          array[nlocalx+1][0][nlocalz+1]=value_p4;
          array[nlocalx+1][nlocaly+1][0]=value_p1;
          array[nlocalx+1][nlocaly+1][nlocalz+1]=value_p3;
       }

    
    



}
MPI_Barrier(new_comm);
/*
if(my3drank==0)
    {
       // std::cout<<"RANK="<<my3drank<<std::endl;
        printf("coord=(%d,%d,%d), rank=%d\n", coords[0],coords[1],coords[2],my3drank);
        int starty, startz, endy,endz;
       if (dims[1]>1)
        {
         starty=0;
         endy =nlocaly+2;
        }
      else
        {
         starty=1;
         endy =nlocaly+1;
         }
     dims[2]> 1?startz=0:startz=1;
     dims[2]> 1?endz=nlocalz+2:endz=nlocalz+1;


        for(int i=0;i<nlocalx+2;i++)
        {
     
            for(int j=0;j<nlocaly+2;j++)
            {
                for(int k=0;k<nlocalz+2;k++)
                {
                    
                    printf("Phi_old[%d][%d][%d]=%lf\n",i,j,k,array[i][j][k]);
                }
            }
        }
  
    } 

*/



 }

void Parallel_CahnHill3D::ProcessRank(MPI_Comm new_comm,int* newcoords)
 {
      MPI_Comm_rank(new_comm, &my3drank);
      MPI_Comm_size(new_comm, &size);
      int coords[3];
      MPI_Cart_coords(new_comm,my3drank,3,coords);
      newcoords[0]=coords[0]+1;
   
   
   
 }


void Parallel_CahnHill3D::ComputeLaplacianBase_p(MPI_Comm new_comm)
    {
      MPI_Comm_rank(new_comm, &my3drank);
      MPI_Comm_size(new_comm, &size);
      int dims[3],coords[3],periods[3];
      MPI_Cart_get (new_comm,3,dims,periods, coords );
       double AP=0.0, BP=0.0,  ATP=0.0, CP=0.0;

      if((dims[1]==1) && (dims[2]==1)) // if number of processor along y axis is 1 ie no spliting
        {
           printf("I am in laplacianBase with Pz=1\n");
        

    // printf("rank=%d,phi(0,0,0)=%lf\n",my3drank,PHI_old_p(0,0,0));
          
          setIndex_p(1,nlocaly);
        int i;
   #pragma omp parallel private(i,AP,BP,CP,ATP) num_threads(N_Threads)
   {
     #pragma omp for  
          for( i=1;i<=nlocalx;i++)
            {
             for(int j=1;j<=nlocaly;j++)
               {
                for(int k=1;k<=nlocalz;k++)
                  {
                                        
                                            
                    AP = C1*(PHI_old_p(i+1,j,k) + PHI_old_p(i-1,j,k)
                       + PHI_old_p(i,up1x[j],k) + PHI_old_p(i,down1x[j],k)
                       + PHI_old_p(i,j,up1x[k]) + PHI_old_p(i,j,down1x[k]) );

                    BP =C2*(PHI_old_p(i-1,up1x[j],k) +  PHI_old_p(i-1,down1x[j],k)
                      + PHI_old_p(i+1,up1x[j],k) +  PHI_old_p(i+1,down1x[j],k)
                      + PHI_old_p(i,down1x[j],up1x[k]) +  PHI_old_p(i,down1x[j],down1x[k])
                      + PHI_old_p(i,up1x[j],up1x[k]) +  PHI_old_p(i,up1x[j],down1x[k])
                      + PHI_old_p(i-1,j,up1x[k]) +  PHI_old_p(i-1,j,down1x[k])
                      + PHI_old_p(i+1,j,up1x[k]) +  PHI_old_p(i+1,j,down1x[k]));

                   CP=C3*( PHI_old_p(i-1,down1x[j],down1x[k]) + PHI_old_p(i-1,up1x[j],down1x[k])
                      + PHI_old_p(i-1,down1x[j],up1x[k]) + PHI_old_p(i-1,up1x[j],up1x[k])
                      + PHI_old_p(i+1,down1x[j],down1x[k]) + PHI_old_p(i+1,up1x[j], down1x[k])
                      + PHI_old_p(i+1,down1x[j],up1x[k]) + PHI_old_p(i+1,up1x[j],up1x[k]));
                 ATP = AP + BP + CP;
                 gamma_p(i,j,k)= g(PHI_old_p(i,j,k))-PHI_old_p(i,j,k) + D*(ATP-PHI_old_p(i,j,k));
               //  printf("rank=%d: gamma_p(%d,%d,%d)=%lf\n",my3drank,i,j,k,gamma_p(i,j,k));

            
                 }
              }
            }
             printf("rank=%d,phi(0,0,0)=%lf\n",my3drank,PHI_old_p(0,0,0));
   }
     }
  else
  {   
   // printf("C1=%lf, C2=%lf, C3=%lf, D=%lf\n",CahnHill3D::C1,C2,C3,D);
   int i;
   #pragma omp parallel private(i,AP,BP,CP,ATP) num_threads(N_Threads)
   {
     #pragma omp for
      for( i=1;i<=nlocalx;i++)
         {
          for(int j=1;j<=nlocaly;j++)
            {
             for(int k=1;k<=nlocalz;k++)
              {
                AP = C1*(PHI_old_p(i+1,j,k) + PHI_old_p(i-1,j,k)
                       + PHI_old_p(i,j+1,k) + PHI_old_p(i,j-1,k)
                       + PHI_old_p(i,j,k+1) + PHI_old_p(i,j,k-1) );

                BP =C2*(PHI_old_p(i-1,j+1,k) +  PHI_old_p(i-1,j-1,k)
                      + PHI_old_p(i+1,j+1,k) +  PHI_old_p(i+1,j-1,k)
                      + PHI_old_p(i,j-1,k+1) +  PHI_old_p(i,j-1,k-1)
                      + PHI_old_p(i,j+1,k+1) +  PHI_old_p(i,j+1,k-1)
                      + PHI_old_p(i-1,j,k+1) +  PHI_old_p(i-1,j,k-1)
                      + PHI_old_p(i+1,j,k+1) +  PHI_old_p(i+1,j,k-1));

                CP=C3*( PHI_old_p(i-1,j-1,k-1) + PHI_old_p(i-1,j+1,k-1)
                      + PHI_old_p(i-1,j-1,k+1) + PHI_old_p(i-1,j+1,k+1)
                      + PHI_old_p(i+1,j-1,k-1) + PHI_old_p(i+1,j+1, k-1)
                      + PHI_old_p(i+1,j-1,k+1) + PHI_old_p(i+1,j+1,k+1));
                 ATP = AP + BP + CP;
                 gamma_p(i,j,k)= g(PHI_old_p(i,j,k))-PHI_old_p(i,j,k) + D*(ATP-PHI_old_p(i,j,k));
             /*    if(my3drank==0)
                  {
                    printf("gamma_p(%d,%d,%d)=%lf\n",i,j,k,gamma_p(i,j,k));
                    
                  }*/


                 }
              }
            }

   }
  }

    }

void Parallel_CahnHill3D::setSecond_laplacian(MPI_Comm new_comm)
    {
       MPI_Comm_rank(new_comm, &my3drank);    
      MPI_Comm_size(new_comm, &size);        
      double AP=0.0, BP=0.0,  CP=0.0;
      int dims[3],coords[3], periods[3];
      MPI_Cart_get (new_comm,3,dims,periods, coords );

  if((dims[2]==1) && (dims[1]==1))
     {
         printf("In second laplacian Pz=1\n");
        
           setIndex_p(1,nlocaly);
           int i;
      #pragma omp parallel private(i,AP,BP,CP) num_threads(N_Threads)
        {
          #pragma omp for
            for( i=1;i<=nlocalx;i++)            
              {                                   
               for(int j=1;j<=nlocaly;j++)        
                 {                                
                  for(int k=1;k<=nlocalz;k++)     
                     { 
                     
                  
                        AP = C1*(gamma_p(i+1,j,k)     + gamma_p(i-1,j,k)
                               + gamma_p(i,up1x[j],k) + gamma_p(i,down1x[j],k)
                               + gamma_p(i,j,up1x[k]) + gamma_p(i,j,down1x[k]) );

                        BP =C2*(gamma_p(i-1,up1x[j],k)       +  gamma_p(i-1,down1x[j],k)
                              + gamma_p(i+1,up1x[j],k)       +  gamma_p(i+1,down1x[j],k)
                              + gamma_p(i,down1x[j],up1x[k]) +  gamma_p(i,down1x[j],down1x[k])
                              + gamma_p(i,up1x[j],up1x[k])   +  gamma_p(i,up1x[j],down1x[k])
                              + gamma_p(i-1,j,up1x[k])       +  gamma_p(i-1,j,down1x[k])
                              + gamma_p(i+1,j,up1x[k])       +  gamma_p(i+1,j,down1x[k]));

                      CP =C3*( gamma_p(i-1,down1x[j],down1x[k]) + gamma_p(i-1,up1x[j],down1x[k])
                             + gamma_p(i-1,down1x[j],up1x[k])   + gamma_p(i-1,up1x[j],up1x[k])
                             + gamma_p(i+1,down1x[j],down1x[k]) + gamma_p(i+1,up1x[j], down1x[k])
                             + gamma_p(i+1,down1x[j],up1x[k])   + gamma_p(i+1,up1x[j],up1x[k]));
                       
                       Laplacian2_p(i,j,k)=AP + BP + CP;
                       
                       
                     }
                 }
               }
        }
       }

else
   {   
     int i;
   #pragma omp parallel private(i,AP,BP,CP) num_threads(N_Threads)
    {
     #pragma omp for                                      
      for( i=1;i<=nlocalx;i++)            
         {                                   
          for(int j=1;j<=nlocaly;j++)        
            {                                
             for(int k=1;k<=nlocalz;k++)     
              {                              
                AP = C1*(gamma_p(i+1,j,k) + gamma_p(i-1,j,k)
                       + gamma_p(i,j+1,k) + gamma_p(i,j-1,k)
                       + gamma_p(i,j,k+1) + gamma_p(i,j,k-1) );
                                             
                BP =C2*(gamma_p(i-1,j+1,k) +  gamma_p(i-1,j-1,k)
                      + gamma_p(i+1,j+1,k) +  gamma_p(i+1,j-1,k)
                      + gamma_p(i,j-1,k+1) +  gamma_p(i,j-1,k-1)
                      + gamma_p(i,j+1,k+1) +  gamma_p(i,j+1,k-1)
                      + gamma_p(i-1,j,k+1) +  gamma_p(i-1,j,k-1)
                      + gamma_p(i+1,j,k+1) +  gamma_p(i+1,j,k-1));
                                             
                CP=C3*( gamma_p(i-1,j-1,k-1) + gamma_p(i-1,j+1,k-1)
                      + gamma_p(i-1,j-1,k+1) + gamma_p(i-1,j+1,k+1)
                      + gamma_p(i+1,j-1,k-1) + gamma_p(i+1,j+1, k-1)
                      + gamma_p(i+1,j-1,k+1) + gamma_p(i+1,j+1,k+1));

                Laplacian2_p(i,j,k)=AP + BP + CP;
            }
          }
        }
    }
   }
   }

void Parallel_CahnHill3D::WriteToFile(MPI_Comm new_comm)
    {  

        MPI_Comm_rank(new_comm, &rank);
      MPI_Comm_size(new_comm, &size);

      char fname[200];
      sprintf(fname,"solution%d.dat",rank);
      FILE* fp;
      fp = fopen(fname,"w");
      for(int i=1;i<nlocalx+1;i++)
         {
           for(int j =1; j<nlocaly+1;j++)
             {
               for(int k=1;k<nlocalz+1;k++)
                 {
                  fprintf(fp,"%20.6f", PHI_p[i][j][k]);
                 }

             }
          // fprintf(fp,"\n");

        }
     fclose(fp);


 }


void Parallel_CahnHill3D::WriteToFile_MPI(MPI_Comm new_comm)
  {
         MPI_Status status; 
        MPI_File     fh;  
        MPI_Datatype filetype;
       int dims[3],coords[3], periods[3], start_indices[3], localsizes[3];
       //int globalsizes[3],localsizes[3],memsizes[3];
       //globalsizes[0]=Nx;
       //globalsizes[1]=Ny;
       //globalsizes[2]=Nz;
       //localsizes[0]=nlocalx; localsizes[1]=nlocaly; localsizes[2]=nlocalz;
       localsizes[0]=(Nx/Procx)+Nx%Procx; localsizes[1]=(Ny/Procy) + Ny%Procy; localsizes[2]=(Nz/Procz)+Nz%Procz;
      int  array_of_gsizes[3], array_of_distribs[3],array_of_dargs[3], array_of_psizes[3];
      for(int i=0;i<3;i++)
      {
          array_of_gsizes[i] = 128;
          array_of_distribs[i] = MPI_DISTRIBUTE_BLOCK;
          array_of_dargs[i] = MPI_DISTRIBUTE_DFLT_DARG;
       //   array_of_psizes[i] = 2;
          
      }
     array_of_psizes[0] =Procx; array_of_psizes[1] =Procy; array_of_psizes[2] =Procz;
      
      
       MPI_Cart_get (new_comm,3,dims,periods, coords );
       MPI_Comm_rank(new_comm, &my3drank);
       MPI_Comm_size(new_comm, &size);
       MPI_Cart_coords(new_comm,my3drank,3,coords);

       start_indices[0]=coords[0]*localsizes[0];
       start_indices[1]=coords[1]*localsizes[1];
       start_indices[2]=coords[2]*localsizes[2];
       printf("rank=%d,r=%d,start_indices=(%d,%d,%d)\n",my3drank,Nx%Procx,start_indices[0],start_indices[1],start_indices[2]);     
      
      
     // MPI_Type_create_subarray(3, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);
      
      
     MPI_Type_create_darray(size, my3drank, 3,array_of_gsizes, array_of_distribs, array_of_dargs, array_of_psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);
      
       MPI_Type_commit(&filetype);
       
        char outputfilename[]="datafile_mpi";
        char filememview[]="native";     
       MPI_File_open(new_comm, outputfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &fh); 
       MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, filememview,MPI_INFO_NULL);
     //  MPI_File_write(fh, PHI_p, nlocalx*nlocaly*nlocalz, MPI_DOUBLE, &status);
       MPI_File_write(fh, PHI_local_result, nlocalx*nlocaly*nlocalz, MPI_DOUBLE, &status);
      

       MPI_File_close(&fh);
       MPI_Barrier(new_comm);
       MPI_Type_free(&filetype);
        
  }
void Parallel_CahnHill3D::ReadFile_MPI(MPI_Comm new_comm)
  {
        MPI_Status status;
        MPI_File     fh;
        MPI_Datatype filetype;
        
        double *PHI_read_in=new double [Nx*Ny*Nz];      
        
       // new_comm= CreatCartesianTopology();
        MPI_Comm_rank(new_comm, &my3drank);
        MPI_Comm_size(new_comm, &size);
      
        //MPI_get_file_size
        /*
       int  array_of_gsizes[3], array_of_distribs[3],array_of_dargs[3], array_of_psizes[3];
       for(int i=0;i<3;i++)
          {
            array_of_gsizes[i] = 8;
            array_of_distribs[i] = MPI_DISTRIBUTE_BLOCK;
            array_of_dargs[i] = MPI_DISTRIBUTE_DFLT_DARG;
       //   array_of_psizes[i] = 2;
          
          }
       array_of_psizes[0] =Procx; array_of_psizes[1] =Procy; array_of_psizes[2] =Procz;*/
      int dims[3],coords[3], periods[3], start_indices[3], localsizes[3];
       int globalsizes[3];//,memsizes[3];
       globalsizes[0]=Nx;
       globalsizes[1]=Ny;
       globalsizes[2]=Nz;
       //localsizes[0]=nlocalx; localsizes[1]=nlocaly; localsizes[2]=nlocalz;
       localsizes[0]=(Nx/Procx)+Nx%Procx; localsizes[1]=(Ny/Procy) + Ny%Procy; localsizes[2]=(Nz/Procz)+Nz%Procz;
       MPI_Cart_get (new_comm,3,dims,periods, coords );
       MPI_Comm_rank(new_comm, &my3drank);
       MPI_Comm_size(new_comm, &size);
       MPI_Cart_coords(new_comm,my3drank,3,coords);

       start_indices[0]=coords[0]*localsizes[0];
       start_indices[1]=coords[1]*localsizes[1];
       start_indices[2]=coords[2]*localsizes[2];
       //printf("rank=%d,r=%d,start_indices=(%d,%d,%d)\n",my3drank,Nx%Procx,start_indices[0],start_indices[1],start_indices[2]);     
      
      
    
       
//       MPI_Type_create_darray(size, my3drank, 3,array_of_gsizes, array_of_distribs, array_of_dargs, array_of_psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_create_subarray(3, globalsizes, localsizes, start_indices,MPI_ORDER_C, MPI_DOUBLE, &filetype);

      
       MPI_Type_commit(&filetype);
     
      
        
         char inputfilename[]="datafile_mpi";//"output_0.dat";
         char filememview[]="native";

        //offset = (my3drank/size)*Ny*nlocalx * sizeof(double)  + (rank % 4) * n_local_cols * size
        MPI_File_open(new_comm, inputfilename, MPI_MODE_RDONLY,MPI_INFO_NULL, &fh);
       // MPI_File_seek( fh, 0, MPI_SEEK_SET );
       printf("READIN\n"); 
        MPI_File_set_view(fh,0,MPI_DOUBLE,filetype,filememview,MPI_INFO_NULL);
        MPI_File_read_all( fh, PHI_read_in, nlocalx*nlocaly*nlocalz, filetype, &status );
        
        MPI_File_close(&fh);
        MPI_Barrier(new_comm);
        MPI_Type_free(&filetype);
      
        for(int i=1;i<=nlocalx;i++)
          {
           for(int j=1;j<=nlocaly;j++) 
              {
               for(int k=1;k<=nlocalz;k++)
                 {
                   PHI_old_p[i][j][k]=PHI_read_in[((i-1)*nlocaly + (j-1))*nlocalz + k-1];
                //   if(my3drank==0)
                 //  {
              //     printf("PHI_old_p[%d][%d][%d]=%lf\n",i,j,k,PHI_old_p[i][j][k]);
                 //  }
                 }
              }
          }


        //FreeCommunicator( new_comm);
        delete [] PHI_read_in;


  }

void Parallel_CahnHill3D::ReadFile(std::ifstream& infile)
 {

     
            MPI_Comm new_comm;
            new_comm= CreatCartesianTopology();

            //MPI_Status status;
            MPI_Comm_rank(new_comm, &rank);
            MPI_Comm_size(new_comm, &size);
            tag = 201;
            int coords[3];
            double *u1, ***u2;
              u1= new double[Nx*Ny*Nz];
              //std::memset(u1, 0, Nx*Ny*Nz*sizeof(double));
             u2 = new double**[Nx];
             for(int i=0;i<Nx; i++)
             {
               u2[i]=new double*[Ny];
               for(int j=0;j<Ny;j++)
                  {
                    u2[i][j]  =  new double[Nz];  
                  }
             }
             
              for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {
                      for(int k=0;k<Nz;k++)
                        {
                          infile>>u2[i][j][k];
                        }
                      }

                    }
                
            if(rank==0)
             {
      
              for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {
                      for(int k=0;k<Nz;k++)
                        {
                          infile>>u1[(i*Ny + j)*Nz +k];
                        }
                      }

                    }
             }
        MPI_Bcast( u1, Nx*Ny*Nz, MPI_DOUBLE, 0, new_comm); 

          MPI_Barrier(new_comm);
 
          MPI_Cart_coords(new_comm,rank,3,coords);
           
             
             for(int i=1;i<=nlocalx;i++)
                {

                   for(int j=1;j<=nlocaly;j++)
                       {
                        for(int k=1;k<=nlocalz;k++)
                         {
                       
                           //PHI_old_p[i][j][k]=u1[ ((rank*nlocalx+i-1)*nlocaly +j-1)*nlocalz+ k-1];

                          // PHI_old_p[i][j][k]=u1[((coords[0]*nlocalx + i-1)*nlocaly + coords[1]*nlocaly+j-1)*nlocalz + coords[2]*nlocalz + k-1];
                          PHI_old_p[i][j][k]=u2[coords[0]*nlocalx + i-1][coords[1]*nlocaly + j-1][coords[2]*nlocalz + k-1];

                          // if(rank==4)
                          //  {
                            //  printf("rank=%d, phi_old_p[%d][%d][%d]=%lf\n",rank,i,j,k,PHI_old_p[i][j][k]);
                           // }
                         }
                       }

               }

        for(int i=0;i<Nx;i++)
           {
            for(int j=0;j<Ny;j++)
               { 
                delete[] u2[i][j];
          
                } 
            delete[] u2[i];
           }  
    
         delete []u1;
         delete []u2;
         FreeCommunicator( new_comm);

}
void Parallel_CahnHill3D:: setIndex_p(int start, int end)
    {
        for(int s=start; s<=end; s++)
           {
                up1x[s]=s+1;
                down1x[s]=s-1;
            }
            down1x[start]=end;
            up1x[end]=start;
 
   } 
