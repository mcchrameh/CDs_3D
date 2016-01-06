#ifndef _CAHNHILL_H
#define _CAHNHILL_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "mpi.h"



class CahnHill3D
{ 
  protected:
       int nx, ny,nz;
       double delta_t, M, delta_x, delta_y,delta_z,u,b,k, C1,C2,C3;
       double tau,A,F,v,D,dxx,B,r,Max_time_iteration;
       
       double ***PHI,***PHI_old,***gamma,***Laplacian2,*Q, *x,*y,*z;
       int *down1x, *down2x, *up1x, *up2x,*down1x_2, *up1x_2;
       
      
      
       void UpdateSolution();
       void initialCondition();
       void setLaplacianBase();
       void setIndex(int start,int end);
       void setIndex2(int start,int end);

       double g(double phi);
       void SetSecondLaplacian2();
       void FiniteDifferenceScheme();
      
      // void Read_input_parameters(int *integers, double *rationals);
      
                                                               
  public:
        CahnHill3D();
        CahnHill3D(int Nx,int Ny,int Nz);
       virtual ~CahnHill3D();
       void  display();
       std::string make_output_filename(int index)
           {  
             std::ostringstream ss;
            ss << "output_" << index << ".dat"; 
             return ss.str(); 
            }
       void save_vtk(double ***Q1, int Nx, int Ny, int Nz);
       void result(int count);
       void Solve();
       void Solver2();
       void ReadFilein(std::ifstream& infile);

       friend double **alloc_2d_int(int rows,int cols);
 };
     
      

class Parallel_CahnHill3D//: public  CahnHill3D
{
   protected:
            double delta_t, M, delta_x, delta_y,delta_z,u,b,k, C1,C2,C3;
            double tau,A,F,v,D,dxx,B,r,Max_time_iteration;
            int  nlocalx,nlocaly,nlocalz, remainder, N, Nx,Ny,Nz, Procx,Procy,Procz;
            int  rank, size, tag,low, ROW,COL,my3drank, right_nbr, left_nbr, up_nbr,down_nbr, inner_nbr,outer_nbr;
            int  offset,startrow_x,endrow_x, startcol_y,endcol_y, startdepth_z, enddepth_z;
            double ***PHI_p,***PHI_old_p,***gamma_p,***Laplacian2_p, *PHI_local_result, *PHI_global_result;
            //double **,*u_firstrow_z,*u_lastrow_y,*u_firstrow_y,*u_lastrow_x,*u_firstrow_x;
            double  *matrix_right_oxz, *matrix_left_oxz,*matrix_north_oyz,*matrix_south_oyz,*matrix_inner_oyx, *matrix_outer_oyx,
                    *array_left_oxz,*array_right_oxz,*array_north_oyz, *array_south_oyz, *array_inner_oxy,*array_outer_oxy, *array_x_left,
                     *array_x_right, value_p1,value_p2,value_p3,value_p4, value1_p1,value1_p2,value1_p3,value1_p4;
           // double *Boundary_top1, *Boundary_bottom1, *array_gamma_top,*array_gamma_bottom, *u_store,*bc,*U;
          //  int  *right1x,*left1x;

           int *down1x, *down2x, *up1x, *up2x,*down1x_2, *up1x_2;
   public:

        // Parallel_CahnHill3D();
         Parallel_CahnHill3D();
        ~Parallel_CahnHill3D();
         void  Initialize_parallel(MPI_Comm new_comm);
         void  ExchangeData(MPI_Comm new_comm, double ***array);

         void  ComputeLaplacianBase_p(MPI_Comm new_comm);
         double g(double phi);
         void  setSecond_laplacian(MPI_Comm new_comm);
         void  FiniteDifferenceScheme_p(MPI_Comm new_comm);
         void  UpdateSolution_p(MPI_Comm new_comm);
         void  WriteToFile_MPI(MPI_Comm new_comm);
         void  ProcessRank(MPI_Comm new_comm, int* newcoords);
         void Read_input_parameters(int *integers, double *rationals);
         void  parallel_solver();
         void  WriteToFile(MPI_Comm new_comm);
         void  ReadFile(std::ifstream& infile);
         void  ReadFile_MPI(MPI_Comm new_comm);
         MPI_Comm  CreatCartesianTopology();
         void  FreeCommunicator(MPI_Comm new_comm);
         void setIndex_p(int start, int end);
    

  

};

#endif
