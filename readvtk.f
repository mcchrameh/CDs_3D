   !!   PROGRAM rEAD(P, NX,NY,NZ)
      SUBROUTINE rEAD(P, NX,NY,NZ)
      IMPLICIT NONE
      INTEGER I,J,K,NX,NY,NZ,N3
    !!  PARAMETER (NX=32)
    !!  PARAMETER (NY=32)
    !!  PARAMETER (NZ=32)
      PARAMETER (N3=NX*NY*NZ)
      DOUBLE PRECISION P(NX,NY,NZ)
      character (LEN = 40) :: name1,name2,name3,name4,name5
      character (LEN = 40) :: name6,name7,name8,name9,name10  
      OPEN(unit=2,FILE="output_0.dat")
      OPEN(40,FILE='file.vtk')
     
	DO I = 1,NX
	     DO J = 1,NY
	          DO K = 1,NZ
                PRINT*,I,J,K
        !             READ(2,*) P(I,J,K)
      ENDDO
      ENDDO
      ENDDO
        
      name1='# '//'vtk '//'DataFile '//'Version '//'3.0'
      name2='vtkfile'
      name3='ASCII'
      name4='DATASET '//'STRUCTURED_POINTS'
      name5='DIMENSIONS'
      name6='ORIGIN'
      name7='SPACING'
      name8='POINT_DATA' 
      name9='SCALARS '//'scalars '//'float'
      name10='LOOKUP_TABLE '//'default'
      
!      WRITE(40,*)'# vtk DataFile Version 3.0'
!      WRITE(40,*)'vtkfile'
!      WRITE(40,*)'ASCII'
!      WRITE(40,*)'DATASET STRUCTURED_POINTS'
!      WRITE(40,*)'DIMENSIONS',NX,Ny,NZ
!      WRITE(40,*)'ORIGIN',0,0,0
!      WRITE(40,*)'SPACING',1,1,1
!      WRITE(40,*)'POINT_DATA',N3 
!      WRITE(40,*)'SCALARS scalars float 1'
!      WRITE(40,*)'LOOKUP_TABLE default'
      WRITE(40,*) name1
      WRITE(40,*) name2
      WRITE(40,*) name3
      WRITE(40,*) name4
      WRITE(40,*) name5,NX,Ny,NZ
      WRITE(40,*) name6,0,0,0
      WRITE(40,*) name7,1,1,1
      WRITE(40,*) name8,N3
      WRITE(40,*) name9,1
      WRITE(40,*) name10
  
     
	DO I = 1,NX
	     DO J = 1,NY
!	          DO K = 1,NZ
           
      write(40,255) (P(K,J,I),K=1,NZ)
                    PRINT*, K,I,J
      enddo     
      write (40,255)      
      ENDDO
255   FORMAT (1024f10.6)
      close(2)
      close(40)
      end program read  
