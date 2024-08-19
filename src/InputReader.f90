!===============================================================================
      subroutine Read_Input_and_Allocate
!===============================================================================
      use mainvar

     !------------------------------------------------
      read(5,*) Fname
      read(5,*) QC_interface
      read(5,*) QC_route
      read(5,*) nproc, mem
      read(5,*) charge
      read(5,*) multi(1), multi(2)
      read(5,*) Method
      read(5,*) Basis
      read(5,*) OptMethod
      read(5,*) StpLim
      read(5,*) maxstep
      read(5,*) ConvPar

     !Threshold paramteres for convergence
      if(trim(ConvPar) == "New_Conv_Params") then
         read(5,*) TDE, TDXMax, TDXRMS, TGMax, TGRMS
      else
        !Default values same as in Gaussian (except TDE)

        !energy difference
         TDE=5.d-5
        !max and averaged norm of step (old X --> new X) 
         TDXMax=4.d-3; TDXRMS=2.5d-3
        !max and averaged norm of gradient
         TGMax=7.d-4;  TGRMS=5.d-4 
      endif
     !------------------------------------------------

     !------------------------------------------------
      open(7,file="Initial_Geom.xyz")

      read(7,*) natom
      ndim  = 3*natom
      call allocate_mainvars

      read(7,*)

      do i=1,natom
         jx=3*i-2
         jy=3*i-1
         jz=3*i
         read(7,*) Atom(i),X_1(jx),X_1(jy),X_1(jz)
      enddo
      close(7)
     !------------------------------------------------

      end subroutine
!===============================================================================
