      Program Optimizer
      use mainvar
      implicit none

      integer :: i,istep, jx,jy,jz
      integer :: Conv


     !Reading X_1 initial geom and other info for QC calculation
      call Read_Input_and_Allocate

      call Create_QC_interface

      call Initialize_Optimizer

     !-------------------------------------------------------------------
     ! Optimization of minimum energy crossing-point
     !-------------------------------------------------------------------
      do istep=1,maxstep
        
        !Calculate energies and forces at the new postion X_2 = X_3
         call Calc_Ene_Force(natom,istep,Fname,Atom,X_2,Ea,Eb,Ga,Gb)


        !Compute the new Effective Gradient:
        !Geff_2 = fac * DeltaE * PerpG  +  ParG
         call Effective_Gradient(ndim,Ea,Eb,Ga,Gb,ParG,PerpG,Geff_2)

        !---------------------------------------------------------------------
        !In all methods iH_2 and X_3 are the outputs, everyting else is input
        !---------------------------------------------------------------------
         if(istep == 1) then
          !We need one more previous step to start any of the optimizers 
           call Simple_SteepestDescent(ndim,X_2,Geff_2,X_3)
           iH_2 = iH_1
         else   
           select case (trim(OptMethod))
           case("BFGS_old") 
             call BFGS_J_UpdateX(ndim,StpLim,X_1,X_2,Geff_1,Geff_2,iH_1,iH_2,X_3)

           case("BFGS")
             call BFGS_new_Update(ndim,StpLim,X_1,X_2,Geff_1,Geff_2,iH_1,iH_2,X_3)

           case("BBGrad1")
             call BB_TwoPointGrad(1,ndim,StpLim,X_1, X_2, Geff_1, Geff_2, X_3)

           case("BBGrad2")
             call BB_TwoPointGrad(2,ndim, StpLim, X_1, X_2, Geff_1, Geff_2, X_3)

           case("Broyden")
             call DFP_Update(ndim,StpLim,X_1,X_2,Geff_1,Geff_2,iH_1,iH_2,X_3)

           case("DFP")
             call DFP_Update(ndim,StpLim,X_1,X_2,Geff_1,Geff_2,iH_1,iH_2,X_3)

           case("SR1")
             call SR1_Update(ndim,StpLim,X_1,X_2,Geff_1,Geff_2,iH_1,iH_2,X_3)

           case default
             stop "Invalid Optimzer! You can use: BBGrad1, BBGrad1, &
                                              BFGS, DFP, SR1, Broyden"
           end select
         endif
        !-----------------------------------------------------------------------


         call TestConv_and_Print(istep,Conv)


        !Saving old coords
         X_1 = X_2;     X_2 = X_3


        !Saving old effective gradient and inverse Hessian
         Geff_1 = Geff_2;   iH_1 = iH_2

         if(Conv == 1) exit

      enddo !end of optimizition cycle
     !----------------------------------------------------------------

      write(6,*)

      if(istep < maxstep) then
        write(6,*) "!!!!!Convergence!!!!!"
        write(6,*) " End of Optimization"
      else
        write(6,*) "Sorry, NO Convergence :("
      endif



      end program


