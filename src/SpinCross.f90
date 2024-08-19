!===================================================================================
      module mainvar
      implicit none

     !Variables for the QC input
      integer :: nproc,mem,charge
      integer, dimension(2) :: multi
      character(len=100) :: Method, Basis, Fname, QC_interface
      character(len=100) :: OptMethod 
      character(len=150) :: QC_route

     !All other general variables 
      integer :: natom,ndim,maxstep
      real*8  :: Ea, Eb
      character(len=20) :: ConvPar, StpLim

      real*8,  allocatable :: iH_1(:,:), iH_2(:,:), ParG(:), PerpG(:)
      real*8,  allocatable :: Ga(:), Gb(:), Geff_1(:), Geff_2(:)
      real*8,  allocatable :: X_1(:), X_2(:), X_3(:)

      character(len=2), allocatable :: Atom(:)

     !Convergence parameters
     !real*8, parameter :: TGMax=7.d-4,TGRMS=5.d-4
     !real*8, parameter :: TDE=5.d-5,TDXMax=4.d-3,TDXRMS=2.5d-3
      real*8 :: TGMax,TGRMS, TDE,TDXMax,TDXRMS

      real*8, parameter :: tokJmol=2625.5d0, tokcalmol=627.5d0


      contains 
        !-----------------------------------
         subroutine allocate_mainvars
         implicit none

         allocate(iH_1(ndim,ndim)) 
         allocate(iH_2(ndim,ndim))
         allocate(ParG(ndim))     
         allocate(PerpG(ndim))     
         allocate(Ga(ndim))        
         allocate(Gb(ndim))
         allocate(Geff_1(ndim))
         allocate(Geff_2(ndim))
         allocate(X_1(ndim))
         allocate(X_2(ndim))
         allocate(X_3(ndim))
         allocate(Atom(natom))

         end subroutine
        !-----------------------------------

      end module
!=====================================================================================
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



!===============================================================================
      subroutine Initialize_Optimizer
!===============================================================================
!     Initialize Effective gradient and inverse Hessian
!     iH_1 values of 0.7 correspond to Hessians of ~1.4 Hartree/Angs**2
!     for more details see [subroutine Effective_Gradient]
!===============================================================================
      use mainvar
      implicit none
      integer :: i

     !effective gradient
      Geff_1 = 0.d0

     !inverse hessian
      iH_1 = 0.d0
      forall(i=1:ndim) iH_1(i,i) = 0.7d0

      iH_2 = iH_1

     !Geoms in the 0th and the 1st steps are the same
      X_2 = X_1

      end subroutine
!===============================================================================



!===============================================================================
      subroutine Effective_Gradient(N,Ea,Eb,Ga,Gb,ParG,PerpG,Geff)
!===============================================================================
!  Computes the parallel and perpendicular compenents of the Effective Gradient,
!  As well as the effective gradient itself.
       
!  The Hessian term along "normal" coordinates is empirically about 1.4 Hartree / Angstrom**2
!  Using facPP ~ 140 /Hartree means that along this coordinate, too, the Hessian is about right.
!     facPp chosen so that inverse Hessian ~ diagonal 0.7 Ang*Ang/Hartree
!
!  These factors are only really important for the first step
!  The "difference gradient" is typically ca. 0.075 Hartree/Bohr.
!      i.e. 0.14 Hartree/Angstrom.
!  Assuming this is constant, this gives a Hessian term for the func (Ea-Eb)**2
!     of ca. 0.01 Hartree**2 / Angstrom**2  (0.14**2 / 2)
!===============================================================================
      implicit none
       
      integer n, i
      real*8 :: Ea, Eb, overlap, sqr_dG, scal
      real*8,dimension(N) :: Ga, Gb, Geff, ParG, PerpG,dG
      real*8, parameter :: facPP=140.d0, facP=1.d0

      dG = Ga - Gb 

      overlap = dot_product(Ga,dG)

      sqr_dG = dot_product(dG,dG)

      scal = overlap/sqr_dG

     !perpendicular gradient
      PerpG = (Ea - Eb)*dG

     !paralell gradient
      ParG = Ga - dG*scal

     !effective gradient
      Geff = facPP*PerpG + facP*ParG

      end subroutine
!===============================================================================


!===============================================================================
      subroutine BB_TwoPointGrad(Eq, n, StpLim, X_1, X_2, G_1, G_2, X_3)
!===============================================================================
!     Barizilai-Borwein two-point step size gradient method
!     (IMA Journal of Numerical Analysis, 1988, 8, 141-149)              
!
!     X_new = X_old - alpha*G_old                   
! 
!     where alpha depends on the two previous X_1, X_2 and G_1 and G_2                 
! 
!     Input:
!         n = 3*natom (dimension of the problem)
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!
!     Output:
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      
      integer ::  n, nstep, Eq
      real*8 :: Ovrlp_XG, Xsq, Gsq, alpha
      real*8, dimension(n) :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n) :: dX, dG, ChgeX
      character(len=20) :: StpLim
       
      dG = G_2 - G_1
      dX = X_2 - X_1

      Ovrlp_XG = dot_product(dG,dX) ! <dX|dG>
      Xsq      = dot_product(dX,dX) ! <dX|dX>
      Gsq      = dot_product(dG,dG) ! <dG|dG>

      if(Eq == 1) then
         alpha = Ovrlp_XG / Gsq   !Eq(5) in Ref paper
      elseif(Eq == 2) then
         alpha = Xsq / Ovrlp_XG   !Eq(6) in Ref paper
      else
         stop "Problem in BB_TwoPointGrad routine"
      endif

      ChgeX = -alpha*G_2

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX
       
      end subroutine
!===============================================================================

!===============================================================================
      subroutine Step_Limit(n,ChgeX)
!===============================================================================
!     Scale back step length (ChgeX) if its too big
!===============================================================================
      implicit none
      integer :: n !dimension of the problem n=3*natom
      real*8  :: stpl, lgstst,stpmax
      real*8  :: ChgeX(n)
      real*8, parameter :: STPMX = 0.1d0

      stpmax = n * STPMX

      stpl = norm2(ChgeX)


      if(stpl .gt. stpmax) then
         ChgeX = ChgeX / stpl * stpmax
      endif

      lgstst = maxval(abs(ChgeX))

      if(lgstst .gt. STPMX) then
         ChgeX = ChgeX / lgstst * STPMX
      endif

      end subroutine
!===============================================================================



!===============================================================================
      subroutine Simple_SteepestDescent(n, Xold, Grad, Xnew)
!===============================================================================
!     Simple steepest descent
!     x_new = x_old - alpha*Grad   -->   dX = -beta*Grad
!===============================================================================
      implicit none
      integer :: n !dimension of the problem n=3*natom
      real*8, dimension(n) :: ChgeX, Xold, Xnew, Grad
      real*8, parameter :: beta = 0.7d0

      ChgeX = -beta*Grad

      call Step_Limit(n,ChgeX)

      Xnew = Xold + ChgeX

      end subroutine
!===============================================================================
      


 
!===============================================================================
      subroutine BFGS_new_Update(n,StpLim,X_1,X_2,G_1,G_2,iH_1,iH_2,X_3)
!===============================================================================
!     Input:
!         n = 3*natom (dimension of the problem)
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!         iH_1 : previous inverse Hessian
!
!     Output:
!         iH_2 : updated approximate inverse Hessian
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      integer i, j, k, n
      real*8, dimension(n)   :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n,n) :: iH_1, iH_2
      
      real*8, dimension(n)   :: dX, dG, ChgeX
      real*8, dimension(n,n) :: XX_diad, XG_diad, GX_diad
      real*8, dimension(n,n) :: Eye, Left_mat, Right_mat, DD
      real*8 :: overlap, dG_sq, dX_sq
      character(len=20) :: StpLim

      dG = G_2 - G_1
      dX = X_2 - X_1

      dG_sq = dot_product(dG,dG)
      dX_sq = dot_product(dX,dX)
      overlap = dot_product(dG,dX)

     !compute |dX><dG| (here Dirac's bra-ket notation is used)
      call diad(dX,dG,n,XG_diad) 

     !compute |dG><dX|
      call diad(dG,dX,n,GX_diad) 

     !compute |dX><dX|
      call diad(dX,dX,n,XX_diad) 

     !Unit matrix
      Eye = 0.d0 
      forall(i=1:n) Eye(i,i) = 1.0d0

     !Left_mat = 1 - |dX><dG| / <G|X>
      Left_mat =  Eye - XG_diad/overlap

     !Right_mat = 1 - |dG><dX| / <G|X>
      Right_mat =  Eye - GX_diad/overlap

      DD = matmul(iH_1,Right_mat)

     !New updated inverse Hessian
      iH_2 = matmul(Left_mat, DD) + XX_diad/overlap


     !***************************
      ChgeX = -matmul(iH_2,G_2)
     !***************************

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX
       
      end subroutine
!===============================================================================

      
!===============================================================================
      subroutine Broyden_Update(n,StpLim,X_1,X_2,G_1,G_2,iH_1,iH_2,X_3)
!===============================================================================
!     Input:
!         n = 3*natom (dimension of the problem)
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!         iH_1 : previous inverse Hessian
!
!     Output:
!         iH_2 : updated approximate inverse Hessian
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      integer i, j, k, n
      real*8, dimension(n)   :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n,n) :: iH_1, iH_2
      
      real*8, dimension(n)   :: dX, dG, ChgeX, iH_dG, YY
      real*8, dimension(n,n) :: YYX_diad, DD 
      real*8 :: overlap
      character(len=20) :: StpLim

      dG = G_2 - G_1
      dX = X_2 - X_1

      iH_dG = matmul(iH_1,dG)

      YY = dX - iH_dG 

      !<dX|iH_1|dG>
      overlap = dot_product(dX,iH_dG)

     !compute |dX - iH_dG >< dX|
      call diad(YY,dX,n,YYX_diad) 

      DD = matmul(YYX_diad,iH_1)

     !New updated inverse Hessian
      iH_2 = iH_1 + matmul(DD,iH_1) / overlap 


     !Change in coordiante
      ChgeX = -matmul(iH_2,G_2)
       

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX
       
      end subroutine
!===============================================================================


!===============================================================================
      subroutine DFP_Update(n,StpLim,X_1,X_2,G_1,G_2,iH_1,iH_2,X_3)
!===============================================================================
!     Input:
!         n = 3*natom (dimension of the problem)
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!         iH_1 : previous inverse Hessian
!
!     Output:
!         iH_2 : updated approximate inverse Hessian
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      integer i, j, k, n
      real*8, dimension(n)   :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n,n) :: iH_1, iH_2
      
      real*8, dimension(n)   :: dX, dG, ChgeX, iH_dG
      real*8, dimension(n,n) :: XX_diad, GG_diad, DD, GGH 
      real*8 :: overlap, ghg
      character(len=20) :: StpLim

      dG = G_2 - G_1
      dX = X_2 - X_1

      iH_dG = matmul(iH_1,dG)

      overlap = dot_product(dX,dG)

      iH_dG = matmul(iH_1,dG)

      ghg = dot_product(dG,iH_dG)

     !compute |dX><dX|
      call diad(dX,dX,n,XX_diad) 

     !compute |dG><dG|
      call diad(dG,dG,n,GG_diad) 

      GGH = matmul(GG_diad,iH_1)

      DD = matmul(iH_1,GGH)

     !New updated inverse Hessian
      iH_2 = iH_1 + XX_diad/overlap - DD/ghg


     !Change in coordiante
      ChgeX = -matmul(iH_2,G_2)
       

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX
       
      end subroutine
!===============================================================================


!===============================================================================
      subroutine SR1_Update(n,StpLim,X_1,X_2,G_1,G_2,iH_1,iH_2,X_3)
!===============================================================================
!     Input:
!         n = 3*natom (dimension of the problem)
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!         iH_1 : previous inverse Hessian
!
!     Output:
!         iH_2 : updated approximate inverse Hessian
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      integer i, j, k, n
      real*8, dimension(n)   :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n,n) :: iH_1, iH_2
      
      real*8, dimension(n)   :: dX, dG, ChgeX, iH_dG, YY
      real*8, dimension(n,n) :: YY_diad, DD 
      real*8 :: overlap
      character(len=20) :: StpLim

      dG = G_2 - G_1
      dX = X_2 - X_1

      iH_dG = matmul(iH_1,dG)

      YY = dX - iH_dG 

      overlap = dot_product(YY,dG)

     !compute |YY >< YY|
      call diad(YY,dX,n,YY_diad) 

     !New updated inverse Hessian
      iH_2 = iH_1 + YY_diad / overlap 

     !Change in coordiante
      ChgeX = -matmul(iH_2,G_2)
       

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX
       
      end subroutine
!===============================================================================



!===============================================================================
      subroutine diad(a0,b0,nmax,c)
!--------------------------------------------------
!     compute outer/diadic/tensor product of two vector
!     c = |a><b|
!--------------------------------------------------
      implicit none 
      integer :: nmax, i, j
      real*8  :: c(nmax,nmax)
      real*8  :: a0(nmax),b0(nmax)

       do i=1,nmax
         do j=1,nmax
            c(i,j) = a0(i)*b0(j)
         enddo
       enddo

      end subroutine
!===============================================================================





!===============================================================================
      subroutine BFGS_J_UpdateX(n,StpLim,X_1,X_2,G_1,G_2,Hi_1,Hi_2,X_3)
!===============================================================================
!     Jeremy Harvey's version of BFGS from the old MECP code
!     and refreshed a little bit
!===============================================================================
!     Input:
!         n = 3*natom (dimension of the problem)
!         nstep : iteration step
!         X_1 and X_2 : previous two coordinates
!         G_1 and G_2 : previous two effective gradients
!         Hi_1 : previous inverse Hessian
!
!     Output:
!         Hi_2 : updated approximate inverse Hessian
!         X_3 : the new updated coordinate towards the MECP
!===============================================================================
      implicit none
      
      ! Specially Adapted BFGS routine from Numerical Recipes
      
      integer i, j, k, n
      real*8, dimension(n) :: X_1, X_2, X_3, G_1, G_2
      real*8, dimension(n) :: DelG, HDelG, ChgeX, DelX, w
      real*8, dimension(n,n) :: Hi_1, Hi_2
      real*8 :: fac, fad, fae, dum1, dum2
      character(len=20) :: StpLim

      DelG = G_2 - G_1
      DelX = X_2 - X_1

      HDelG = matmul(Hi_1,DelG)

      fac = dot_product(DelG,DelX)
      fae = dot_product(DelG,HDelG)

      fac = 1.d0/fac
      fad = 1.d0/fae

      w = fac*DelX - fad*HDelG 

      do i=1,n
         do j=1,n
              dum1 = fac*delx(i)*delx(j)
              dum2 = fad*HDelG(i)*HDelG(j) + fae*w(i)*w(j)

              Hi_2(i,j) = Hi_1(i,j) + dum1 - dum2
         enddo
      enddo


      ChgeX = -matmul(Hi_2,G_2)

      if(StpLim == "StepLim_ON") call Step_Limit(n,ChgeX)

      X_3 = X_2 + ChgeX

      end subroutine
!===============================================================================
 

!===============================================================================
      subroutine TestConv_and_Print(istep, Conv)
!===============================================================================
!     Checks convergence, and updates report file
!     There are 5 criteria for testing convergence
!     They are the same as in Gaussian (Except DeltaE):
!     Av.DeltaX, Max.DeltaX, Av.Grad., Max.Grad., DeltaE
!===============================================================================
      use mainvar
      implicit none
      integer :: istep, Conv
      real*8 :: DE, DXMax, Gmax, GRMS, DXRMS
      real*8, dimension(3*natom) :: DeltaX
      character(len=100) :: form2

      DE = Ea - Eb

      DeltaX = X_3 - X_2

     !Max value of Delta X and the Effective gradient
      DXMax  = maxval(abs(DeltaX))
      Gmax   = maxval(abs(Geff_2))

     !Average norm of the Delta X and the Effective gradient
      GRMS   = norm2(Geff_2)/sqrt(float(3*natom))
      DXRMS  = norm2(DeltaX)/sqrt(float(3*natom))

      Conv = 0

      if((GMax <= TGMax) .and. (GRMS <= TGRMS) .and. (DXMax <= TDXMax) &
                  .and. (DXRMS <= TDXRMS) .and. (abs(DE) <= TDE)) then
        !**********
         Conv = 1
        !**********
      endif


     !----------------------------------------------------------------------
     ! Print Log info into a file and the screen as well
     !----------------------------------------------------------------------
      if(istep == 1) then
        open(18,file="Traj_Opt.xyz")
        open(12,file="Report_File_of_Optimization.log")

        call Print_Header_for_Report(6) 
        call Print_Header_for_Report(12) 
      endif

      form2="(i5,2f20.8,f18.5,2x,4f15.5)"
      form2=trim(form2)

      write(6,form2)  istep,Ea,Eb,DE*tokJmol,DXMax,DXRMS,Gmax,GRMS
      call flush(6)
      write(12,form2) istep,Ea,Eb,DE*tokJmol,DXMax,DXRMS,Gmax,GRMS
      call flush(12)


      Conv = 0
      if((GMax <= TGMax) .and. (GRMS <= TGRMS) .and. (DXMax <= TDXMax) &
                  .and. (DXRMS <= TDXRMS) .and. (DE <= TDE)) then 
        !**********
         Conv = 1
        !**********
      endif
      
     !----------------------------------------------------------------------
     !Print the geometry along the optimization pathway
     !----------------------------------------------------------------------
      call Write_Geometry(18,istep,natom,Atom,Ea,Eb,X_3)

      end subroutine
!===============================================================================


!===============================================================================
      subroutine Print_Header_for_Report(ifile)
!===============================================================================
      use mainvar
      integer :: ifile
      character(len=100) :: form1, form0, form00


      form1="(a5,2a20,a18,2x,4a15)"
      form1=trim(form1)
      form0="(a45,f18.5,2x,4f15.5)"
      form0=trim(form0)
      form00="(a45,a18,2x,4a15)"
      form00=trim(form00)

     !Header for the progress of the calculation
      write(ifile,form0) "Convergence Threshold Values: ", &
                       TDE*tokJmol,TDXMax,TDXRMS,TGMax,TGRMS
      write(ifile,form00)     "                              ", &
                     "  |  ", "  |  ", "  |  ", "  |  ", "  |  "
      write(ifile,form00)     "                              ", &
                    "  V  ", "  V  ", "  V  ", "  V  ", "  V  "
      write(ifile,form1) "Step"," Ea [Multi_A]", " Eb [Multi_B]", &
                     " dE[kJ/mol]" , "Max[dX]", "RMS[dX]",   &
                     "Max[Geff]","RMS[Geff]"


      end subroutine
!===============================================================================



!===============================================================================
      subroutine Write_Geometry(ifile,nstep,natom,Atom,Ea,Eb,X)
!===============================================================================
      implicit none
      integer :: i, jx,jy,jz, natom, ifile, nstep
      character(len=2), dimension(natom) :: Atom
      real*8 :: Ea,Eb
      real*8, dimension(3*natom) ::  X

      write(ifile,*) natom
      write(ifile,*) "nstep: ", nstep, "Energy-A:", Ea, "Energy-B:", Eb 
      do i=1,natom
         jx=3*i-2
         jy=3*i-1
         jz=3*i
         write(ifile,'(a2,5x,3f15.6)') atom(i),X(jx),X(jy),X(jz)
      enddo
      call flush(ifile)

      end subroutine
!===============================================================================



!===============================================================================
      subroutine Calc_Ene_Force(natom,nstep,Fname,Atom,X,Ea,Eb,grad_a,grad_b)
!===============================================================================
      implicit none
      integer :: i,jx,jy,jz,natom,nstep
      character(len=100) :: Fname
      character(len=2), dimension(natom) :: Atom
      real*8 :: Ea,Eb
      real*8, dimension(3*natom) :: grad_a, grad_b, X

      real*8, parameter :: bohr=0.529177d0


     !----------------------------------------------------------
      open(16,file=trim(Fname)//".xyz")
      call Write_Geometry(16,nstep,natom,Atom,99.9d0,99.9d0,X)
      close(16)

      call system("./Run_ForceCalc.sh")

     !----------------------------------------------------------
      open(17,file="force_energy.dat")

     !read the energy and forces from specific interface output
      read(17,*) Ea, Eb
      do i=1,3*natom
         read(17,*) grad_a(i),grad_b(i)
      enddo
      close(17)
     !---------------------------------------------------------

     !Grad in Hartree/Angstrom
      grad_a = grad_a/bohr
      grad_b = grad_b/bohr

      end subroutine
!==================================================================
      






!===============================================================================
      subroutine Create_QC_interface
!===============================================================================
!     Preparation of the corresponding BASH script to run QC
!===============================================================================
      use mainvar
      implicit none
      character(len=3) :: stpctrl

      if(StpLim == "StepLim_ON") then
        stpctrl = "on"
      else
        stpctrl = "off"
      endif


      open(12,file="Run_ForceCalc.sh")        
      open(13,file=Trim(Fname)//".inp")        

     !Create BASH script to run QC and extract energies and gradients
      call Orca_BASH(12)
      call system("chmod +x Run_ForceCalc.sh")

     !Create Orca Input file
      call Orca_Input(13)

      close(12)
      close(13)

      write(6,*)
      write(6,"(a)") "+++++++++++++++++++++++++++++++++++++++++++++++++"
      write(6,"(a)") " Quantum Chem Interface is successfully created"
      write(6,"(a)") "-------------------------------------------------"
      write(6,"(2a)")       " Program:     ",trim(QC_interface)
      write(6,"(2a)")       " Method:      ",trim(Method)
      write(6,"(2a)")       " Basis:       ",trim(Basis)
      write(6,"(a,i0)")     " Nprocs:      ",nproc
      write(6,"(a,i0)")     " Memory[MB]:  ",mem
      write(6,"(a,i0)")     " Charge:      ",charge
      write(6,"(a,i0,a,i0)")" Multi[A,B]:  ",multi(1)," ", multi(2)
      write(6,"(2a)")       " Optimizer:   ",trim(OptMethod)
      write(6,"(2a)")       " Step Limit:  ",trim(stpctrl)
      write(6,"(a,i0)")     " Maxiter:     ",maxstep
      write(6,*)        
      write(6,"(a)")       " Available Optimizers: "
      write(6,"(a)")       " BBGrad1, BBGrad2, BFGS, DFP, SR1, Broyden" 
      write(6,"(a)") "+++++++++++++++++++++++++++++++++++++++++++++++++"
        

      if(trim(ConvPar) == "New_Conv_Params") then
        write(6,*) 
        write(6,"(a)")    " !!! Convergence parameters are changed !!!" 
        write(6,"(a)")    "               Are you sure?               "
        write(6,*)         
        write(6,"(a)")    " Default threshold values for convergence &
                                            (in Hartree/Angstrom): "
        write(6,"(a)")    " (Same as in Gaussian except DeltaE) " 
        write(6,"(a)")        
        write(6,"(a)")    "   Delta[E]  = 5.0E-5  (0.13128 kJ/mol)"
        write(6,"(a)")    "   Max[dX]   = 4.0E-3   RMS[dX]   = 2.5E-3"
        write(6,"(a)")    "   Max[Grad] = 7.0E-4   RMS[Grad] = 5.0E-4"
        write(6,*)
        write(6,"(a)")    " Instead of these you use:"
      endif

      write(6,*)

      end subroutine
!===============================================================================
 
!===============================================================================
      subroutine Orca_Input(ifl)
!===============================================================================
!     nproc    = number of processors
!     mem      = memory in MB
!     Method   = e.g.: Method = "UKS B3LYP engrad", or Method = "MP2 engrad"
!     Basis    = basis set
!     charge   = charge of the chemical moetiy
!     multi    = spin multiplicity vector for the two statates
!     Fname    = name of file
!     ifl      = number of file, where to write
!===============================================================================
      use mainvar
      implicit none
      integer :: ifl

      write(ifl,*) "%MaxCore ", mem
      write(ifl,*) "%pal nprocs ", nproc, " end"
      write(ifl,*) 
      write(ifl,*) "! ",trim(Method)//" ", trim(Basis)
      write(ifl,*) 
      write(ifl,*) "* xyzfile ",charge,multi(1),"  ",trim(Fname)//".xyz"
      write(ifl,*) 
      write(ifl,*) 
      write(ifl,*) 
      write(ifl,*) "$new_job" 
      write(ifl,*)
      write(ifl,*) "! ",trim(Method)//" ", trim(Basis)
      write(ifl,*)
      write(ifl,*) "* xyzfile ",charge,multi(2),"  ",trim(Fname)//".xyz"

      end subroutine
!===============================================================================



 
!===============================================================================
      subroutine Orca_BASH(ifl)
!===============================================================================
!     Fname   = name of file
!     QC_rout = pathway (exact link) for Orca, e.g.: /home/$USER/Orca-5.03/orca
!===============================================================================
      use mainvar
      implicit none
      integer :: ifl

      write(ifl,*) "#!/bin/bash"

      write(ifl,*) "fname=",Trim(Fname)

      write(ifl,*) "orca_link=",trim(QC_route)
      write(ifl,"(a13,i0)") "natom3_plus1=",(3*natom)+1
      write(ifl,"(a13,i0)") "natom3_plus2=",(3*natom)+2

      write(ifl,*)

      write(ifl,*) "rm -rf orca_tmp"

      write(ifl,*) "mkdir orca_tmp"

      write(ifl,*) "cp $fname.xyz $fname.inp orca_tmp"

      write(ifl,*) "cd orca_tmp"

      write(ifl,*)

      write(ifl,*) "$orca_link  $fname.inp > $fname.out"

      write(ifl,*)

      write(ifl,*) 'grep -A2 "The current total energy" ', &
                       "$fname.engrad|tail -1 > a1.txt"
      write(ifl,*) 'grep -A$natom3_plus2 "The current gradient" ', &
                        "$fname.engrad|tail -$natom3_plus1 > a2.txt"

      write(ifl,*) "cat a1.txt a2.txt > force1.dat"

      write(ifl,*)

      write(ifl,*) 'grep -A2 "The current total energy" ', &
                      "${fname}_job2.engrad|tail -1 > b1.txt"
      write(ifl,*) 'grep -A$natom3_plus2 "The current gradient" ', &
                  "${fname}_job2.engrad|tail -$natom3_plus1 > b2.txt"

      write(ifl,*)

      write(ifl,*) "cat b1.txt b2.txt > force2.dat"

      write(ifl,*) "paste force1.dat force2.dat > force_energy.dat"

      write(ifl,*)

      write(ifl,*) "rm -rf a1.txt a2.txt b1.txt b2.txt"

      write(ifl,*) "cp force_energy.dat ../"

      write(ifl,*) "cd .."

      write(ifl,*) "rm -rf orca_tmp"


      end subroutine
!===============================================================================

      
