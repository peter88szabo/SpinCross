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
 
