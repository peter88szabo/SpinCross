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

