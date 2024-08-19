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
