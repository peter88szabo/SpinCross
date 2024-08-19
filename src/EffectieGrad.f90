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



      
