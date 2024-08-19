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
 
