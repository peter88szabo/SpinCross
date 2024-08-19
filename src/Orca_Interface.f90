
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

      
