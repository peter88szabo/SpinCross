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

