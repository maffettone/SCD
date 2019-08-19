module sidechain_decoration

  ! Module to contain the side chain decoration and manipulation
  ! MODIFIED 25/08/2017 by PMM
  use types
  use maths
  use class_proteins
  use random
  use lib_dir

  implicit none

  private
  public :: decorate, calc_chi

  type fragment
     character(len=3) :: res_label                          !label for residue of interest
     type(atom),allocatable :: atoms(:)                     !atom list for sidechain fragment
     type(atom) :: Ca                                       !Carbon alpha for reference origin
  end type fragment

  type chi_rel
     character(len=4) :: chain(4,5)                         !chain of atom labels for dihedral
     character(len=4) :: axis(2,5)                          !Axis for dihedral rotation
     logical :: chi_on(5)                                   !Logical for chi num relative to res
  end type chi_rel

contains

  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: decorate(protein)
  ! SUBROUTINE: calc_chi(peptide_chain,chi_relations)
  ! SUBROUTINE: find_clash(peptide_chain,adjacency_matrix,vanderwalls_symbols,vdw_radii)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! SUBROUTINE: import_fragments(fragments)
  ! SUBROUTINE: add_sc(peptide_chain)
  ! SUBROUTINE: set_chi(peptide_chain,residue_index, chi_index,chi_value,chi_relations)
  ! SUBROUTINE: build_kDOP(peptide_chain,adjacency_mat,chi_relations)
  ! SUBROUTINE: relax_sc(peptide_chain,adjacency_matrix)
  ! SUBROUTINE: chi_relations(peptide_chain, chi_rel)
  ! FUNCTION: calcE(peptide_chain,adjacency_matrix,vanderwalls_symbols,vdw_radii)
  ! ----------------------------------------------------------------------------

  subroutine decorate(prot)
    ! SUBROUTINE: decorate(protein)
    ! Calling routine to add sidechains to protein if not already,
    ! or to adjust sidechains according to kDOPs
    type(protein), intent(inout) :: prot                    !Protein
    integer,allocatable :: adj(:,:)                         !Adjacency matrix
    integer, allocatable :: focus(:)                        !List of focus indicies
    type(chi_rel),allocatable :: chi_rels(:)                !Chi relations aligned with pep
    logical :: success                                      !Check whether successful MC
    integer :: maxmoves                                     !Maximum MC moves for relaxation
    integer :: max_loop=10                                  !Maximum power of 2, times 100 sweeps
    integer :: i

    success = .false.
    call init_random_seed()
    allocate(chi_rels(size(prot%pep)))
    call chi_relations(prot%pep,chi_rels)
    if (prot%bbonly) then
       call add_sc(prot%pep)
    end if
    call calc_chi(prot%pep,chi_rels)
    call build_kDOP(prot%pep,adj,chi_rels)
    maxmoves = 100*prot%n_res
    write(*,*)
    do i=1,max_loop
       call relax_sc(prot%pep,adj,chi_rels,maxmoves,success,focus)
       if(success) then
          exit
       else
          maxmoves=2*maxmoves
          write(*,*)'Number of moves in run: ', maxmoves
       end if
    end do
    if (.not. success) then
       maxmoves = 100*prot%n_res
       call focus_relax_sc(prot%pep,adj,chi_rels,maxmoves,success,focus)
    end if
  end subroutine decorate
  ! -----------------------------------------------------------------------------

  subroutine relax_sc(pep,adj,chi_rels,maxmoves,success,focus)
    ! SUBROUTINE: relax_sc(polypeptide_chain,adjacency_matrix)
    ! Instead of an MC, this permuts all of the possible options for rotomers (g+,g-,t)
    ! and stops when a sensible configuration in met
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    type(residue),intent(inout) :: pep(:)                   !Peptide chain
    type(chi_rel),intent(in) :: chi_rels(:)                 !Chi relations to pass to set_chi()
    integer, intent(in) :: maxmoves                         !Number of MC moves
    logical,optional,intent(out) :: success                 !Check whether completed at 0 energy
    integer, allocatable,intent(inout) :: focus(:)          !List of focus indicies
    real(dp) :: rotamers(3)                                 !Set of angles for rotomers to check
    real(dp),allocatable :: vdwr(:)                         !van der waals radii
    character(len=2),allocatable :: vdwsym(:)               !van der waals symbols
    real(dp) :: r                                           !random number
    real(dp) :: total_E                                     !Total Energy
    real(dp) :: dE                                          !Change in energy
    real(dp) :: chi_prev                                    !Previous chi, saved for reversion
    real(dp) :: kT                                          !Temperature
    real(dp) :: dT                                          !Change in Temperature
    real(dp) :: fudge=1.4                                   !Radius fudge factor
    integer :: res_idx                                      !Residue index
    integer :: chi_idx                                      !Dihedral index
    integer :: info
    integer :: unit
    character(len=80) :: errmsg
    integer :: i,k

    !------------------------------------ PARAMETER SET----------------------------------------
    rotamers = [-60._dp,180._dp,60._dp]
    call init_random_seed()
    !------------------------------------ PARAMETER SET----------------------------------------
    ! Pulling in cov radii instead van der Waals radii
    info=0
    open(newunit=unit, file=trim(cov_lib), status='old',IOSTAT = info,IOMSG=errmsg)
    if(info /=0)then
       write(*,*)errmsg
       stop 'Failure in scd.f90'
    end if
    allocate(vdwr(10),vdwsym(10))
    read(unit,*);read(unit,*)
    do i=1,size(vdwr)
       read(unit,*)vdwsym(i),vdwr(i)
    end do
    vdwr = vdwr*fudge
    close(unit)
    total_E = calcE(pep,adj,vdwsym,vdwr)
    kT = total_E*2.0_dp
    dT = ((10.**-6)/kT)**(1./maxmoves)
    do i=1,maxmoves
       do
          ! Only bothering with residues which have sidechains
          call random_number(r)
          res_idx = ceiling(r*size(pep))
          if(size(pep(res_idx)%chi)==0)cycle
          call random_number(r)
          chi_idx = ceiling(r*size(pep(res_idx)%chi))
          exit
       end do
       call random_number(r)
       k = ceiling(r*size(rotamers))
       chi_prev = pep(res_idx)%chi(chi_idx)
       dE = localE(pep,adj,res_idx,vdwsym,vdwr)
       call set_chi(pep,res_idx,chi_idx,rotamers(k),chi_rels)
       dE = localE(pep,adj,res_idx,vdwsym,vdwr) -dE
       call random_number(r)
       if (r < exp(-1._dp*dE/kT)) then
          total_E = total_E + dE
          if (total_E<10.**(-3)) then
             write(*,*)i
             exit
          end if
       else
          call set_chi(pep,res_idx,chi_idx,chi_prev,chi_rels)
       end if
       kT = kT * dT
    end do

    if (total_E<10.**(-3) .and. present(success)) then
       success = .true.
    else if(present(success)) then
       success = .false.
       call find_clash(pep,adj,vdwsym,vdwr)
       call find_focus(pep,adj,vdwsym,vdwr,focus)
    end if
  end subroutine relax_sc
  ! -----------------------------------------------------------------------------


  subroutine focus_relax_sc(pep,adj,chi_rels,maxmoves,success,focus)
    ! SUBROUTINE: relax_sc(polypeptide_chain,adjacency_matrix)
    ! Instead of an MC, this permuts all of the possible options for rotomers (g+,g-,t)
    ! and stops when a sensible configuration in met
    ! ONLY RUNS ON INPUT FOCUS SET
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    type(residue),intent(inout) :: pep(:)                   !Peptide chain
    type(chi_rel),intent(in) :: chi_rels(:)                 !Chi relations to pass to set_chi()
    integer, intent(in) :: maxmoves                         !Number of MC moves
    integer, intent(in) :: focus(:)                         !Residue indicies to focus on
    logical,optional,intent(out) :: success                 !Check whether completed at 0 energy
    real(dp) :: rotamers(3)                                 !Set of angles for rotomers to check
    real(dp),allocatable :: vdwr(:)                         !van der waals radii
    character(len=2),allocatable :: vdwsym(:)               !van der waals symbols
    real(dp) :: r                                           !random number
    real(dp) :: total_E                                     !Total Energy
    real(dp) :: dE                                          !Change in energy
    real(dp) :: chi_prev                                    !Previous chi, saved for reversion
    real(dp) :: kT                                          !Temperature
    real(dp) :: dT                                          !Change in Temperature
    real(dp) :: fudge=1.4                                   !Radius fudge factor
    integer :: res_idx                                      !Residue index
    integer :: chi_idx                                      !Dihedral index
    integer :: info
    integer :: unit
    character(len=80) :: errmsg
    integer :: i,k

    !------------------------------------ PARAMETER SET----------------------------------------
    rotamers = [-60._dp,180._dp,60._dp]
    call init_random_seed()
    !------------------------------------ PARAMETER SET----------------------------------------
    ! Pulling in cov radii instead van der Waals radii
    info=0
    open(newunit=unit, file=trim(cov_lib), status='old',IOSTAT = info,IOMSG=errmsg)
    if(info /=0)then
       write(*,*)errmsg
       stop 'Failure in scd.f90'
    end if
    allocate(vdwr(10),vdwsym(10))
    read(unit,*);read(unit,*)
    do i=1,size(vdwr)
       read(unit,*)vdwsym(i),vdwr(i)
    end do
    vdwr = vdwr*fudge
    close(unit)
    total_E = calcE(pep,adj,vdwsym,vdwr)
    kT = total_E*10.0_dp
    dT = ((10.**-6)/kT)**(1./maxmoves)
    do i=1,maxmoves
       call random_number(r)
       res_idx = ceiling(r*size(focus))
       res_idx = focus(res_idx)
       chi_idx = ceiling(r*size(pep(res_idx)%chi))
       call random_number(r)
       k = ceiling(r*size(rotamers))
       chi_prev = pep(res_idx)%chi(chi_idx)
       dE = localE(pep,adj,res_idx,vdwsym,vdwr)
       call set_chi(pep,res_idx,chi_idx,rotamers(k),chi_rels)
       dE = localE(pep,adj,res_idx,vdwsym,vdwr) -dE
       call random_number(r)
       if (r < exp(-1._dp*dE/kT)) then
          total_E = total_E + dE
          if (total_E<10.**(-3)) then
             write(*,*)i
             exit
          end if
       else
          call set_chi(pep,res_idx,chi_idx,chi_prev,chi_rels)
       end if
       kT = kT * dT
    end do

    if (total_E<10.**(-3) .and. present(success)) then
       success = .true.
    else if(present(success)) then
       success = .false.
       call find_clash(pep,adj,vdwsym,vdwr)
    end if
  end subroutine focus_relax_sc
  ! -----------------------------------------------------------------------------

  function calcE(pep,adj,vdwsym,vdwr) result(r)
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    type(residue),intent(in) :: pep(:)                      !Peptide chain
    real(dp) :: vdwr(:)                                     !van der waals radii
    character(len=2) :: vdwsym(:)                           !van der waals symbols
    real(dp) :: r                                           !return value of total E
    real(dp) :: r1,r2                                       !Radii of interest
    integer :: i,j,k,l,m

    r = 0._dp
    do i=1,size(adj,1)
       do j=i,size(adj,2)
          if(adj(j,i) > 0) then
             do k=1,size(pep(j)%backbone)
                do l=1,size(vdwr)
                   if(trim(pep(j)%backbone(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%backbone(k)%coords)**2) < (r1+r2)**2) then
                      r = r+1.0
                   end if
                end do
             end do
             do k=1,size(pep(j)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(j)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%backbone)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%backbone(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%backbone(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      r = r+1.0
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      r = r+1.0
                   end if
                end do
             end do
          end if
       end do
    end do
  end function calcE
  ! -----------------------------------------------------------------------------

  subroutine find_clash(pep,adj,vdwsym,vdwr)
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    type(residue),intent(in) :: pep(:)                      !Peptide chain
    real(dp) :: vdwr(:)                                     !van der waals radii
    character(len=2) :: vdwsym(:)                           !van der waals symbols
    real(dp) :: E                                           !Energy term
    real(dp) :: r1,r2                                       !Radii of interest
    integer :: i,j,k,l,m


    write(*,"('WARNING: Decorated protein may have clashes. Check carefully.')")

    do i=1,size(adj,1)
       do j=i,size(adj,2)
          if(adj(j,i) > 0) then
             do k=1,size(pep(j)%backbone)
                do l=1,size(vdwr)
                   if(trim(pep(j)%backbone(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%backbone(k)%coords)**2) < (r1+r2)**2) then
                      write(*,"('Sidechain: ',I3,' ',A3,' ',A4,'   Backbone: ',I3,' ',A3,' ',A4)") &
                           i,pep(i)%label,pep(i)%sidechain(m)%label,j,pep(j)%label,pep(j)%backbone(k)%label
                      write(*,"('Minimum Distance: ',F5.3,'   Actual Distance: ',F5.3)") &
                           (r1+r2),sqrt(sum((pep(i)%sidechain(m)%coords-pep(j)%backbone(k)%coords)**2))
                   end if
                end do
             end do
             do k=1,size(pep(j)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(j)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%backbone)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%backbone(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%backbone(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      write(*,"('Backbone: ',I3,' ',A3,' ',A4,'   Sidechain: ',I3,' ',A3,' ',A4)") &
                           i,pep(i)%label,pep(i)%backbone(m)%label,j,pep(j)%label,pep(j)%sidechain(k)%label
                      write(*,"('Minimum Distance: ',F5.3,'   Actual Distance: ',F5.3)") &
                           (r1+r2),sqrt(sum((pep(i)%backbone(m)%coords-pep(j)%sidechain(k)%coords)**2))
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      write(*,"('Sidechain: ',I3,' ',A3,' ',A4,'   Sidechain: ',I3,' ',A3,' ',A4)") &
                           i,pep(i)%label,pep(i)%sidechain(m)%label,j,pep(j)%label,pep(j)%sidechain(k)%label
                      write(*,"('Minimum Distance: ',F5.3,'   Actual Distance: ',F5.3)") &
                           (r1+r2),sqrt(sum((pep(i)%sidechain(m)%coords-pep(j)%sidechain(k)%coords)**2))
                   end if
                end do
             end do
          end if
       end do
    end do
    write(*,*)
  end subroutine find_clash
  ! -----------------------------------------------------------------------------

  subroutine find_focus(pep,adj,vdwsym,vdwr,focus)
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    type(residue),intent(in) :: pep(:)                      !Peptide chain
    integer, allocatable,intent(inout) :: focus(:)          !List of focus indicies
    real(dp) :: vdwr(:)                                     !van der waals radii
    character(len=2) :: vdwsym(:)                           !van der waals symbols
    real(dp) :: E                                           !Energy term
    real(dp) :: r1,r2                                       !Radii of interest
    integer :: i,j,k,l,m
    integer :: count

    if (allocated(focus)) then
       deallocate(focus)
    end if
    count=0
    do i=1,size(adj,1)
       do j=i,size(adj,2)
          if(adj(j,i) > 0) then
             do k=1,size(pep(j)%backbone)
                do l=1,size(vdwr)
                   if(trim(pep(j)%backbone(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%backbone(k)%coords)**2) < (r1+r2)**2) then
                      count=count+2
                   end if
                end do
             end do
             do k=1,size(pep(j)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(j)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%backbone)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%backbone(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%backbone(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      count=count+2
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      count=count+2
                   end if
                end do
             end do
          end if
       end do
    end do

    allocate(focus(count))
    count=1
    do i=1,size(adj,1)
       do j=i,size(adj,2)
          if(adj(j,i) > 0) then
             do k=1,size(pep(j)%backbone)
                do l=1,size(vdwr)
                   if(trim(pep(j)%backbone(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%backbone(k)%coords)**2) < (r1+r2)**2) then
                      focus(count) = i
                      count=count+1
                      focus(count) =j
                      count=count+1
                   end if
                end do
             end do
             do k=1,size(pep(j)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(j)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r1 = vdwr(l)
                      exit
                   end if
                end do
                do m=1,size(pep(i)%backbone)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%backbone(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%backbone(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      focus(count) = i
                      count=count+1
                      focus(count) =j
                      count=count+1
                   end if
                end do
                do m=1,size(pep(i)%sidechain)
                   do l=1,size(vdwr)
                      if(trim(pep(i)%sidechain(m)%sym)==trim(vdwsym(l))) then
                         r2 = vdwr(l)
                         exit
                      end if
                   end do
                   if(sum((pep(i)%sidechain(m)%coords-pep(j)%sidechain(k)%coords)**2) < (r1+r2)**2) then
                      focus(count) = i
                      count=count+1
                      focus(count) =j
                      count=count+1
                   end if
                end do
             end do
          end if
       end do
    end do
  end subroutine find_focus
  ! -----------------------------------------------------------------------------

  
  function localE(pep,adj,idx,vdwsym,vdwr) result(r)
    integer, intent(in) :: adj(:,:)                         !Adjacency matrix
    integer, intent(in) :: idx                              !residue of interest
    type(residue),intent(in) :: pep(:)                      !Peptide chain
    real(dp) :: vdwr(:)                                     !van der waals radii
    character(len=2) :: vdwsym(:)                           !van der waals symbols
    real(dp) :: r                                           !return value of total E
    real(dp) :: r1,r2                                       !Radii of interest
    integer :: i,j,k,l

    r = 0._dp

    do i=1,size(adj,2)
       if(adj(i,idx)>0) then
          do j=1,size(pep(idx)%backbone)
             do l=1,size(vdwr)
                if(trim(pep(idx)%backbone(j)%sym)==trim(vdwsym(l))) then
                   r1 = vdwr(l)
                   exit
                end if
             end do
             do k=1,size(pep(i)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(i)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r2 = vdwr(l)
                      exit
                   end if
                end do
                if(sum((pep(i)%sidechain(k)%coords-pep(idx)%backbone(j)%coords)**2) < (r1+r2)**2) then
                   r = r+1.0
                end if
             end do
          end do
          do j=1,size(pep(idx)%sidechain)
             do l=1,size(vdwr)
                if(trim(pep(idx)%sidechain(j)%sym)==trim(vdwsym(l))) then
                   r1 = vdwr(l)
                   exit
                end if
             end do
             do k=1,size(pep(i)%backbone)
                do l=1,size(vdwr)
                   if(trim(pep(i)%backbone(k)%sym)==trim(vdwsym(l))) then
                      r2 = vdwr(l)
                      exit
                   end if
                end do
                if(sum((pep(i)%backbone(k)%coords-pep(idx)%sidechain(j)%coords)**2) < (r1+r2)**2) then
                   r = r+1.0
                end if
             end do
             do k=1,size(pep(i)%sidechain)
                do l=1,size(vdwr)
                   if(trim(pep(i)%sidechain(k)%sym)==trim(vdwsym(l))) then
                      r2 = vdwr(l)
                      exit
                   end if
                end do
                if(sum((pep(i)%sidechain(k)%coords-pep(idx)%sidechain(j)%coords)**2) < (r1+r2)**2) then
                   r = r+1.0
                end if
             end do
          end do
       end if
    end do
  end function localE
  ! -----------------------------------------------------------------------------

  
  subroutine build_kDOP(pep,adj,chi_rels)
    ! SUBROUTINE: build_kDOP(peptide_chain,adjacency_mat,chi_relations)
    ! Makes 7 dimensional oriented polytopes which encompass all rotomer permutations
    ! of a sidechain, and then calculates adjacency matrix for inclusion in repulsive calc
    integer,allocatable,intent(inout) :: adj(:,:)           !Adjacency matrix
    type(residue),intent(inout) :: pep(:)                   !Peptide chain
    type(chi_rel),intent(in) :: chi_rels(:)                 !Chi relations to pass to set_chi()
    real(dp),allocatable :: kDOP(:,:,:)                     !kDOP for each sidechain
    real(dp) :: rotamers(3)                                 !Set of angles to check for rotomers
    real(dp),allocatable :: vdwr(:)                         !van der waals radii
    character(len=2),allocatable :: vdwsym(:)               !van der waals symbols
    real(dp) :: basis(3,7)                                  !basis vectors for kDOP
    real(dp) :: xmin,xmax                                   !Dimensional minimum and maximum
    integer :: info
    integer :: unit
    character(len=80) :: errmsg
    integer :: i,j,j1,j2,j3,j4,j5,k,l

    rotamers = [-60._dp,180._dp,60._dp]
    allocate(adj(size(pep),size(pep)))
    allocate(kDOP(2,7,size(pep)))

    ! VDW radii for larger than needed kDOP
    info=0
    open(newunit=unit, file=trim(vdw_lib), status='old',IOSTAT = info,IOMSG=errmsg)
    if(info /=0)then
       write(*,*)errmsg
       stop 'Failure in scd.f90'
    end if
    allocate(vdwr(10),vdwsym(10))
    read(unit,*);read(unit,*)
    do i=1,size(vdwr)
       read(unit,*)vdwsym(i),vdwr(i)
    end do
    close(unit)



!     ! SANITY CHECK TEST
!     do k=1,size(pep)
!        if (trim(pep(k)%label) =='PRO')then
!           call set_chi(pep,k,1,0._dp,chi_rels)
!           call set_chi(pep,k,2,0._dp,chi_rels)
!           open(newunit=unit, file='/Users/alggroup/Desktop/testing.pdb', status='replace',IOSTAT = info,IOMSG=errmsg)
! 15        format(A6,I5,T13,A4,A1,A3,T22,A1,I4,A1,T31,3F8.3,2F6.2,T77,A2,A2)
!           j=1
!           do i=1,size(pep(k)%backbone)
!              write(unit,15)'ATOM  ',j,pep(k)%backbone(i)%label,' ',pep(k)%label,'A',1,' ',&
!                   pep(k)%backbone(i)%coords,1._dp,0._dp,pep(k)%backbone(i)%sym,'  '
!              j=j+1
!           end do
!           do i=1,size(pep(k)%sidechain)
!              write(unit,15)'ATOM  ',j,pep(k)%sidechain(i)%label,' ',pep(k)%label,'A',1,' ',&
!                   pep(k)%sidechain(i)%coords,1._dp,0._dp,pep(k)%sidechain(i)%sym,'  '
!              j=j+1
!           end do
!           stop
!        end if
!     end do
!     ! SANITY CHECK TEST



    
    ! Setting up normal basis vectors for 7-DOP
    basis(:,1) = [0._dp,0._dp,1._dp]
    basis(:,2) = [0._dp,1._dp,0._dp]
    basis(:,3) = [1._dp,0._dp,0._dp]
    basis(:,4) = sqrt(2._dp)*(basis(:,1) + basis(:,2) + basis(:,3))
    basis(:,5) = sqrt(2._dp)*(basis(:,1) - basis(:,2) + basis(:,3))
    basis(:,6) = sqrt(2._dp)*(basis(:,1) + basis(:,2) - basis(:,3))
    basis(:,7) = sqrt(2._dp)*(basis(:,1) - basis(:,2) - basis(:,3))
    do i=1,7
       basis(:,i) = basis(:,i)/norm2(basis(:,i))
    end do

    ! Calculating starting with backbone for all
    ! Updates min/max atom by atom
    do i=1,size(pep)
       if (pep(i)%heteroatom) cycle
       do j=1,size(pep(i)%backbone)
          ! Finding correct radius
          do k=1,size(vdwr)
             if(trim(pep(i)%backbone(j)%sym)==trim(vdwsym(k))) then
                l=k
                exit
             end if
          end do
          ! Projecting radius and adjusting kDOP
          do k=1,size(basis,2)
             xmin = dot_product(basis(:,k),pep(i)%backbone(j)%coords) - vdwr(l)
             xmax = xmin+2*vdwr(l)
             if(j==1) then
                kDOP(1,k,i) = xmin
                kDOP(2,k,i) = xmax
                cycle
             end if
             if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
             if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
          end do
       end do
    end do

    do i=1,size(pep)
       select case(size(pep(i)%chi))
       case(0)
          cycle
       case(1)
          do j1=1,size(rotamers)
             call set_chi(pep,i,1,rotamers(j1),chi_rels)
             do j=1,size(pep(i)%sidechain)
                do k=1,size(vdwr)
                   if(trim(pep(i)%sidechain(j)%sym)==trim(vdwsym(k))) then
                      l=k
                      exit
                   end if
                end do
                do k=1,size(basis,2)
                   xmin = dot_product(basis(:,k),pep(i)%sidechain(j)%coords) - vdwr(l)
                   xmax = xmin+2*vdwr(l)
                   if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
                   if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
                end do
             end do
          end do
       case(2)
          do j1=1,size(rotamers)
             call set_chi(pep,i,1,rotamers(j1),chi_rels)
             do j2 =1,size(rotamers)
                call set_chi(pep,i,2,rotamers(j2),chi_rels)
                do j=1,size(pep(i)%sidechain)
                   do k=1,size(vdwr)
                      if(trim(pep(i)%sidechain(j)%sym)==trim(vdwsym(k))) then
                         l=k
                         exit
                      end if
                   end do
                   do k=1,size(basis,2)
                      xmin = dot_product(basis(:,k),pep(i)%sidechain(j)%coords) - vdwr(l)
                      xmax = xmin+2*vdwr(l)
                      if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
                      if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
                   end do
                end do
             end do
          end do
       case(3)
          do j1=1,size(rotamers)
             call set_chi(pep,i,1,rotamers(j1),chi_rels)
             do j2 =1,size(rotamers)
                call set_chi(pep,i,2,rotamers(j2),chi_rels)
                do j3 =1,size(rotamers)
                   call set_chi(pep,i,3,rotamers(j3),chi_rels)
                   do j=1,size(pep(i)%sidechain)
                      do k=1,size(vdwr)
                         if(trim(pep(i)%sidechain(j)%sym)==trim(vdwsym(k))) then
                            l=k
                            exit
                         end if
                      end do
                      do k=1,size(basis,2)
                         xmin = dot_product(basis(:,k),pep(i)%sidechain(j)%coords) - vdwr(l)
                         xmax = xmin+2*vdwr(l)
                         if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
                         if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
                      end do
                   end do
                end do
             end do
          end do
       case(4)
          do j1=1,size(rotamers)
             call set_chi(pep,i,1,rotamers(j1),chi_rels)
             do j2 =1,size(rotamers)
                call set_chi(pep,i,2,rotamers(j2),chi_rels)
                do j3 =1,size(rotamers)
                   call set_chi(pep,i,3,rotamers(j3),chi_rels)
                   do j4 =1,size(rotamers)
                      call set_chi(pep,i,4,rotamers(j4),chi_rels)
                      do j=1,size(pep(i)%sidechain)
                         do k=1,size(vdwr)
                            if(trim(pep(i)%sidechain(j)%sym)==trim(vdwsym(k))) then
                               l=k
                               exit
                            end if
                         end do
                         do k=1,size(basis,2)
                            xmin = dot_product(basis(:,k),pep(i)%sidechain(j)%coords) - vdwr(l)
                            xmax = xmin+2*vdwr(l)
                            if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
                            if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
                         end do
                      end do
                   end do
                end do
             end do
          end do
       case(5)
          do j1=1,size(rotamers)
             call set_chi(pep,i,1,rotamers(j1),chi_rels)
             do j2 =1,size(rotamers)
                call set_chi(pep,i,2,rotamers(j2),chi_rels)
                do j3 =1,size(rotamers)
                   call set_chi(pep,i,3,rotamers(j3),chi_rels)
                   do j4 =1,size(rotamers)
                      call set_chi(pep,i,4,rotamers(j4),chi_rels)
                      do j5 =1,size(rotamers)
                         call set_chi(pep,i,5,rotamers(j5),chi_rels)
                         do j=1,size(pep(i)%sidechain)
                            do k=1,size(vdwr)
                               if(trim(pep(i)%sidechain(j)%sym)==trim(vdwsym(k))) then
                                  l=k
                                  exit
                               end if
                            end do
                            do k=1,size(basis,2)
                               xmin = dot_product(basis(:,k),pep(i)%sidechain(j)%coords) - vdwr(l)
                               xmax = xmin+2*vdwr(l)
                               if(xmin < kDOP(1,k,i)) kDOP(1,k,i) = xmin
                               if(xmax > kDOP(2,k,i)) kDOP(2,k,i) = xmax
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end select
    end do

    ! Buildind adjacency matrix
    adj = 1
    do i=1,size(pep)
       do j=1,size(pep)
          do k=1,size(kDOP,2)
             if(kDOP(1,k,i)>kDOP(2,k,j) .or. kDOP(2,k,i)<kDOP(1,k,j)) then
                adj(j,i) = 0
                exit
             end if
          end do
       end do
    end do

    do i=1,size(adj,2)
       adj(i,i)=0
    end do
    ! SANITY CHECK TEST
    ! do i=1,size(pep)
    !    write(*,'(27i3)') adj(:,i)
    ! end do
    ! open(newunit=unit, file='/Users/alggroup/Desktop/testing.pdb', status='replace',IOSTAT = info,IOMSG=errmsg)
    ! 15  format(A6,I5,T13,A4,A1,A3,T22,A1,I4,A1,T31,3F8.3,2F6.2,T77,A2,A2)
    ! j=1
    ! k=5
    ! write(*,*)pep(k)%label
    ! do i=1,size(pep(k)%backbone)
    !    write(unit,15)'ATOM  ',j,'CA',' ',pep(k)%label,'A',1,' ',pep(k)%backbone(i)%coords,1._dp,0._dp,'C ','  '
    !    j=j+1
    ! end do
    ! do i=1,size(pep(k)%sidechain)
    !    write(unit,15)'ATOM  ',j,'CA',' ',pep(k)%label,'A',1,' ',pep(k)%sidechain(i)%coords,1._dp,0._dp,'C ','  '
    !    j=j+1
    ! end do
    ! do i=1,size(basis,2)
    !    write(unit,15)'ATOM  ',j,'O ',' ','XXX','A',1,' ',basis(:,i)*kDOP(1,i,k),1._dp,0._dp,'O ','  '
    !    j=j+1
    ! end do
    ! do i=1,size(basis,2)
    !    write(unit,15)'ATOM  ',j,'N ',' ','XXX','A',1,' ',basis(:,i)*kDOP(2,i,k),1._dp,0._dp,'N ','  '
    !    j=j+1
    ! end do
    ! stop
    ! SANITY CHECK TEST
  end subroutine build_kDOP
  ! -----------------------------------------------------------------------------

  subroutine calc_chi(pep,chi_rels)
    ! SUBROUTINE: calc_chi(peptide_chain,chi_relations)
    type(residue),intent(inout) :: pep(:)                   !Polypeptide chain
    type(chi_rel),intent(in) :: chi_rels(:)                 !Chi relations
    real(dp) :: chain(3,4)                                  !Chain of xyz coordinates
    real(dp) ::v1(3),v2(3),v3(3)                            !Set of unit normals
    real(dp) :: x,y,mat(3,3)                                !Matrix and numbers for calc
    integer :: i,j,k,l                                      !Dummy ints

    chain = 0._dp
    ! Locates the coordinates for the relevant backbone atoms
    do i=1,size(pep)
       if (pep(i)%heteroatom) then
          allocate(pep(i)%chi(0))
          cycle
       end if
       k=0
       do j=1,size(chi_rels(i)%chi_on)
          if(chi_rels(i)%chi_on(j)) k=k+1
       end do
       allocate(pep(i)%chi(k))
       if (k==0) cycle
       ! Calculates the chi values for each respective protein
       pep(i)%chi = 0._dp
       do j=1,size(pep(i)%chi)
          chainLoop: do k=1,size(chain,2)
             do l=1,size(pep(i)%backbone)
                if(trim(pep(i)%backbone(l)%label)==trim(chi_rels(i)%chain(k,j))) then
                   chain(:,k) = pep(i)%backbone(l)%coords
                   cycle chainLoop
                end if
             end do
             do l=1,size(pep(i)%sidechain)
                if(trim(pep(i)%sidechain(l)%label)==trim(chi_rels(i)%chain(k,j))) then
                   chain(:,k) = pep(i)%sidechain(l)%coords
                   cycle chainLoop
                end if
             end do
          end do chainLoop

          ! Calculates plane normals, relation to axis of rotation for sign of angle
          ! Completes calculation for only relevant chi, leaves others set to 0.0
          v1=chain(:,2)-chain(:,1); v2=chain(:,3)-chain(:,2); v3=chain(:,4)-chain(:,3)
          v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
          mat(:,1) = v2
          mat(:,2) = cross(v1,v2)
          mat(:,3) = cross(v2,v3)
          x = dot_product(mat(:,2),mat(:,3))
          y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
          pep(i)%chi(j) = atan2(y,x)*180._dp/pi
       end do
    end do

    ! SANITY CHECK TEST
    ! do i=1,size(pep)
    !    write(*,*) pep(i)%label,pep(i)%chi
    ! end do
    ! stop
    ! SANITY CHECK TEST
  end subroutine calc_chi
  ! -----------------------------------------------------------------------------

  subroutine set_chi(pep,res_idx,chi_idx,chi_val,chi_rels)
    ! SUBROUTINE: set_chi(peptide_chain,residue_index, chi_index,chi_value,chi_relations)
    ! Depends on all of the atoms dependent on an axis of rotation following the
    ! axis members in the sidechain atom list. This is the case for the add_sc from the lib
    type(residue),intent(inout) :: pep(:)                   !Polypeptide chain
    integer,intent(in) :: res_idx                           !Residue index
    integer,intent(in) :: chi_idx                           !X1,X2,X3....
    real(dp),intent(in) :: chi_val                          !new value for angle in degrees
    type(chi_rel),intent(in) :: chi_rels(:)                 !Chi relations
    real(dp) :: d_angle                                     !Change in angle
    real(dp),dimension(3) :: origin,axis,coords,vector      !Defined for changes
    real(dp) :: rot(3,3)                                    !Rotation matrix
    integer :: i,k

    d_angle = (chi_val - pep(res_idx)%chi(chi_idx))*pi/180._dp

    ! Locating axis of rotation
    do i=1,size(pep(res_idx)%backbone)
       if(trim(pep(res_idx)%backbone(i)%label)==trim(chi_rels(res_idx)%axis(1,chi_idx))) then
          origin = pep(res_idx)%backbone(i)%coords
       end if
    end do
    do i=1,size(pep(res_idx)%sidechain)
       if(trim(pep(res_idx)%sidechain(i)%label)==trim(chi_rels(res_idx)%axis(1,chi_idx))) then
          origin = pep(res_idx)%sidechain(i)%coords
       else if(trim(pep(res_idx)%sidechain(i)%label)==trim(chi_rels(res_idx)%axis(2,chi_idx))) then
          axis = pep(res_idx)%sidechain(i)%coords-origin
          k=i
       end if
    end do

    call RotMatrix(d_angle,axis,rot)

    do i=k,size(pep(res_idx)%sidechain)
       vector = pep(res_idx)%sidechain(i)%coords - origin
       vector = matmul(rot,vector)
       coords = origin+vector
       pep(res_idx)%sidechain(i)%coords = coords
    end do
    pep(res_idx)%chi(chi_idx) = chi_val

  end subroutine set_chi
  ! -----------------------------------------------------------------------------
  subroutine add_sc(pep)
    ! SUBROUTINE: add_sc(peptide_chain)
    ! Collects the fragments of sidechains from a pdb file then places them with respec to the alpha carbon
    type(residue), intent(inout) :: pep(:)                  !Polypeptide chain
    type(fragment),allocatable :: fragments(:)              !Fragments of sidechains to add
    real(dp) :: r1(3),r2(3),r3(3)                           !Vector to next point
    real(dp) :: rotv(3)                                     !Rotation vector
    real(dp) :: Rmat(3,3)                                   !Rotation matrix
    real(dp) :: angle                                       !Change of angle
    real(dp) :: coords(3)                                   !Coordinates temp
    character(len=4) :: atom_label                          !Atom label
    integer :: i,j,k                                        !Dummy ints

    call import_fragments(fragments)

    ! First finds the unit vector which points along the outward vertex of a
    ! Tetrahedron for a right handed structure (r2)
    ! Also dictates L-enantiomer
    do i=1,size(pep)
       atom_label = 'CA'
       do j=1,size(pep(i)%backbone)
          if(trim(pep(i)%backbone(j)%label)==trim(atom_label)) then
             coords = pep(i)%backbone(j)%coords
             exit
          end if
       end do
       atom_label = 'N'
       do j=1,size(pep(i)%backbone)
          if(trim(pep(i)%backbone(j)%label)==trim(atom_label)) then
             r1 = pep(i)%backbone(j)%coords - coords
             exit
          end if
       end do
       atom_label = 'C'
       do j=1,size(pep(i)%backbone)
          if(trim(pep(i)%backbone(j)%label)==trim(atom_label)) then
             r3 = pep(i)%backbone(j)%coords - coords
             exit
          end if
       end do
       r1 = r1/norm2(r1)
       r3 = r3/norm2(r3)
       r2 = r1
       angle = -0.5*acos(dot_product(r3,r1))
       rotv =  cross(r3,r1)
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)

       rotv = cross(r2,rotv)
       angle = -1*(180.-109.5/2.)*pi/180
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)

       ! Translate to align alpha carbons
       do j=1,size(fragments)
          if(trim(fragments(j)%res_label) == trim(pep(i)%label)) then
             r1 = coords - fragments(j)%CA%coords
             if(allocated(pep(i)%sidechain)) deallocate(pep(i)%sidechain)
             allocate(pep(i)%sidechain(size(fragments(j)%atoms)))
             do k=1,size(fragments(j)%atoms)
                pep(i)%sidechain = fragments(j)%atoms
             end do
             do k=1,size(pep(i)%sidechain)
                pep(i)%sidechain(k)%coords = pep(i)%sidechain(k)%coords+r1
             end do
          else if(trim(fragments(j)%res_label) == 'HIL' .and. trim(pep(i)%label)=='HIS') then
             r1 = coords - fragments(j)%CA%coords
             if(allocated(pep(i)%sidechain)) deallocate(pep(i)%sidechain)
             allocate(pep(i)%sidechain(size(fragments(j)%atoms)))
             do k=1,size(fragments(j)%atoms)
                pep(i)%sidechain = fragments(j)%atoms
             end do
             do k=1,size(pep(i)%sidechain)
                pep(i)%sidechain(k)%coords = pep(i)%sidechain(k)%coords+r1
             end do
          end if
       end do

       ! Rotates the sidechain to align with the designated CA-CB vector
       do j=1,size(pep(i)%sidechain)
          if(trim(pep(i)%sidechain(j)%label) == 'CB') then
             r1 = pep(i)%sidechain(j)%coords - coords
             exit
          end if
       end do
       r1 = r1/norm2(r1)
       rotv  = cross(r1,r2)
       angle = acos(dot_product(r1,r2))
       call RotMatrix(angle,rotv,Rmat)

       do j=1,size(pep(i)%sidechain)
          r1 = pep(i)%sidechain(j)%coords - coords
          r1 = matmul(Rmat,r1)
          pep(i)%sidechain(j)%coords = coords + r1
       end do
    end do

  end subroutine add_sc
  ! -----------------------------------------------------------------------------

  subroutine import_fragments(fragments)
    ! Uses the libraries to build fragments for all types of amino acids
    type(fragment),allocatable,intent(inout) :: fragments(:)!Fragments of sidechains to add
    integer :: unit                                         !File units
    character(len=100) :: line                              !Read line
    character(len=100) :: targ                              !Target in search
    character(len=4) :: atom_label                          !Atom label
    character(len=3) :: res_label                           !Residue label
    integer :: i,j,k,l                                      !Dummy 
    integer :: info                                         !Error
    integer ::nFrag                                         !Number of fragments
    character(len=80) :: errmsg                             !Error msg


    info=0
    open(newunit=unit, file=trim(geom_lib), status='old',IOSTAT = info,IOMSG=errmsg)
    if(info /=0)then
       write(*,*)errmsg
       stop 'Failure in scd.f90'
    end if

    nFrag=1
    read(unit,'(A17,A4)')line,targ
    do
       read(unit,'(A17,A4)',iostat=info,iomsg=errmsg)line,atom_label
       if(info /=0)then
          exit
       end if
       if(trim(atom_label)/=trim(targ)) then
          targ=atom_label
          nFrag=nFrag+1
       end if
    end do
    allocate(fragments(nFrag))
    rewind(unit)
    read(unit,'(A17,A3)',iostat=info,iomsg=errmsg)line,fragments(1)%res_label
    do i=1,size(fragments)-1
       do
          read(unit,'(A17,A3)',iostat=info,iomsg=errmsg)line,res_label
          if(info /=0)then
             exit
          end if
          if(trim(res_label)/=trim(fragments(i)%res_label)) then
             fragments(i+1)%res_label=res_label
             exit
          end if
       end do
    end do
    rewind(unit)
    
    do i=1,size(fragments)
       ! Building fragments from appropriate bond lengths and angles
       ! First read to count atoms, then read to assign
       j=0;k=0;
       do
          read(unit,'(A12,A4,A1,A3)',iostat=info)line,atom_label,line,res_label
          k=k+1
          if(trim(res_label)/=trim(fragments(i)%res_label) .or. info/=0) exit
          if (.not. bb_list(atom_label)) then
             j=j+1
          end if
       end do
       do l=1,k
          backspace(unit)
       end do
       allocate(fragments(i)%atoms(j))
       j=1
       do
          read(unit,'(A12,A4,A1,A3,A11)',advance='no',iostat=info)line,atom_label,line,res_label,line
          if(trim(res_label)/=trim(fragments(i)%res_label).or. info/=0) exit
          if (.not. bb_list(atom_label)) then
             fragments(i)%atoms(j)%label=atom_label
             fragments(i)%atoms(j)%sym = atom_label(1:1)
             fragments(i)%atoms(j)%idx=0
             read(unit,*)fragments(i)%atoms(j)%coords(:)
             j=j+1
          else if(trim(atom_label) == 'CA') then
             fragments(i)%CA%label=atom_label
             fragments(i)%CA%sym = atom_label(1:1)
             fragments(i)%CA%idx=0
             read(unit,*)fragments(i)%CA%coords(:)
          else
             read(unit,*)
          end if
       end do
       backspace(unit)
    end do
  end subroutine import_fragments
  ! -----------------------------------------------------------------------------

  subroutine chi_relations(pep,chi_rels)
    type(chi_rel), intent(inout) :: chi_rels(:)
    type(residue), intent(in) :: pep(:)
    character(len=10) :: targ
    character(len=100) :: line
    character(len=80) :: errmsg                             !Error msg
    integer :: info
    integer :: unit
    integer ::i,j

    open(newunit=unit, file=trim(chi_lib), status='old',IOSTAT = info,IOMSG=errmsg)
    if(info /=0)then
       write(*,*)errmsg
       stop 'Failure in scd.f90'
    end if

    do i=1,size(pep)
       targ = 'CHI1'
       do
          read(unit,*) line
          if (trim(line) == trim(targ)) exit
       end do
       targ  = pep(i)%label
       chi_rels(i)%chi_on =.false.
       outer: do j=1,5
          do
             read(unit,'(A8)',advance='no',iostat=info) line
             if(line(1:len_trim(targ)) == trim(targ)) then
                chi_rels(i)%chi_on(j) = .true.
                read(unit,*) chi_rels(i)%axis(:,j),chi_rels(i)%chain(:,j)
                exit
             else if (line(1:3) == 'END')then
                exit outer
             else
                read(unit,*)
             end if
          end do
       end do outer
       rewind(unit)
    end do
  end subroutine chi_relations
  ! -----------------------------------------------------------------------------
end module sidechain_decoration
