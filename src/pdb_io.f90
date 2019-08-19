! LAST EDIT: Phil Maffettone 2017_01_06
! Updated to include reading sequence
module pdb_io

  use types
  use maths
  use class_proteins

  implicit none

  private :: sequence2geometry
contains

  ! --------------------Public Contents-----------------------------------------
  ! read_pdb(filename,protein,[backbone_only,info_integer,error_message])
  ! write_pdb(filename,protein,[backbone_only])
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! sequence2geometry()
  ! ----------------------------------------------------------------------------

  subroutine read_pdb(filename,prot,bbonly,info,errmsg)
    ! Subroutine to read in pdb file, to peptide pep.
    ! Uses atom class different from that used in phonon calculations (no mass or number)
    ! Residue naming appropriate to backbone in PDB (1CFC.pdb was used as reference)
    ! Backbone candidates are in class_proteins function
    ! Defaults to a backbone only Protein object

    character(len=*),intent(in) :: filename                 !Filename
    type(protein),intent(inout) :: prot                     !Protein
    type(residue), allocatable :: pep(:)                    !Polypeptide
    logical, intent(in), optional :: bbonly                 !Backbone only
    character(len=4) :: atom_label                          !Atom_label
    character(len=3),allocatable :: res_label(:)            !All AA residue labels
    character(len=2) :: sym                                 !Atomic Symbol
    character(len=80) :: line,scanfor                       !Buffer and search
    integer :: i,j,k                                        !Dummy integers
    integer :: unit                                         !File I/O unit
    integer :: n_res                                        !Number of residues
    integer :: n_atoms                                      !Number of atoms
    integer :: seq_n                                        !Sequence number
    integer :: n_lines                                      !Number of lines
    type(atom),allocatable :: atom_list(:)                  !Temporary list of atoms
    real(dp) :: coords(3)                                   !Atomic Coordinates
    integer, optional,intent(out) :: info                   !Flag, 0 for Success
    integer :: info_loc                                     !Local info flag
    character(len=80), optional, intent(out) :: errmsg      !Error Message
    character(len=80) :: errmsg_loc                         !Local error message
    
    info_loc = 0
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Reading primary protein structure
    do
       read(unit,'(A)') scanfor
       if (scanfor(1:6) == 'SEQRES') exit
       if (scanfor(1:10) == 'REMARK 465') then
          i=465
          errmsg = '(Missing residues) - please check SEQRES and ATOM records accordingly'
       end if
    end do
    backspace(unit)
    read(unit,"(A13,i4)")line,n_res
    allocate(pep(n_res),res_label(n_res))
    backspace(unit)
    backspace(unit)
    k=1
    do i=1,(n_res/13+1)
       read(unit,*)
       read(unit,'(A19)',advance='no')line
       do j=1,13
          read(unit,'(A3,A1)',advance='no')res_label(k),line
          if(k==n_res) then
             exit
          else
             k=k+1
          end if
       end do
    end do

    do
       read(unit,'(A6)', IOSTAT = info_loc)scanfor
       if(info_loc<0) then      !EOF only sequence data avail
          errmsg_loc = 'End of input pdb file reached prior to structural information'
       end if
       if(scanfor(1:5) == 'MODEL')exit
       if(scanfor == 'ATOM  ' .or. scanfor == 'HETATM') then
          backspace(unit); exit;
       end if
    end do

    do i=1,n_res
       ! Count the atoms in the residue
       n_atoms=0; n_lines=0
       ! First pull sequence number
       do
          read(unit,"(A6,A16,i4)")scanfor,line,seq_n
          if (scanfor == 'ATOM  ' .or. scanfor == 'HETATM') then
             exit
          end if
       end do
       backspace(unit)
       do
          read(unit,"(A6,A16,i4)")scanfor,line,k
          n_lines=n_lines+1
          if(scanfor == 'TER   ')exit
          ! Ignores other information such as UISO
          if (scanfor == 'ATOM  ' .or. scanfor == 'HETATM') then
             if (k/=seq_n) then
                exit
             else
                n_atoms=n_atoms+1
             end if
          end if
       end do
       allocate(atom_list(n_atoms))
       do j=1,n_lines
          backspace(unit)
       end do
       ! Build a list of atoms then pass to the residue constructor
       do j=1,n_atoms
          do
             read(unit,"(A6)",advance = 'no') scanfor
             if(scanfor == 'ATOM  ' .or. scanfor == 'HETATM') then
                read(unit,"(A6,A4,A1,A3,A10,3f8.3,A22,A2)")&
                     line,atom_label,line,res_label(i),line,coords,line,sym
                atom_list(j) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
                exit
             else
                read(unit,*)line
             end if
          end do
       end do
       if(present(bbonly)) then
          pep(i) = residue_(res_label(i),atom_list,bbonly)
       else
          pep(i) = residue_(res_label(i),atom_list)
       end if
       deallocate(atom_list)
    end do
    close(unit)

    call calc_dihedral(pep)
    if(present(bbonly)) then
       prot = protein_(pep,bbonly)
    else
       prot = protein_(pep, bbonly = .true.)
    end if

  end subroutine read_pdb
  ! --------------------------------------------------------------------------------------------

  subroutine write_pdb(filename,prot,bb_only)
    ! Writes simplified .PDB file from a protein structure and a file name

    character(len=*),intent(in) :: filename                  !Filename
    type(protein),intent(in):: prot                          !Protein object
    character(len=3) :: res_label                            !Residue label
    character(len=4) :: atom_label                           !Atom label
    character(len=2) :: sym                                  !Atomic symbol
    real(dp) :: coords(3)                                    !Atomic coordinates
    integer :: n_res,n_atoms                                 !Total number of each
    integer :: i,j,k                                         !Dummy integers
    integer :: unit                                          !File I/O unit
    logical, optional :: bb_only                             !Backbone only

    ! TODO: YET TO INCLUDE Proper HEADER AND REMARKS
    n_res = prot%n_res

    open(newunit=unit,file=filename,status='replace')

    ! Writing Header section
    write(unit,"(A6,T11,A40,A9,T63,A4)")'HEADER', 'NULL','DATE', 'ID'
    write(unit,"(A6,T9,A2,A70)")'TITLE','  ',adjustl(trim(filename))
    k=1
    do i=1,(n_res/13+1)
       write(unit,"(A6,T8,I3,T12,A1,T14,I4,A2)",advance='no') 'SEQRES',i,'A',n_res,'  '
       do j=1,13
          res_label = get_label(prot%pep(k))
          write(unit,"(A3,A1)",advance='no')res_label, ' '
          if(k==n_res) then
             exit
          else
             k=k+1
          end if
       end do
       write(unit,*)
    end do

15  format(A6,I5,T13,A4,A1,A3,T22,A1,I4,A1,T31,3F8.3,2F6.2,T77,A2,A2)
    k=0
    do i=1,n_res
       n_atoms = size(prot%pep(i)%backbone)
       res_label = get_label(prot%pep(i))
       do j=1,n_atoms
          k=k+1
          atom_label = get_label(prot%pep(i)%backbone(j))
          coords = get_coords(prot%pep(i)%backbone(j))
          sym = adjustr(get_sym(prot%pep(i)%backbone(j)))
          write(unit,15)'ATOM  ',k,atom_label,' ',res_label,'A',i,' ',coords,1._dp,0._dp,sym,'  '
       end do
       if(present(bb_only)) then
          if(bb_only) cycle
       end if
       n_atoms = size(prot%pep(i)%sidechain)
       if(n_atoms<=0)cycle
       do j=1,n_atoms
          k=k+1
          atom_label = get_label(prot%pep(i)%sidechain(j))
          coords = get_coords(prot%pep(i)%sidechain(j))
          sym = adjustr(get_sym(prot%pep(i)%sidechain(j)))
          write(unit,15)'ATOM  ',k,atom_label,' ',res_label,'A',i,' ',coords,1._dp,0._dp,sym,'  '
       end do
    end do
    write(unit,"(A6,I5,T13,A4,A1,A3,T22,A1,I4,T80)") 'TER   ',k,'    ',' ',res_label,'A',i
    write(unit,"(A)") 'END'
    close(unit)

  end subroutine write_pdb
  ! -----------------------------------------------------------------------------

  subroutine sequence2geometry(pep,res_label,n_res)
    ! Local subroutine generates geometry from sequence information only
    type(residue), allocatable,intent(inout) :: pep(:)      !Polypeptide
    character(len=4) :: atom_label                          !Atom_label
    character(len=3),allocatable,intent(inout) :: res_label(:) !All AA residue labels
    integer, intent(inout) :: n_res                         !Numbrer of residues
    character(len=2) :: sym                                 !Atomic Symbol
    type(atom),allocatable :: atom_list(:)                  !Temporary list of atoms
    real(dp) :: r1(3),r2(3)                                 !Vector to next point
    real(dp) :: rotv(3)                                     !Rotation vector
    real(dp) :: Rmat(3,3)                                   !Rotation matrix
    real(dp) :: angle                                       !Change of angle
    real(dp) :: coords(3)                                   !Atomic Coordinates
    integer :: i,j                                          !Dummy Integers
    integer :: n_atoms                                      !Number of atoms
    logical, allocatable :: amino(:)                        !List if het atoms or not proper AA
    logical :: het                                          !Marking heteroatoms moving forward

    ! Removing potential heteroatoms
    allocate(amino(size(res_label)))
    do i=1,size(res_label)
       amino(i) = amino_list(res_label(i))
    end do
    het= .false.
    j=0
    do i=1,size(amino)
       if (amino(i)) then
          if(het) then
             write(*,*) 'Poor sequence input as PDB. Improper amino in middle of chain.'
             write(*,*) i, res_label(i)
             stop
          end if
          j=j+1
       else
          het = .true.
       end if
    end do
    deallocate(pep)
    n_res = j
    allocate(pep(n_res))

    ! The chain will be built straight and two dimensional
    ! then use a rotation matrix about the z axis to bend
    ! Values from Engh_1991
    n_atoms = 4
    allocate(atom_list(n_atoms))
    rotv = [0., 0., 1.]
    coords = 0.
    r1 = [1., 0., 0.]
    sym= 'N '
    atom_label = sym
    atom_list(1) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
    r2 = 1.458*r1
    angle = pi - (121.7*pi/180.)
    call RotMatrix(angle,rotv,Rmat)
    r2 = matmul(Rmat,r2)
    coords = coords+r2
    sym = 'C '
    atom_label = 'CA'
    atom_list(2) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
    r1 = r2/norm2(r2)
    r2 = 1.525*r1
    angle = pi - (111.0*pi/180.)
    call RotMatrix(angle,rotv,Rmat)
    r2 = matmul(Rmat,r2)
    coords = coords+r2
    sym = 'C '
    atom_label = sym
    atom_list(3) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
    r1 = r2/norm2(r2)
    r2 = 1.231*r1
    angle = pi - (120.1*pi/180.)
    call RotMatrix(angle,rotv,Rmat)
    r2 = matmul(Rmat,r2)
    coords = coords+r2
    sym = 'O '
    atom_label = sym
    atom_list(4) =  atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
    coords = coords - r2                                    !Bring coords back to C
    pep(1) = residue_(res_label(1),atom_list)

    do i=2,n_res
       r2 = 1.329*r1
       angle = (117.2*pi/180.) - pi                         !Opposite dir of C=O
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)
       coords = coords+r2
       sym = 'N '
       atom_label = sym
       atom_list(1) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
       r1 = r2/norm2(r2)
       r2 = 1.458*r1
       angle = pi - (121.7*pi/180.)
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)
       coords = coords+r2
       sym = 'C '
       atom_label = 'CA'
       atom_list(2) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
       r1 = r2/norm2(r2)
       r2 = 1.525*r1
       angle = pi - (111.0*pi/180.)
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)
       coords = coords+r2
       sym = 'C '
       atom_label = sym
       atom_list(3) = atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
       r1 = r2/norm2(r2)
       r2 = 1.231*r1
       angle = pi - (120.1*pi/180.)
       call RotMatrix(angle,rotv,Rmat)
       r2 = matmul(Rmat,r2)
       coords = coords+r2
       sym = 'O '
       atom_label = sym
       atom_list(4) =  atom_(label= adjustl(atom_label),sym=adjustl(sym),coords=coords)
       coords = coords - r2                                    !Bring coords back to C
       pep(i) = residue_(res_label(i),atom_list)
    end do

  end subroutine sequence2geometry
end module pdb_io
