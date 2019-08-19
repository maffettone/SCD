! LAST EDIT: Phil Maffettone 2017-01-06
module class_proteins

  ! This class is for the construction of relevant proteins. Proteins will be constructed
  ! as arrays of amino acids called "residues." Each residue will contain the backbone atoms
  ! and dihedral angles, as well as a sidechain class of atoms and dihedrals.

  use types
  use maths

  implicit none

  private
  public :: residue,residue_,atom,atom_,change_coords,get_label, get_sym,&
       get_coords,change_dihedral,get_dihedral, protein, protein_,calc_dihedral,&
       amino_list, bb_list

  interface get_label
     module procedure atom_label
     module procedure residue_label
  end interface get_label


  type atom
     character(len=4) :: label                              !User atom label
     character(len=2) :: sym                                !Atomic symbol
     real(dp) :: coords(3)                                  !Atomic coordinates
     integer :: idx                                         !Index for running matrix of coordinates
  end type atom

  type residue
     character(len=3) :: label                              !Standard amino acid label
     real(dp) :: phi,psi,w                                  !Angles in degrees
     real(dp),allocatable :: chi(:)                         !Angles in degrees
     type(atom),allocatable :: sidechain(:),backbone(:)     !Split arrays of atoms
     logical :: heteroatom                                  !For non residue components of protein
  end type residue

  type protein
     type(residue), allocatable :: pep(:)                   !List of residues (peptide)
     real(dp), allocatable :: coordinates(:,:)              !Running matrix of coordinates (3,n_atoms)
     logical, allocatable :: ca_mask(:)                     !Demarcates alpha carbons in coordinate list
     logical, allocatable :: c_mask(:)                      !Demarcates any carbons in coordinate list
     logical, allocatable :: o_mask(:)                      !Demarcates any oxygens in coordinate list
     logical, allocatable :: n_mask(:)                      !Demarcates any nitrogens in coordinate list
     logical, allocatable :: h_mask(:)                      !Demarcates any hydrogens in coordinate list
     logical, allocatable :: s_mask(:)                      !Demarcates any sulfurs in coordinate list
     character(len=2) :: species(5)                         !List of relevant symbols
     integer :: spec_count(5)                               !Count of atoms of each species
     integer :: n_res                                       !Number of residues
     logical :: bbonly                                      !Backbone only
  end type protein

contains

  ! --------------------Public Contents-----------------------------------------
  ! FUNCTION: protein_(residue_list,[bbonly])
  ! FUNCTION: residue_(label,atom_list)
  ! FUNCTION: atom_(label,symbol,coordinates,[idx])
  ! SUBROUTINE: change_coords(atom,new_coordinates)
  ! FUNCTION: get_label(atom or residue)
  ! FUNCTION: get_sym(atom)
  ! FUNCTION: get_coords(atom)
  ! FUNCTION: get_dihedral(residue,angle_type)
  ! SUBROUTINE: change_dihedral(residue, change_in_angle, angle_type)
  ! SUBROUTINE: calc_dihedral(residue_list)
  ! FUNCTION: bb_list(atom_label)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! FUNCTION: count_spec(residue_list,atomic_symbol,bbonly)
  ! ----------------------------------------------------------------------------

  function atom_(label,sym,coords,idx) result(r)
    ! Atom Constructor
    character(len=4), intent(in) :: label                   !Atom label
    character(len=2), optional, intent(in) :: sym           !Atomic Symbol
    real(dp) :: coords(3)                                   !Atomic coordinates
    integer, optional :: idx                                !Reference index
    integer :: i                                            !Dummy
    type(atom) :: r                                         !Return atom
    if(present(sym)) then
       r%sym = sym
    else
       r%sym = label(1:2)
       write(*,"(A)") 'No element was given, assigned by first letter of label'
    end if
    r%label = trim(adjustl(label))

    ! Symbol check and removing capital letter
    if(.not. alphabetical(r%sym(1:1)) .or. .not. alphabetical(r%sym(2:2))) then
       i=1
       do
          if( alphabetical(label(i:i))) then
             exit
          else
             i = i+1
          end if
       end do
       r%sym = label(i:i+1)
    end if
    if(.not.(ichar(r%sym(2:2))>=iachar("a") .and. ichar(r%sym(2:2))<=iachar("z"))) then
       r%sym(2:2) = ' '
    end if
    
    r%coords = coords
    if (present(idx)) then
       r%idx = idx
    else
       r%idx = 0
    end if
  end function atom_
  ! ----------------------------------------------------------------------------

  function atom_label(x) result(r)
    type(atom),intent(in) :: x
    character(len=4) :: r
    r = x%label
  end function atom_label
  ! ----------------------------------------------------------------------------

  function get_sym(x) result(r)
    type(atom),intent(in) :: x
    character(len=2) :: r
    r = x%sym
  end function get_sym
  ! ----------------------------------------------------------------------------

  function get_coords(x) result(r)
    type(atom),intent(in) :: x
    real(dp) :: r(3)
    r = x%coords
  end function get_coords
  ! ----------------------------------------------------------------------------

  subroutine change_coords(a,x,prot)
    type(atom), intent(inout) :: a                          !Atom
    real(dp), intent(in) :: x(3)                            !New Coordinates
    type(protein), intent(inout) :: prot                    !Protein object for complete coords
    a%coords = x
    prot%coordinates(:,a%idx) = x
  end subroutine change_coords
  ! ----------------------------------------------------------------------------

  function residue_(label,atoms,bbonly) result(r)
    ! FUNCTION:  residue_(label,atom_list,[backbone_only])
    ! Initializes residue by allocating apropriate size, and populating with list of atoms
    ! If no atoms are given, all set to XXXX for addition later on
    ! Terminal residues have a OH group or NH3 (instead of NH) in the backbone
    ! Backbone only is assumed
    type(atom), intent(in) :: atoms(:)                      !List of atoms
    character(len=3), intent(inout) :: label                !Residue Label
    logical, intent(in), optional :: bbonly                 !Backbone Only
    character(len=4) :: atom_label                          !Atom Label
    type(residue) :: r                                      !Returned residue
    integer :: n                                            !Number of atoms
    integer :: i,j,k                                        !Dummy integers
    integer :: n_bb                                         !Number of atoms in backbone
    r%phi=0._dp;r%psi=0._dp;r%w=0._dp;
    r%heteroatom = .false.
    call upper_case(label)
    r%label = label
    n = size(atoms)

    ! Counting number in the backbone
    n_bb=0
    do i=1,n
       atom_label = get_label(atoms(i))
       if(bb_list(atom_label)) n_bb=n_bb+1
    end do

    if(present(bbonly) .and. .not. bbonly) then
       allocate(r%sidechain(n-n_bb))
    else
       if(allocated(r%sidechain)) deallocate(r%sidechain)
    end if
    allocate(r%backbone(n_bb))
    if(.not. amino_list(label)) r%heteroatom = .true.

    ! Adding atoms
    j=1
    k=1
    do i=1,n
       atom_label = get_label(atoms(i))
       if (bb_list(atom_label)) then
          r%backbone(j) = atoms(i)
          j=j+1
       elseif(present(bbonly) .and. .not. bbonly) then
          r%sidechain(k) = atoms(i)
          k=k+1
       end if
    end do
  end function residue_
  ! -----------------------------------------------------------------------------


  function residue_label(res) result(r)
    type(residue),intent(in) :: res
    character(len=3) :: r
    r = res%label
  end function residue_label
  ! ------------------------------------------------------------------------------

  function protein_(pep,bbonly) result(r)
    type(protein) :: r                                      !Constructed protein
    type(residue), intent(in) :: pep(:)                     !List of residues
    logical, optional :: bbonly                             !Backbone only
    integer :: n_atoms                                      !Number of atoms
    integer :: i,j,k                                        !Dummy integer

    allocate(r%pep(size(pep)))
    !call calc_dihedral(r%pep)
    r%pep = pep
    r%n_res = size(pep)
    r%species(1) = 'C '
    r%species(2) = 'O '
    r%species(3) = 'N '
    r%species(4) = 'H '
    r%species(5) = 'S '
    ! Backbone only refers only to the calculations the protein is used for.
    ! Sidechain atoms may appear in the residue, however they will be disregarded for
    ! calculation purposes.
    if (present(bbonly)) then
       r%bbonly = bbonly
    else
       r%bbonly = .true.
    end if

    n_atoms = 0
    do i=1,5
       r%spec_count(i) = count_spec(pep,r%species(i),bbonly)
       n_atoms = n_atoms + r%spec_count(i)
    end do

    allocate(r%coordinates(3,n_atoms))
    allocate(r%ca_mask(n_atoms),r%c_mask(n_atoms),r%o_mask(n_atoms),r%n_mask(n_atoms),&
         r%h_mask(n_atoms), r%s_mask(n_atoms))
    r%coordinates = 0._dp
    r%ca_mask = .false.
    r%c_mask = .false.
    r%o_mask = .false.
    r%n_mask = .false.
    r%h_mask = .false.
    r%s_mask = .false.

    k=0
    if(r%bbonly) then
       do i=1,size(pep)
          do j=1,size(pep(i)%backbone)
             k=k+1
             r%pep(i)%backbone(j)%idx = k
             r%coordinates(:,k) = pep(i)%backbone(j)%coords
             if(trim(pep(i)%backbone(j)%label) == 'CA')r%ca_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='C ') r%c_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='O ') r%o_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='N ') r%n_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='H ') r%h_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='S ') r%s_mask(k) = .true.
          end do
       end do
       do i=1,size(pep)
          allocate(r%pep(i)%sidechain(0))
       end do
    else
        do i=1,size(pep)
          do j=1,size(pep(i)%backbone)
             k=k+1
             r%pep(i)%backbone(j)%idx = k
             r%coordinates(:,k) = pep(i)%backbone(j)%coords
             if(trim(pep(i)%backbone(j)%label) == 'CA')r%ca_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='C ') r%c_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='O ') r%o_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='N ') r%n_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='H ') r%h_mask(k) = .true.
             if(pep(i)%backbone(j)%sym =='S ') r%s_mask(k) = .true.
          end do
          do j=1,size(pep(i)%sidechain)
             k=k+1
             r%pep(i)%sidechain(j)%idx = k
             r%coordinates(:,k) = pep(i)%sidechain(j)%coords
             if(pep(i)%sidechain(j)%sym =='C ') r%c_mask(k) = .true.
             if(pep(i)%sidechain(j)%sym =='O ') r%o_mask(k) = .true.
             if(pep(i)%sidechain(j)%sym =='N ') r%n_mask(k) = .true.
             if(pep(i)%sidechain(j)%sym =='H ') r%h_mask(k) = .true.
             if(pep(i)%sidechain(j)%sym =='S ') r%s_mask(k) = .true.
          end do
       end do
    end if
  end function protein_
  !-----------------------------------------------------------------------------

  subroutine calc_dihedral(pep)
    ! Depends on backbone being named with N-CA-C-
    ! For initial calculation of, w, phi, psi
    ! The chain runs CA(i-1)--C(i-1)--N(i)--CA(i)--C(i)--N(i+1)
    type(residue), intent(inout) :: pep(:)              !Input protein chain
    real(dp) :: chain(3,6)                              !Chain of xyz coordinates
    real(dp) ::v1(3),v2(3),v3(3)                        !Set of unit normals
    real(dp) :: x,y,mat(3,3)                            !Matrix and numbers for calc
    integer :: i,j


    ! Locates the coordinates for the relevant backbone atoms
    do i=2,size(pep)-1
       if(pep(i+1)%heteroatom)exit                          !Only works for heteroatoms after chain
       do j=1, size(pep(i-1)%backbone)
          if(trim(pep(i-1)%backbone(j)%label) == 'CA') then
             chain(:,1) = pep(i-1)%backbone(j)%coords
          end if
          if(trim(pep(i-1)%backbone(j)%label) == 'C') then
             chain(:,2) = pep(i-1)%backbone(j)%coords
          end if
       end do
       do j=1, size(pep(i)%backbone)
          if(trim(pep(i)%backbone(j)%label) == 'N') then
             chain(:,3) = pep(i)%backbone(j)%coords
          end if
          if(trim(pep(i)%backbone(j)%label) == 'CA') then
             chain(:,4) = pep(i)%backbone(j)%coords
          end if
          if(trim(pep(i)%backbone(j)%label) == 'C') then
             chain(:,5) = pep(i)%backbone(j)%coords
          end if
       end do
       do j=1,size(pep(i+1)%backbone)
          if(trim(pep(i+1)%backbone(j)%label) == 'N') then
             chain(:,6) = pep(i+1)%backbone(j)%coords
          end if
       end do

       ! Calculates plane normals, relation to axis of rotation for sign of angle
       v1=chain(:,2)-chain(:,1); v2=chain(:,3)-chain(:,2); v3=chain(:,4)-chain(:,3)
       v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
       mat(:,1) = v2
       mat(:,2) = cross(v1,v2)
       mat(:,3) = cross(v2,v3)
       x = dot_product(mat(:,2),mat(:,3))
       y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
       pep(i)%w = atan2(y,x)*180._dp/pi
       v1=chain(:,3)-chain(:,2); v2=chain(:,4)-chain(:,3); v3=chain(:,5)-chain(:,4)
       v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
       mat(:,1) = v2
       mat(:,2) = cross(v1,v2)
       mat(:,3) = cross(v2,v3)
       x = dot_product(mat(:,2),mat(:,3))
       y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
       pep(i)%phi = atan2(y,x)*180._dp/pi
       v1=chain(:,4)-chain(:,3); v2=chain(:,5)-chain(:,4); v3=chain(:,6)-chain(:,5)
       v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
       mat(:,1) = v2
       mat(:,2) = cross(v1,v2)
       mat(:,3) = cross(v2,v3)
       x = dot_product(mat(:,2),mat(:,3))
       y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
       pep(i)%psi = atan2(y,x)*180._dp/pi
    end do

    ! Repeats the above for terminal amino acids (first and last before heteroatom)
    ! No phi or omega for N terminus
    ! No psi for C terminus
    do j=1, size(pep(i-1)%backbone)
       if(trim(pep(i-1)%backbone(j)%label) == 'CA') then
          chain(:,1) = pep(i-1)%backbone(j)%coords
       end if
       if(trim(pep(i-1)%backbone(j)%label) == 'C') then
          chain(:,2) = pep(i-1)%backbone(j)%coords
       end if
    end do
    do j=1, size(pep(i)%backbone)
       if(trim(pep(i)%backbone(j)%label) == 'N') then
          chain(:,3) = pep(i)%backbone(j)%coords
       end if
       if(trim(pep(i)%backbone(j)%label) == 'CA') then
          chain(:,4) = pep(i)%backbone(j)%coords
       end if
       if(trim(pep(i)%backbone(j)%label) == 'C') then
          chain(:,5) = pep(i)%backbone(j)%coords
       end if
    end do
    v1=chain(:,2)-chain(:,1); v2=chain(:,3)-chain(:,2); v3=chain(:,4)-chain(:,3)
    v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
    mat(:,1) = v2
    mat(:,2) = cross(v1,v2)
    mat(:,3) = cross(v2,v3)
    x = dot_product(mat(:,2),mat(:,3))
    y=dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
    pep(i)%w = atan2(y,x)*180._dp/pi
    v1=chain(:,3)-chain(:,2); v2=chain(:,4)-chain(:,3); v3=chain(:,5)-chain(:,4)
    v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
    mat(:,1) = v2
    mat(:,2) = cross(v1,v2)
    mat(:,3) = cross(v2,v3)
    x = dot_product(mat(:,2),mat(:,3))
    y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
    pep(i)%phi = atan2(y,x)*180._dp/pi
    pep(i)%psi=0._dp

    i=1
    do j=1, size(pep(i)%backbone)
       if(trim(pep(i)%backbone(j)%label) == 'N') then
          chain(:,3) = pep(i)%backbone(j)%coords
       end if
       if(trim(pep(i)%backbone(j)%label) == 'CA') then
          chain(:,4) = pep(i)%backbone(j)%coords
       end if
       if(trim(pep(i)%backbone(j)%label) == 'C') then
          chain(:,5) = pep(i)%backbone(j)%coords
       end if
    end do
    do j=1,size(pep(i+1)%backbone)
       if(trim(pep(i+1)%backbone(j)%label) == 'N') then
          chain(:,6) = pep(i+1)%backbone(j)%coords
       end if
    end do
    v1=chain(:,4)-chain(:,3); v2=chain(:,5)-chain(:,4); v3=chain(:,6)-chain(:,5)
    v1=v1/norm2(v1); v2=v2/norm2(v2);v3=v3/norm2(v3)
    mat(:,1) = v2
    mat(:,2) = cross(v1,v2)
    mat(:,3) = cross(v2,v3)
    x = dot_product(mat(:,2),mat(:,3))
    y = dot_product(mat(:,1),cross(mat(:,2),mat(:,3)))
    pep(i)%psi = atan2(y,x)*180._dp/pi
    pep(i)%phi = 0._dp
    pep(i)%w = 0._dp
  end subroutine calc_dihedral
  ! ----------------------------------------------------------------------------

  subroutine change_dihedral(res,d_angle, option)
    ! Changes the dihedral by a value of delta
    type(residue), intent(inout) :: res             !Single residue of a protein
    character(len=*), intent(in) :: option          !Option for angle type
    character(len=10) :: angle_type                 !Internal angle type
    real(dp),intent(in) ::  d_angle                 !Change in angle (radians)

    angle_type = option
    call upper_case(angle_type)
    select case(trim(angle_type))
    case('OMEGA')
       res%w = res%w + d_angle*180._dp/pi
       if (res%w>180.) then
          res%w = res%w - 360._dp
       else if(res%w<-180.) then
          res%w = res%w + 360._dp
       end if
    case('PSI')
       res%psi = res%psi + d_angle*180._dp/pi
       if (res%psi>180.) then
          res%psi = res%psi - 360._dp
       else if(res%psi<-180.) then
          res%psi = res%psi + 360._dp
       end if
    case('PHI')
       res%phi = res%phi + d_angle*180._dp/pi
       if (res%phi>180.) then
          res%phi =  res%phi - 360._dp
       else if(res%phi<-180.) then
          res%phi = res%phi +360._dp
       end if
    end select
  end subroutine change_dihedral
  ! -----------------------------------------------------------------------------

  function get_dihedral(res,option) result(r)
    type(residue), intent(inout) :: res              !Single residue of a protein
    character(len=*), intent(in) :: option           !Option for angle type
    real(dp) :: r                                    !Dihedral angle (degrees)
    character(len=10) :: angle_type                  !Local angle type
    angle_type = option
    call upper_case(angle_type)
    select case(trim(angle_type))
    case('OMEGA')
       r=res%w
    case('PHI')
       r=res%phi
    case('PSI')
       r=res%psi
    case default
       r=0.0_dp
    end select
  end function get_dihedral
  ! -----------------------------------------------------------------------------

  logical function bb_list(label) result(r)
    ! Checks whether atom label is within know list of backbone labels
    character(len=4), intent(in) :: label
    r = .false.
    if(label=='N   ')r=.true.
    if(label=='CA  ')r=.true.
    if(label=='C   ')r=.true.
    if(label=='O   ')r=.true.
    if(label=='HA  ')r=.true.
    if(label=='HA1 ')r=.true.
    if(label=='HA2 ')r=.true.
    if(label=='HA3 ')r=.true.
    if(label=='H1  ')r=.true.
    if(label=='H2  ')r=.true.
    if(label=='H3  ')r=.true.
    if(label=='OXT ')r=.true.
    if(label=='H   ')r=.true.
  end function bb_list
  ! -----------------------------------------------------------------------------

  function count_spec(pep,sym,bbonly) result(r)
    ! Counts the number of atoms of a certain type in a protein (N_a)
    character(len=2), intent(in) :: sym                     !Atomic symbol
    type(residue), intent(in) :: pep(:)                     !Protein
    logical, optional :: bbonly                             !For backbone only calculations
    integer :: i,j                                          !Dummy integers
    integer :: r                                            !Return count

    r=0
    if(present(bbonly) .and. bbonly) then
       do i=1,size(pep)
          do j=1,size(pep(i)%backbone)
             if (get_sym(pep(i)%backbone(j)) == sym) then
                r = r + 1
             end if
          end do
       end do
    else
       do i=1,size(pep)
          do j=1,size(pep(i)%backbone)
             if (get_sym(pep(i)%backbone(j)) == sym) then
                r = r + 1
             end if
          end do
          do j=1,size(pep(i)%sidechain)
             if (get_sym(pep(i)%sidechain(j)) == sym) then
                r = r+ 1
             end if
          end do
       end do
    end if
  end function count_spec
  !-----------------------------------------------------------------------------

  function sym2num(sym) result(r)
    ! Takes in atomic symbol and sends back the corresponding index
    integer :: r
    character(len=2), intent(in) :: sym

    select case(sym)
    case('C ')
       r=1
    case('O ')
       r=2
    case('N ')
       r=3
    case('H ')
       r=4
    case('S ')
       r=5
    end select
  end function sym2num
  !-----------------------------------------------------------------------------

  function amino_list(label) result(r)
    ! Checks if label is technical amino acid for other funcitonality
    ! http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/
    character(len=3), intent(in) :: label
    logical :: r
    r=.false.
    select case(label)
    case('ALA')
       r = .true.
    case('ARG')
       r = .true.
    case('ASN')
       r= .true.
    case('ASP')
       r= .true.
    case('CYS')
       r= .true.
    case('GLN')
       r = .true.
    case('GLU')
       r = .true.
    case('GLY')
       r = .true.
    case('HIS')
       r = .true.
    case('ILE')
       r= .true.
    case('LEU')
       r= .true.
    case('LYS')
       r= .true.
    case('MET')
       r= .true.
    case('PHE')
       r= .true.
    case('PRO')
       r= .true.
    case('SER')
       r= .true.
    case('THR')
       r= .true.
    case('TRP')
       r= .true.
    case('TYR')
       r= .true.
    case('VAL')
       r= .true.
    case default
       r=.false.
    end select
  end function amino_list
  
  ! CURRENTLY UNUSED
  ! function aa_size(label) result(r)
  !   ! switches in for total amino acid residual size, after peptide bond formation
  !   ! Unit formula minus 3 (H2O)
  !   ! http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/
  !   character(len=3), intent(in) :: label
  !   integer :: r
  !   r=0
  !   select case(label)
  !   case('ALA')
  !      r = 10
  !   case('ARG')
  !      r = 23
  !   case('ASN')
  !      r= 14
  !   case('ASP')
  !      r= 13
  !   case('CYS')
  !      r= 11
  !   case('GLN')
  !      r = 17
  !   case('GLU')
  !      r = 16
  !   case('GLY')
  !      r = 7
  !   case('HIS')
  !      r = 17
  !   case('ILE')
  !      r= 19
  !   case('LEU')
  !      r= 19
  !   case('LYS')
  !      r= 21
  !   case('MET')
  !      r= 17
  !   case('PHE')
  !      r= 20
  !   case('PRO')
  !      r= 14
  !   case('SER')
  !      r= 10
  !   case('THR')
  !      r= 14
  !   case('TRP')
  !      r= 24
  !   case('TYR')
  !      r= 21
  !   case('VAL')
  !      r= 16
  !   case default
  !      write(*,*) "Inappropriate amino acid or residue label detected!"
  !   end select
  ! end function aa_size

  ! OBSOLETE WITH PUBLIC BAKCBONE AND SIDECHAIN
  ! function residue_size(res,opt) result(r)
  !   type(residue),intent(in) :: res
  !   character(len=4), intent(in), optional :: opt
  !   character(len=4) :: option
  !   integer :: r
  !   if(present(opt))then
  !      option=opt
  !      call upper_case(option)
  !      if(option == 'SIDE') then
  !         r = size(res%sidechain)
  !      else if (option == 'BACK') then
  !         r = size(res%backbone)
  !      else
  !         write(*,"(A)") 'Incorrect option given to funciton residue_size. Total size returned.'
  !      end if
  !   else
  !      r=size(res%backbone)+size(res%sidechain)
  !   end if
  ! end function residue_size
  ! ------------------------------------------------------------------------------


end module class_proteins
