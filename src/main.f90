program main

  ! The following will 
  use lib_dir
  use types
  use maths
  use class_proteins
  use sidechain_decoration
  use pdb_io

  implicit none

  type(protein) :: prot                                     !Protein
  character(len=100) :: out_fname                           !Output filename
  character(len=100) :: in_fname                            !Input filename
  character(len=100) :: errmsg                              !Error messages
  character(len=100) :: contact_error                       !Contact error message
  integer :: info                                           !Info integer
  integer :: i                                              !Dummy integer

  info=0
  contact_error='Now seems like a good time to complain to Phil. phillip[dot]maffettone[at]chem[dot]ox[dot]ax[dot]uk'
  
  i = command_argument_count()
  if (i<1) then
     write(*,*)'Please start the program by writing ./SCD [input file name]'
     write(*,*)contact_error
     stop
  else if (i>=2) then
     call get_command_argument(1,in_fname)
     call get_command_argument(2,out_fname)
  else
     call get_command_argument(1,in_fname)
     out_fname = trim(install_dir)//'SCDOutput/output.pdb'
     call system('mkdir -p '//trim(install_dir)//'SCDOutput')
  end if

  call read_pdb(in_fname,prot,bbonly=.true.,info=info,errmsg=errmsg)
  if(info/=0) then
     write(*,*)errmsg
     write(*,*)in_fname
     write(*,*)contact_error
     stop
  end if

  call decorate(prot)

  call write_pdb(out_fname,prot,.false.)

  write(*,*)
  write(*,"('Decorated protein located at: ',A80)")out_fname
  write(*,*)
end program main
