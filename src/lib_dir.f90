module lib_dir
  ! Saved locations of libraries and dependents
  ! Needs to be updated for each installation

  implicit none

  public

  ! Installation directory
  character(len=100),parameter,public :: install_dir='/Users/alggroup/Documents/Development/SCD/'

  ! Geometry library from AMBER parm94 PDB file
  character(len=100),parameter,public :: geom_lib= trim(install_dir)//'lib/protein.amber94.edit.pdb'

  ! Using table from CCP14 for dihedral library
  character(len=100),parameter,public :: chi_lib= trim(install_dir)//'lib/sc_dih.lib'

  ! Van der Waals radii are implemented in a sensible ratio
  character(len=100),parameter,public :: vdw_lib= trim(install_dir)//'lib/vdwr.lib'

  ! Covalent radii are implemented in place of Van der Waals with a fudje
  character(len=100),parameter,public :: cov_lib= trim(install_dir)//'lib/covr.lib'


end module lib_dir
