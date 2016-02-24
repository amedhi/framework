module latticelibrary
!---------------------------------------------------------------------
! Version    : 2.0.0
! Description: Contains definiations of lattice.
!
! Author     : Amal Medhi <amedhi@physics.iisc.ernet.in>
! Creation   : September, 2010
! Revision   : October, 2011
!---------------------------------------------------------------------
use constants
use errstat
use utils
use taskparm
!
implicit none
!
private

integer(int32), parameter, public :: SQUARE = 1
integer(int32), parameter, public :: TWOSITE_SQUARE = 2
integer(int32), parameter, public :: BILAYER_SQUARE = 3
integer(int32), parameter, public :: HONEYCOMB = 10
integer(int32), parameter, public :: SW_HONEYCOMB = 11
integer(int32), parameter, public :: CHAIN = 15
integer(int32), parameter, public :: CHAIN_BCC = 16

integer(int32), parameter :: max_natom = 10
integer(int32), parameter :: max_nbond = 20
integer(int32), parameter :: max_nn = 4

! module wide error status
integer(int32), private :: minfo = ok

! 'atom' definition
type, public ::  atom_t
  integer(int32) :: id
  integer(int32) :: type
  integer(int32) :: atomid
  integer(int32) :: offset(3)
  real(dp) :: coord(3)
  real(dp) :: cellcoord(3)
end type atom_t

!'bond' definition
type, public ::  bond_t
  integer(int32) :: neighb
  integer(int32) :: type
  type (atom_t) :: source
  type (atom_t) :: target
  real(dp) :: vector(3)
end type bond_t

!unit cell definition
type, public :: unitcell_t
  integer(int32) :: natom
  integer(int32) :: nbond
  integer(int32) :: max_atom_type_val
  integer(int32) :: max_bond_type_val
  integer(int32) :: max_neighb
  integer(int32) :: nbond_total
  type (atom_t), allocatable :: atom(:)
  type (bond_t), allocatable :: bond(:)
contains
  procedure, private :: create
  procedure, private :: add_atom
  procedure, private :: add_bond
end type unitcell_t

!lattice 
!basis vector
type, public :: basisvector_t
  real(dp) :: a1(3), a2(3), a3(3)
end type basisvector_t

! each dimension
type, public :: spatialdim_t
  integer(int32) :: size
  character(len=12) :: boundary
  character(len=12) :: periodicity
end type spatialdim_t

! lattice
type, public :: lattice_t
  type (unitcell_t) :: unitcell
  type (spatialdim_t) :: dim1
  type (spatialdim_t) :: dim2
  type (spatialdim_t) :: dim3
  type (basisvector_t) :: basis_vector
  !integer(int32) :: max_neighb
  integer(int32) :: id
  integer(int32) :: ldim
  character(len=30) :: name
contains
  procedure, public :: init => init_sub
  procedure, private :: construct
end type lattice_t

contains
!---------------------------------------------------------------------
!
subroutine init_sub(lattice, parms, info)
  !Purpose: 
  !Arguments:
  !1.
  implicit none

  ! calling parameters
  class (lattice_t) :: lattice
  type (taskparm_t), intent(in) :: parms
  integer(int32), intent(out) :: info

  ! local varibales
  integer(int32) :: ldim, lid, n
  integer(int32) :: lsize1, lsize2, lsize3, max_neighb
  character(len=30) :: lattice_name, lattice_long_name
  character(len=12) :: boundary1, boundary2, boundary3
  real(dp) :: aa, bb, cc, x, y
  real(dp) :: av1(3), av2(3), uv1(3), uv2(3), uv3(3), coord(3)
  type (basisvector_t) :: basis_vector
  type (spatialdim_t) :: dim1, dim2, dim3
  type (unitcell_t) :: unitcell

  ! lattice constrants
  aa = 1.0_dp; bb = aa; cc = aa;

  ! read lattice attributes
  lattice_name = parms%set_value("lattice", "NULL", info)
  lattice_long_name = to_upper(lattice_name)
  lsize1 = parms%set_value("lsize1", 1, info)
  lsize2 = parms%set_value("lsize2", 1, info)
  lsize3 = parms%set_value("lsize3", 1, info)

  ! lattcie short name
  n = index(lattice_long_name, "(NN")
  if (n > 0) then
    lattice_name = lattice_long_name(1:n-1)
  else
    lattice_name = lattice_long_name
  endif

  ! max neighb
  max_neighb = 1
  n = len_trim(lattice_long_name)
  if (index(lattice_long_name, "(NN0)") > 0) then
    max_neighb = 0
  else if (index(lattice_long_name, "(NN1)") > 0) then
    max_neighb = 1
  else if (index(lattice_long_name, "(NN2)") > 0) then
    max_neighb = 2
  else if (index(lattice_long_name, "(NN3)") > 0) then
    max_neighb = 3
  else if (index(lattice_long_name, "(NN") > 0) then
    call errmsg(error, "latticelibrary::init", "undefined NN specification in '"//trim(lattice_long_name)//"'")
    info=error; return
  endif

  ! boundary conditions
  boundary1 = parms%set_value("bc1", "periodic", info)
  boundary2 = parms%set_value("bc2", "periodic", info)
  boundary3 = parms%set_value("bc3", "periodic", info)

  ! checks
  if (lsize1<1 .or. lsize2<1 .or. lsize3<1) then 
    call errmsg(error, "latticelibrary::init", "invaild lattice size")
    info = error; return
  end if

  select case(to_upper(trim(boundary1)))
  case ("PERIODIC")
    boundary1 = "periodic"
  case ("ANTIPERIODIC")
    boundary1 = "antiperiodic"
  case ("TILTPERIODIC")
    boundary1 = "tiltperiodic"
  case ("OPEN")
    boundary1 = "open"
  case default
    call errmsg(error, "latticelibrary::define", "Unrecognized boundary condition, &
                        '"//trim(boundary1)//"'")
    info=error; return
  end select
  select case(to_upper(trim(boundary2)))
  case ("PERIODIC")
    boundary2 = "periodic"
  case ("ANTIPERIODIC")
    boundary2 = "antiperiodic"
  case ("TILTPERIODIC")
    boundary2 = "tiltperiodic"
  case ("OPEN")
    boundary2 = "open"
  case default
    call errmsg(error, "latticelibrary::define", "Unrecognized boundary condition, &
                        '"//trim(boundary2)//"'")
    info=error; return
  end select
  select case(to_upper(trim(boundary3)))
  case ("PERIODIC")
    boundary3 = "periodic"
  case ("ANTIPERIODIC")
    boundary3 = "antiperiodic"
  case ("TILTPERIODIC")
    boundary3 = "tiltperiodic"
  case ("OPEN")
    boundary3 = "open"
  case default
    call errmsg(error, "latticelibrary::define", "Unrecognized boundary condition, &
                        '"//trim(boundary3)//"'")
    info=error; return
  end select
  ! in case of unit size
  if (lsize1 == 1) boundary1 = "open"
  if (lsize2 == 1) boundary2 = "open"
  if (lsize3 == 1) boundary3 = "open"

  ! allocate
  if (allocated(unitcell%atom)) deallocate(unitcell%atom)
  allocate (unitcell%atom(max_natom), stat=info)
  call catch(info, alloc_t, "define_lattice")

  if (allocated(unitcell%bond)) deallocate(unitcell%bond)
  allocate (unitcell%bond(max_nbond), stat=info)
  call catch(info, alloc_t, "define_lattice")

  ! define lattices
  minfo = ok
  !------------SWITCH OFF ALL EXTRA 2ND NEIGHBOUR LINKS-----------
  select case(to_upper(lattice_name))
  case ("SQUARE")
    !lattice id & dimension
    lid = SQUARE
    ldim = 2

    ! basis vectors
    basis_vector%a1 = [aa,     0.0_dp, 0.0_dp]
    basis_vector%a2 = [0.0_dp, aa,     0.0_dp]
    basis_vector%a3 = [0.0_dp, 0.0_dp, 0.0_dp]

    ! extent
    dim1%size = lsize1
    dim2%size = lsize2
    dim3%size = 1
    dim1%boundary = boundary1
    dim2%boundary = boundary2
    dim3%boundary = "open"

    ! create an empty unitcell 
    n = 2*max_neighb
    call unitcell%create(atoms=1, bonds=n, max_neighb=max_neighb)

    ! add atoms
    call unitcell%add_atom(type=1, at_site=1, coord=[0.0_dp, 0.0_dp, 0.0_dp])

    ! add bonds
    if (max_neighb >= 1) then
      call unitcell%add_bond(type=1, ngb=1, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[1,0,0])
      call unitcell%add_bond(type=2, ngb=1, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[0,1,0])
    endif

    if (max_neighb >= 2) then
      call unitcell%add_bond(type=3, ngb=2, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[1,1,0])
      call unitcell%add_bond(type=4, ngb=2, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[-1,1,0])
    endif

    if (max_neighb >= 3) then
      call unitcell%add_bond(type=5, ngb=3, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[2,0,0])
      call unitcell%add_bond(type=6, ngb=3, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[0,2,0])
    endif
    !-----SQUARE------

  case ("TWOSITE SQUARE")
    !lattice id & dimension
    lid = TWOSITE_SQUARE
    ldim = 2

    ! basis vectors
    x = sqrt(2.0_dp)*aa
    basis_vector%a1 = [x,      0.0_dp, 0.0_dp]
    basis_vector%a2 = [0.0_dp, x,     0.0_dp]
    basis_vector%a3 = [0.0_dp, 0.0_dp, 0.0_dp]

    ! extent
    dim1%size = lsize1
    dim2%size = lsize2
    dim3%size = 1
    dim1%boundary = boundary1
    dim2%boundary = boundary2
    dim3%boundary = "open"

    ! create an empty unitcell 
    if (max_neighb /= 1) then
      call errmsg(error, "latticelibrary:init", "this lattice is not defined for max_neighb", max_neighb)
      info = error; return
    endif
    call unitcell%create(atoms=2, bonds=4, max_neighb=1)

    ! add atoms
    x = 0.5_dp*sqrt(2.0_dp)*aa
    call unitcell%add_atom(type=1, at_site=1, coord=[0.0_dp, 0.0_dp, 0.0_dp])
    call unitcell%add_atom(type=2, at_site=2, coord=[x, x, 0.0_dp])

    ! add bonds
    call unitcell%add_bond(type=1, ngb=1, src=1, src_offset=[0,0,0], tgt=2, tgt_offset=[0,0,0])
    call unitcell%add_bond(type=2, ngb=1, src=2, src_offset=[0,0,0], tgt=1, tgt_offset=[1,0,0])
    call unitcell%add_bond(type=3, ngb=1, src=2, src_offset=[0,0,0], tgt=1, tgt_offset=[0,1,0])
    call unitcell%add_bond(type=4, ngb=1, src=2, src_offset=[0,0,0], tgt=1, tgt_offset=[1,1,0])
    !-----TWOSITE_SQUARE------

  case ("HONEYCOMB")
    ! lattice id & dimension
    lid = HONEYCOMB
    ldim = 2

    ! basis vectors
    ! aa = A-A or B-B atom distance
    x = 0.5_dp*aa; y = -0.5_dp*sqrt(3.0_dp)*aa 
    basis_vector%a1 = [x, y, 0.0_dp]
    x = 0.5_dp*aa; y = +0.5_dp*sqrt(3.0_dp)*aa 
    basis_vector%a2 = [x, y, 0.0_dp]
    basis_vector%a3 = [0.0_dp, 0.0_dp, 0.0_dp]

    ! extent
    dim1%size = lsize1
    dim2%size = lsize2
    dim3%size = 1
    dim1%boundary = boundary1
    dim2%boundary = boundary2
    dim3%boundary = "open"

    ! create an empty unitcell 
    select case(max_neighb)
    case (0); n = 0; case (1); n = 3; case (2); n = 9;
    case default
      call errmsg(error, "latticelibrary:init", "this lattice is not defined for max_neighb", max_neighb)
      info = error; return
    end select
    call unitcell%create(atoms=2, bonds=n, max_neighb=max_neighb)

    ! add atoms
    call unitcell%add_atom(type=1, at_site=1, coord=[0.0_dp, 0.0_dp, 0.0_dp])
    y = aa/sqrt(3.0_dp);
    call unitcell%add_atom(type=2, at_site=2, coord=[0.0_dp, y, 0.0_dp])
      
    ! add bonds
    if (max_neighb >= 1) then
      ! 1st A-B bond
      call unitcell%add_bond(type=1, ngb=1, src=1, src_offset=[0,0,0], tgt=2, tgt_offset=[0,0,0])
      ! 2nd A-B bond
      call unitcell%add_bond(type=2, ngb=1, src=1, src_offset=[0,0,0], tgt=2, tgt_offset=[1,0,0])
      ! 3rd A-B bond
      call unitcell%add_bond(type=3, ngb=1, src=2, src_offset=[0,0,0], tgt=1, tgt_offset=[0,1,0])
    endif

    if (max_neighb >= 2) then
      ! 1st A-A bond
      call unitcell%add_bond(type=4, ngb=2, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[1,0,0])
      ! 2nd A-A bond
      call unitcell%add_bond(type=5, ngb=2, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[0,1,0])
      ! 3rd A-A bond
      call unitcell%add_bond(type=6, ngb=2, src=1, src_offset=[0,0,0], tgt=1, tgt_offset=[1,1,0])
      ! 1st B-B bond
      call unitcell%add_bond(type=7, ngb=2, src=2, src_offset=[0,0,0], tgt=2, tgt_offset=[1,0,0])
      ! 2nd B-B bond
      call unitcell%add_bond(type=8, ngb=2, src=2, src_offset=[0,0,0], tgt=2, tgt_offset=[0,1,0])
      ! 3rd B-B bond
      call unitcell%add_bond(type=9, ngb=2, src=2, src_offset=[0,0,0], tgt=2, tgt_offset=[1,1,0])
    endif
    !-----HONEYCOMB------

  case ("SW_HONEYCOMB")
    ! lattice id & dimension
    lid = SW_HONEYCOMB
    ldim = 2

    ! define a few vectors
    ! aa = A-A or B-B atom distance in HONEYCOMB
    x = aa; y = sqrt(3.0_dp)*aa ! double the values for HONEYCOMB
    av1 = [x, -y, 0.0_dp]
    av2 = [x, +y, 0.0_dp]

    ! the 3 nn vectors in HONEYCOMB
    x = 0.0_dp; y = aa/sqrt(3.0_dp)
    uv1 = [x, y, 0.0_dp]
    x = 0.5_dp*aa;  y = (0.5_dp * aa)/sqrt(3.0_dp)
    uv2 = [+x, -y, 0.0_dp]
    uv3 = [-x, -y, 0.0_dp]

    ! basis vectors
    basis_vector%a1 = av1
    basis_vector%a2 = av2
    basis_vector%a3 = [0.0_dp, 0.0_dp, 0.0_dp]

    ! extent
    dim1%size = lsize1
    dim2%size = lsize2
    dim3%size = 1
    dim1%boundary = boundary1
    dim2%boundary = boundary2
    dim3%boundary = "open"

    ! create an empty unitcell
    select case(max_neighb)
    case (0); n = 0; case (1); n = 12; 
    case default
      call errmsg(error, "latticelibrary:init", "this lattice is not defined for max_neighb", max_neighb)
      info = error; return
    end select
    call unitcell%create(atoms=8, bonds=n, max_neighb=max_neighb)

    ! add atoms
    coord = [0.0_dp, 0.0_dp, 0.0_dp]
    call unitcell%add_atom(type=1, at_site=1, coord=coord)
    coord = uv1
    call unitcell%add_atom(type=2, at_site=2, coord=coord)
    coord = 0.5_dp * av1
    call unitcell%add_atom(type=3, at_site=3, coord=coord)
    coord = 0.5_dp * av1 + uv1
    call unitcell%add_atom(type=4, at_site=4, coord=coord)
    coord = 0.5_dp * av2 
    call unitcell%add_atom(type=5, at_site=5, coord=coord)
    coord = 0.5_dp * av2 + uv1
    call unitcell%add_atom(type=6, at_site=6, coord=coord)
    x = (0.5_dp * aa)/sqrt(3.0_dp)
    coord = 0.5_dp * (av1 + av2) + [-x, +x, 0.0_dp]
    call unitcell%add_atom(type=7, at_site=7, coord=coord)
    coord = 0.5_dp * (av1 + av2) + uv1 + [+x, -x, 0.0_dp]
    call unitcell%add_atom(type=8, at_site=8, coord=coord)

    ! add bonds
    if (max_neighb >= 1) then
      !------------intra-cell bonds------------
      call unitcell%add_bond(type=1, ngb=1, src=1, src_offset=[0,0,0], tgt=2, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=2, ngb=1, src=3, src_offset=[0,0,0], tgt=4, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=3, ngb=1, src=5, src_offset=[0,0,0], tgt=6, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=4, ngb=1, src=7, src_offset=[0,0,0], tgt=8, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=5, ngb=1, src=1, src_offset=[0,0,0], tgt=4, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=6, ngb=1, src=2, src_offset=[0,0,0], tgt=5, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=7, ngb=1, src=4, src_offset=[0,0,0], tgt=7, tgt_offset=[0,0,0])
      call unitcell%add_bond(type=8, ngb=1, src=5, src_offset=[0,0,0], tgt=7, tgt_offset=[0,0,0])
      !------------inter-cell bonds------------
      call unitcell%add_bond(type=9, ngb=1, src=3, src_offset=[0,0,0], tgt=2, tgt_offset=[1,0,0])
      call unitcell%add_bond(type=10, ngb=1, src=6, src_offset=[0,0,0], tgt=1, tgt_offset=[0,1,0])
      call unitcell%add_bond(type=11, ngb=1, src=8, src_offset=[0,0,0], tgt=6, tgt_offset=[1,0,0])
      call unitcell%add_bond(type=12, ngb=1, src=8, src_offset=[0,0,0], tgt=3, tgt_offset=[0,1,0])
    endif
    !-----SW HONEYCOMB------
  case default
    call errmsg(error, "latticelibrary:init", "lattice "//trim(lattice_long_name)//" not found" )
    info = error; return
  end select
  if (minfo /= ok) return

  ! construct
  call lattice%construct(lattice_long_name, lid, ldim, basis_vector, dim1, dim2, dim3, unitcell, info)

  return
end subroutine init_sub
!---------------------------------------------------------------------
!
subroutine create(this, atoms, bonds, max_neighb) 
  ! Purpose: 
  ! Arguments:
  !1.
  implicit none

  ! calling parameters
  class (unitcell_t) :: this
  integer(int32), intent(in) :: atoms
  integer(int32), intent(in) :: bonds
  integer(int32), intent(in) :: max_neighb

  ! checks
  if (atoms < 1) then
    call errmsg(alert, "latticelibrary::unitcell%create", "atom number must be > 0")
    minfo = error; return
  endif
  if (bonds < 1) then
    call errmsg(alert, "latticelibrary::unitcell%create", "bond number must be > 0")
    minfo = error; return
  endif
  if (max_neighb < 1) then
    call errmsg(alert, "latticelibrary::unitcell%create", "max_neighb must be > 0")
    minfo = error; return
  endif

  ! allocate
  if (allocated(this%atom)) deallocate(this%atom)
  allocate (this%atom(atoms), stat=minfo)
  call catch(minfo, alloc_t, "latticelibrary::unitcell%create")

  if (allocated(this%bond)) deallocate(this%bond)
  allocate (this%bond(bonds), stat=minfo)
  call catch(minfo, alloc_t, "latticelibrary::unitcell%create")

  ! at this time, it is a empty unit cell
  this%natom = 0
  this%nbond = 0
  this%max_atom_type_val = 0
  this%max_bond_type_val = 0
  this%max_neighb = max_neighb
  ! TEMPORARY MEASURE
  this%nbond_total = bonds

  return
end subroutine create
!---------------------------------------------------------------------
!
subroutine add_atom(this, type, at_site, coord) 
  ! Purpose: 
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (unitcell_t) :: this
  integer(int32), intent(in) :: type
  integer(int32), intent(in) :: at_site
  real(dp), intent(in) :: coord(:)

  ! local variables
  integer(int32) :: n

  ! checks
  if (type < 1) then
    call errmsg(alert, "latticelibrary::unitcell%add_atom", "only 'type value > 0' allowed")
    minfo = error; return
  endif
  if (this%natom == size(this%atom)) then
    call errmsg(error, "latticelibrary:unitcell%add_atom", "can't add, this atom exceeds the limit set")
    minfo = error; return
  endif

  if (this%max_atom_type_val < type) this%max_atom_type_val = type
  this%natom = this%natom + 1
  n = this%natom
  this%atom(n)%id = n
  this%atom(n)%type = type
  this%atom(n)%atomid = at_site
  this%atom(n)%coord(1) = coord(1)
  this%atom(n)%coord(2) = coord(2)
  this%atom(n)%coord(3) = coord(3)
  ! the above asignment is done this way because of a 'ifort bug' 
  this%atom(n)%offset = [0,0,0]
  this%atom(n)%cellcoord = [0.0_dp,0.0_dp,0.0_dp]

  return
end subroutine add_atom
!---------------------------------------------------------------------
!
subroutine add_bond(this, type, ngb, src, src_offset, tgt, tgt_offset) 
  ! Purpose: 
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (unitcell_t) :: this
  integer(int32), intent(in) :: type
  integer(int32), intent(in) :: ngb
  integer(int32), intent(in) :: src
  integer(int32), intent(in) :: src_offset(:)
  integer(int32), intent(in) :: tgt
  integer(int32), intent(in) :: tgt_offset(:)

  ! local variables
  integer(int32) :: n

  ! checks
  if (type < 1) then
    call errmsg(alert, "latticelibrary::unitcell%add_bond", "only 'type value > 0' allowed")
    minfo = error; return
  endif
  if (ngb > this%max_neighb) then
    call errmsg(error, "latticelibrary:unitcell%add_bond", "ngb exceeds max_neighb value set")
    minfo = error; return
  endif
  if (src > this%natom) then
    call errmsg(error, "latticelibrary:unitcell%add_bond", "source for this bond does not exist")
    minfo = error; return
  endif
  if (tgt > this%natom) then
    call errmsg(error, "latticelibrary:unitcell%add_bond", "target for this bond does not exist")
    minfo = error; return
  endif
  if (this%nbond == size(this%bond)) then
    call errmsg(error, "latticelibrary:unitcell%add_bond", "can't add, this bond exceeds the limit set")
    minfo = error; return
  endif

  if (this%max_bond_type_val < type) this%max_bond_type_val = type
  this%nbond = this%nbond + 1
  n = this%nbond
  this%bond(n)%neighb = ngb
  this%bond(n)%type = type
  this%bond(n)%source = this%atom(src)
  this%bond(n)%source%offset(1) = src_offset(1)
  this%bond(n)%source%offset(2) = src_offset(2)
  this%bond(n)%source%offset(3) = src_offset(3)

  this%bond(n)%target = this%atom(tgt)
  this%bond(n)%target%offset(1) = tgt_offset(1)
  this%bond(n)%target%offset(2) = tgt_offset(2)
  this%bond(n)%target%offset(3) = tgt_offset(3)
  ! the above asignment is done this way because of a 'ifort bug'

  !write(*,*) this%bond(n)%target%offset, tgt_offset
  !read(*,*)

  return
end subroutine add_bond
!---------------------------------------------------------------------
!
function make_atom(id, type, atomid, coord) result(atom)
  ! Purpose: 
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  integer(int32), intent(in) :: id
  integer(int32), intent(in) :: type
  integer(int32), intent(in) :: atomid
  real(dp), intent(in) :: coord(:)
  ! result
  type (atom_t) :: atom

  atom%id = id
  atom%type = type
  atom%atomid = atomid
  atom%coord = coord
  atom%offset = [0,0,0]
  atom%cellcoord = [0.0_dp,0.0_dp,0.0_dp]

return
end function make_atom
!---------------------------------------------------------------------
!
function link_atom(atom, offset) result(linked_atom)
  ! Purpose: 
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  type (atom_t), intent(in) :: atom
  integer(int32), intent(in) :: offset(:)
  ! result
  type (atom_t) :: linked_atom

  linked_atom = atom
  linked_atom%offset = offset

return
end function link_atom
!---------------------------------------------------------------------
!
subroutine construct_old(lattice, name, lid, ldim, basis_vector, dim1, dim2, dim3, unitcell) 
  ! Purpose: 
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (lattice_t) :: lattice
  character(len=*) :: name
  integer(int32), intent(in) :: lid
  integer(int32), intent(in) :: ldim
  type (basisvector_t), intent(in) :: basis_vector
  type (spatialdim_t), intent(in) :: dim1
  type (spatialdim_t), intent(in) :: dim2
  type (spatialdim_t), intent(in) :: dim3
  type (unitcell_t), intent(in) :: unitcell

  lattice%name = name
  lattice%id = lid
  lattice%ldim = ldim
  !lattice%max_neighb = max_neighb
  lattice%basis_vector = basis_vector
  lattice%dim1 = dim1
  lattice%dim2 = dim2
  lattice%dim3 = dim3
  lattice%unitcell = unitcell

  ! incorporate antiperiodic & tilted periodic boundary condition
  lattice%dim1%periodicity = dim1%boundary
  lattice%dim2%periodicity = dim2%boundary
  lattice%dim3%periodicity = dim3%boundary
  if (dim1%boundary == "antiperiodic") lattice%dim1%boundary = "periodic"
  if (dim2%boundary == "antiperiodic") lattice%dim2%boundary = "periodic"
  if (dim3%boundary == "antiperiodic") lattice%dim3%boundary = "periodic"

  if (dim1%boundary == "tiltperiodic") lattice%dim1%boundary = "periodic"
  if (dim2%boundary == "tiltperiodic") lattice%dim2%boundary = "periodic"
  if (dim3%boundary == "tiltperiodic") lattice%dim3%boundary = "periodic"

return
end subroutine construct_old
!---------------------------------------------------------------------
!
subroutine construct(lattice, name, lid, ldim, basis_vector, dim1, dim2, dim3, unitcell, info) 
  ! Purpose: 
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (lattice_t) :: lattice
  character(len=*) :: name
  integer(int32), intent(in) :: lid
  integer(int32), intent(in) :: ldim
  type (basisvector_t), intent(in) :: basis_vector
  type (spatialdim_t), intent(in) :: dim1
  type (spatialdim_t), intent(in) :: dim2
  type (spatialdim_t), intent(in) :: dim3
  type (unitcell_t), intent(in) :: unitcell
  integer(int32), intent(out) :: info

  ! local variables
  integer(int32) :: i
  integer(int32), allocatable :: iwork(:)

  lattice%name = name
  lattice%id = lid
  lattice%ldim = ldim
  !lattice%max_neighb = max_neighb
  lattice%basis_vector = basis_vector
  lattice%dim1 = dim1
  lattice%dim2 = dim2
  lattice%dim3 = dim3

  ! incorporate antiperiodic & tilted periodic boundary condition
  lattice%dim1%periodicity = dim1%boundary
  lattice%dim2%periodicity = dim2%boundary
  lattice%dim3%periodicity = dim3%boundary
  if (dim1%boundary == "antiperiodic") lattice%dim1%boundary = "periodic"
  if (dim2%boundary == "antiperiodic") lattice%dim2%boundary = "periodic"
  if (dim3%boundary == "antiperiodic") lattice%dim3%boundary = "periodic"

  if (dim1%boundary == "tiltperiodic") lattice%dim1%boundary = "periodic"
  if (dim2%boundary == "tiltperiodic") lattice%dim2%boundary = "periodic"
  if (dim3%boundary == "tiltperiodic") lattice%dim3%boundary = "periodic"

  ! check & add unitcell 
  if (unitcell%natom /= size(unitcell%atom)) then
    call errmsg(error, "latticelibrary::construct", "you have not added all the atoms")
    info = error; return
  endif
  if (unitcell%nbond /= size(unitcell%bond)) then
    call errmsg(error, "latticelibrary::construct", "you have not added all the bonds")
    info = error; return
  endif
  ! see if atom type values are in contiguous range
  if (allocated(iwork)) deallocate(iwork)
  allocate(iwork(unitcell%max_atom_type_val))
  iwork(:) = 0
  forall (i=1:unitcell%natom) iwork(unitcell%atom(i)%type) = i
  do i = 1, unitcell%max_atom_type_val
    if (iwork(i) == 0) then
      call errmsg(error, "latticelibrary::construct", "atom 'type values' are not in contiguous range")
      info = error; return
    endif
  enddo
  ! see if bond type values are in contiguous range
  if (allocated(iwork)) deallocate(iwork)
  allocate(iwork(unitcell%max_bond_type_val))
  iwork(:) = 0
  forall (i=1:unitcell%nbond) iwork(unitcell%bond(i)%type) = i
  do i = 1, unitcell%max_bond_type_val
    if (iwork(i) == 0) then
      call errmsg(error, "latticelibrary::construct", "bond 'type values' are not in contiguous range")
      info = error; return
    endif
  enddo
  ! add unitcell
  lattice%unitcell = unitcell

  info = ok
return
end subroutine construct
!---------------------------------------------------------------------
!
end module latticelibrary
