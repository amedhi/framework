module latticegraph
!-----------------------------------------------------------------------
! Version    : 2.0.0
! Description: Provides definition of lattice. Creates graph of a 
!              given lattice. It transforms the input lattice by 
!              extending the underlying unitcell to respect lattice 
!              translational symmetry. Stores the transformed lattice &
!              make it public along with the graph. All subsequent
!              operations on graph elements is defined with repect to
!              this transformed lattice and NOT the original one.
!
! Author     : Amal Medhi <amedhi@physics.iisc.ernet.in>
! Creation   : September, 2010
! Revision   : October, 2011
!-----------------------------------------------------------------------
use constants
use errstat
use latticelibrary
!
implicit none
! 
private

! lattice graph definition
! vertex 
type, private ::  vertex_t
  integer(int32) :: id
  integer(int32) :: uid
  integer(int32) :: type
  integer(int32) :: atomid
  integer(int32) :: class
  integer(int32) :: ndirlink
  integer(int32), allocatable :: link_id(:)
  integer(int32) :: cellpos(3)
  integer(int32) :: n1neighb
  integer(int32) :: n2neighb
  integer(int32), allocatable :: nb_list1(:)
  integer(int32), allocatable :: nb_list2(:)
  real(dp) :: coord(3)
  real(dp) :: cellcoord(3)
end type vertex_t

!edge
type, private ::  edge_t
  integer(int32) :: id
  integer(int32) :: type
  integer(int32) :: class
  integer(int32) :: phase
  integer(int32) :: src_id
  integer(int32) :: tgt_id
  real(dp) :: vector(3)
end type edge_t

type, private ::  reflect_t
  integer(int32) :: at
  integer(int32) :: offset
end type reflect_t

!graph
type, public :: latticegraph_t
  type (vertex_t), allocatable, public :: vertex(:)
  type (edge_t), allocatable, public :: edge(:)
  integer(int32), public :: nvertex=0
  integer(int32), public :: nedge=0
  integer(int32), public :: nvclass=0
  integer(int32), public :: max_vertex_type = 0
  integer(int32), public :: max_edge_type = 0
  integer(int32), public :: netype=0
  type (lattice_t), public :: glattice 
  logical, public :: tilted_lattice = .false.

  private	  
  integer(int32) :: neclass=0
  integer(int32) :: npcell=0
  integer(int32) :: ncell=0
  integer(int32) :: atomid_range=0 !range of 'site id' per unitcell
  integer(int32) :: tilted_shift=0
  integer(int32) :: tilted_size1, tilted_size2
  integer(int32), allocatable :: tilted_bottom(:)
  integer(int32), allocatable :: tilted_top(:)

  !Implement rotattional symmetry -- rotate the 'lattice' clockwise
  !about the origin wrt an axis
  integer(int32) :: xrotate(3,3)
  integer(int32) :: yrotate(3,3)
  integer(int32) :: zrotate(3,3)

  !Implement reflection symmetry -- reflect a 'site' about
  !about an axis
  type (reflect_t) :: xreflect
  type (reflect_t) :: yreflect
  type (reflect_t) :: zreflect

  procedure (get_cell_number), pointer :: cell_number_fn => null()
contains
  procedure, public :: construct => construct_sub
  procedure, public :: translate_vertex => translate_vertex_fun
  procedure, public :: get_next_bravindex => get_next_bravindex_fun
  procedure, public :: print_unitgraph => print_unitgraph_sub
  procedure, public :: print_graph => print_graph_sub
  procedure, public :: get_class_siteid_map => get_class_siteid_map_sub
  procedure, public :: symmetry_reduce => symmetry_reduce_sub
  !procedure, public :: construct_nnlist => construct_nnlist_sub
  procedure, private :: check_validity
  procedure, private :: construct_glattice
  procedure, private :: get_cell_number
  procedure, private :: cell_number_fun_tiltbc
  procedure, private :: bravindex_tilted_lattice
  procedure, private :: get_vertex_id
  procedure, private :: translate_cell
  procedure, private :: boundary_wrap
  procedure, private :: boundary_wrap_tiltbc
  procedure, private :: get_vertex_class
  procedure, private :: get_edge_class
  procedure, private, nopass :: get_rotmatrix
  procedure, private, nopass :: serialize
end type latticegraph_t

contains
!
!---------------------------------------------------------------------
!
subroutine construct_sub(this, inlattice, info) 
  ! Purpose: First transform the input lattice by extending the 
  ! unitcell & other attributes to respect translational symmetry. 
  ! Creates graph of the transformed lattice.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t) :: this
  type (lattice_t), intent(in) :: inlattice
  integer(int32), intent(out) :: info

  ! local variables
  integer(int32) :: i, j, m, n
  integer(int32) :: id, src_id, tgt_id, natom, nbond, nbond_total, nvert, nedge
  integer(int32) :: s_phase, t_phase
  integer(int32) :: bravindex(3), cellpos(3), offset(3)
  integer(int32), allocatable :: vertex_class(:)
  integer(int32), allocatable :: edge_class(:)
  real(dp) :: coord(3)
  type(unitcell_t) :: unitcell
  type (vertex_t), allocatable :: vertex(:)
  type (edge_t), allocatable :: edge(:)

  ! check if vaild 'unitcell'
  call this%check_validity(inlattice, info)
  if (info /= ok) return
  
  ! store 'lattice' info, reorder attributes
  call this%construct_glattice( inlattice )

  ! short names
  natom = this%glattice%unitcell%natom
  nbond = this%glattice%unitcell%nbond
  nbond_total = this%glattice%unitcell%nbond_total
  unitcell = this%glattice%unitcell

  ! construct the vertex 'class'-id s
  if (allocated(vertex_class)) deallocate(vertex_class)
  allocate (vertex_class(natom), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )
  call this%get_vertex_class(unitcell, vertex_class, this%nvclass)

  ! construct the edge 'class'-id s
  if (allocated(edge_class)) deallocate(edge_class)
  allocate (edge_class(nbond), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )
  call this%get_edge_class(unitcell, edge_class, this%neclass)

  ! vertex construction: 1st stage
  nvert = (this%ncell * natom)
  if (allocated(vertex)) deallocate(vertex)
  allocate (vertex(nvert), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )

  bravindex = [1,1,1]
  cells: do i = 1, this%ncell
    ! translate cell
    offset = bravindex - [1,1,1]
    unitcell = this%translate_cell(this%glattice%unitcell, offset)

    ! add vertices
    atoms: do m = 1, natom
      n = this%get_vertex_id( unitcell%atom(m)%id, unitcell%atom(m)%offset )
      vertex(n)%id = n
      vertex(n)%uid = unitcell%atom(m)%id
      vertex(n)%type = unitcell%atom(m)%type
      vertex(n)%atomid = unitcell%atom(m)%atomid
      vertex(n)%class = vertex_class(m)
      vertex(n)%cellpos = unitcell%atom(m)%offset
      vertex(n)%coord  = unitcell%atom(m)%coord 
      vertex(n)%cellcoord  = unitcell%atom(m)%cellcoord 
    end do atoms

    ! next bravais index
    bravindex = this%get_next_bravindex(bravindex)
  end do cells

  ! edge construction: 1st stage
  nedge = (this%ncell * nbond) 
  if (allocated(edge)) deallocate(edge)
  allocate (edge(nedge), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )

  ! maximum no of link per vertex can be 'nbond'
  do i = 1, nvert
    if (allocated(vertex(i)%link_id)) deallocate(vertex(i)%link_id)
    allocate (vertex(i)%link_id(nbond), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )
    ! init
    vertex(i)%ndirlink = 0

    !neighbors of each vertex (allocate with max value initially)
    if (allocated(vertex(i)%nb_list1)) deallocate(vertex(i)%nb_list1)
    allocate (vertex(i)%nb_list1(2*nbond_total), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )

    if (allocated(vertex(i)%nb_list2)) deallocate(vertex(i)%nb_list2)
    allocate (vertex(i)%nb_list2(2*nbond_total), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )

    vertex(i)%n1neighb = 0
    vertex(i)%n2neighb = 0
  end do

  bravindex = [1,1,1]
  unitcell = this%glattice%unitcell
  n = 0
  cells_again: do i = 1, this%ncell
    bonds: do m = 1, nbond_total
      ! source vertex
      cellpos = bravindex + unitcell%bond(m)%source%offset
      cellpos = this%boundary_wrap(cellpos, s_phase)
      src_id = this%get_vertex_id( unitcell%bond(m)%source%id, cellpos )
      if (src_id == 0) cycle !link skuttled due to open boundary

      ! target vertex
      cellpos = bravindex + unitcell%bond(m)%target%offset
      offset = cellpos - [1,1,1] !before boundary wrap, save this cell offset
      cellpos = this%boundary_wrap(cellpos, t_phase)
      tgt_id = this%get_vertex_id( unitcell%bond(m)%target%id, cellpos )
      if (tgt_id == 0) cycle !-do-

      ! link exists, store it
      if (m <= nbond) then
        n = n + 1
        edge(n)%id = n
        edge(n)%type = unitcell%bond(m)%type
        edge(n)%src_id = src_id
        edge(n)%tgt_id = tgt_id
        edge(n)%class = edge_class(m)
        ! If the egde is wrapped across boundary odd number of times,
        ! then phase = -1
        edge(n)%phase = s_phase*t_phase
        ! it's a directed link to the 'source vertex'
        vertex(src_id)%ndirlink = vertex(src_id)%ndirlink + 1
        vertex(src_id)%link_id( vertex(src_id)%ndirlink ) = n

        ! edge vector
        id = unitcell%bond(m)%target%id
        coord = offset(1)*this%glattice%basis_vector%a1 &
              + offset(2)*this%glattice%basis_vector%a2 &
              + offset(3)*this%glattice%basis_vector%a3
        coord = coord + unitcell%atom(id)%coord 
        edge(n)%vector = coord - vertex(src_id)%coord 
      endif

      ! neighbour list
      if (unitcell%bond(m)%neighb == 1) then
        vertex(src_id)%n1neighb = vertex(src_id)%n1neighb + 1
	      j = vertex(src_id)%n1neighb
	      vertex(src_id)%nb_list1(j) = tgt_id

        vertex(tgt_id)%n1neighb = vertex(tgt_id)%n1neighb + 1
	      j = vertex(tgt_id)%n1neighb
	      vertex(tgt_id)%nb_list1(j) = src_id

      else if (unitcell%bond(m)%neighb == 2) then
        vertex(src_id)%n2neighb = vertex(src_id)%n2neighb + 1
	      j = vertex(src_id)%n2neighb
	      vertex(src_id)%nb_list2(j) = tgt_id

        vertex(tgt_id)%n2neighb = vertex(tgt_id)%n2neighb + 1
	      j = vertex(tgt_id)%n2neighb
	      vertex(tgt_id)%nb_list2(j) = src_id
      else
        call errmsg(alert, "graph::construct_graph", "found neighbor range > 2, ignored")
      endif
    end do bonds

    ! next bravais index
    bravindex = this%get_next_bravindex(bravindex)
  end do cells_again
  ! actual number of edge
  nedge = n

  ! neighbor list
  !call this%construct_nnlist()

  ! graph
  if (allocated(this%vertex)) deallocate(this%vertex)
  allocate (this%vertex(nvert), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )
  do i = 1, nvert
    if (allocated(this%vertex(i)%link_id)) deallocate(this%vertex(i)%link_id)
    allocate (this%vertex(i)%link_id( vertex(i)%ndirlink ), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )

    if (allocated(this%vertex(i)%nb_list1)) deallocate(this%vertex(i)%nb_list1)
    allocate (this%vertex(i)%nb_list1( vertex(i)%n1neighb ), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )

    if (allocated(this%vertex(i)%nb_list2)) deallocate(this%vertex(i)%nb_list2)
    allocate (this%vertex(i)%nb_list2( vertex(i)%n2neighb ), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_graph' )
  end do
  if (allocated(this%edge)) deallocate(this%edge)
  allocate (this%edge(nedge), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_graph' )
  this%nvertex = nvert
  this%nedge = nedge
  forall (i = 1:nvert) this%vertex(i) = vertex(i)
  forall (i = 1:nedge) this%edge(i) = edge(i)

  !check
  !if ( .true. ) then
  if ( .false. ) then
    call this%print_unitgraph()
    call this%print_graph()
    stop
  end if

  !check neighbour list
  if ( .false. ) then
    do i = 1, this%nvertex
      write(*,'(a,i5,1x,a)') "vertex = ", i, "NN:"
      do j = 1, this%vertex(i)%n1neighb
        write(*,'(2x,i3)') this%vertex(i)%nb_list1(j)
      enddo
      write(*,'(a,i5,1x,a)') "vertex = ", i, "NNN:"
      do j = 1, this%vertex(i)%n1neighb
        write(*,'(2x,i3)') this%vertex(i)%nb_list2(j)
      enddo
      write(*,*)
    enddo
    stop
  endif

  return
end subroutine construct_sub
!---------------------------------------------------------------------
!
!subroutine construct_nnlist_sub(this) 
  ! Purpose: Construct nearest neighbor list
  !
  ! Arguments:
  ! 1.
!  implicit none

  ! calling parameters
!  class (latticegraph_t) :: this

  ! local variables
!  return
!end subroutine construct_nnlist_sub
!---------------------------------------------------------------------
!
subroutine get_class_siteid_map_sub(this, map) 
  ! Purpose: Vertex class id to site id mapping
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(out) :: map(:)

  ! local variables
  integer(int32) :: i, cid

  if (size(map) /= this%nvclass) then
    call errmsg(error, "graph::get_class_siteid_map", "the map array size &
                should be equal to 'nvclass'")
  end if
  do i = 1, this%glattice%unitcell%natom
    cid = this%vertex(i)%class
    map(cid) = this%vertex(i)%atomid
  end do

  return
end subroutine get_class_siteid_map_sub
!---------------------------------------------------------------------
!
subroutine symmetry_reduce_sub(this, i, j, ri, rj) 
  ! Purpose: Given a bond specifiled by its connecting vertices i, j, the 
  ! subroutine outputs a bond which related to the input bond by lattice
  ! translational, rototational and reflection symmetry, but with the lowest
  ! values of the two connecting vertex id-s.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: i
  integer(int32), intent(in) :: j
  integer(int32), intent(out) :: ri
  integer(int32), intent(out) :: rj

  !local variables
  integer(int32) :: si, sj
  integer(int32) :: n, offset(3), intvec(3)
  logical :: reduced

  if (i<1 .or. i>this%nvertex) return
  if (j<1 .or. j>this%nvertex) return
  if (this%vertex(i)%id < this%vertex(j)%id) then
    si = this%vertex(i)%id; sj = this%vertex(j)%id;
  else
    si = this%vertex(j)%id; sj = this%vertex(i)%id;
  endif

  !reduce by translation
  ri = this%vertex(si)%uid
  offset = this%vertex(ri)%cellpos - this%vertex(si)%cellpos
  rj = this%translate_vertex(sj, offset)

  !write(*,*) i, j
  !(further 'reduction' done if vertices have same uid)
  if (this%vertex(ri)%uid == this%vertex(rj)%uid) then
    !-------reduce by rotation-----
    !(only  (0,x) to (x,0) of vertices with same uid)
    reduced = .false.
    intvec = this%vertex(rj)%cellpos - this%vertex(ri)%cellpos
    if (intvec(2)==0 .and. intvec(3)/=0) then
      intvec(1:3) = matmul(this%yrotate(1:3,1:3), intvec(1:3)) 
      reduced = .true.
    endif
    if (intvec(1)==0 .and. intvec(2)/=0) then
      intvec(1:3) = matmul(this%zrotate(1:3,1:3), intvec(1:3)) 
      reduced = .true.
    endif
    if (reduced) then
      offset(1:3) = intvec(1:3) + [1,1,1]
      rj = this%get_vertex_id(this%vertex(rj)%uid, offset)
    endif

    !-------reduce by reflection-----
    reduced = .false.
    intvec = this%vertex(rj)%cellpos - this%vertex(ri)%cellpos
    if (intvec(1) > this%xreflect%at) then
      n = intvec(1) - this%xreflect%at
      intvec(1) = this%xreflect%at - n - this%xreflect%offset
      reduced = .true.
    endif
    if (intvec(2) > this%yreflect%at) then
      n = intvec(2) - this%yreflect%at
      intvec(2) = this%yreflect%at - n - this%yreflect%offset
      reduced = .true.
    endif
    if (intvec(3) > this%zreflect%at) then
      n = intvec(3) - this%zreflect%at
      intvec(3) = this%zreflect%at - n - this%zreflect%offset
      reduced = .true.
    endif
    if (reduced) then
      offset(1:3) = intvec(1:3) + [1,1,1]
      rj = this%get_vertex_id(this%vertex(rj)%uid, offset)
    endif
  endif
  !write(*,*) ri, rj
  !read(*,*)

  return
end subroutine symmetry_reduce_sub
!---------------------------------------------------------------------
!
function translate_vertex_fun(this, v, intvec) result(translated_v)
  ! Purpose: Tranlates input vertex by a bravais index vector according 
  ! to the underlying lattice translational symmetry. Returns zero if 
  ! the translation takes the vertex out of an 'open' boundary.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: v
  integer(int32), intent(in) :: intvec(3)
  ! result
  integer(int32) :: translated_v

  ! local variables
  integer(int32) :: cellpos(3)

  translated_v = 0
  if (v<1 .or. v>this%nvertex) return
  cellpos = this%vertex(v)%cellpos + intvec
  cellpos = this%boundary_wrap( cellpos )
  translated_v = this%get_vertex_id( this%vertex(v)%uid, cellpos )

  return
end function translate_vertex_fun
!---------------------------------------------------------------------
!
subroutine check_validity(this, lattice, info) 
  ! Purpose: Suggestive.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t) :: this
  type (lattice_t), intent(in) :: lattice
  integer(int32), intent(out) :: info

  ! local variables
  integer(int32) :: i, j, id

  ! unitcell atom 'id'-s must be all distinct 
  do i = 1, lattice%unitcell%natom
    id = lattice%unitcell%atom(i)%id
    do j = 1, lattice%unitcell%natom
      if (i == j) cycle
	if (id == lattice%unitcell%atom(j)%id) then
	  call errmsg(error, 'graph::constructor_sub', "distinct atoms can not have identical 'id'" )
	  info = 1; return
      end if
    end do
  end do

  ! if 'tilted pbc'
  if (trim(lattice%dim3%periodicity) == "tiltperiodic") then
    call errmsg(error, 'graph::check_validity', "tilted PBC along 3rd &
         dimension is not allowed currently" )
    info = 1; return
  end if

  if (trim(lattice%dim1%periodicity) == "tiltperiodic") then
    if (trim(lattice%dim2%periodicity) /= "tiltperiodic") then
	call errmsg(error, 'graph::check_validity', "if dimension 1 has tilted PBC, &
           dimension 2 must have it" )
      info = 1; return
    end if
    if (mod(lattice%dim1%size, 2)==0 .or. mod(lattice%dim2%size, 2)==0) then
	call errmsg(error, 'graph::check_validity', "In case tilted PBC, &
           lattice size should be odd" )
      info = 1; return
    end if
    !if (mod(lattice%dim1%size)<3 .or. mod(lattice%dim2%size)<3) then
!	call errmsg(error, 'graph::check_validity', "In case tilted PBC, &
!             lattice dimension should be at least 3")
!        info = 1; return
!      end if
    this%tilted_lattice = .true.
  end if
  
  info = 0
  return
end subroutine check_validity
!---------------------------------------------------------------------
!
subroutine construct_glattice(this, inlattice) 
  ! Purpose: Construct the extended lattice for the graph. 
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t) :: this
  type (lattice_t), intent(in) :: inlattice

  ! local variables
  integer(int32) :: i, m, n, nt, id, src_id, tgt_id
  integer(int32) :: natom, nbond, nbond_total
  integer(int32), allocatable :: iwork(:)
  integer(int32) :: bravindex(3), cellpos(3), offset(3), src_offset(3), tgt_offset(3)
  real(dp) :: x, y, bvec(3), xvec(3), a(3,3)
  type (atom_t), allocatable :: atom(:)
  type (bond_t), allocatable :: bond(:)
  type (basisvector_t) :: basis_vector
  type (unitcell_t) :: unitcell
  logical :: tilted_lattice

  ! Even if there is "tilted bc", switch it off temporarily
  tilted_lattice = this%tilted_lattice
  this%tilted_lattice = .false.

  ! the given unitcell
  unitcell = inlattice%unitcell

  ! relabel atom 'id'-s in continuous order (if already not so), update bond connectivity info
  do m = 1, unitcell%natom
    unitcell%atom(m)%id = m
    id = inlattice%unitcell%atom(m)%id
    do n = 1, unitcell%nbond_total
      if (id == inlattice%unitcell%bond(n)%source%id) unitcell%bond(n)%source%id = m
      if (id == inlattice%unitcell%bond(n)%target%id) unitcell%bond(n)%target%id = m
    end do
  end do

  ! serialize atom 'atomid'-s, if already not so
  ! order
  if (allocated(iwork)) deallocate(iwork)
  allocate( iwork(inlattice%unitcell%natom), stat=estat )
  call catch( estat, alloc_t, 'graph::construct_glattice' )
  forall (i = 1:unitcell%natom) iwork(i) =  unitcell%atom(i)%atomid
  call this%serialize( iwork )
  forall (i = 1:unitcell%natom) unitcell%atom(i)%atomid = iwork(i) 
  ! update background info
  this%atomid_range = unitcell%atom(iwork(unitcell%natom))%atomid

  ! Construct 'extended' unitcell (= one_bravais_lattice_point + basis).
  ! 'basis' atoms are constructed from given unitcell atoms & boundary conditions.
  !
  ! If the resulting lattice is 1 dimensional, rotate the lattice so that the 
  ! basis vector (stored in 'bvec') is along x.
  !
  ! First define the lattice to have only the 'extended' unitcell.
  this%glattice = inlattice
  this%glattice%unitcell = unitcell
  this%glattice%ldim = 3

  if (trim(inlattice%dim1%boundary) == "periodic") then
    this%glattice%dim1%size = 1
    bvec = inlattice%basis_vector%a1
  else 
    this%glattice%ldim = this%glattice%ldim - 1
  end if

  if (trim(inlattice%dim2%boundary) == "periodic") then
    bvec = inlattice%basis_vector%a2
    this%glattice%dim2%size = 1
  else 
    this%glattice%ldim = this%glattice%ldim - 1
  end if

  if (trim(inlattice%dim3%boundary) == "periodic") then
    bvec = inlattice%basis_vector%a3
    this%glattice%dim3%size = 1
  else 
    this%glattice%ldim = this%glattice%ldim - 1
  end if

  ! if 1 dimensional, rotate the lattice to align 'bvec' along x
  basis_vector = inlattice%basis_vector
  if (this%glattice%ldim == 1) then
    x = sqrt(dot_product(bvec,bvec))
    xvec = [1.0_dp, 0.0_dp, 0.0_dp]
    a = this%get_rotmatrix(bvec, xvec)
    ! redefine basis vectors
    if (trim(this%glattice%dim1%boundary) == "periodic") then
      basis_vector%a1 = [x, 0.0_dp, 0.0_dp]
      basis_vector%a2 = matmul(a, basis_vector%a2)
      basis_vector%a3 = matmul(a, basis_vector%a3)
    end if
    if (trim(this%glattice%dim2%boundary) == "periodic") then
      basis_vector%a1 = matmul(a, basis_vector%a1)
      basis_vector%a2 = [x, 0.0_dp, 0.0_dp]
      basis_vector%a3 = matmul(a, basis_vector%a3)
    end if
    if (trim(this%glattice%dim3%boundary) == "periodic") then
      basis_vector%a1 = matmul(a, basis_vector%a1)
      basis_vector%a2 = matmul(a, basis_vector%a2)
      basis_vector%a3 = [x, 0.0_dp, 0.0_dp]
    end if
    ! rotate atom coordinates
    do m = 1, unitcell%natom  
      xvec = this%glattice%unitcell%atom(m)%coord
      this%glattice%unitcell%atom(m)%coord = matmul(a, xvec)
    end do
    ! update new basis vector
    this%glattice%basis_vector = basis_vector
  end if

  ! update other private info
  this%npcell = (this%glattice%dim1%size * this%glattice%dim2%size)
  this%ncell = (this%npcell * this%glattice%dim3%size)

  ! now construct the extended unitcell
  natom = this%ncell*this%glattice%unitcell%natom
  nbond = this%ncell*this%glattice%unitcell%nbond
  nbond_total = this%ncell*this%glattice%unitcell%nbond_total

  ! allocate
  if (allocated(atom)) deallocate(atom)
  allocate (atom(natom), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_glattice' )

  if (allocated(bond)) deallocate(bond)
  allocate (bond(nbond_total), stat=estat)
  call catch( estat, alloc_t, 'graph::construct_glattice' )

  ! collect the atoms
  ! set all atom 'offset' to be 1
  forall (i = 1:unitcell%natom) this%glattice%unitcell%atom(i)%offset = [1,1,1]
  bravindex = [1,1,1]
  cells: do i = 1, this%ncell
    ! translate cell
    offset = bravindex - [1,1,1]
    unitcell = this%translate_cell( this%glattice%unitcell, offset )
    ! add atoms
    atoms: do m = 1, unitcell%natom
      n = this%get_vertex_id( unitcell%atom(m)%id, unitcell%atom(m)%offset )
      atom(n) = unitcell%atom(m)
	atom(n)%id = n
	atom(n)%offset = [1,1,1] !all atoms are in the extended unitcell
    end do atoms
    ! next bravais index
    bravindex = this%get_next_bravindex(bravindex)
  end do cells

  ! connect the bonds
  bravindex = [1,1,1]
  unitcell = this%glattice%unitcell
  n = 0; nt = 0
  cells_again: do i = 1, this%ncell
    bonds: do m = 1, unitcell%nbond_total
      ! source atom
      src_offset = unitcell%bond(m)%source%offset
      cellpos = bravindex + src_offset
      if (cellpos(1)>=1 .and. cellpos(1)<=this%glattice%dim1%size) src_offset(1) = 0
      if (cellpos(2)>=1 .and. cellpos(2)<=this%glattice%dim2%size) src_offset(2) = 0
      if (cellpos(3)>=1 .and. cellpos(3)<=this%glattice%dim3%size) src_offset(3) = 0
      cellpos = this%boundary_wrap( cellpos )
      src_id = this%get_vertex_id( unitcell%bond(m)%source%id, cellpos )
      if (src_id == 0) cycle !bond skuttled due to open boundary

      ! target atom
      tgt_offset = unitcell%bond(m)%target%offset
      cellpos = bravindex + tgt_offset
      if (cellpos(1)>=1 .and. cellpos(1)<=this%glattice%dim1%size) tgt_offset(1) = 0
      if (cellpos(2)>=1 .and. cellpos(2)<=this%glattice%dim2%size) tgt_offset(2) = 0
      if (cellpos(3)>=1 .and. cellpos(3)<=this%glattice%dim3%size) tgt_offset(3) = 0
      cellpos = this%boundary_wrap( cellpos )
      tgt_id = this%get_vertex_id( unitcell%bond(m)%target%id, cellpos )
      if (tgt_id == 0) cycle !-do-

      ! bond exists, add it
      if (m <= unitcell%nbond) n = n + 1
      nt = nt + 1
      bond(nt) = unitcell%bond(m)
      bond(nt)%source = atom(src_id)
      bond(nt)%target = atom(tgt_id)
      bond(nt)%source%offset = src_offset
      bond(nt)%target%offset = tgt_offset
    end do bonds

    ! next bravais index
    bravindex = this%get_next_bravindex(bravindex)
  end do cells_again
  ! actual number of bonds
  nbond = n
  nbond_total = nt

  ! finally the desired lattice, with extended unitcell
  ! basis vectors
  this%glattice%basis_vector = basis_vector
  ! name
  if (this%glattice%ldim == 1) this%glattice%name = "1D Lattice"
  if (this%glattice%ldim == 2) this%glattice%name = "2D Lattice"
  if (this%glattice%ldim == 3) this%glattice%name = "3D Lattice"
  ! extent
  this%glattice%dim1 = inlattice%dim1
  this%glattice%dim2 = inlattice%dim2
  this%glattice%dim3 = inlattice%dim3
  if (trim(this%glattice%dim1%boundary) /=  "periodic") then
    this%glattice%dim1%size = 1
    this%glattice%basis_vector%a1 = [0.0_dp, 0.0_dp, 0.0_dp]
  end if
  if (trim(this%glattice%dim2%boundary) /=  "periodic") then
    this%glattice%dim2%size = 1
    this%glattice%basis_vector%a2 = [0.0_dp, 0.0_dp, 0.0_dp]
  end if
  if (trim(this%glattice%dim3%boundary) /=  "periodic") then
    this%glattice%dim3%size = 1
    this%glattice%basis_vector%a3 = [0.0_dp, 0.0_dp, 0.0_dp]
  end if

  ! unitcell
  deallocate( this%glattice%unitcell%atom )
  deallocate( this%glattice%unitcell%bond )

  if (allocated(this%glattice%unitcell%atom)) deallocate(this%glattice%unitcell%atom)
  allocate( this%glattice%unitcell%atom(natom), stat=estat )
  call catch( estat, alloc_t, 'graph::construct_glattice' )

  if (allocated(this%glattice%unitcell%bond)) deallocate(this%glattice%unitcell%bond)
  allocate( this%glattice%unitcell%bond(nbond_total), stat=estat )
  call catch( estat, alloc_t, 'graph::construct_glattice' )

  this%glattice%unitcell%natom = natom
  this%glattice%unitcell%nbond = nbond
  this%glattice%unitcell%nbond_total = nbond_total
  forall (m = 1:natom) this%glattice%unitcell%atom(m) = atom(m)
  forall (m = 1:nbond_total) this%glattice%unitcell%bond(m) = bond(m)

  ! background info
  this%npcell = (this%glattice%dim1%size * this%glattice%dim2%size)
  this%ncell = (this%npcell * this%glattice%dim3%size)
  if (allocated(iwork)) deallocate(iwork)
  allocate( iwork(natom), stat=estat )
  call catch( estat, alloc_t, 'graph::construct_glattice' )
  forall (m = 1:natom) iwork(m) = atom(m)%atomid
  this%atomid_range = maxval(iwork)

  !max site 'type' value
  this%max_vertex_type = 0
  do m = 1, this%glattice%unitcell%natom
    if (this%glattice%unitcell%atom(m)%type > this%max_vertex_type) then
      this%max_vertex_type = this%glattice%unitcell%atom(m)%type
    endif
  enddo

  !max bond 'type' value
  this%max_edge_type = 0
  do m = 1, this%glattice%unitcell%nbond
    if (this%glattice%unitcell%bond(m)%type > this%max_edge_type) then
      this%max_edge_type = this%glattice%unitcell%bond(m)%type
    endif
  enddo


  ! In case of "tilted boundary conditions", we need customized 
  ! attributes.
  this%tilted_lattice = tilted_lattice
  if (this%tilted_lattice) then
    this%tilted_shift = 1
    this%tilted_size1 = (this%glattice%dim1%size+1)/2
    this%tilted_size2 = this%glattice%dim2%size+1
    this%npcell = this%npcell + 1
    this%ncell = this%ncell + this%glattice%dim3%size

    if (allocated(this%tilted_bottom)) deallocate(this%tilted_bottom)
    allocate (this%tilted_bottom(this%glattice%dim1%size), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_glattice' )

    if (allocated(this%tilted_top)) deallocate(this%tilted_top)
    allocate (this%tilted_top(this%glattice%dim1%size), stat=estat)
    call catch( estat, alloc_t, 'graph::construct_glattice' )

    m = (this%glattice%dim1%size+1)/2
    do i = 1, m
      this%tilted_bottom(i) = 1
      this%tilted_top(i) = this%glattice%dim2%size
    end do
    do i = m+1, this%glattice%dim1%size
      this%tilted_bottom(i) = 2
      this%tilted_top(i) = this%glattice%dim2%size+1
    end do
    this%tilted_top(m) = this%glattice%dim2%size+1

  else
    this%tilted_shift = 0
  end if


  !------Rotational Symmetry------
  !About x-axis
  !----default----
  this%xrotate(1,1:3) = [1,0,0]
  this%xrotate(2,1:3) = [0,1,0]
  this%xrotate(3,1:3) = [0,0,1]
  if (.not. this%tilted_lattice) then !--do not know what, in this case
    if (trim(this%glattice%dim2%periodicity)=="periodic" .and. this%glattice%dim3%periodicity=="periodic") then
      x = dot_product(this%glattice%basis_vector%a2(1:3), this%glattice%basis_vector%a2(1:3))
      y = dot_product(this%glattice%basis_vector%a3(1:3), this%glattice%basis_vector%a3(1:3))
      if (abs(x-y)<1.0E-10_dp) then
        this%xrotate(2,1:3) = [0, 0,+1]
        this%xrotate(3,1:3) = [0,-1, 0]
      endif
    endif
  endif

  !About y-axis
  !----default----
  this%yrotate(1,1:3) = [1,0,0]
  this%yrotate(2,1:3) = [0,1,0]
  this%yrotate(3,1:3) = [0,0,1]
  if (.not. this%tilted_lattice) then
    if (trim(this%glattice%dim1%periodicity)=="periodic" .and. this%glattice%dim3%periodicity=="periodic") then
      x = dot_product(this%glattice%basis_vector%a1(1:3), this%glattice%basis_vector%a1(1:3))
      y = dot_product(this%glattice%basis_vector%a3(1:3), this%glattice%basis_vector%a3(1:3))
      if (abs(x-y)<1.0E-10_dp) then
        this%yrotate(1,1:3) = [ 0, 0,+1]
        this%yrotate(3,1:3) = [-1, 0, 0]
      endif
    endif
  endif

  !About z-axis
  !----default----
  this%zrotate(1,1:3) = [1,0,0]
  this%zrotate(2,1:3) = [0,1,0]
  this%zrotate(3,1:3) = [0,0,1]
  if (.not. this%tilted_lattice) then 
    if (trim(this%glattice%dim1%periodicity)=="periodic" .and. this%glattice%dim2%periodicity=="periodic") then
      x = dot_product(this%glattice%basis_vector%a1(1:3), this%glattice%basis_vector%a1(1:3))
      y = dot_product(this%glattice%basis_vector%a2(1:3), this%glattice%basis_vector%a2(1:3))
      if (abs(x-y)<1.0E-10_dp) then
        this%zrotate(1,1:3) = [0, +1, 0]
        this%zrotate(2,1:3) = [-1, 0, 0]
      endif
    endif
  endif
  !write(*,*) this%zrotate(1,1:3)
  !write(*,*) this%zrotate(2,1:3)
  !write(*,*) this%zrotate(3,1:3)
  !stop

  !Reflection operation
  !About a plane _|_ to 'cell x-axis'
  !---default---
  this%xreflect%at = 0
  this%xreflect%offset = 0
  if (.not. this%tilted_lattice) then 
    if (trim(this%glattice%dim1%boundary)=="periodic") then
      if (mod(this%glattice%dim1%size,2) == 0) then
        this%xreflect%at = this%glattice%dim1%size/2
	this%xreflect%offset = 0
      else
        this%xreflect%at = (this%glattice%dim1%size+1)/2
	this%xreflect%offset = 1
      endif
    endif
  endif

  !About a plane _|_ to 'cell y-axis'
  this%yreflect%at = 0
  this%yreflect%offset = 0
  if (.not. this%tilted_lattice) then 
    if (trim(this%glattice%dim2%boundary)=="periodic") then
      if (mod(this%glattice%dim2%size,2) == 0) then
        this%yreflect%at = this%glattice%dim2%size/2
	this%yreflect%offset = 0
      else
        this%yreflect%at = (this%glattice%dim2%size+1)/2
	this%yreflect%offset = 1
      endif
    endif
  endif

  !About a plane _|_ to 'cell z-axis'
  this%zreflect%at = 0
  this%zreflect%offset = 0
  if (.not. this%tilted_lattice) then 
    if (trim(this%glattice%dim3%boundary)=="periodic") then
      if (mod(this%glattice%dim3%size,2) == 0) then
        this%zreflect%at = this%glattice%dim3%size/2
	this%zreflect%offset = 0
      else
        this%zreflect%at = (this%glattice%dim3%size+1)/2
	this%zreflect%offset = 1
      endif
    endif
  endif

  !print info
  !check: if (.true.) then
  check: if (.false.) then
    write (*,'(/,1x,2a)') "graph::construct_glattice:", "'glattice' info"
    write (*,'(1x,a)')  "---------------------------------------------------"
    write (*,'(1x,2a)') 'name = ', trim(this%glattice%name)
    write (*,'(1x,a,i3)') 'lattice id = ', this%glattice%id
    write (*,'(1x,a,i3)') 'dimension = ', this%glattice%ldim
    write (*,'(1x,a,3f6.3)') 'basisvec_a1 = ', this%glattice%basis_vector%a1
    write (*,'(1x,a,3f6.3)') 'basisvec_a2 = ', this%glattice%basis_vector%a2
    write (*,'(1x,a,3f6.3)') 'basvisec_a3 = ', this%glattice%basis_vector%a3
    write (*,'(1x,a,i4,a,i4,a,i4)') 'lattice size = ', this%glattice%dim1%size, ' x', &
          this%glattice%dim2%size, ' x', this%glattice%dim3%size
    write (*,'(1x,6a)') 'boundary = ', trim(this%glattice%dim1%boundary), '-', &
          trim(this%glattice%dim2%boundary), '-', trim(this%glattice%dim3%boundary)
    write (*,'(1x,a)') 'atoms:'
    do m = 1, natom
      write (*,'(1x,a,i3,a,a,i3,a,i3,a,i3,a,3f6.2,a,3i3,a)') "atom =", m, ": ", &
            "id =", this%glattice%unitcell%atom(m)%id, &
            ", type =", this%glattice%unitcell%atom(m)%type, &
	      ", atomid =", this%glattice%unitcell%atom(m)%atomid, &
	      ", coord = (", this%glattice%unitcell%atom(m)%coord, &
	      " ), offset = (", this%glattice%unitcell%atom(m)%offset, " )"
    end do
    write (*,'(1x,a)') 'bonds:'
    do m = 1, nbond
      write (*,'(1x,a,i3,a,a,i3,a,i3,a,3i3,a,i3,a,3i3,a)') "bond =", m, ": ", &
            "type =", this%glattice%unitcell%bond(m)%type, &
            ", source_id =", this%glattice%unitcell%bond(m)%source%id, &
	      ", offset = (", this%glattice%unitcell%bond(m)%source%offset, &
	      " ), target_id =", this%glattice%unitcell%bond(m)%target%id, &
	      ", offset = (", this%glattice%unitcell%bond(m)%target%offset, " )"
    end do
  end if check

  return
end subroutine construct_glattice
!---------------------------------------------------------------------
!
function get_next_bravindex_fun(this, current_index) result(next_index)
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: current_index(3)
  ! result
  integer(int32) :: next_index(3)

  if (this%tilted_lattice) then
    next_index = this%bravindex_tilted_lattice(current_index)
    return
  end if

  if (minval(current_index)<1 &
        .or. current_index(1)>this%glattice%dim1%size &
	  .or. current_index(2)>this%glattice%dim2%size &
	  .or. current_index(3)>this%glattice%dim3%size) then
    call errmsg(alert, "graph::get_next_bravindex", "current index out of range")
    return
  end if
  next_index(1) = current_index(1)+1
  next_index(2) = current_index(2)
  next_index(3) = current_index(3)
  if (next_index(1) > this%glattice%dim1%size) then
    next_index(1) = 1
    next_index(2) = current_index(2)+1
    if (next_index(2) > this%glattice%dim2%size) then
      next_index(2) = 1
      next_index(3) = current_index(3)+1
      if (next_index(3) > this%glattice%dim3%size) then
	  next_index = [0,0,0]
      end if
    end if
  end if

  return
end function get_next_bravindex_fun
!---------------------------------------------------------------------
!
function bravindex_tilted_lattice(this, current_index) result(next_index)
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: current_index(3)
  ! result
  integer(int32) :: next_index(3)

  ! tilted_size1 = (lattice%dim1%size+1)/2
  ! tilted_size2 = (lattice%dim2%size+1)
  if (minval(current_index)<1 &
        .or. current_index(1)>this%glattice%dim1%size &
	  .or. current_index(2)>this%tilted_size2 &
	  .or. current_index(3)>this%glattice%dim3%size) then
    call errmsg(alert, "graph::bravindex_tilted_lattice", "current index out of range")
    return
  end if

  next_index(1) = current_index(1)+1
  next_index(2) = current_index(2)
  next_index(3) = current_index(3)
  if (current_index(2) == 1) then
    if (next_index(1) > this%tilted_size1) then
      next_index(1) = 1
      next_index(2) = current_index(2)+1
      if (next_index(2) > this%tilted_size2) then
        next_index(2) = 1
        next_index(3) = current_index(3)+1
        if (next_index(3) > this%glattice%dim3%size) then
	    next_index = [0,0,0]
        end if
      end if
    end if
  else if (current_index(2) == this%glattice%dim2%size) then
    if (next_index(1) > this%glattice%dim1%size) then
      next_index(1) = this%tilted_size1
      next_index(2) = this%tilted_size2
      if (next_index(2) > this%tilted_size2) then
        next_index(1) = 1
        next_index(2) = 1
        next_index(3) = current_index(3)+1
        if (next_index(3) > this%glattice%dim3%size) then
	    next_index = [0,0,0]
        end if
      end if
    end if
  else
    if (next_index(1) > this%glattice%dim1%size) then
      next_index(1) = 1
      next_index(2) = current_index(2)+1
      if (next_index(2) > this%tilted_size2) then
        next_index(1) = 1
        next_index(2) = 1
        next_index(3) = current_index(3)+1
        if (next_index(3) > this%glattice%dim3%size) then
	    next_index = [0,0,0]
        end if
      end if
    end if
  end if

  return
end function bravindex_tilted_lattice
!---------------------------------------------------------------------
!
function translate_cell(this, incell, bravindex) result(newcell)
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  type (unitcell_t), intent(in) :: incell
  integer(int32), intent(in) :: bravindex(3)
  ! result
  type (unitcell_t) :: newcell

  ! list of local variables
  integer(int32) :: m, id_offset, cl_offset
  integer(int32) :: cellpos(3)
  real(dp) :: R(3)

  ! copy
  newcell = incell

  ! modify
  cellpos = incell%atom(1)%offset + bravindex 
  !cl_offset = this%get_cell_number(cellpos) - 1
  cl_offset = bravindex(1) + this%glattice%dim1%size * bravindex(2) &
           + this%npcell * bravindex(3)
  ! Taking care of tilted lattice
  if (this%tilted_lattice) then
    if (incell%atom(1)%offset(2)<=1) then
      if (cellpos(2)>1 .and. cellpos(2)<this%tilted_size2) then
        cl_offset = cl_offset - (this%tilted_size1-1)
	    else if (cellpos(2)>=this%tilted_size2) then
        cl_offset = cl_offset - 2*(this%tilted_size1-1)
      end if
    else 
      if (cellpos(2)>=this%tilted_size2) then
        cl_offset = cl_offset - (this%tilted_size1-1)
	    end if
    end if
  end if
  id_offset = cl_offset * this%atomid_range

  R = (bravindex(1))*this%glattice%basis_vector%a1 &
    + (bravindex(2))*this%glattice%basis_vector%a2 &
    + (bravindex(3))*this%glattice%basis_vector%a3
  do m = 1, incell%natom
    newcell%atom(m)%atomid = incell%atom(m)%atomid +  id_offset
    newcell%atom(m)%offset = cellpos
    newcell%atom(m)%coord = incell%atom(m)%coord + R
    newcell%atom(m)%cellcoord = R
  end do

  return
end function translate_cell
!---------------------------------------------------------------------
!
function get_cell_number(this, cellpos) result(nth_cell)
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: cellpos(3)
  ! result
  integer(int32) :: nth_cell

  if (this%tilted_lattice) then 
    nth_cell = this%cell_number_fun_tiltbc(cellpos)
    return
  end if

  if (cellpos(1)<1 .or. cellpos(2)<1 .or. cellpos(3)<1) then
    nth_cell = 0; return
  end if
  nth_cell = cellpos(1) + this%glattice%dim1%size * (cellpos(2)-1) &
           + this%npcell * (cellpos(3) - 1)

  return
end function get_cell_number
!---------------------------------------------------------------------
!
function cell_number_fun_tiltbc(this, cellpos) result(nth_cell)
  ! Purpose: 
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: cellpos(3)
  ! result
  integer(int32) :: nth_cell

  if (cellpos(1)<1 .or. cellpos(2)<1 .or. cellpos(3)<1) then
    nth_cell = 0; return
  end if

  if (cellpos(2) == 1) then 
    nth_cell = cellpos(1) + this%glattice%dim1%size * (cellpos(2)-1) &
             + this%npcell * (cellpos(3)-1)
  else if (cellpos(2)>1 .and. cellpos(2)<this%tilted_size2) then
    nth_cell = cellpos(1) + this%glattice%dim1%size * (cellpos(2)-1) &
             + this%npcell * (cellpos(3)-1) - (this%tilted_size1-1)
  else
    nth_cell = cellpos(1) + this%glattice%dim1%size * (cellpos(2)-1) &
             + this%npcell * (cellpos(3)-1) - 2*(this%tilted_size1-1)
  end if

  return
end function cell_number_fun_tiltbc
!---------------------------------------------------------------------
!
function get_vertex_id(this, uid, cellpos) result(id)
  ! Purpose: Suggestive.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: uid
  integer(int32), intent(in) :: cellpos(3)
  ! result
  integer(int32) :: id, cellno

  cellno = this%get_cell_number(cellpos)
  if (cellno < 1) then
    id = 0; return
  end if
  id = (cellno-1)*this%glattice%unitcell%natom + uid

  return
end function get_vertex_id
!---------------------------------------------------------------------
!
function boundary_wrap(this, cellpos, phase) result(newpos)
  ! Purpose: Unitcell coordinates get folded if boundary condition is
  ! periodic.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: cellpos(3)
  integer(int32), intent(out), optional :: phase
  integer(int32) :: newpos(3)

  ! local variables
  integer(int32) :: n

  if (this%tilted_lattice) then
    if (present(phase)) then
      newpos = this%boundary_wrap_tiltbc(cellpos, phase)
    else
      newpos = this%boundary_wrap_tiltbc(cellpos)
    end if
    return
  end if

  newpos = cellpos
  ! to wrap, temporarlity number sites from 0 to N-1, take modulo N, then increament 1
  if (newpos(1)<1 .or. newpos(1)>this%glattice%dim1%size) then
    if (trim(this%glattice%dim1%boundary) == "periodic") then
      newpos(1) = newpos(1) - 1
      newpos(1) = modulo(newpos(1), this%glattice%dim1%size)
      newpos(1) = newpos(1) + 1
    else
      newpos(1) = 0
    end if
  end if

  if (newpos(2)<1 .or. newpos(2)>this%glattice%dim2%size) then
    if (trim(this%glattice%dim2%boundary) == "periodic") then
      newpos(2) = newpos(2) - 1
      newpos(2) = modulo(newpos(2), this%glattice%dim2%size)
      newpos(2) = newpos(2) + 1
    else
      newpos(2) = 0
    end if
  end if

  if (newpos(3)<1 .or. newpos(3)>this%glattice%dim3%size) then
    if (trim(this%glattice%dim3%boundary) == "periodic") then
      newpos(3) = newpos(3) - 1
      newpos(3) = modulo(newpos(3), this%glattice%dim3%size)
      newpos(3) = newpos(3) + 1
      if (present(phase)) phase = -1
    else
      newpos(3) = 0
    end if
  end if

  ! Determine phase: Phase is negative, if the edge gets wrapped across 
  ! 'antiperiodic' boundary odd number of times
  if_phase: if (present(phase)) then
    phase = 1
    if (trim(this%glattice%dim1%periodicity) == "antiperiodic") then
      if (cellpos(1)<1) then
	  n = abs(cellpos(1))/this%glattice%dim1%size + 1
	  if (mod(n,2) /= 0) phase = -phase
      else if (cellpos(1)>this%glattice%dim1%size) then
	  n = (cellpos(1)-1)/this%glattice%dim1%size
	  if (mod(n,2) /= 0) phase = -phase
	endif
    endif

    if (trim(this%glattice%dim2%periodicity) == "antiperiodic") then
      if (cellpos(2)<1) then
	  n = abs(cellpos(2))/this%glattice%dim2%size + 1
	  if (mod(n,2) /= 0) phase = -phase
      else if (cellpos(2)>this%glattice%dim2%size) then
	  n = (cellpos(2)-1)/this%glattice%dim2%size
	  if (mod(n,2) /= 0) phase = -phase
	endif
    endif

    if (trim(this%glattice%dim3%periodicity) == "antiperiodic") then
      if (cellpos(3)<1) then
	  n = abs(cellpos(3))/this%glattice%dim3%size + 1
	  if (mod(n,2) /= 0) phase = -phase
      else if (cellpos(3)>this%glattice%dim3%size) then
	  n = (cellpos(3)-1)/this%glattice%dim3%size
	  if (mod(n,2) /= 0) phase = -phase
	endif
    endif
  endif if_phase

  return
end function boundary_wrap
!---------------------------------------------------------------------
!
function boundary_wrap_tiltbc(this, cellpos, phase) result(newpos)
  ! Purpose: Unitcell coordinates get folded if boundary condition is
  ! periodic.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in) :: cellpos(3)
  integer(int32), intent(out), optional :: phase
  integer(int32) :: newpos(3)

  ! local variables
  integer(int32) :: n, n1, n2, n3

  !write(*,*) "cellpos=", cellpos
  n1 = cellpos(1)
  n2 = cellpos(2)
  n3 = cellpos(3)

  ! First round, wrap the cordinate-1
  if (n1 < 1) then
    n1 = n1 - 1
    n1 = modulo(n1, this%glattice%dim1%size)
    n1 = n1 + 1
    ! the tilted shift
    n2 = n2 + this%tilted_shift
  else if (n1>this%glattice%dim1%size) then
    ! do same as before except for negative 'tilted_shift'
    n1 = n1 - 1
    n1 = modulo(n1, this%glattice%dim1%size)
    n1 = n1 + 1
    ! the tilted shift
    n2 = n2 - this%tilted_shift
  end if

  ! 2nd round, wrap cordinate-2. This necessitates wrapping cordinate-2 as
  ! well as these two are interwinned.
  do while (n1<1 .or. n1>this%glattice%dim1%size .or. n2<this%tilted_bottom(n1) &
            .or. n2>this%tilted_top(n1))
    ! to wrap, temporarlity number sites from 0 to N-1, take modulo N, then increament 1
    if (n2 < this%tilted_bottom(n1)) then
	if (n1 <= (this%glattice%dim1%size+1)/2) then
        n2 = n2 - 1
        n2 = modulo(n2, this%glattice%dim2%size)
        n2 = n2 + 1
	else
        n2 = n2 - 2
        n2 = modulo(n2, this%glattice%dim2%size)
        n2 = n2 + 2
	end if
        ! the tilted shift
      n1 = n1 - this%tilted_shift
    else if (n2>this%tilted_top(n1)) then
      ! do same as before except for negative 'tilted_shift'
      n2 = n2 - 1
      n2 = modulo(n2, this%glattice%dim2%size)
      n2 = n2 + 1
      ! the tilted shift
      n1 = n1 + this%tilted_shift
    end if

    if (n1 < 1) then
      n1 = n1 - 1
      n1 = modulo(n1, this%glattice%dim1%size)
      n1 = n1 + 1
      ! the tilted shift
      n2 = n2 + this%tilted_shift
    else if (n1>this%glattice%dim1%size) then
      ! do same as before except for negative 'tilted_shift'
      n1 = n1 - 1
      n1 = modulo(n1, this%glattice%dim1%size)
      n1 = n1 + 1
      ! the tilted shift
      n2 = n2 - this%tilted_shift
    end if
  end do !while

  if (n3<1 .or. n3>this%glattice%dim3%size) then
    if (trim(this%glattice%dim3%boundary) == "periodic") then
      n3 = n3 - 1
      n3 = modulo(n3, this%glattice%dim3%size)
      n3 = n3 + 1
    else
      n3 = 0
    end if
  end if
  newpos = [n1,n2,n3]
  !write(*,*) "nespos=", newpos; read(*,*)

  ! Determine phase: Phase is negative, if the edge gets wrapped across 
  ! 'antiperiodic' boundary odd number of times
  if_phase: if (present(phase)) then
    phase = 1
    if (trim(this%glattice%dim3%periodicity) == "antiperiodic") then
      if (cellpos(3)<1) then
	  n = abs(cellpos(3))/this%glattice%dim3%size + 1
	  if (mod(n,2) /= 0) phase = -phase
      else if (cellpos(3)>this%glattice%dim3%size) then
	  n = (cellpos(3)-1)/this%glattice%dim3%size
	  if (mod(n,2) /= 0) phase = -phase
	endif
    endif
  endif if_phase

  return
end function boundary_wrap_tiltbc
!---------------------------------------------------------------------
!
function get_rotmatrix(r, rp) result(a)
  ! Purpose: Calculates rotation matrix which would rotate vector 'r' 
  ! to align it to common origin vector 'rp'. Rotation axis is 
  ! along 'rp x r'. For the formula used to calculate the matrix,
  ! see H. Goldstein, Classical Mechanics, Eq. 4-92 through 4-96.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  real(dp), intent(in) :: r(3)
  real(dp), intent(in) :: rp(3)
  real(dp) :: a(3,3)

  ! local variables
  real(dp) :: e0, e1, e2, e3, phi
  real(dp) :: n(3)

  ! angle between "r" and "r'"
  e1 = sqrt(dot_product(r, r))
  e2 = sqrt(dot_product(rp, rp))
  phi = acos(dot_product(r, rp)/(e1*e2))
  if (abs(phi) < dp_tol) then
    a(:,:)=0.0_dp
    a(1,1)=1.0_dp; a(2,2)=1.0_dp; a(3,3)=1.0_dp;
    return
  end if

  ! unit vector perpendicular to 'rp' and 'r'
  ! \hat{n} = |\vec{r'} x \vec{r}|
  n(1) = rp(2)*r(3) - rp(3)*r(2)
  n(2) = rp(3)*r(1) - rp(1)*r(3)
  n(3) = rp(1)*r(2) - rp(2)*r(1)
  n = n/sqrt(dot_product(n, n))

  ! the 'e' parameters
  phi = phi*0.50_dp
  e0 = cos(phi)
  e1 = n(1)*sin(phi)
  e2 = n(2)*sin(phi)
  e3 = n(3)*sin(phi)

  ! rotation matrix
  a(1,1) = e0*e0 + e1*e1 - e2*e2 - e3*e3
  a(1,2) = 2.0_dp*(e1*e2 + e0*e3)
  a(1,3) = 2.0_dp*(e1*e3 - e0*e2)

  a(2,1) = 2.0_dp*(e1*e2 - e0*e3)
  a(2,2) = e0*e0 - e1*e1 + e2*e2 - e3*e3
  a(2,3) = 2.0_dp*(e2*e3 + e0*e1)

  a(3,1) = 2.0_dp*(e1*e3 + e0*e2)
  a(3,2) = 2.0_dp*(e2*e3 - e0*e1)
  a(3,3) = e0*e0 - e1*e1 - e2*e2 + e3*e3

  return
end function get_rotmatrix
!---------------------------------------------------------------------
!
subroutine get_vertex_class(this, unitcell, class_id, nvclass) 
  ! Purpose: Atoms in a unitcell with a unique combination of 
  ! attributes- {atom_site, type} gets a unique 'class_id'.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  type (unitcell_t), intent(in) :: unitcell
  integer(int32), intent(out) :: class_id(:)
  integer(int32), intent(out) :: nvclass

  ! list of local variables
  integer(int32) :: i, j, k
  integer(int32) :: natom, min_typ, max_typ, min_sid, max_sid
  integer(int32) :: iwork1(unitcell%natom)
  integer(int32), allocatable :: iwork2(:,:)

  natom = unitcell%natom
  !determine ranges
  !site range
  forall (i = 1:natom) iwork1(i) = unitcell%atom(i)%atomid
  min_sid = minval(iwork1)
  max_sid = maxval(iwork1)
  !type range
  forall (i = 1:natom) iwork1(i) = unitcell%atom(i)%type
  min_typ = minval(iwork1)
  max_typ = maxval(iwork1)

  !workspace
  if (allocated(iwork2)) deallocate(iwork2)
  allocate( iwork2(min_sid:max_sid, min_typ:max_typ), stat=estat )
  call catch( estat, alloc_t, "graph::get_vertex_class" )

  !assign unique class_id for unique {atomid, type} pair
  k = 1
  do i = min_sid, max_sid
    do j = min_typ, max_typ
      iwork2(i,j) = k
	k = k + 1
    end do
  end do
  forall (i = 1:natom)
    class_id(i) = iwork2( unitcell%atom(i)%atomid, unitcell%atom(i)%type )
  end forall

  ! serialize
  call this%serialize( class_id )

  ! no of distinct classes
  nvclass = maxval( class_id )

  return
end subroutine get_vertex_class
!---------------------------------------------------------------------
!
subroutine get_edge_class(this, unitcell, class_id, nclass) 
  ! Purpose: Bonds in a unitcell with a unique combination of 
  ! attributes- {src_site, tgt_site, bond_type} gets a unique 
  ! 'class id'.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  type (unitcell_t), intent(in) :: unitcell
  integer(int32), intent(out) :: class_id(:)
  integer(int32), intent(out) :: nclass

  ! list of local variables
  integer(int32) :: i, j, k, n
  integer(int32) :: id, natom, nbond, min_typ, max_typ, min_sid, max_sid
  integer(int32) :: iwork1(unitcell%natom)
  integer(int32) :: iwork2(unitcell%nbond)
  integer(int32), allocatable :: iwork3(:,:,:)

  natom = unitcell%natom
  ! determine ranges
  ! vertex class range
  forall (i = 1:natom) iwork1(i) = unitcell%atom(i)%atomid
  min_sid = minval(iwork1)
  max_sid = maxval(iwork1)
  ! type range
  nbond = unitcell%nbond
  forall (i = 1:nbond) iwork2(i) = unitcell%bond(i)%type
  min_typ = minval(iwork2)
  max_typ = maxval(iwork2)

  ! workspace
  if (allocated(iwork3)) deallocate(iwork3)
  allocate( iwork3(min_sid:max_sid, min_sid:max_sid, min_typ:max_typ), stat=estat )
  call catch( estat, alloc_t, "graph::get_edge_class" )

  ! assign unique id for unique {src_site, tgt_site, bond_type} combination
  n = 1
  do i = min_sid, max_sid
    do j = min_sid, i
      do k = min_typ, max_typ
	  ! 'source-target' & 'target-source' combination give idential class
        iwork3(i,j,k) = n
        iwork3(j,i,k) = n
	  n = n + 1
      end do
    end do
  end do
  do n = 1, nbond
    id = unitcell%bond(n)%source%id
    i = unitcell%atom(id)%atomid
    id = unitcell%bond(n)%target%id
    j = unitcell%atom(id)%atomid
    k = unitcell%bond(n)%type
    class_id(n) = iwork3( i, j, k )
  end do

  ! serialize
  call this%serialize( class_id )

  ! no of distinct classes
  nclass = maxval( class_id )

  return
end subroutine get_edge_class
!---------------------------------------------------------------------
!
subroutine get_edge_class_alt(this, unitcell, vclass_id, eclass_id, neclass) 
  ! Purpose: Bonds with a unique attribute combination- 
  ! {type, src_class, tgt_class} gets a unique 'class id'.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  type (unitcell_t), intent(in) :: unitcell
  integer(int32), intent(in) :: vclass_id(:)
  integer(int32), intent(out) :: eclass_id(:)
  integer(int32), intent(out) :: neclass

  ! list of local variables
  integer(int32) :: i, j, k, n
  integer(int32) :: atom, natom, nbond, min_typ, max_typ, min_cls, max_cls
  integer(int32), allocatable :: iwork(:,:,:)

  natom = unitcell%natom
  ! determine ranges
  ! vertex class range
  forall (i = 1:natom) eclass_id(i) = vclass_id(i)
  min_cls = minval(vclass_id)
  max_cls = maxval(vclass_id)
  ! type range
  nbond = unitcell%nbond
  forall (i = 1:nbond) eclass_id(i) = unitcell%bond(i)%type
  min_typ = minval(eclass_id)
  max_typ = maxval(eclass_id)

  ! workspace
  if (allocated(iwork)) deallocate(iwork)
  allocate( iwork(min_cls:max_cls, min_cls:max_cls, min_typ:max_typ), stat=estat )
  call catch( estat, alloc_t, "graph::get_edge_class" )

  ! assign unique id for unique {src_class, tgt_class, type} combination
  n = 1
  do i = min_cls, max_cls
    do j = min_cls, i
      do k = min_typ, max_typ
	  ! 'source-target' & 'target-source' is identified
        iwork(i,j,k) = n
        iwork(j,i,k) = n
	  n = n + 1
      end do
    end do
  end do
  do n = 1, nbond
    atom = unitcell%bond(n)%source%id
    !i = vclass_id(atom)
    i = unitcell%atom(atom)%atomid
    atom = unitcell%bond(n)%target%id
    !j = vclass_id(atom)
    j = unitcell%atom(atom)%atomid
    k = unitcell%bond(n)%type
    eclass_id(n) = iwork( i, j, k )
  end do

  ! serialize
  call this%serialize( eclass_id )

  ! no of distinct classes
  neclass = maxval( eclass_id )

  return
end subroutine get_edge_class_alt
!---------------------------------------------------------------------
!
subroutine serialize(iarray) 
  ! Purpose: Integer input array elements are redefined so that in an
  ! ascending ordered sequence of the elements, successive 
  ! elements, if differ, differ by 1. The lowest element of the
  ! array is assigned value 1.
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  integer(int32), intent(inout) :: iarray(:)

  ! local variables
  integer(int32) :: i, j, k, n, ival, ndim
  integer(int32) :: idx(size(iarray))
  integer(int32) :: iwork(size(iarray))

  ndim = size(iarray)
  forall (i = 1:ndim) idx(i) = i
  iwork = iarray
  ! order
  do i = 1, ndim
    ival = iwork(idx(i))
    k = i
    do j = i, ndim
      if (ival > iwork(idx(j))) then
	  ival = iwork(idx(j))
	  k = j 
	end if
    end do
    n = idx(k)
    idx(k) = idx(i)
    idx(i) = n
  end do
  ! serialize
  iarray(idx(1)) = 1
  do i = 2, ndim
    if (iwork(idx(i)) > iwork(idx(i-1))) then
      iarray(idx(i)) = iarray(idx(i-1)) + 1
    else
      iarray(idx(i)) = iarray(idx(i-1)) 
    end if
  end do

  return
end subroutine serialize
!---------------------------------------------------------------------
!
subroutine print_unitgraph_sub(this, ofunit) 
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in), optional :: ofunit

  ! local variables
  integer(int32) :: of = stdout 
  integer(int32) :: i, sv, tv

  if (present(ofunit)) of = ofunit
  write (of,'(/,a)') '# graph:: graph of a unitcell'
  write (of,'(2a)')  '#', repeat('-', 79)
  write (of,'(a)') '# vertices:'
  do i = 1, this%glattice%unitcell%natom
    write (of,'(a,i4,a,a,i3,a,i3,a,i4,a,3f6.2,a)') '# id=', this%vertex(i)%id, ':', &
          ' type=', this%vertex(i)%type, ', class=', this%vertex(i)%class, &
	    ', atomid=', this%vertex(i)%atomid, ', coord=(', this%vertex(i)%coord, ')'
  end do
  write (of,'(/,a)') '# edges:'
  do i = 1, this%glattice%unitcell%nbond
    write (*,'(a,i4,a,a,i3,a,i3,a,i4,a,i4,a,3f6.2,a)') "# id=", this%edge(i)%id, ":", &
          " type=", this%edge(i)%type, ", class=", this%edge(i)%class, &
          ", src=", this%edge(i)%src_id, &
	    ", tgt=", this%edge(i)%tgt_id, &
	    ", vector= (", this%edge(i)%vector, ")"
  end do
  write (of,'(/,2a)')  '#', repeat('-', 79)
  write (of,'(a,18x,a,/)') '# plot data: edge source coordinate \\ target coordinate', 'bond id'
  do i = 1, this%glattice%unitcell%nbond
    sv = this%edge(i)%src_id
    tv = this%edge(i)%tgt_id
    write (of,'(1x,3ES23.15, i6)') this%vertex(sv)%coord(1), this%vertex(sv)%coord(2), &
          this%vertex(sv)%coord(3), this%edge(i)%id
    write (of,'(1x,3ES23.15,//)') this%vertex(tv)%coord(1), this%vertex(tv)%coord(2), &
          this%vertex(tv)%coord(3)
  end do

  return
end subroutine print_unitgraph_sub
!---------------------------------------------------------------------
!
subroutine print_graph_sub(this, ofunit) 
  ! Purpose: Suggestive
  !
  ! Arguments:
  ! 1.
  implicit none

  ! calling parameters
  class (latticegraph_t), intent(in) :: this
  integer(int32), intent(in), optional :: ofunit

  ! local variables
  integer(int32) :: of = stdout 
  integer(int32) :: i, j, sv, tv

  if (present(ofunit)) of = ofunit
  write (of,'(/,a)') '# graph:: lattice graph'
  write (of,'(2a)')  '#', repeat('-', 76)
  write (of,'(a,18x,a,/)') '# plot data: edge source coordinate \\ target coordinate', 'bond id'
  do i = 1, this%nedge
    sv = this%edge(i)%src_id
    tv = this%edge(i)%tgt_id
    write (of,'(1x,3ES23.15, i6)') this%vertex(sv)%coord(1), this%vertex(sv)%coord(2), &
          this%vertex(sv)%coord(3), this%edge(i)%id
    write (of,'(1x,3ES23.15,//)') this%vertex(tv)%coord(1), this%vertex(tv)%coord(2), &
          this%vertex(tv)%coord(3)
  end do
  write (of,'(a)') '# vertices'
  write (of,'(2a)')  '#', repeat('-', 90)
  do i = 1, this%nvertex
    write (of,'(a,i6,a,a,i3,a,i3,a,i6,a,3f6.2,a,3i3,a)') '# id=', this%vertex(i)%id, ':', &
          ' type=', this%vertex(i)%type, ', class=', this%vertex(i)%class, &
  	    ', atomid=', this%vertex(i)%atomid, ', coord=(', this%vertex(i)%coord, &
          '), cell=(', this%vertex(i)%cellpos, ')'
  end do
  write (of,'(/,a)') '# vertex links'
  write (of,'(2a)')  '#', repeat('-', 76)
  do i = 1, this%nvertex
    write (of,'(a,i6,a)') '# vertex =', this%vertex(i)%id, ':'
    do j = 1, this%vertex(i)%ndirlink
      write (of,'(a,8x,a,i2,a,i6)') '#', 'link-', j, ' =', this%vertex(i)%link_id(j)
    end do
  end do

  return
end subroutine print_graph_sub
!---------------------------------------------------------------------
!
end module latticegraph
