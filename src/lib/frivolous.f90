!< FriVolous, Finite Volume block-structured Fortran abstract class.
module frivolous
!< FriVolous, Finite Volume block-structured Fortran abstract class.

!-----------------------------------------------------------------------------------------------------------------------------------
use frivolous_kinds
use frivolous_adt_conservative, only : conservative
use vecfor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public :: finite_volume_block
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: finite_volume_block
  !< Finite Volume block-structured Fortran abstract class.
  !<
  !< Object-Oriented designed class for Finite Volume block-structured numerical computations: it allows the easy handling of
  !< metrics data for the robust and efficient computation of numerical spatial operators in the framework of Finite Volume Methods
  !< (FVM). It is based on a simple yet powerful Abstract Data Type (ADT) that is overloaded with useful methods for handling the
  !< back-end operations necessary to compute a numerical spatial operator independently by the actual state variable being
  !< numerically integrated.
  !<
  !< Let us assume that the fluid domain \(D\) is decomposed in \(N_b\) structured blocks \(D^b\), each subdivided in
  !< \(N_i \times N_j \times N_k\) disjoint hexahedrons \(D_{ijk}^b\) such that \(\bigcup D_{ijk}^b = D^b\).
  !< The Finite Volume abstract class is designed to aid the computations of the spatial operators into each block:
  !< $$
  !<\frac{\partial}{{\partial t}}\int\limits_{V_{ijk}} {\overrightarrow U dV}  =
  !<-\sum\limits_{s = 1}^6 {\int\limits_{S_s} {\left(\overline{\overline {F}}\right) \cdot \overrightarrow n dS}} +
  !< \int\limits_{V_{ijk}} {\overrightarrow {{Q}} dV}\label{eq:rans-cons-num}
  !< $$
  !< where \(S_s\) is the \(s^{th}\) face of the finite volume \(D_{ijk}\) whose measure is \(V_{ijk}\).
  !<
  !< A structured block is composed of hexahedron finite volumes with quadrilateral faces using the
  !< following internal numeration for nodes and faces:
  !<```
  !< /|\Z
  !<  |                            F(4)         _ F(6)
  !<  |                            /|\          /!
  !<  |                        7    |          /    8
  !<  |                         *------------------*
  !<  |                        /|   |        /    /|
  !<  |                       / |   |       /    / |
  !<  |                      /  |   |      /    /  |
  !<  |                     /   |   |     /    /   |
  !<  |                    /    |   |    +    /    |
  !<  |                   /     |   |        /     |
  !<  |                  /      |   +       /      |
  !<  |                 /      3|          /       |4
  !<  |                /        * --------/--------*
  !<  |      F(1)<----/----+   /         /        /
  !<  |              *------------------*    +-------->F(2)
  !<  |             5|       /          |6      /
  !<  |              |      /           |      /
  !<  |              |     /        +   |     /
  !<  |              |    /         |   |    /
  !<  |              |   /      +   |   |   /
  !<  |              |  /      /    |   |  /
  !<  |              | /      /     |   | /
  !<  |              |/      /      |   |/
  !<  |              *------------------*
  !<  |             1      /        |    2
  !<  |                   /        \|/
  !<  |   _ Y           |/_       F(3)
  !<  |   /|         F(5)
  !<  |  /
  !<  | /
  !<  |/                                                    X
  !<  O----------------------------------------------------->
  !<```
  !< Each hexadron cells is faces-connected to its neighboring, thus the cells build a structured block with implicit
  !< connectivity, e.g. in 2D space a FriVolous block could be as the following:
  !<```
  !<                 _ J
  !<                 /|                          _____
  !<               5+ ...*----*----*----*----*...     |
  !<               /    /    /    /    /    /         |
  !<              /    /    /    /    /    /          |
  !<            4+ ...*----*----*----*----*...        |
  !<            /    /    /    /    /    /            |
  !<           /    /    /    /    /    /             |
  !<         3+ ...*----*----*----*----*...           |  Structured block of 4x4 Finite Volumes
  !<         /    /    / FV /    /    /               |
  !<        /    /    /    /    /    /                |
  !<      2+ ...*----*----*----*----*...              |
  !<      /    /    /    /    /    /                  |
  !<     /    /    /    /    /    /                   |
  !<   1+ ...*----*----*----*----*...                 |
  !<   /     .    .    .    .    .                    |
  !<  /      .    .    .    .    .               _____
  !< O-------+----+----+----+----+-------------> I
  !<         1    2    3    4    5
  !<```
  !<The nodes of cells are not required to be on the Cartesian coordinates, thus allowing a general
  !< curvilinear mesh: the are 3 implicit coordinate lines, *i*, *j* and *k* that are not required to be orthogonal.
  private
  ! block dimensions
  integer(I_P)              :: Ng(1:6)=[0, 0, 0, 0, 0, 0] !< Number of ghost cells (e.g. used for imposing boundary conditions).
  integer(I_P)              :: Ni=0                       !< Number of cells along the i-th implicit coordinate.
  integer(I_P)              :: Nj=0                       !< Number of cells along the j-th implicit coordinate.
  integer(I_P)              :: Nk=0                       !< Number of cells along the k-th implicit coordinate.
  ! dynamic array members
  type(vector), allocatable :: node(:,:,:)                !< Nodes coordinates.
  type(vector), allocatable :: FNi(:,:,:)                 !< Faces normal along i-th implicit coordinate. NormL2 is the face area.
  type(vector), allocatable :: FNj(:,:,:)                 !< Faces normal along j-th implicit coordinate. NormL2 is the face area.
  type(vector), allocatable :: FNk(:,:,:)                 !< Faces normal along k-th implicit coordinate. NormL2 is the face area.
  real(R_P),    allocatable :: FAi(:,:,:)                 !< Faces area along i-th implicit coordinate. Optionally allocated.
  real(R_P),    allocatable :: FAj(:,:,:)                 !< Faces area along i-th implicit coordinate. Optionally allocated.
  real(R_P),    allocatable :: FAk(:,:,:)                 !< Faces area along i-th implicit coordinate. Optionally allocated.
  real(R_P),    allocatable :: volume(:,:,:)              !< Cells volume.
  ! scalar members
  type(vector)              :: emin                       !< Coordinates of minimum abscissa of the block.
  type(vector)              :: emax                       !< Coordinates of maximum abscissa of the block.
  logical                   :: store_face_area=.false.    !< Activate (separate) faces area storing.
  logical                   :: cartesian=.false.          !< Flag for checking if the block is Cartesian.
  logical                   :: nullify_x=.false.          !< Nullify X direction (2D yz, 1D y or z domain).
  logical                   :: nullify_y=.false.          !< Nullify Y direction (2D xy, 1D x or y domain).
  logical                   :: nullify_z=.false.          !< Nullify Z direction (2D xy, 1D x or y domain).
  contains
    private
    procedure, pass(self), public :: init             !< Init block.
    procedure, pass(self), public :: destroy          !< Destroy block.
    procedure, pass(self), public :: linspace         !< Create a Cartesian block with linearly spaced nodes.
    procedure, pass(self), public :: metrics          !< Compute block metrics.
    procedure, pass(self), public :: node2center      !< Compute cell centers coordinates from cell nodes.
    procedure, pass(self), public :: interpolate2node !< Interpolate cell-centered variable to nodes.
    procedure, pass(self), public :: residuals        !< Compute finite-volume residuals.
    procedure, pass(self), public :: residuals_fast   !< Compute finite-volume residuals. Fast version withoud dot product.
    procedure, pass(self), public :: is_cartesian     !< Return .true. is the block is Cartesian, .false. otherwise.
    final                         :: finalize         !< Finalize block.
    ! operators
    procedure, pass(lhs), public :: assign_block      !< Block = Block.
    ! operators overloading
    generic, public :: assignment(=) => assign_block  !< Overloading = assignament.
endtype finite_volume_block
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public methods
  subroutine init(self, store_face_area, cartesian, nullify_x, nullify_y, nullify_z, emin, emax, Ng, node, Ni, Nj, Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Init block.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(INOUT) :: self            !< Block.
  logical,      optional,     intent(IN)    :: store_face_area !< Activate (separate) faces area storing.
  logical,      optional,     intent(IN)    :: cartesian       !< Flag for checking if the block is Cartesian.
  logical,      optional,     intent(IN)    :: nullify_x       !< Nullify X direction (2D yz, 1D y or z domain).
  logical,      optional,     intent(IN)    :: nullify_y       !< Nullify Y direction (2D xy, 1D x or y domain).
  logical,      optional,     intent(IN)    :: nullify_z       !< Nullify Z direction (2D xy, 1D x or y domain).
  type(vector), optional,     intent(IN)    :: emin            !< Coordinates of minimum abscissa of the block.
  type(vector), optional,     intent(IN)    :: emax            !< Coordinates of maximum abscissa of the block.
  integer(I_P), optional,     intent(IN)    :: Ng(1:)          !< Number of ghost cells.
  integer(I_P), optional,     intent(IN)    :: Ni              !< Number of cells along the i-th implicit coordinate.
  integer(I_P), optional,     intent(IN)    :: Nj              !< Number of cells along the j-th implicit coordinate.
  integer(I_P), optional,     intent(IN)    :: Nk              !< Number of cells along the k-th implicit coordinate.
  type(vector), optional,     intent(IN)    :: node(0:,0:,0:)  !< Nodes coordinates.
  integer(I_P)                              :: Ng_(1:6)        !< Number of ghost cells, dummy varibale.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (((.not.present(Ni)).and.(.not.present(Nj)).and.(.not.present(Nk))) .and. (.not.present(node))) then
    ! raise an error: one between node array or i-j-k dimensions must be provided
    return
  endif

  call self%destroy

  ! set scalar members fisrt
  if (present(emin           )) self%emin            = emin
  if (present(emax           )) self%emax            = emax
  if (present(store_face_area)) self%store_face_area = store_face_area
  if (present(cartesian      )) self%cartesian       = cartesian
  if (present(nullify_x      )) self%nullify_x       = nullify_x
  if (present(nullify_y      )) self%nullify_y       = nullify_y
  if (present(nullify_z      )) self%nullify_z       = nullify_z

  ! set block dimensions
  Ng_ = 0 ; if (present(Ng)) Ng_ = Ng ; self%Ng = Ng_
  if (present(node)) then
    ! infer dimensions from node array
    self%Ni = size(node, dim=1) - self%Ng(1) - self%Ng(2) - 1
    self%Nj = size(node, dim=2) - self%Ng(3) - self%Ng(4) - 1
    self%Nk = size(node, dim=3) - self%Ng(5) - self%Ng(6) - 1
  else
    ! dimensions are explicitly passed
    self%Ni = Ni
    self%Nj = Nj
    self%Nk = Nk
  endif

  ! allocate dynamic array members
  allocate(self%node  (0 - Ng_(1):Ni + Ng_(2), 0 - Ng_(3):Nj + Ng_(4), 0 - Ng_(5):Nk + Ng_(6)))
  allocate(self%FNi   (0 - Ng_(1):Ni + Ng_(2), 1 - Ng_(3):Nj + Ng_(4), 1 - Ng_(5):Nk + Ng_(6)))
  allocate(self%FNj   (1 - Ng_(1):Ni + Ng_(2), 0 - Ng_(3):Nj + Ng_(4), 1 - Ng_(5):Nk + Ng_(6)))
  allocate(self%FNk   (1 - Ng_(1):Ni + Ng_(2), 1 - Ng_(3):Nj + Ng_(4), 0 - Ng_(5):Nk + Ng_(6)))
  allocate(self%volume(1 - Ng_(1):Ni + Ng_(2), 1 - Ng_(3):Nj + Ng_(4), 1 - Ng_(5):Nk + Ng_(6)))
  if (self%store_face_area) then
    allocate(self%FAi (0 - Ng_(1):Ni + Ng_(2), 1 - Ng_(3):Nj + Ng_(4), 1 - Ng_(5):Nk + Ng_(6)))
    allocate(self%FAj (1 - Ng_(1):Ni + Ng_(2), 0 - Ng_(3):Nj + Ng_(4), 1 - Ng_(5):Nk + Ng_(6)))
    allocate(self%FAk (1 - Ng_(1):Ni + Ng_(2), 1 - Ng_(3):Nj + Ng_(4), 0 - Ng_(5):Nk + Ng_(6)))
  endif

  ! set nodes coordinates if passed
  if (present(node)) self%node = node
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy block.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(INOUT) :: self  !< Block.
  ! type(finite_volume_block)                 :: fresh !< Instance of a new, fresh block used for reset scalars to defaults.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%node  )) deallocate(self%node  )
  if (allocated(self%FNi   )) deallocate(self%FNi   )
  if (allocated(self%FNj   )) deallocate(self%FNj   )
  if (allocated(self%FNk   )) deallocate(self%FNk   )
  if (allocated(self%FAi   )) deallocate(self%FAi   )
  if (allocated(self%FAj   )) deallocate(self%FAj   )
  if (allocated(self%FAk   )) deallocate(self%FAk   )
  if (allocated(self%volume)) deallocate(self%volume)
  ! self = fresh
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  elemental subroutine linspace(self, emin, emax, error)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a Cartesian block with linearly spaced nodes.
  !<
  !< @note If the extents (emin, emax) of the block are not passed, the values eventually passed during block initialization are
  !< used.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(INOUT) :: self    !< Block.
  type(vector), optional,     intent(IN)    :: emin    !< Coordinates of minimum abscissa of the block.
  type(vector), optional,     intent(IN)    :: emax    !< Coordinates of maximum abscissa of the block.
  integer(I_P), optional,     intent(OUT)   :: error   !< Error code.
  integer(I_P)                              :: error_  !< Dummy error code.
  type(vector)                              :: delta   !< Diagonal of block bounding-box.
  real(R_P)                                 :: delta_x !< X component of diagonal of block bounding-box.
  real(R_P)                                 :: delta_y !< Y component of diagonal of block bounding-box.
  real(R_P)                                 :: delta_z !< Z component of diagonal of block bounding-box.
  integer(I_P)                              :: i       !< Counter.
  integer(I_P)                              :: j       !< Counter.
  integer(I_P)                              :: k       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  error_ = 0
  self%cartesian = .true.
  if (present(emin)) self%emin = emin
  if (present(emax)) self%emax = emax
  if (self%emin/=self%emax) then
    delta = (emax - emin)/(self%Ni * ex + self%Nj * ey + self%Nk * ez)
    delta_x = delta.dot.ex
    delta_y = delta.dot.ey
    delta_z = delta.dot.ez
    do k=0 - self%Ng(5), self%Nk + self%Ng(6)
      do j=0 - self%Ng(3), self%Nj + self%Ng(4)
        do i=0 - self%Ng(1), self%Ni + self%Ng(2)
          self%node(i, j, k) = emin + (i * delta_x) * ex + (j * delta_y) * ey + (k * delta_z) * ez
        enddo
      enddo
    enddo
  else
    error_ = 1
  endif
  if (present(error)) error = error_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine linspace

  elemental subroutine metrics(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute metrics block.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(INOUT) :: self         !< Block.
  type(vector)                              :: triplet(1:9) !< Dummy vectors.
  real(R_P)                                 :: signi        !< Sign of direction of normals along i-th coordinate.
  real(R_P)                                 :: signj        !< Sign of direction of normals along i-th coordinate.
  real(R_P)                                 :: signk        !< Sign of direction of normals along i-th coordinate.
  integer(I_P)                              :: i            !< Counter.
  integer(I_P)                              :: j            !< Counter.
  integer(I_P)                              :: k            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! check direction of normals
  ! set at the middle of the block
  i = max(1, self%Ni)
  j = max(1, self%Nj)
  k = max(1, self%Nk)
  ! i coordinate
  triplet(1) = self%node(i, j,   k) - self%node(i, j-1, k-1)
  triplet(2) = self%node(i, j-1, k) - self%node(i, j,   k-1)
  triplet(3) = triplet(1).cross.triplet(2)
  triplet(4) = (0.25_R_P * (self%node(i,   j, k) + self%node(i,   j-1, k) + self%node(i,   j, k-1) + self%node(i,   j-1, k-1))) - &
               (0.25_R_P * (self%node(i-1, j, k) + self%node(i-1, j-1, k) + self%node(i-1, j, k-1) + self%node(i-1, j-1, k-1)))
  signi = sign(1._R_P, (triplet(3).dot.triplet(4)))
  ! j coordinate
  triplet(1) = self%node(i, j, k  ) - self%node(i-1, j, k-1)
  triplet(2) = self%node(i, j, k-1) - self%node(i-1, j, k  )
  triplet(3) = triplet(1).cross.triplet(2)
  triplet(4) = (0.25_R_P * (self%node(i, j,   k) + self%node(i-1, j,   k) + self%node(i, j,   k-1) + self%node(i-1, j,   k-1))) - &
               (0.25_R_P * (self%node(i, j-1, k) + self%node(i-1, j-1, k) + self%node(i, j-1, k-1) + self%node(i-1, j-1, k-1)))
  signj = sign(1._R_P, (triplet(3).dot.triplet(4)))
  ! k coordinate
  triplet(1) = self%node(i,   j, k) - self%node(i-1, j-1, k)
  triplet(2) = self%node(i-1, j, k) - self%node(i,   j-1, k)
  triplet(3) = triplet(1).cross.triplet(2)
  triplet(4) = (0.25_R_P * (self%node(i, j, k  ) + self%node(i-1, j, k  ) + self%node(i, j-1, k  ) + self%node(i-1, j-1, k  ))) - &
               (0.25_R_P * (self%node(i, j, k-1) + self%node(i-1, j, k-1) + self%node(i, j-1, k-1) + self%node(i-1, j-1, k-1)))
  signk = sign(1._R_P, (triplet(3).dot.triplet(4)))
  ! compute faces normals
  do k=1, self%Nk
    do j=1, self%Nj
      do i=0, self%Ni
        call self%FNi(i, j, k)%face_normal4(pt1 = self%node(i, j-1, k-1), &
                                            pt2 = self%node(i, j  , k-1), &
                                            pt3 = self%node(i, j  , k  ), &
                                            pt4 = self%node(i, j-1, k  ))
        self%FNi(i, j, k) = self%FNi(i, j, k) * signi
      enddo
    enddo
  enddo
  do k=1, self%Nk
    do j=0, self%Nj
      do i=1, self%Ni
        call self%FNj(i, j, k)%face_normal4(pt1 = self%node(i-1, j, k-1), &
                                            pt2 = self%node(i-1, j, k  ), &
                                            pt3 = self%node(i  , j, k  ), &
                                            pt4 = self%node(i  , j, k-1))
        self%FNj(i, j, k) = self%FNj(i, j, k) * signj
      enddo
    enddo
  enddo
  do k=0, self%Nk
    do j=1, self%Nj
      do i=1, self%Ni
        call self%FNk(i, j, k)%face_normal4(pt1 = self%node(i-1, j-1, k), &
                                            pt2 = self%node(i  , j-1, k), &
                                            pt3 = self%node(i  , j  , k), &
                                            pt4 = self%node(i-1, j  , k))
        self%FNk(i, j, k) = self%FNk(i, j, k) * signk
      enddo
    enddo
  enddo
  ! compute finite volumes
  do k=1, self%Nk
    do j=1, self%Nj
      do i=1, self%Ni
        triplet(1) = self%node(i  , j  , k  ) - self%node(i  , j-1, k-1) + self%node(i-1, j  , k  ) - self%node(i-1, j-1, k-1)
        triplet(2) = self%node(i  , j  , k  ) - self%node(i-1, j  , k-1)
        triplet(3) = self%node(i  , j  , k-1) - self%node(i-1, j-1, k-1)
        triplet(4) = self%node(i-1, j  , k  ) - self%node(i-1, j-1, k-1)
        triplet(5) = self%node(i  , j  , k  ) - self%node(i-1, j  , k-1) + self%node(i  , j-1, k  ) - self%node(i-1, j-1, k-1)
        triplet(6) = self%node(i  , j  , k  ) - self%node(i-1, j-1, k  )
        triplet(7) = self%node(i  , j  , k  ) - self%node(i  , j-1, k-1)
        triplet(8) = self%node(i  , j-1, k  ) - self%node(i-1, j-1, k-1)
        triplet(9) = self%node(i  , j  , k  ) - self%node(i-1, j-1, k  ) + self%node(i  , j  , k-1) - self%node(i-1, j-1, k-1)
        self%volume(i, j, k) = (triplet(1).dot.(triplet(2).cross.triplet(3))) + &
                               (triplet(4).dot.(triplet(5).cross.triplet(6))) + &
                               (triplet(7).dot.(triplet(8).cross.triplet(9)))
      enddo
    enddo
  enddo
  ! nullify normals for 2D or 1D domains
  if (self%nullify_x) then
    self%FNi = (self%FNi.paral.ex) + (0._R_P * ey      ) + (0._R_P * ez      )
    self%FNj = (0._R_P * ex      ) + (self%FNj.paral.ey) + (self%FNj.paral.ez)
    self%FNk = (0._R_P * ex      ) + (self%FNk.paral.ey) + (self%FNk.paral.ez)
  endif
  if (self%nullify_y) then
    self%FNi = (self%FNi.paral.ex) + (0._R_P * ey      ) + (self%FNi.paral.ez)
    self%FNj = (0._R_P * ex      ) + (self%FNj.paral.ey) + (0._R_P * ez      )
    self%FNk = (self%FNk.paral.ex) + (0._R_P * ey      ) + (self%FNk.paral.ez)
  endif
  if (self%nullify_z) then
    self%FNi = (self%FNi.paral.ex) + (self%FNi.paral.ey) + (0._R_P * ez      )
    self%FNj = (self%FNj.paral.ex) + (self%FNj.paral.ey) + (0._R_P * ez      )
    self%FNk = (0._R_P * ex      ) + (0._R_P * ey      ) + (self%FNk.paral.ez)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine metrics

  pure function node2center(self) result(center)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute cell centers coordinates from cell nodes.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(IN) :: self          !< Block.
  type(vector), allocatable              :: center(:,:,:) !< Cell centers coordinates.
  integer(I_P)                           :: i             !< Counter.
  integer(I_P)                           :: j             !< Counter.
  integer(I_P)                           :: k             !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(center(1 - self%Ng(1):self%Ni + self%Ng(2), 1 - self%Ng(3):self%Nj + self%Ng(4), 1 - self%Ng(5):self%Nk + self%Ng(6)))
  do k=1 - self%Ng(5), self%Nk + self%Ng(6)
    do j=1 - self%Ng(3), self%Nj + self%Ng(4)
      do i=1 - self%Ng(1), self%Ni + self%Ng(2)
        center(i, j, k) = (self%node(i,   j,   k  ) + &
                           self%node(i-1, j,   k  ) + &
                           self%node(i  , j-1, k  ) + &
                           self%node(i  , j  , k-1) + &
                           self%node(i-1, j-1, k-1) + &
                           self%node(i  , j-1, k-1) + &
                           self%node(i-1, j  , k-1) + &
                           self%node(i-1, j-1, k  )) * 0.125_R_P
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction node2center

  pure subroutine interpolate2node(self, var_cell, var_node)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate cell-centered variable to nodes.
  !<
  !< @note The interpolation is linear and based on the volume-weights.
  !<
  !< @note Only internal cells are considered, ghost ones are trimmed.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(IN)  :: self                                                  !< Block.
  real(R_P),                  intent(IN)  :: var_cell(1-self%Ng(1):, 1-self%Ng(3):, 1-self%Ng(5):) !< Cell-centered variable.
  real(R_P),                  intent(OUT) :: var_node(0-self%Ng(1):, 0-self%Ng(3):, 0-self%Ng(5):) !< Node-centered variable.
  real(R_P), allocatable                  :: var_cell_framed(:,:,:)                                !< Cell-centered var framed.
  real(R_P), allocatable                  :: volume_framed(:,:,:)                                  !< Volume framed.
  integer(I_P)                            :: i                                                     !< Counter.
  integer(I_P)                            :: j                                                     !< Counter.
  integer(I_P)                            :: k                                                     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! building framed variable and volume
  allocate(var_cell_framed(0:self%Ni+1, 0:self%Nj+1, 0:self%Nk+1)) ; var_cell_framed = 0._R_P
  var_cell_framed(1:self%Ni, 1:self%Nj, 1:self%Nk) = var_cell(1:self%Ni, 1:self%Nj, 1:self%Nk)
  allocate(volume_framed(0:self%Ni+1, 0:self%Nj+1, 0:self%Nk+1)) ; volume_framed = 0._R_P
  volume_framed(1:self%Ni, 1:self%Nj, 1:self%Nk) = self%volume(1:self%Ni, 1:self%Nj, 1:self%Nk)
  ! check frames
  if (self%Ng(1)>0) then
    var_cell_framed(0, 1:self%Nj, 1:self%Nk) = var_cell(   0, 1:self%Nj, 1:self%Nk)
    volume_framed(  0, 1:self%Nj, 1:self%Nk) = self%volume(0, 1:self%Nj, 1:self%Nk)
  endif
  if (self%Ng(2)>0) then
    var_cell_framed(self%Ni+1, 1:self%Nj, 1:self%Nk) = var_cell(   self%Ni+1, 1:self%Nj, 1:self%Nk)
    volume_framed(  self%Ni+1, 1:self%Nj, 1:self%Nk) = self%volume(self%Ni+1, 1:self%Nj, 1:self%Nk)
  endif
  if (self%Ng(3)>0) then
    var_cell_framed(1:self%Ni, 0, 1:self%Nk) = var_cell(   1:self%Ni, 0, 1:self%Nk)
    volume_framed(  1:self%Ni, 0, 1:self%Nk) = self%volume(1:self%Ni, 0, 1:self%Nk)
  endif
  if (self%Ng(4)>0) then
    var_cell_framed(self%Ni, 1:self%Nj+1, 1:self%Nk) = var_cell(   self%Ni, 1:self%Nj+1, 1:self%Nk)
    volume_framed(  self%Ni, 1:self%Nj+1, 1:self%Nk) = self%volume(self%Ni, 1:self%Nj+1, 1:self%Nk)
  endif
  if (self%Ng(5)>0) then
    var_cell_framed(1:self%Ni, 1:self%Nj, 0) = var_cell(   1:self%Ni, 1:self%Nj, 0)
    volume_framed(  1:self%Ni, 1:self%Nj, 0) = self%volume(1:self%Ni, 1:self%Nj, 0)
  endif
  if (self%Ng(6)>0) then
    var_cell_framed(self%Ni, 1:self%Nj, 1:self%Nk+1) = var_cell(   self%Ni, 1:self%Nj, 1:self%Nk+1)
    volume_framed(  self%Ni, 1:self%Nj, 1:self%Nk+1) = self%volume(self%Ni, 1:self%Nj, 1:self%Nk+1)
  endif
  ! interpolate on nodes
  do k=0, self%Nk
    do j=0, self%Nj
      do i=0, self%Ni
        var_node(i, j, k) = (var_cell_framed(i+1, j+1, k+1) * volume_framed(i+1, j+1, k+1) &
                          +  var_cell_framed(i  , j+1, k+1) * volume_framed(i  , j+1, k+1) &
                          +  var_cell_framed(i+1, j  , k+1) * volume_framed(i+1, j  , k+1) &
                          +  var_cell_framed(i  , j  , k+1) * volume_framed(i  , j  , k+1) &
                          +  var_cell_framed(i+1, j+1, k  ) * volume_framed(i+1, j+1, k  ) &
                          +  var_cell_framed(i  , j+1, k  ) * volume_framed(i  , j+1, k  ) &
                          +  var_cell_framed(i+1, j  , k  ) * volume_framed(i+1, j  , k  ) &
                          +  var_cell_framed(i  , j  , k  ) * volume_framed(i  , j  , k  ))/sum(volume_framed(i:i+1, j:j+1, k:k+1))
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate2node

  subroutine residuals(self, Fi, Fj, Fk, res)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute finite-volume residuals.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(IN)    :: self                                             !< Block.
  class(conservative),        intent(IN)    :: Fi(0-self%Ng(1):, 1-self%Ng(3):, 1-self%Ng(5):)  !< Fluxes on *i* interfaces.
  class(conservative),        intent(IN)    :: Fj(1-self%Ng(1):, 0-self%Ng(3):, 1-self%Ng(5):)  !< Fluxes on *j* interfaces.
  class(conservative),        intent(IN)    :: Fk(1-self%Ng(1):, 1-self%Ng(3):, 0-self%Ng(5):)  !< Fluxes on *k* interfaces.
  class(conservative),        intent(INOUT) :: res(1:, 1:, 1:)                                  !< Residuals.
  integer(I_P)                              :: i                                                !< Counter.
  integer(I_P)                              :: j                                                !< Counter.
  integer(I_P)                              :: k                                                !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! select type(Fi)
  ! class is (conservative)
    do k=1, self%Nk
      do j=1, self%Nj
        do i=1, self%Ni
          res(i, j, k) = ((Fi(i-1, j  ,  k  ).dot.self%FNi(i-1, j  , k  )) - (Fi(i, j, k).dot.self%FNi(i, j, k))  &
                       +  (Fj(i  , j-1,  k  ).dot.self%FNj(i  , j-1, k  )) - (Fj(i, j, k).dot.self%FNj(i, j, k))  &
                       +  (Fk(i  , j  ,  k-1).dot.self%FNk(i  , j  , k-1)) - (Fk(i, j, k).dot.self%FNk(i, j, k))) &
                       / self%volume(i, j, k)
        enddo
      enddo
    enddo
  ! endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine residuals

  subroutine residuals_fast(self, Fi, Fj, Fk, res)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute finite-volume residuals.
  !<
  !< @note this version of *residuals* compuation is based on the assumption that fluxes *Fi, Fj, Fk* are already provided as the
  !< *normal (to the interface) component*, thus the dot product can be omitted. As a matter of facts, during the computation of
  !< fluxes the normal and tangential components are generally computed separately: there is no need to perform the dot product
  !< twice, thus, in those cases, this residual procedure should be faster than the more general above one.
  !<
  !< @note this version of *residuals* assumes that *faces area* are stored *separately* into `self%FA(i,j,k)`, however no check is
  !< done (for not hurting performances).
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(IN)    :: self                                             !< Block.
  class(conservative),        intent(IN)    :: Fi(0-self%Ng(1):, 1-self%Ng(3):, 1-self%Ng(5):)  !< Fluxes on *i* interfaces.
  class(conservative),        intent(IN)    :: Fj(1-self%Ng(1):, 0-self%Ng(3):, 1-self%Ng(5):)  !< Fluxes on *j* interfaces.
  class(conservative),        intent(IN)    :: Fk(1-self%Ng(1):, 1-self%Ng(3):, 0-self%Ng(5):)  !< Fluxes on *k* interfaces.
  class(conservative),        intent(INOUT) :: res(1:, 1:, 1:)                                  !< Residuals.
  integer(I_P)                              :: i                                                !< Counter.
  integer(I_P)                              :: j                                                !< Counter.
  integer(I_P)                              :: k                                                !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! select type(Fi)
  ! class is (conservative)
    do k=1, self%Nk
      do j=1, self%Nj
        do i=1, self%Ni
          res(i, j, k) = ((Fi(i-1, j  ,  k  ) * self%FAi(i-1, j  , k  )) - (Fi(i, j, k) * self%FAi(i, j, k))  &
                       +  (Fj(i  , j-1,  k  ) * self%FAj(i  , j-1, k  )) - (Fj(i, j, k) * self%FAj(i, j, k))  &
                       +  (Fk(i  , j  ,  k-1) * self%FAk(i  , j  , k-1)) - (Fk(i, j, k) * self%FAk(i, j, k))) &
                       / self%volume(i, j, k)
        enddo
      enddo
    enddo
  ! endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine residuals_fast

  elemental function is_cartesian(self) result(is_it)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return .true. is the block is Cartesian, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(IN) :: self  !< Block.
  logical                                :: is_it !< Is it Cartesian or not?
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  is_it = self%cartesian
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_cartesian

  ! private methods
  subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy block.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(finite_volume_block), intent(INOUT) :: self !< Block.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  pure subroutine assign_block(lhs, rhs)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Symmetric assignment block = block.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(finite_volume_block), intent(INOUT) :: lhs !< Left hand side.
  type(finite_volume_block),  intent(IN)    :: rhs !< Right hand side.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
                             lhs%Ng              = rhs%Ng
                             lhs%Ni              = rhs%Ni
                             lhs%Nj              = rhs%Nj
                             lhs%Nk              = rhs%Nk
  if (allocated(rhs%node  )) lhs%node            = rhs%node
  if (allocated(rhs%FNi   )) lhs%FNi             = rhs%FNi
  if (allocated(rhs%FNj   )) lhs%FNj             = rhs%FNj
  if (allocated(rhs%FNk   )) lhs%FNk             = rhs%FNk
  if (allocated(rhs%FAi   )) lhs%FAi             = rhs%FAi
  if (allocated(rhs%FAj   )) lhs%FAj             = rhs%FAj
  if (allocated(rhs%FAk   )) lhs%FAk             = rhs%FAk
  if (allocated(rhs%volume)) lhs%volume          = rhs%volume
                             lhs%emin            = rhs%emin
                             lhs%emax            = rhs%emax
                             lhs%store_face_area = rhs%store_face_area
                             lhs%cartesian       = rhs%cartesian
                             lhs%nullify_x       = rhs%nullify_x
                             lhs%nullify_y       = rhs%nullify_y
                             lhs%nullify_z       = rhs%nullify_z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_block
endmodule frivolous
