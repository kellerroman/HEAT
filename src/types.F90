module types
   use const, only : ip,dp
   implicit none
   type :: tBoundaryCondition
      integer :: BC_Type
      integer :: dist(2,2)
      !< Start und Stop der Jeweiligen Radbedigung
      !< (Start=1/Stop=2,Dir1=1/Dir2=2)
      !< i Face: Dir1 = k, Dir2 = j
      !< j Face: Dir1 = k, Dir2 = i
      !< k Face: Dir1 = j, Dir2 = i
      integer :: ncell(2)
      !< Ausdehnung der Randbedigung
      real(kind=dp), allocatable :: data(:,:)

      integer :: CPU_id
      integer :: a(4,3)
      !< Koeffizientenmatrix für den Übertrag der Indize Laufrichtungen
      !< i1 = a1,1 + a2,1*i+a3,1*j,a4,1*k
      !< j1 = a1,2 + a2,2*i+a3,2*j,a4,2*k
      !< k1 = a1,3 + a2,3*i+a3,3*j,a4,3*k

   end type tBoundaryCondition
   type :: tFace
      type(tBoundaryCondition), allocatable :: BC(:)
      integer :: nBC
   end type tFace

   type :: tBlock
      integer :: nPkt(3)
      integer :: nCell(3)
      real(kind=dp), allocatable :: T(:,:,:,:)
      real(kind=dp), allocatable :: a(:,:,:,:)
      real(kind=dp), allocatable :: flux(:,:,:,:)
      real(kind=dp), allocatable :: res(:,:,:,:)
      real(kind=dp), allocatable :: xyz(:,:,:,:)
      real(kind=dp), allocatable :: cell_vol(:,:,:)
!< Cellvolumen invers
      real(kind=dp), allocatable :: schwerpunkt(:,:,:,:)

      real(kind=dp), allocatable :: Edge_area (:,:,:,:)
!< Edge Length (I,J,K, EDGE_DIR)
!< EDGE_DIR: 1 ist in i-Richtung (West&Ost), 2 ist inw j-Richtung (Süd&Nord),...
!<  i+1 West ist Ost, j+1 Süd ist Nord, k+1 Back ist Front

      real(kind=dp), allocatable :: dn(:,:,:,:)
      real(kind=dp), allocatable :: dt(:,:,:)
      real(kind=dp), allocatable :: dn2_dt(:,:,:)
      !< Zellgröße für die berechnung des Zeitschritts 0.5/(1/dx²+1/dy²+1/dz²)
      !< als dn wird max(n) + min(n) verwendet
      real(kind=dp), allocatable :: mat_geo(:,:,:,:)
      !< Matrix Einträge (pos,i,j,k)
      !< pos: -3: k-1
      !<      -2: j-1
      !<      -1: i-1
      !<       0: i,j,k
      !<       1: i+1
      !<       2: j+1
      !<       3: k+1
      type(tFace)    :: face(6)
   !< GIBT BLOCKVERBDINUNG DES BLOCKES AN: (FACE)
   !< FACE: 1 = W, 2 = E, 3 = S, 4 = N, 5 = B, 6 = F
   end type tBlock
end module
