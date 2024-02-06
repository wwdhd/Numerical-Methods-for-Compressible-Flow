   ! Higher-order version  
   ! Reference coordinates
   ! Degrees of freedom are numbered from 0 to MaxNumbDegrees
   ! Euler equations
   ! Property of Cranfield University


  Module Defi

   Use basis
   Use qr

  Implicit None


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  ! mesh related data types and variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  ! declaration of types
   Type Element
     Integer   NumberOfSides
     Integer   ElementType
     Real(8)   Center(2)
     Real(8)   Area,RAD,LIMITER
     Integer   NeighborsId(4)
     Real(8)   VertexCoordinates(2,4)
     Real(8)   Normals(2,4)
     Real(8)   Edges(4)
     Integer   VertexId(4)
     Integer   NeighborsSidesNum(4)
     Integer   SideBCNum(4)
   End type

   ! single stencil, all coordinates & areas are in reference coordinate system
   Type StencilData
     Integer                 iSize    ! number of cells in  the stencil
     Integer, pointer    ::  iCellIDs(:)              !  => NULL()
     Integer, pointer    ::  iElementNumSides(:)      !  => NULL()
     Real(8), pointer    ::  rElementCentres(:,:)     !  => NULL()
     Real(8), pointer    ::  rElementVertexes(:,:,:)  !  => NULL()
     Real(8), pointer    ::  rElementAreas(:)         !  => NULL()
     Real(8), pointer    ::  rAkl(:,:)                !  => NULL()
     Real(8), pointer    ::  rMa(:,:)                 !  => NULL()
     Real(8), pointer    ::  rQmatrix(:,:),  rRmatrix(:,:) !=> Null()
     Real(8), pointer    ::  rInvMatrix(:,:) !=> Null()
     Real(8), pointer    ::  rFinalMatrix(:,:)! => Null()

   End type

   ! collection of stencils for a cell, used for weno
   Type CellStencilData
     Integer    iNumberOfSides
     Integer    iNumberOfStencils
     Real(8)    rArea
     Real(8)    rJacobian(2,2)
     Real(8)    rInverseJacobian(2,2)
     Real(8), pointer    ::  rCellVertexes(:,:)   !  => NULL()
     Real(8), pointer    ::  rIntBasisFunctions(:)!  => NULL()
     Real(8), pointer    ::  rIndicatorMatrix(:,:) ! => NULL()
     Type(StencilData), pointer :: Stencils(:)    !  => NULL()
   End type

   Type VertexNeibData
     Integer    iNumberOfNeib
     Integer     iNeibIds(10) 
   End type


   ! total number of quadrilaterals
   Integer NElementMax 
   ! total number of vertises
   Integer NVerMax  

   Type(VertexNeibData), allocatable :: VertexNeib(:)

   Real(8)  CFL

   Type(Element), allocatable :: MeshElements(:)
   Real, Allocatable   :: Vertexes(:,:)
   ! first index -  stencil number
   ! 0 - central stencil
   ! 1-3 -  sectorial stencils corresponding to sides
   Type(CellStencilData), allocatable :: CellStencils(:)


   Integer :: NoNeig = -1 
   Integer     NTriMax, NQuadMax

  ! b.c. conditions data
   Type BoundaryValues
	Real(8) nw(4),Tw(4),uw(4),vw(4)
   End type

   !////////////// boundary conditions data ///////////
!   Type(BoundaryValues), allocatable :: MeshBCValues(:)

   Type BCZone
     character(25)  :: SurfaceName =''
     character(25)  :: BcName =''
     integer BCType
     Real(8) rho
     Real(8) u
     Real(8) v
     Real(8) T,p
   ! energy accomodaton coefficient
     Real(8) AlphaE
   End type 


!   Integer, parameter :: MaxNumberBCSets = 10
!   Type(BcZone)  BCTable(MaxNumberBCSets)
   Type(BcZone), allocatable ::  BCTable(:)

   !-----------  Boundary types
   Integer, parameter :: BcTypeReflective  = -10
   Integer, parameter :: BcTypeFreeStream  = -30
   Integer, parameter :: BcTypeOutFlow     = -40
   Integer, parameter :: BcTypePeriodic    = -100
   Integer, parameter :: SYMMETRY    = -60

   Integer NBSets

   ! ///////   Solver variables ///////

   Real(8)    StencilSizeMultiplayer 

   Real(8), allocatable ::  IntBasisFunction(:,:)
   Integer, allocatable   :: LocalElementNumSides(:,:,:)
   Integer :: MaxNumberOfStencils = 4
    Integer TheoNumDeg,PolynomialOrder
   Integer   MaxNumbDegrees
   Integer VertexToSide3(3,3)
   Integer VertexToSide4(4,3)
   Real(8) XPer, YPer,Xmax,Ymax,Xmin,Ymin   

   ! ///////   Solution related variables ///////
   Real(8), allocatable :: CSV(:,:,:), PV(:,:,:), Flux(:,:,:)
   Real(8), allocatable :: ExtValues(:,:,:,:),ExtVertexValues(:,:,:)
   Real(8), allocatable :: DegreesOfFreedom(:,:,:,:)
   Real(8), allocatable :: ExtrapolatedDegreesOfFreedom(:,:,:,:) ! use either two or four-point quadratures

   ! /////// other variables //////////
   Real(8) GauWei(4),GauPoints(4)
   Integer GauPointsNum,LIMITER
   Real(8)  dt, time,t_
   Real(8), parameter :: pi= 3.141592653589793
   integer it
   Real(8) Errors(0:1)
   Real(8) g1,g2,g3,g4,g5,g6,g7,g8
   Real(8) :: gamma = 1.4

   Type CandidateCell
     Integer   ID
     Real(8)   distance 
    End type

   ! --------- variables for calculaiting the execution time
   Real(8) ReconstructionTime,FluxesTime,TotalTime,LSQTime,ExtrapTime,UpdateTime
   Real(8) t0,t1


   Integer SpatialOrder, SchemeType, IniCondType
   Integer FluxType
   Integer, parameter:: Rusanovtype=1,HllcFluxType=2
   
   Integer DataOutPutFreq, MovieOutPutFreq

 End module