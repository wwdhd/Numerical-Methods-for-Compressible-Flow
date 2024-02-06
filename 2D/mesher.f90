   ! Last update: 06/08/2009

  Module Mesher

   Use basis
   Use qr
   Use Defi


  Implicit None
  
  Contains
  
    
! *************************************************************************************************
  Subroutine  InitMesh
   Integer I,k,k_up,L,j
   Real(8)  Xa(2),Xb(2),Xc(2),Xd(2),dl(2),scalar,xco(3),yco(3),normal1(2),normal2(2),s
   Real(8)  x1(2),x2(2),x3(2)
   Real(8) area1,area2
   integer  v1,v2,v3,v4,voi ! ,nbsets
   character a
   character(50) name
   integer NEntry ! number of data records in boundary-condition set
   Integer bcnum,bctype
   Integer status,file_id
   Real(8) nx,ny,Dist


   VertexToSide3(1,1) = 1;   VertexToSide3(1,2) = 1;   VertexToSide3(1,3) = 2
   VertexToSide3(2,1) = 2;   VertexToSide3(2,2) = 2;   VertexToSide3(2,3) = 3
   VertexToSide3(3,1) = 3;   VertexToSide3(3,2) = 3;   VertexToSide3(3,3) = 1


   VertexToSide4 = 0
!   VertexToSide(1,1) = 1;  
    VertexToSide4(1,2) = 1;   VertexToSide4(1,3) = 2
!   VertexToSide(2,1) = 2;  
    VertexToSide4(2,2) = 2;   VertexToSide4(2,3) = 3
!   VertexToSide(3,1) = 3; 
     VertexToSide4(3,2) = 3;   VertexToSide4(3,3) = 4
!   VertexToSide(4,1) = 4; 
     VertexToSide4(4,2) = 4;   VertexToSide4(4,3) = 1


   open(111,file='mesh.ini')
   ! skip first lines
   do k=1,6
   read(111,*) ! a; this is to avoid the problems with a blank line under Vista
 !  print*,a
   enddo


   ! now read data
   read(111,*) NVerMax, NElementMax,k,NBSets

   print*,' NvertexMax=', nvermax

  ! now allocate the arrays of vertexes  and cells

   allocate(Vertexes(2,NVerMax))
   allocate(MeshElements(NElementMax))
   allocate(IntBasisFunction(0:MaxNumbDegrees,NElementMax))

!   Type VertexNeibData
!     Integer    iNumberOfNeib
!     Integer     iNeibIds(10) 
!   End type

    allocate(VertexNeib(NVerMax))
    do i=1,NVerMax
     VertexNeib(i)%iNumberOfNeib = 0
     VertexNeib(i)%iNeibIds = 0
    enddo



   ! array containing a collection of stencils for the cell
   allocate(CellStencils(NElementmax))

  ! allocate the array with b.c. records
  allocate(BCTable(nbsets))


   read(111,*) a
   read(111,*) a
   Do i=1,NVerMax
    Read(111,*) k,Vertexes(1,i),Vertexes(2,i)
   Enddo

   print*,' done reading vertexes'
   

   read(111,*) a
   read(111,*) a

   !1 = Edge  !2 = Quadrilateral  !3 = Triangle  !4 = Brick  !5 = Wedge (Prism)  !6 = Tetrahedron  !7 = Pyramid

   NQuadMax =0
   NTriMax =0

   Do i=1,NElementMax
     Read(111,'(4x,i4)',advance='no') k
     Read(111,'(2x,i2)',advance='no') MeshElements(i)%ElementType
     Read(111,'(i2)',advance='no')    MeshElements(i)%NumberOfSides  ! seems to work, but still weird

    select case(MeshElements(i)%ElementType)
	 case(3)

	      NTriMax = NTriMax + 1
     
          Read(111,*)  MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2), MeshElements(i)%VertexId(3)

          Xa = Vertexes(:,MeshElements(i)%VertexID(1))
          Xb = Vertexes(:,MeshElements(i)%VertexID(2))
          Xc = Vertexes(:,MeshElements(i)%VertexID(3))

          MeshElements(i)%VertexCoordinates(:,1) = xa
          MeshElements(i)%VertexCoordinates(:,2) = xb
          MeshElements(i)%VertexCoordinates(:,3) = xc

          MeshElements(i)%Area = abs(0.5*(XA(1)*XB(2) - XB(1)*XA(2)  + XB(1)*XC(2) - XC(1)*XB(2)+ XC(1)*XA(2) - XA(1)*XC(2)))
          MeshElements(i)%Center(1) = (XA(1) + XB(1)  + XC(1))/3
          MeshElements(i)%Center(2) = (XA(2) + XB(2)  + XC(2))/3      
         

	 case(2)

           NQuadMax = NQuadMax + 1

           Read(111,*)  MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2),&
                        MeshElements(i)%VertexId(3), MeshElements(i)%VertexId(4)

           Xa = Vertexes(:,MeshElements(i)%VertexID(1))
           Xb = Vertexes(:,MeshElements(i)%VertexID(2))
           Xc = Vertexes(:,MeshElements(i)%VertexID(3))
           Xd = Vertexes(:,MeshElements(i)%VertexID(4))

           MeshElements(i)%VertexCoordinates(:,1) = xa
           MeshElements(i)%VertexCoordinates(:,2) = xb
           MeshElements(i)%VertexCoordinates(:,3) = xc
           MeshElements(i)%VertexCoordinates(:,4) = xd

        ! to compute the area, split the quadrilateral into 2 triangles:
     !  xa,xb,xc  and xc,xd,xa
          x1 = xa;          x2 = xb;          x3 = xc
          Area1 = abs(0.5*(X1(1)*X2(2) - X2(1)*X1(2)  + X2(1)*X3(2) - X3(1)*X2(2)+ X3(1)*X1(2) - X1(1)*X3(2)))
	
          x1 = xc;        x2 = xd;        x3 = xa
          Area2 = abs(0.5*(X1(1)*X2(2) - X2(1)*X1(2)  + X2(1)*X3(2) - X3(1)*X2(2)+ X3(1)*X1(2) - X1(1)*X3(2)))

          MeshElements(i)%Area = Area1 + Area2
          MeshElements(i)%Center(1) = (XA(1) + XC(1))/2
          MeshElements(i)%Center(2) = (XA(2) + XC(2))/2
    end select

    ! compute side lengths and normal vectors
   ! L is the edge number
   
    Do L=1,MeshElements(i)%NumberOfSides

	 select case(MeshElements(i)%NumberOfSides)
	 case(3)
           Xa = Vertexes(:,MeshElements(i)%VertexID(VertexToSide3(L,2)))
           Xb = Vertexes(:,MeshElements(i)%VertexID(VertexToSide3(L,3)))
	 case(4)
           Xa = Vertexes(:,MeshElements(i)%VertexID(VertexToSide4(L,2)))
           Xb = Vertexes(:,MeshElements(i)%VertexID(VertexToSide4(L,3)))
     case default
      print*,' error in mesh setup'
      read*           
	 end select

     ! Store the edge size
     Dist =distance(xa,xb)
     MeshElements(i)%Edges(L) = Dist

     ! outward normal vector
     nx =    (xb(2) - xa(2)) /Dist
     ny =  - (xb(1) - xa(1))/Dist
     MeshElements(i)%Normals(1,L)  = nx
     MeshElements(i)%Normals(2,L)  = ny     
   Enddo ! L

  ! prestore info of vertex neibs
    Do L=1,MeshElements(i)%NumberOfSides    
      ! retrieve the number of the vertex
      k =  MeshElements(i)%VertexId(L)
      VertexNeib(k)%iNumberOfNeib = VertexNeib(k)%iNumberOfNeib + 1
      j = VertexNeib(k)%iNumberOfNeib
      VertexNeib(k)%iNeibIds(j) = i
    Enddo ! L

   Enddo ! i=1,NElementMax
   
   open(555,file='vertexneib.txt')
   do i=1,NverMax
   write(555,*)i,VertexNeib(i)%iNumberOfNeib,VertexNeib(i)%iNeibIds(1:VertexNeib(i)%iNumberOfNeib)
   enddo
   close(555)

  ! we have now finished reading the nodes & triangles. now we start reading b.c. sets   

   print*,' now read the boundary table & conditions'
   print*,' nbsets =',Nbsets

   open(222,file='bc_set.ini')
      do k=1,NBSets
        read(222,*) Bctable(k)%BcType
        read(222,*) Bctable(k)%rho
        read(222,*) Bctable(k)%u
        read(222,*) Bctable(k)%v
        read(222,*) Bctable(k)%p

	Select case(Bctable(k)%BcType)
	 case(BcTypePeriodic)
       Bctable(k)%bcname = 'periodic'
	 case(bctypereflective)
       Bctable(k)%bcname = 'reflective'
	 case(bctypefreestream)
       Bctable(k)%bcname ='freestream'
	 case(bctypeoutflow)
       Bctable(k)%bcname ='outflow'
        case(SYMMETRY)
       Bctable(k)%bcname ='SYMMETRY'
     case default
       print*, ' Unknown boundary type ! '
       print*, ' Press any key to stop the code'
       read*
       stop
	End select
      enddo
   close(222)

!    read(111,*) a
! 
! 
! ! skip the fluid part of the mesh
!    read(111,*) a ! element group
!    read(111,*) a ! group
!    read(111,*) k ! '0'
!    
!    if (mod(nelementmax,10) .eq. 0) then
!     k_up = NelementMax/10   
!    else    
!     k_up = NelementMax/10  +1
!    endif    
!    Do k=1,k_up
!     read(111,*)
!    Enddo
!    read(111,*) a
! 
!    
!    read(111,*) a
   print*,' loop over all boundary conditions sets'
   
   Do k=1,NBSets
!     read(111,*) a
!    name = ''
    read(111,*) Name,voi,NEntry
	Bctable(k)%SurfaceName = Name
	Do j=1,NEntry
	  ! read the element number, type (here always = 3) and side number
	  read(111,*) i,voi,L
	  ! now assign the corresponding side of the triangle i the number of the boundary condition
	  MeshElements(i)%SideBCNum(L) = k
	Enddo
	! advance by character line
! 	 read(111,*) a
   Enddo

   close(111)

   print*,' NTriMax     = ',NTrimax
   print*,' NQuadMax    = ',Nquadmax
   print*,' NElementMax = ',NElementMax


   ! compute Xmax,Ymax
   xmax = 0.
   xmin = 10.
   ymax = 0.
   ymin = 10.
  
   Do i=1,NVerMax
    xmax = max(xmax,vertexes(1,i))
    xmin = min(xmin,vertexes(1,i))
    ymax = max(ymax,vertexes(2,i))
    ymin = min(ymin,vertexes(2,i))
   Enddo

   Xper = Xmax-Xmin
   Yper = Ymax-Ymin


  ! Now I have to find a neibour for each side of the element  
  ! Loop over all elements
  
  
   print*,' look for neighbours'

  ! Now I have to find a neibour for each side of the element
  
    file_id = 555
  
    Open(file_id,file='neib.ini',IOSTAT=STATUS,STATUS='OLD')
    IF (STATUS .EQ. 0 ) THEN
     ! read neigs from the file
     print*,' read neibs from the file'
     Do i=1,NElementMax
      select case(MeshElements(i)%NumberOfSides)
      case(3)
        read(file_id,*)MeshElements(i)%NeighborsID(1:3)
      case(4)
        read(file_id,*)MeshElements(i)%NeighborsID(1:4)
      end select   
     Enddo         
     close(file_id)
     
    Else
    
      open(file_id,file='neib.ini')   
   
   ! first assign no neighbour to everyone
   
   Do i=1,NElementMax
     MeshElements(i)%NeighborsId(:)  = NoNeig
   Enddo
  
   
   Do i=1,NElementMax
   
    if (mod(i,nelementmax/10) .eq. 0) print*,' done ',i*1.0/nelementmax,' percent'
    
    
    ! loop over edges
    Do L=1,MeshElements(i)%NumberOfSides
  

    ! if there is already a neiboug assigned, skip the edge
    if (MeshElements(i)%NeighborsID(L) .ne.  NoNeig) goto 101    

    Normal1 = MeshElements(i)%Normals(:,L)	 

	! now I have to find the local vertex number which is EdgeToSide(L,1),EdgeToSide(L,2)
	 select case(MeshElements(i)%NumberOfSides)
	  case(3)
           v1 = MeshElements(i)%VertexId(VertexToSide3(L,2))
           v2 = MeshElements(i)%VertexId(VertexToSide3(L,3))
	  case(4)
           v1 = MeshElements(i)%VertexId(VertexToSide4(L,2))
           v2 = MeshElements(i)%VertexId(VertexToSide4(L,3))
	  end select

	 ! first part of the loop
  	 Do k=i+1,NElementMax
	  If (k .eq. i) cycle
	  ! loop over the edges of this other triangle
	  Do j=1,MeshElements(k)%NumberOfSides
	    select case(MeshElements(k)%NumberOfSides)
		 case(3) 
                   v3 = MeshElements(k)%VertexId(VertexToSide3(j,2))
                   v4 = MeshElements(k)%VertexId(VertexToSide3(j,3))
		 case(4)
                   v3 = MeshElements(k)%VertexId(VertexToSide4(j,2))
                   v4 = MeshElements(k)%VertexId(VertexToSide4(j,3))
 		end select

	    ! now check if the coincide
	    if ( max(abs(v1-v4), abs(v2-v3))  .eq. 0   .or.  &
	         max(abs(v1-v3), abs(v2-v4))  .eq. 0) then
	      ! the edges coincide
	      ! assign the neiboubour for the size L of cell i
	 	  MeshElements(i)%NeighborsID(L) = k		 
	      ! assign the neiboubour for the size j of cell k
	 	  MeshElements(k)%NeighborsID(j) = i		 
 		 Goto 101
 	    endif

	  Enddo ! j
	Enddo ! k 

	 ! second part of the loop
  	 Do k=i-1,1,-1
	  If (k .eq. i) cycle
	  ! loop over the edges of this other triangle
	  Do j=1,MeshElements(k)%NumberOfSides
	    select case(MeshElements(k)%NumberOfSides)
		 case(3) 
                   v3 = MeshElements(k)%VertexId(VertexToSide3(j,2))
                   v4 = MeshElements(k)%VertexId(VertexToSide3(j,3))
		 case(4)
                   v3 = MeshElements(k)%VertexId(VertexToSide4(j,2))
                   v4 = MeshElements(k)%VertexId(VertexToSide4(j,3))
 		end select

	    ! now check if the coincide
	    if ( max(abs(v1-v4), abs(v2-v3))  .eq. 0   .or.  &
	         max(abs(v1-v3), abs(v2-v4))  .eq. 0) then
	      ! the edges coincide
	      ! assign the neiboubour for the size L of cell i
	 	  MeshElements(i)%NeighborsID(L) = k		 
	      ! assign the neiboubour for the size j of cell k
	 	  MeshElements(k)%NeighborsID(j) = i		 
 		 Goto 101
 	    endif

	  Enddo ! j
	Enddo ! k 

     101 continue

    Enddo ! do L=1,NumberOfSides
    
    select case(MeshElements(i)%NumberOfSides)
     case(3)
      write(file_id,*)MeshElements(i)%NeighborsID(1:3)
     case(4)
      write(file_id,*)MeshElements(i)%NeighborsID(1:4)
    end select      
    
   L = 0
  Enddo ! do i= 1,NElementMax
  close(file_id)
  Endif ! neibs


   print*,' write demo.res'

    open(1,file='demores.dat')
    write(1,*) 'TITLE = "Demo""'
    write(1,*)' VARIABLES = "x" "y" "rho"'
    write(1,*)' ZONE N=',NVerMax,' E=   ',NELEMENTMAX,'  F=FEPOINT  ET=QUADRILATERAL'
    Do i=1,NVerMax
      write(1,*), Vertexes(1,i), Vertexes(2,i), exp(- 5*(Vertexes(1,i)-0.5)**2 - 5*(Vertexes(2,i)-0.5)**2)
    Enddo
   
    Do i=1,NElementMax
	if (MeshElements(i)%NumberOfSides .eq. 3) then
	   write(1,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2), MeshElements(i)%VertexId(3), MeshElements(i)%VertexId(1)
	else
	   write(1,*) MeshElements(i)%VertexId(1), MeshElements(i)%VertexId(2), MeshElements(i)%VertexId(3), MeshElements(i)%VertexId(4)
	endif
    Enddo	
   close(1)


  print*,' assign periodic boundaries'


  ! Now I have to look at cells without the neigbough and assign periodic boundaries
  Do i=1,NElementMax

   ! loop over edges
   Do L=1,MeshElements(i)%NumberOfSides

!    if ( MeshElements(i)%NeighborsID(L) .eq. NoNeig .and.  &
!	     MeshElements(i)%SideBCNum(L) .eq. BcTypePeriodic)  then

    ! first check that the side does not have a neigbour
     k = MeshElements(i)%NeighborsID(L)
    if (k  .ne. NoNeig)  cycle
   
    ! now check that this is periodic condition
      ! retrieve the number of the boundary
      BcNum = MeshElements(i)%SideBCNum(L)
      ! retrieve bc type of the boundary
      BcType  =Bctable(BcNum)%BCType

	  If (bctype .ne. BcTypePeriodic) cycle      

     normal1 = MeshElements(i)%Normals(:,L)   
	
	 ! now I have to find the local vertex number which is EdgeToSide(L,1),EdgeToSide(L,2) minus pereodinc thing
	 select case(MeshElements(i)%NumberOfSides)
	 case(3)
            Xa = Vertexes(:,MeshElements(i)%VertexID(VertexToSide3(L,2)))
            Xb = Vertexes(:,MeshElements(i)%VertexID(VertexToSide3(L,3)))
	 case(4)
            Xa = Vertexes(:,MeshElements(i)%VertexID(VertexToSide4(L,2)))
            Xb = Vertexes(:,MeshElements(i)%VertexID(VertexToSide4(L,3)))
	 end select

	 ! now check if it is horizontal or vertical
	 If (abs(MeshElements(i)%Normals(2,L)) .le. 1e-5) then
	   ! vertical edge
            xa(1) = xa(1) - XPer*MeshElements(i)%Normals(1,L)
            xb(1) = xb(1) - XPer*MeshElements(i)%Normals(1,L)
	 else
	   ! horizontal edge
            xa(2) = xa(2) - YPer*MeshElements(i)%Normals(2,L)
            xb(2) = xb(2) - YPer*MeshElements(i)%Normals(2,L)
	 endif

    ! now search for the oppozite 
	 ! Loop over all other triangles and search for the same edge
  	 Do k=1,NElementMax

	  If (k .eq. i) cycle

	  ! loop over the edges of this other triangle
	  Do j=1,MeshElements(k)%NumberOfSides

	   select case(MeshElements(k)%NumberOfSides)
	    case(3)
             Xc = Vertexes(:,MeshElements(k)%VertexID(VertexToSide3(j,2)))
             Xd = Vertexes(:,MeshElements(k)%VertexID(VertexToSide3(j,3)))
	    case(4)
             Xc = Vertexes(:,MeshElements(k)%VertexID(VertexToSide4(j,2)))
             Xd = Vertexes(:,MeshElements(k)%VertexID(VertexToSide4(j,3)))
	   end select

	    ! now check if the coincide
	    if ( max(distance(xa,xd),distance(xb,xc)) .le. 1e-5) then
	      ! the edges coincide
	 	 MeshElements(i)%NeighborsID(L) = k		 
 		 Goto 102
 	    endif

	    if (distance(xa,xc) .le. 1e-5 .and.  distance(xb,xd) .le. 1e-5) then
	     ! the edges coincide
		 MeshElements(i)%NeighborsID(L) = k		 
	       goto 102
	    endif
	  Enddo ! j=1,NumberOfSides
	Enddo ! k=1,NElementMax

    102 continue
   Enddo ! do L=1,NumberOfSides
   L = 0
   Enddo ! do i= 1,NElementMax



   Print*,' Done reading the mesh'

  End subroutine


 ! ****************************
  Real(8) function distance(x,y)
   real(8) x(2), y(2)
    distance = sqrt(  (x(1) - y(1))**2 + (x(2) - y(2))**2 )
  End function


   ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function  CalcIntBasisFunction(i)
    integer i,k
     real(8) v1(2),v2(2),v3(2),v4(2)
     real(8),dimension(0:MaxNumbDegrees) :: CalcIntBasisFunction
     
     CalcIntBasisFunction(:) = 0.
     select case(MeshElements(i)%NumberOfSides)
        case(3)
           v1 = CellStencils(i)%rCellVertexes(1,:)
           v2 = CellStencils(i)%rCellVertexes(2,:)
           v3 = CellStencils(i)%rCellVertexes(3,:)
           v4  = 0. 
        case(4)
           v1 = CellStencils(i)%rCellVertexes(1,:)
           v2 = CellStencils(i)%rCellVertexes(2,:)
           v3 = CellStencils(i)%rCellVertexes(3,:)
           v4 = CellStencils(i)%rCellVertexes(4,:)
     end select
     Do k=1, MaxNumbDegrees
      ! here use the fact  that integrals are so far zero
       CalcIntBasisFunction(k) = ComputeBasisFunctionIntegralsOverElement(v1,v2,v3,v4,MeshElements(i)%NumberOfSides,k,i)
      Enddo
   End function



  ! ------------------------------------
  Subroutine ComputeNeibSidesID
   Integer i,L,K,M,norm,NeibID
   Real vel, UL, UR,US,xa(2),xb(2),xs(2),xs2(2),U,Du,dl(2),scalar
   Real normal1(2), normal2(2),s,signvel

   Do i=1,NElementMax
     ! loop over edges
      Do L=1,MeshElements(i)%NumberOfSides
         select case(MeshElements(i)%NumberOfSides)
             case(3)
               Xa = MeshElements(i)%VertexCoordinates(:,VertexToSide3(L,2))
               Xb = MeshElements(i)%VertexCoordinates(:,VertexToSide3(L,3))
             case(4)
               Xa = MeshElements(i)%VertexCoordinates(:,VertexToSide4(L,2))
               Xb = MeshElements(i)%VertexCoordinates(:,VertexToSide4(L,3))
          end select
          xs = (xa+xb)/2

          ! now find the corresponding face of the neigbour
          ! instead of normal vector use the side center
          NeibID = MeshElements(i)%NeighborsId(L)
          If (neibid <0) goto 101
          do k=1,MeshElements(NeibID)%NumberOfSides
            select case(MeshElements(NeibID)%NumberOfSides)
             case(3)
               Xa = MeshElements(NeibID)%VertexCoordinates(:,VertexToSide3(k,2))
               Xb = MeshElements(NeibID)%VertexCoordinates(:,VertexToSide3(k,3))
             case(4)
               Xa = MeshElements(NeibID)%VertexCoordinates(:,VertexToSide4(k,2))
               Xb = MeshElements(NeibID)%VertexCoordinates(:,VertexToSide4(k,3))
            end select
            xs2 = (xa+xb)/2

	        if (abs(xs2(1) - xs(1)) .ge. XPer/2)     xs2(1) = xs2(1) + XPer*sign(1.D0,xs(1) - xs2(1))
  	        if (abs(xs2(2) - xs(2)) .ge. YPer/2)     xs2(2) = xs2(2) + YPer*sign(1.D0,xs(2) - xs2(2))

            s = distance(xs,xs2)
            if (s .le. 1e-5) then
                MeshElements(i)%NeighborsSidesNum(L) = k
                goto 101
   	        endif	    
	     enddo
         ! nothing is found
         print*,' Error!!',' i=',i,' L=',L
         read*

         101 continue

      Enddo ! L=1,3

   Enddo  ! i = 1, NElementMax

  End subroutine




 ! *************************************************************

  Real(8) function TriangleArea(xa,xb,xc)
    real(8) xa(2),xb(2),xc(2)
     TriangleArea = abs(0.5*(XA(1)*XB(2) - XB(1)*XA(2)  + XB(1)*XC(2) - XC(1)*XB(2)+ XC(1)*XA(2) - XA(1)*XC(2)))
   End function


 ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ! compute the integrals of basis functions over the given triangle triangle
  Real(8) Function ComputeBasisFunctionIntegralsOverElement(v1,v2,v3,v4,NumSides,L,i)
    integer L,StMax,i,NumSides
	Real(8) x,y,v1(2),v2(2),v3(2),v4(2),s1,s2
  
     Select case(NumSides)
      case(3)

	   ComputeBasisFunctionIntegralsOverElement =  &
	          ComputeBasisFunctionIntegralsOverTria(v1,v2,v3,L,i)*TriangleArea(v1,v2,v3)

      case(4)

       s1 =  ComputeBasisFunctionIntegralsOverTria(v1,v2,v3,L,i)*TriangleArea(v1,v2,v3) 
       s2 =  ComputeBasisFunctionIntegralsOverTria(v3,v4,v1,L,i)*TriangleArea(v3,v4,v1)
	   ComputeBasisFunctionIntegralsOverElement = s1 + s2

     End select
  End function


  Subroutine Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,v1,v2,v3)
    Real(8)  qpoints(2,12), qweights(12),v1(2),v2(2),v3(2)
!    Real(8) :: qpoints3_coeff(3,3) = (/0.5D0,0.5D0,0.D0,0.5D0,0.5D0,0.D0/)
    Integer NumberOfPoints
     
    Select case(Spatialorder)


      case(1)
        ! this quadrature is exact for parabolas
        NumberOfPoints = 1
        qweights(1:NumberOfPoints)   = 1.D0
        qpoints(:,1) = (v1 + v2 + v3)/3

       case(2,3,4)
         ! this quadrature is exact for parabolas
 
         NumberOfPoints = 3
 
         qweights(1:NumberOfPoints)   = 1.D0/3
         qpoints(:,1) = 0.5D0*v1 + 0.5D0*v2 + 0.0D0*v3
         qpoints(:,2) = 0.5D0*v1 + 0.0D0*v2 + 0.5D0*v3
         qpoints(:,3) = 0.0D0*v1 + 0.5D0*v2 + 0.5D0*v3 
 

       case(5,6)
 
       ! this quadrature is accurate for polynomials of up to fifth order

        NumberOfPoints = 7

        ! quadrature points
        ! central point
        qpoints(:,1) = (v1+v2+v3)/3.
        qweights(1)   = 0.225

	    ! first group
        qweights(2:4)   = 0.125939180544827
        qpoints(:,2) = 0.797426985353087*v1 + 0.101286507323456*v2 + 0.101286507323456*v3
        qpoints(:,3) = 0.101286507323456*v1 + 0.797426985353087*v2 + 0.101286507323456*v3
        qpoints(:,4) = 0.101286507323456*v1 + 0.101286507323456*v2 + 0.797426985353087*v3

	    ! second  group
        qweights(5:7) = 0.132394152788506
        qpoints(:,5)  = 0.059715871789770*v1 + 0.470142064105115*v2 + 0.470142064105115*v3
        qpoints(:,6)  = 0.470142064105115*v1 + 0.059715871789770*v2 + 0.470142064105115*v3
        qpoints(:,7)  = 0.470142064105115*v1 + 0.470142064105115*v2 + 0.059715871789770*v3

     case default

       ! this quadrature is accurate for polynomials of up to seventh order

        NumberOfPoints = 12

        ! quadrature points
	    ! first group
        qweights(1:3)   = 0.050844906370207
        qpoints(:,1) = 0.873821971016996*v1 + 0.063089014491502*v2 + 0.063089014491502*v3
        qpoints(:,2) = 0.873821971016996*v2 + 0.063089014491502*v1 + 0.063089014491502*v3
        qpoints(:,3) = 0.873821971016996*v3 + 0.063089014491502*v1 + 0.063089014491502*v2

	    ! second group
        qweights(4:6)   = 0.116786275726379
        qpoints(:,4) = 0.501426509658179*v1 + 0.249286745170910*v2 + 0.249286745170910*v3
        qpoints(:,5) = 0.501426509658179*v2 + 0.249286745170910*v1 + 0.249286745170910*v3
        qpoints(:,6) = 0.501426509658179*v3 + 0.249286745170910*v1 + 0.249286745170910*v2

	    ! third group
        qweights(7:12)   = 0.082851075618374
	    ! (a,b,c)
        qpoints(:,7) = 0.636502499121399*v1 + 0.310352451033785*v2 + 0.053145049844816*v3
	    ! (a,c,b)
        qpoints(:,8) = 0.636502499121399*v1 + 0.310352451033785*v3 + 0.053145049844816*v2
	    ! (b,a,c)
        qpoints(:,9) = 0.636502499121399*v2 + 0.310352451033785*v1 + 0.053145049844816*v3
	    ! (b,c,a)
        qpoints(:,10) = 0.636502499121399*v2 + 0.310352451033785*v3 + 0.053145049844816*v1
    	! (c,a,b)
        qpoints(:,11) = 0.636502499121399*v3 + 0.310352451033785*v1 + 0.053145049844816*v2
	    ! (c,b,a)
        qpoints(:,12) = 0.636502499121399*v3 + 0.310352451033785*v2 + 0.053145049844816*v1

     End select
  
   End subroutine



 ! *******************************************************************
 ! compute the integrals of basis functions over the given triangle triangle
    Real(8) Function ComputeBasisFunctionIntegralsOverTria(v1,v2,v3,L,i)
      integer k,L,StMax,i,NumberOfPoints
      Real(8) x,y,v1(2),v2(2),v3(2),s
      Real(8) qpoints(2,12), qweights(12), int
 
      CALL Set2DIntegrationRule(qpoints,qweights,NumberOfPoints,v1,v2,v3)
      ! now compute the integrals
      Int = 0.
      Do k=1, NumberOfPoints 
            Int = Int + BasisFunctions(qpoints(1,k),qpoints(2,k),L,i)*qweights(k)
      Enddo
	  ComputeBasisFunctionIntegralsOverTria = Int
  End function


 ! sort the arrays
 Subroutine BubbleSort(N,List)
  integer N,i,j
  Type(CandidateCell) List(1:N), p

  do i=1, N
   do j=1, N-i
    if (List(j)%distance > List(j+1)%distance) then
      p= List(j)
      List(j) = List(j+1)
      List(J+1) = p
    endif
   enddo ! j
  enddo  ! i
 End subroutine

 ! ----------------------
 Subroutine ComputeJacobians(v1,v2,v3,Jacobian,InverseJacobian,Determ)
   real(8) v1(2),v2(2),v3(2),determ,Jacobian(2,2),InverseJacobian(2,2)

      Jacobian(1,1) = V2(1) - V1(1); 	 Jacobian(1,2) = V3(1) - V1(1)
	  Jacobian(2,1) = V2(2) - V1(2); 	 Jacobian(2,2) = V3(2) - V1(2)
      Determ = Jacobian(1,1)*Jacobian(2,2) - Jacobian(1,2)*Jacobian(2,1)
	  InverseJacobian(1,1) = V3(2) - V1(2);  InverseJacobian(1,2) = -(V3(1) - V1(1))
	  InverseJacobian(2,1) = -(V2(2) - V1(2));   InverseJacobian(2,2) = V2(1) - V1(1)
      InverseJacobian(:,:) = InverseJacobian(:,:)/Determ

 End subroutine


  ! ! define the transformation for cell i
  Subroutine DefineTransformationVertexes(i,v1,v2,v3)
   integer i
   Real(8) v1(2),v2(2),v3(2)

   Select case(MeshElements(i)%NumberOfSides)
     case(3)
          v1  = MeshElements(i)%VertexCoordinates(:,1) !  LocalElementVertexes(0,i,0,1,:)
          v2  = MeshElements(i)%VertexCoordinates(:,2) !  LocalElementVertexes(0,i,0,2,:)
          v3  = MeshElements(i)%VertexCoordinates(:,3) !  LocalElementVertexes(0,i,0,3,:)
     case(4)
          v1  = MeshElements(i)%VertexCoordinates(:,1) !  LocalElementVertexes(0,i,0,1,:)
          v2  = MeshElements(i)%VertexCoordinates(:,2) !  LocalElementVertexes(0,i,0,2,:)
          v3  = MeshElements(i)%VertexCoordinates(:,4) !  LocalElementVertexes(0,i,0,4,:)
    End select
 
  End subroutine

 ! *****************************************************************************
  Real(8) function BasisFunctions(x,y,number,i)
   integer number,i
   real(8) x,y,s
!   BasisFunctions = UnmodifiedBasisFunctions(X,Y,number) - IntBasisFunction(number,i)/LocalElementArea(0,i,0)
   BasisFunctions = UnmodifiedBasisFunctions(X,Y,number) - &
       CellStencils(i)%rIntBasisFunctions(number)/CellStencils(i)%Stencils(0)%rElementAreas(0)

  end function
 

  
  End module
