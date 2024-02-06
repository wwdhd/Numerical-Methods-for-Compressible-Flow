
 Module QR

 Implicit None

  Contains
 
  real(8) function l2norm(x,n)
   integer, intent(in) :: n
   real(8) x(n)
   integer i
 
    l2norm = 0.d0
    do i=1,n
     l2norm = l2norm + x(i)**2
    enddo
    l2norm = dsqrt(l2norm)
  end function
  
  
  ! construct a householder vector v that 
  ! annihilates all but the first component of x
  Subroutine house(x,v,n)
    integer, intent(in) :: n
    real(8) x(n),v(n)

    v = x
    v(1) = x(1) + sign(1.d0,x(1))*l2norm(x,n)
    
  End subroutine


  ! construct a householder reflection matrix
  ! from a Householder vector v 
  Subroutine ComputeHouseMatrix(P,v,n)
    integer, intent(in) :: n
    real(8) P(n,n),v(n),vnorm
    integer i,j

    P = 0.D0
    do i=1,n
      P(i,i) = 1.d0
    enddo
   
    vnorm = l2norm(v,N)
    v = v/vnorm
    
    do i=1,n
    do j=1,n
      p(i,j) = p(i,j) - 2*v(i)*v(j)
    enddo
    enddo
    
  End subroutine

 !%%%%%%%%%%%%%%%%%%
  Subroutine TransposeMatrix(A,AT,N)
    integer, intent(in) :: n
    real(8), intent(in) :: A(N,N)
    real(8), intent(out) ::  AT(N,N)
    integer i,j
    do i=1,n
    do j=1,n
      at (i,j) = a(j,i)
    enddo
    enddo
  End subroutine


  ! QR decomposition
!  Subroutine  QRDecomposition(A,Q0,R,N)
  Subroutine  QRDecomposition(A,Q,R,InvertedQR,N)
    integer, intent(in) :: n
    real(8), intent(in) :: A(N,N)
    real(8), intent(out) ::  Q(N,N),R(N,N),InvertedQR(N,N)
    real(8)  Identity(N,N),V(N),P(N,N),InvR(N,N),QT(N,N)
!    real(8) test(N)
    integer i,l,pdim
    real(8), allocatable :: v1(:),x(:)
 
    Identity = 0.D0
    do i=1,n
      Identity(i,i) = 1.d0
    enddo

    Q = Identity
    R = A


    do l=1,n
     ! allocate vector and reflection matrix
     pdim = n-l+1
     allocate(v1(pdim),x(pdim)) 
     x = R(L:N,L)
     ! compute the partial vector    
     CALL House(x,v1,pdim)
     v = 0.
     v(L:N) = v1
     ! compute the reflection matrix
     CALL ComputeHouseMatrix(P,v,N)
     ! construct the Q(l) matrix
     R = MATMUL(P,R)
     Q = MATMUL(Q,P)
     deallocate(x,v1)
    enddo


    101 format( 15(2x,E9.2))
!     print*,' final matrix Q'
!     do i=1,n
!       write(*,101) Q(i,:)
!     enddo
!     print*,' final matrix R'
!     do i=1,n
!       write(*,101) R(i,:)
!     enddo
!	 read*


    ! just in case - make 
!    do i=1,N
!  	 do j=1,i-1
!	   R(i,j) =0.
!	 enddo
!   enddo

!     print*,' final matrix R'
!     do i=1,n
!       write(*,101) R(i,:)
!     enddo
!	 read*

!     A = A - MATMUL(Q,R)
!     print*,' Error in decomposition:'
!     do i=1,n
!       write(*,101) A(i,:)
!     enddo

     CALL Invert(R,InvR,N)

      ! transpose Q
     CALL TransposeMatrix(Q,QT,N)

!   final inverted R^(-1)*Q^(-1)

   InvertedQR = matmul(InvR,QT)

 
  End subroutine

 
 Subroutine Invert(R,InvR,N)
  implicit real (t)
  integer N
  real(8), intent(in) :: R(N,N)
  real(8), intent(out):: InvR(N,N)
!  Real(8) InvR2(N,N)

  integer i,j,k

 
  InvR = 0.d0


  do i=N,1,-1
    invr(i,i) = 1./r(i,i)
    do j=i+1,N
	  invr(i,j) = 0.
	  do k= 1,j-1
	   invr(i,j) = invr(i,j) - r(k,j)*invr(i,k)
	  enddo
 	  invr(i,j) =invr(i,j) /r(j,j)
	enddo
  enddo

!  InvR2= matmul(R,Invr)
!
!  101 format( 15(2x,F6.2))
!  print*,' final matrix InvR2'
!  do i=1,n
!       write(*,101) InvR2(i,:)
!  enddo
!  read*

  End subroutine

 End module