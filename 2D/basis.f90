  ! last release: 20/07/2008
  
  Module Basis

  Implicit None


  Contains

 ! *****************************************************************************
  Real(8) function UnmodifiedBasisFunctions(x,y,number)
   integer number,i
   real(8) x,y,s


   s = 0.
  
   select case(number)
    case(0)
     S = 1.D0
    case(1)
     S = x
    case(2)
     S = y
    case(3)
      S = x**2
    case(4)
      S = x*y
    case(5)
      S = y**2
   case(6)
       S = x**3
    case(7)
       S = x**2*y
    case(8)
       S = x*y**2
    case(9)
       S = y**3
    case(10)
        S = x**4
    case(11)
        S = x**3*y
    case(12)
        S = x**2*y**2
    case(13)
        S = x*y**3
    case(14)
        S = y**4
    case(15)
       S = x**5
    case(16)
        S = x**4*y
    case(17)
        S = x**3*y**2
    case(18)
        S = x**2*y**3
    case(19)
        S = x*y**4
    case(20)
        S = y**5
    case(21)
        S = x**6
    case(22)
        S = x**5*y
    case(23)
        S = x**4*y**2
    case(24)
        S = x**3*y**3
    case(25)
        S = x**2*y**4
    case(26)
        S = x*y**5
    case(27)
        S = y**6
    end select

    UnmodifiedBasisFunctions = S
  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsX(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(1)
      t0 = 1
     case(3)
      t0 = 2*x
     case(4)
      t0 = y
     case(6)
      t0 = 3*x**2
     case(7)
      t0 = 2*x*y
     case(8)
      t0 = y**2
     case(10)
      t0 = 4*x**3
     case(11)
      t0 = 3*x**2*y
     case(12)
      t0 = y**3
     case(15) 
      t0 = 5*x**4
     case(16) 
      t0 = 4*x**3*y
     case(17) 
      t0 = 3*x**2*y**2
     case(18) 
      t0 = 2*x*y**3
     case(19) 
      t0 = y**4
     case(21) 
      t0 = 6*x**5
     case(22) 
      t0 = 5*x**4*y
     case(23) 
      t0 = 4*x**3*y**2
     case(24) 
      t0 = 3*x**2*y**3
     case(25)
      t0 = 2*x*y**4
     case(26) 
      t0 = y**5
    End select

    BasisFunctionsX = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)

     case(2) 
      t0 = 1
     case(4) 
      t0 = x
     case(5) 
      t0 = 2*y
     case(7) 
      t0 = x**2
     case(8) 
      t0 = 2*x*y
     case(9) 
      t0 = 3*y**2
     case(11) 
      t0 = x**3
     case(12) 
      t0 = 3*x*y**2
     case(14) 
      t0 = 4*y**3
     case(16) 
      t0 = x**4
     case(17) 
      t0 = 2*x**3*y
     case(18) 
      t0 = 3*x**2*y**2
     case(19) 
      t0 = 4*x*y**3
     case(20) 
      t0 = 5*y**4
     case(22) 
      t0 = x**5
     case(23) 
      t0 = 2*x**4*y
     case(24) 
      t0 = 3*x**3*y**2
     case(25) 
      t0 = 4*x**2*y**3
     case(26) 
      t0 = 5*x*y**4
     case(27) 
      t0 = 6*y**5
    End select

    BasisFunctionsY= t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXX(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   Select case(number)
     case(  3  ) 
      t0 = 2
     case(  6  ) 
      t0 = 6*x
     case(  7  ) 
      t0 = 2*y
     case(  10  ) 
      t0 = 12*x**2
     case(  11  ) 
      t0 = 6*x*y
     case(  15  ) 
      t0 = 20*x**3
     case(  16  ) 
      t0 = 12*x**2*y
     case(  17  ) 
      t0 = 6*x*y**2
     case(  18  ) 
      t0 = 2*y**3
     case(  21  ) 
      t0 = 30*x**4
     case(  22  ) 
      t0 = 20*x**3*y
     case(  23  ) 
      t0 = 12*x**2*y**2
     case(  24  ) 
      t0 = 6*x*y**3
     case(  25  ) 
      t0 = 2*y**4
   End select

    BasisFunctionsXX = t0

  end function



 ! *****************************************************************************
  Real(8) function BasisFunctionsYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0

   select case(number)

    case(5) 
      t0 = 2
    case(8) 
      t0 = 2*x
    case(9) 
      t0 = 6*y
    case(12) 
      t0 = 6*x*y
    case(14) 
      t0 = 12*y**2
    case(17) 
      t0 = 2*x**3
    case(18) 
      t0 = 6*x**2*y
    case(19) 
      t0 = 12*x*y**2
    case(20) 
      t0 = 20*y**3
    case(23) 
      t0 = 2*x**4
    case(24) 
      t0 = 6*x**3*y
    case(25) 
      t0 = 12*x**2*y**2
    case(26) 
      t0 = 20*x*y**3
    case(27) 
      t0 = 30*y**4

   end select
   BasisFunctionsYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXY(x,y,number)
   integer number
   real(8) x,y,t0


   t0 = 0.

   select case(number)

    case( 4 ) 
      t0 = 1
    case( 7 ) 
      t0 = 2*x
    case( 8 ) 
      t0 = 2*y
    case( 11 ) 
      t0 = 3*x**2
    case( 12 ) 
      t0 = 3*y**2
    case( 16 ) 
      t0 = 4*x**3
    case( 17 ) 
      t0 = 6*x**2*y
    case( 18 ) 
      t0 = 6*x*y**2
    case( 19 ) 
      t0 = 4*y**3
    case( 22 ) 
      t0 = 5*x**4
    case( 23 ) 
      t0 = 8*x**3*y
    case( 24 ) 
      t0 = 9*x**2*y**2
    case( 25 ) 
      t0 = 8*x*y**3
    case( 26 ) 
      t0 = 5*y**4
    end select

   BasisFunctionsXY = t0

  end function

! %%%%%%%%%%%%%%%%
 ! *****************************************************************************
  Real(8) function BasisFunctionsXXX(x,y,number)
   integer number
   real(8) x,y,t0


   t0 = 0.

   select case(number)
    case( 6 ) 
      t0 = 6
    case( 10 ) 
      t0 = 24*x
    case( 11 ) 
      t0 = 6*y
    case( 15 ) 
      t0 = 60*x**2
    case( 16 ) 
      t0 = 24*x*y
    case( 17 ) 
      t0 = 6*y**2
    case( 21 ) 
      t0 = 120*x**3
    case( 22 ) 
      t0 = 60*x**2*y
    case( 23 ) 
      t0 = 24*x*y**2
    case( 24 ) 
      t0 = 6*y**3
   end select

    BasisFunctionsXXX = t0
  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXY(x,y,number)
   integer number
   real(8) x,y,t0


   t0 = 0.
   select case(number)
     case(7) 
      t0 = 2
     case(11) 
      t0 = 6*x
     case(16) 
      t0 = 12*x**2
     case(17) 
      t0 = 12*x*y
     case(18) 
      t0 = 6*y**2
     case(22) 
      t0 = 20*x**3
     case(23) 
      t0 = 24*x**2*y
     case(24) 
      t0 = 18*x*y**2
     case(25) 
      t0 = 8*y**3
   end select
   BasisFunctionsXXY = t0

  end function



 ! *****************************************************************************
  Real(8) function BasisFunctionsXYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0  = 0.
   select case(number)

     case( 8  ) 
      t0 = 2
     case( 12  ) 
      t0 = 6*y
     case( 17  ) 
      t0 = 6*x**2
     case( 18  ) 
      t0 = 12*x*y
     case( 19  ) 
      t0 = 12*y**2
     case( 23  ) 
      t0 = 8*x**3
     case( 24  ) 
      t0 = 18*x**2*y
     case( 25  ) 
       t0 = 24*x*y**2
     case( 26  ) 
      t0 = 20*y**3
   end select

   BasisFunctionsXYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(  9  ) 
      t0 = 6
     case(  12  ) 
      t0 = 6*x
     case(  14  ) 
      t0 = 24*y
     case(  18  ) 
      t0 = 6*x**2
     case(  19  ) 
      t0 = 24*x*y
     case(  20  ) 
      t0 = 60*y**2
     case(  24  ) 
      t0 = 6*x**3
     case(  25  ) 
      t0 = 24*x**2*y
     case(  26  ) 
      t0 = 60*x*y**2
     case(  27  ) 
      t0 = 120*y**3
   end select

   BasisFunctionsYYY = t0

  end function



 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXX(x,y,number)
   integer number
   real(8) x,y,t0

    t0 = 0.
    select case(number)
      case(10)
       t0 = 24
      case(15)
       t0 = 120*x
      case(16)
       t0 = 24*y
      case(21)
       t0 = 360*x**2
      case(22)
       t0 = 120*x*y
      case(23)
       t0 = 24*y**2
    End select

    BasisFunctionsXXXX = t0
  end function




 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.

   select case(number)
     case(16)
      t0 = 24*x
     case(17)
      t0 = 12*y
     case(22)
      t0 = 60*x**2
     case(23)
      t0 = 48*x*y
     case(24)
      t0 = 18*y**2
   end select

   BasisFunctionsXXXY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(17)
      t0 = 12*x
     case(18)
      t0 = 12*y
     case(23)
      t0 = 24*x**2
     case(24)
      t0 = 36*x*y
     case(25)
      t0 = 24*y**2
   end select

   BasisFunctionsXXYY = t0

  end function

 ! *****************************************************************************
  Real(8) function BasisFunctionsXYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(12)
      t0 = 6
    case(18)
      t0 = 12*x
    case(19)
      t0 = 24*y
    case(24)
      t0 = 18*x**2
    case(25)
      t0 = 48*x*y
    case(26)
      t0 = 60*y**2
   end select

   BasisFunctionsXYYY = t0

  end function

 ! *****************************************************************************
  Real(8) function BasisFunctionsYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(14)
      t0 = 24
    case(19)
      t0 = 24*x
    case(20)
      t0 = 120*y
    case(25)
      t0 = 24*x**2
    case(26)
      t0 = 120*x*y
    case(27)
      t0 = 360*y**2
  end select

   BasisFunctionsYYYY = t0

  end function


! for  6th order WENO

 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXXX(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(15)
      t0 = 120
    case(21)
      t0 = 720*x
    case(22)
      t0 = 120*y
   end select

   BasisFunctionsXXXXX = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXXY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(16)
      t0 = 24
     case(22)
      t0 = 120*x
     case(23)
      t0 = 48*y
   end select

   BasisFunctionsXXXXY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(17)
      t0 = 12
     case(23)
      t0 = 48*x
     case(24)
      t0 = 36*y
    end select

   BasisFunctionsXXXYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(18)
      t0 = 12
    case(24)
      t0 = 36*x
    case(25)
      t0 = 48*y
   end select

   BasisFunctionsXXYYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(19)
      t0 = 24
    case(25)
      t0 = 48*x
    case(26)
      t0 = 120*y
    end select
   BasisFunctionsXYYYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsYYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(20)
      t0 = 120
    case(26)
      t0 = 120*x
   case(27)
      t0 = 720*y
   end select

   BasisFunctionsYYYYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXXXX(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(21)
      t0 = 720
   end select

   BasisFunctionsXXXXXX = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXXXY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(22)
      t0 = 120
   end select

   BasisFunctionsXXXXXY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXXYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(23)
      t0 = 48
   end select

   BasisFunctionsXXXXYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXXXYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(24)
      t0 = 36
   end select

   BasisFunctionsXXXYYY = t0

  end function



 ! *****************************************************************************
  Real(8) function BasisFunctionsXXYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
    case(25)
      t0 = 48

   end select

   BasisFunctionsXXYYYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsXYYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(26)
      t0 = 120
   end select

   BasisFunctionsXYYYYY = t0

  end function


 ! *****************************************************************************
  Real(8) function BasisFunctionsYYYYYY(x,y,number)
   integer number
   real(8) x,y,t0

   t0 = 0.
   select case(number)
     case(27)
      t0 = 720
   end select

   BasisFunctionsYYYYYY = t0

  end function

  End module
