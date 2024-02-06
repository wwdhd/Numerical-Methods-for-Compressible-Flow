
   ! Last update: 06/08/2009
   ! Higher-order version  
   ! Reference coordinates
   ! Degrees of freedom are numbered from 0 to MaxNumbDegrees
   ! Euler equations

   Use basis
   Use qr
   Use defi
   Use solver

   Implicit none
   
  ! %%%%%%%%%%%%%%%%%%  MAIN BODY OF THE CODE %%%%%%%%%%%%%%%%%%%%%%

  ! Commense the code
  Call INITALL
  CALL WriteInfo 
   

  print*,' call output'
  CALL OUTPUT

  ! commense time stepping
  print*
  Print*,' Commense time marching'

  ReconstructionTime = 0.
  FluxesTime = 0.
  LSQTime = 0.
  ExtrapTime = 0.
  UpdateTime = 0.

  call cpu_time(t0)

  do 
   Call ComputeTimeStep(dt)
   dt = min(dt,time-t_)  
!   print*,' dt=',dt
!   if (spatialorder .eq. 1) then
!     Call  FirstOrder
!   else     
     Call  ThirdOrderTVD
!   endif     
   t_=t_+dt
   it = it+1
   if ( mod(it, 10) .eq. 0)  then
     write(*,'(2x,a,i,3(a,1x,e11.4))')'it =',it,'    t =',t_,'    t/Time =',t_/time
   endif
   
   if ( mod(it, DataOutPutFreq) .eq. 0) then
     print*,'writing into output files'
     call   OUTPUT
   Endif     

   if ( mod(it, MovieOutPutFreq) .eq. 0) then
     Call   OutputTecplot
   Endif     
   
   if ( mod(it, 200) .eq. 0) then 
      CALL   WriteRestart
      Call   WriteInfo
   endif      

   if (abs(t_ -time)/time .le. 1d-5) then
     goto 333
   endif
 Enddo

 333 continue

 call cpu_time(t1)

 print*,' t = :',t_
 print*,' number of time steps:',it
 call  WriteInfo
 if (inicondtype .eq. 1)  CALL ComputeConvergenceError 
 print*,'writing into the file'
 call   OUTPUT
 Call   OutputTecplot
 print*,'done'
 
 ! close movie file
 CLOSE(111)


 TotalTime = t1-t0
 print*
 print*
 open(1,file='timereport.txt')
 write(*,'(a,f6.2,a)')'   total time to run:',Totaltime, 'sec'
 write(*,'(a,f6.2,a)')'   reconstruction time:',100*Reconstructiontime/TotalTime, ' %'
 write(*,'(a,f6.2,a)')'   fluxes time:',100*FluxesTime/TotalTime, ' %'
 write(*,'(a,f6.2,a)')'   update time:',100*UpdateTime/TotalTime, ' %'
 print*
 write(*,'(a,f6.2,a)')'   LSQ time:',100*LSQTime/Reconstructiontime, ' % of reconstruction'
 write(*,'(a,f6.2,a)')'   Extrapolation time:',100*ExtrapTime/Reconstructiontime, ' % of reconstruction'

 write(1,'(a,f6.2,a)')'   total time to run:',Totaltime, 'sec'
 write(1,'(a,f6.2,a)')'   reconstruction time:',100*Reconstructiontime/TotalTime, ' %'
 write(1,'(a,f6.2,a)')'   fluxes time:',FluxesTime/TotalTime, ' %'
 write(1,'(a,f6.2,a)')'   update time:',UpdateTime/TotalTime, ' %'
 print*
 write(1,'(a,f6.2,a)')'   LSQ time:',100*LSQTime/Reconstructiontime, ' % of reconsturction'
 write(1,'(a,f6.2,a)')'   Extrapolation time:',100*ExtrapTime/Reconstructiontime, ' % of reconsturction'
 close(1)

 
  End
	
  
	
	
