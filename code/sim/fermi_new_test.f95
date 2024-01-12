program feb_prog

  implicit none

  integer(kind=4) :: i, j,N,M,L,mc,mh,mw,clock,Nc_tot,g,k,H,nc_start,nw_start
  integer(kind=4) :: N_player, death,birth, victimizer, nh,nc,nw,game
  real(kind=4) :: w, alphac,alphaw,betah,betas,betac
  real(kind=4) ::  gamma,P, A,t, epsilon,delta,mu,P_tot,newP
  integer(kind=4) :: errflag,control,dt
  integer(kind=4), allocatable, dimension (:) :: player, role
  real(kind=4), allocatable, dimension (:) ::    payoff
  integer, allocatable, dimension (:) :: seed
  integer size


!------------------- 1=honest, 2=organized criminal, 3=lone wolf ---------------------------------!

open (unit = 13, file = "input")

read (13,*) nc_start
read (13,*) nw_start

read (13,*) alphac
read (13,*) alphaw
read (13,*) betas
read (13,*) betah
read (13,*) betac
read (13,*) gamma
read (13,*) epsilon
read (13,*) delta
read (13,*) mu


read (13,*) N
read (13,*) M
read (13,*) g
read (13,*) game

close(13)

P_tot=0.
newP=0.
dt=0

  allocate(player(N), stat=errflag)
  allocate(payoff(N), stat=errflag)
  allocate(role(N), stat=errflag)

  CALL RANDOM_SEED(size = size)
  ALLOCATE(seed(size))
  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  open(unit=12, file='output.dat')

do H=1,game,1

  dt=0

  nW=nw_start
  nC=nc_start
  nH=N-nW-nC
  role(1:nH)=1
  role(nH+1:nH+nC)=2
  role(nC+nH+1:N)=3

  do i=1,N,1
    player(i)=i
  end do
  payoff=0.

  nh=0
  nc=0
  nw=0
  do i=1,N,1
    if(role(i)==1) nh=nh+1
    if(role(i)==2) nc=nc+1
    if(role(i)==3) nw=nw+1
  end do
  Nc_tot=nc
  L=0
  print *, 'dt=', dt,H,L, 'H=',nh,'OC=', nc,'W=',nw


  write(12,*)  H, L,dt, nh,nc,nw



do while(dt<=50000) !L=1,30000,1
  payoff=0.
  do k=1,g,1 !-------------------------------------------------single round g game ---------!

!---------------randomization vector player --------------!
    do j=N,1,-1
      call random_number(A)
      t = player(1+FLOOR(j*A))
      player(1+FLOOR(j*A)) = player(j)
      player(j) = t
    end do

!----------------- make play all the sample --------------!
    do j=1,N/M,1  !------ do N/M times groups
     !--------- single group of M=10
      victimizer = j*M-M+1

      !---- select the 1' position as acting character, the other are victim----!


        if(role(player(victimizer))==3) then !------------------------wolf---!
          nh=0
          nc=0
          nw=0
          do i=2,M,1
            if(role(player(j*M-M+i))==1) nh=nh+1
            if(role(player(j*M-M+i))==2) nc=nc+1
            if(role(player(j*M-M+i))==3) nw=nw+1
          end do
          P=1.+delta*(nc/M-1)
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer))= payoff(player(victimizer))+(M-1.)*alphaw
            do i=2,M,1
              payoff(player(j*M-M+i)) =payoff(player(j*M-M+i)) -alphaw
            end do


          !--------------- investigation -----!
          !----- by the state ---------------------------------------------------!
          P = 1./M
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer)) =payoff(player(victimizer)) -betas
          end if
          !----- by the honest ---------------------------------------------------!
          P = 1.*nh/(M*M)
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer)) =payoff(player(victimizer)) -betah
          end if
        !----- by the criminal ----------------------------------------------------!
          P = 1.*nc/(M*M)
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer)) =payoff(player(victimizer)) -betac
          end if
        end if
        end if


        if(role(player(victimizer))==2) then    !-----------criminal ------!
          nh=0
          nc=1
          nw=0
          do i=2,M,1
            if(role(player(j*M-M+i))==1) then
               nh=nh+1
               payoff(player(j*M-M+i)) =payoff(player(j*M-M+i)) -alphac
            end if
            if(role(player(j*M-M+i))==2) nc=nc+1
            if(role(player(j*M-M+i))==3) then
              nw=nw+1
              payoff(player(j*M-M+i)) =payoff(player(j*M-M+i)) -alphac
            end if
          end do
          payoff(player(victimizer)) =payoff(player(victimizer))+ 1.*(M-nc)*alphac/nc
          do i=2,M,1
            if(role(player(j*M-M+i))==2) payoff(player(j*M-M+i))=payoff(player(j*M-M+i))+1.*(M-nc)*alphac/nc
          end do
          !--------------- investigation -----!
          !---------------- by the state ------------------------------------!
          !---------------- by the state ------------------------------------!
          P = 1.*nc/M
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer)) =payoff(player(victimizer)) -betas
            do i=2,M,1
               if(role(player(j*M-M+i))==2) payoff(player(j*M-M+i))=payoff(player(j*M-M+i))- betas*gamma
            end do
          end if
          !----- by the honest ---------------------------------------------------!
          P = 1.*nc*nh/(M*M)
          call random_number(A)
          if(A<P) then
            payoff(player(victimizer)) =payoff(player(victimizer)) -betah
            do i=2,M,1
               if(role(player(j*M-M+i))==2) payoff(player(j*M-M+i))=payoff(player(j*M-M+i))- gamma*betah
            end do
          end if
        end if
      end do
  end do !----------------------------------------------------end of part 1-2 ------------------!


  payoff=payoff/(1.*g)
  control=0
  do while(control==0)
    dt=dt+1
!--------------- evolution -----------------------------------------------------------------------!
    call random_number(A)
    death =1+FLOOR(N*A)
  !----------------------mutation----------------------------!
    call random_number(A)
    if(A<mu) then
      call random_number(A)
      role(player(death)) = 1+FLOOR(3*A)
      control=control+1
   else
      call random_number(A)
      birth =1+FLOOR(N*A)
      if(role(player(death))/=role(player(birth))) then
        call random_number(A)
        P= 1./(1.+exp(-epsilon*(payoff(player(birth))-payoff(player(death)))))
        if(A<P) then
           role(player(death)) = role(player(birth))
           control=control+1
         end if
      end if
    end if

  end do
!--------------------------------------------end DB process ------------------------------------------------!


  !------------------- write on file --------------------------------------------!
      nh=0
      nc=0
      nw=0
      do i=1,N,1
        if(role(i)==1) nh=nh+1
        if(role(i)==2) nc=nc+1
        if(role(i)==3) nw=nw+1
      end do
      Nc_tot=nc
      write(12,*) H, L, dt,nh,nc,nw

end do ! end of game-------------------------------------------------------------------------------------------!
  print *, dt, H, L, 'H=',nh,'OC=', nc,'W=',nw
end do


!------------------- write on file --------------------------------------------!
!    nh=0
!    nc=0
!    nw=0
!    do i=1,N,1
!      if(role(i)==1) nh=nh+1
!      if(role(i)==2) nc=nc+1
!      if(role(i)==3) nw=nw+1
!    end do
!    Nc_tot=nc
!    write(12,*) L, nh,nc,nw

!    print *, j, 'H=',nh,'OC=', nc,'W=',nw, role(victimizer),Nc_tot, criminal_payoff
  close(12)



  deallocate(player,payoff,role, seed)




end program
