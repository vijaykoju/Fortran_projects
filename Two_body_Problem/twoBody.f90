program twoBody
implicit none

integer, parameter :: dk = selected_real_kind(15,307)
integer :: numOfObj ! number of objects
integer :: AllocateStatus, DeAllocateStatus
real (kind=dk), dimension(:,:), allocatable :: ObjProp ! properties of the objects
real (kind=dk), dimension(:), allocatable :: x, y, vx, vy, m, fx, fy, pe, ke, time
real (kind=dk), dimension(:,:), allocatable :: pv1, xp, yp, energy
real (kind=dk) :: dx, dy, dvx, dvy,ax, ay, drSqrd, dr, df, dfx, dfy, peTot, keTot
integer :: tFinal, nSteps
real :: dt, ee
integer :: i, j, k, n ! looping variables
integer ::  A ! variable for command line argument parsing
character (len=10) :: option1 ! variable for command line argument parsing

A = iargc() ! read input from the command line
if (A<1) then ! if no argument given, print the info and stop the program
  write(*,*)
  write(*,'(a)') "Please provide the number of objects you would like to simulate."
  write(*,'(a)') "usage:   ./twoBody option1"
  write(*,'(a)') "         option1 = integer"
  write(*,*)
  stop  
else
  call getarg(1,option1) ! grab the first command line argument and store it in option1
  read(option1,*) numOfObj ! convert string(option1) to integer(numOfObj)
end if

!numOfObj = 10
allocate(ObjProp(numOfObj,5),x(numOfObj),y(numOfObj),vx(numOfObj),vy(numOfObj),STAT=AllocateStatus)
allocate(m(numOfObj),fx(numOfObj),fy(numOfObj),pe(numOfObj),ke(numOfObj), STAT=AllocateStatus) ! memory allocation
if (AllocateStatus/=0) stop ! not enough memory to allocate

ObjProp = dataMatrix("finalData.dat",numOfObj) ! object property matrix
x = ObjProp(1:numOfObj,1) ! displacement in x-direction (m)
y = ObjProp(1:numOfObj,2) ! displacement in y-direction (m)
vx = ObjProp(1:numOfObj,3) ! velocity in x-direction (m/s)
vy = ObjProp(1:numOfObj,4) ! velocity in y-direction (m/s)
m = ObjProp(1:numOfObj,5) ! mass 

! set the final times
tFinal = 1
dt = 0.001
nSteps = tFinal/dt
ee = 0.05 ! softening length used to prevent close approach problems

allocate(pv1(nsteps,6),xp(nsteps,numOfObj),yp(nsteps,numOfObj),time(nsteps),energy(nsteps,3), STAT=AllocateStatus)
if (AllocateStatus/=0) stop
do k = 1,nsteps ! loop over the time range

  do n = 1,numOfObj ! loop over for each source object
    fx(n) = 0 ! initialize force in x-direction to 0 (N)
    fy(n) = 0 ! initialize force in y-direction to 0 (N)
    pe(n) = 0 ! initialize potential energy (J)
    ke(n) = 0 ! initialize kinetic energy (J)
  end do

  do i = 1,numOfObj-1 ! loop over for each source object
    do j = i+1,numOfObj ! loop over the other objects
      dx = x(i) - x(j) ! small change in x-displacement (m)
      dy = y(i) - y(j) ! small change in y-displacement (m)
      drSqrd = dx*dx+dy*dy+ee*ee ! sq. distance between source and other object (m^2)
      dr = sqrt(drSqrd) ! distance between source and other object (m)
      df = m(i)*m(j)/drSqrd ! total force on source due to other object (N)
      dfx = df*dx/dr ! x-component of the total force (N)
      dfy = df*dy/dr ! y-component of the total force (N)

      ! update the target and source
      fx(i) = fx(i) - dfx ! (N)
      fy(i) = fy(i) - dfy ! (N)
      fx(j) = fx(j) + dfx ! (N)
      fy(j) = fy(j) + dfy ! (N)

      ! update potential energy
      pe(i) = pe(i) - m(i)*m(j)/dr/2 ! (J)
      pe(j) = pe(j) - m(i)*m(j)/dr/2 ! (J)
    end do
  end do

  do i = 1,numOfObj
    ke(i) = 0.5*m(i)*(vx(i)*vx(i) + vy(i)*vy(i)) ! calculate kinetic energy (J)
    ax = fx(i)/m(i) ! acceleration in x-direction (m/s^2)
    ay = fy(i)/m(i) ! acceleration in y-direction (m/s^2)
    dvx = ax*dt ! small change in velocity in x-direction (m/s)
    dvy = ay*dt ! small change in velocity in y-direction (m/s)
    
    ! update the velocities
    vx(i) = vx(i) + dvx ! (m/s)
    vy(i) = vy(i) + dvy ! (m/s)

    ! change in positions
    dx = vx(i)*dt ! (m)
    dy = vy(i)*dt ! (m)

    ! update the position
    x(i) = x(i) + dx ! (m)
    y(i) = y(i) + dy ! (m)
  end do
  
  pv1(k,1) = x(1)
  pv1(k,2) = y(1)
  pv1(k,3) = vx(1)
  pv1(k,4) = vy(1)
  pv1(k,5) = fx(1)
  pv1(k,6) = fy(1)
  
  xp(k,:) = x(:)
  yp(k,:) = y(:)

  peTot = sum(pe)
  keTot = sum(ke)

  time(k) = k*dt
  energy(k,1) = peTot
  energy(k,2) = keTot
  energy(k,3) = peTot+keTot
end do
!print*, peTot
open(unit=2, file = "positionData.dat")
do k = 1,nsteps
  do i = 1,numOfObj
    write(2,'(f10.6,f10.6$)') xp(k,i), yp(k,i)
  end do
  write(2,*)
end do
close(3)

open(unit=3, file = "energyData.dat")
do k = 1,nsteps
  write(3,*) k, energy(k,1), energy(k,2), energy(k,3)
end do
close(3)
deallocate(ObjProp,x,y,vx,vy,m,fx,fy,pe,ke,pv1,xp,yp,time,energy, STAT=DeAllocateStatus)
!call system("gnuplot -persist gnu.gp")
!call system("rm positionData.dat energyData.dat")

contains
! function for reading data from a file
function dataMatrix(fileName,rows)
  integer, parameter :: dk = selected_real_kind(15, 307)
  character(*)::fileName
  integer:: lines, rows
  real (kind=dk), dimension(rows,5):: dataMatrix

  open(unit=3, file=fileName) ! open the file for reading
  do lines = 1,rows ! iterate through each line
    read(3,*,end=11) dataMatrix(lines,1), dataMatrix(lines,2), dataMatrix(lines,3), dataMatrix(lines,4), dataMatrix(lines,5)
  end do
  11 close(3) ! close the file
end function dataMatrix

end program twoBody
