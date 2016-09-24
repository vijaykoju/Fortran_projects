!Author: *Filmon H. Gebreyesus
!	 *Vijay Koju
!COMS6100-Assignment #12 (Steady State Temperature Distribution)
!December 06, 2012


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program steadyStateParallel
implicit none
!Declaration for variables
integer, parameter:: sgle = selected_real_kind(6,37)
integer::t,n=768,m=768,i,j,iteration=0,tempabove50 = 0
real (kind = sgle), dimension(768,768):: tmp,new_tmp
integer (kind = sgle):: count, count_rate, count_max
real (kind = sgle):: tol
!character (len = 4):: fileNum
!character (len = 25):: fileName

call system_clock(count,count_rate,count_max)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                    initial and fixed conditions                        !!!!!   
!$OMP PARALLEL DO  !Parallelize
do j = 1,m 
	do i = 1,n
		tmp(i,j) = 50 ! initial all the grid to 50 degrees
	end do
end do
!$OMP END PARALLEL DO
new_tmp = tmp


tmp(200,500) = 100 ! temperature at grid (200,500)

!$OMP PARALLEL DO
do i = 1,n
	tmp(i,1) = 0 ! temperature at the bottom boundary
	tmp(i,m) = 0 ! temperature at the top boundary
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do j = 1,m
	tmp(n,j) = 0 ! temperature at the right boundary
	tmp(1,j) = 100 ! temperature at the left boundary
end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do j = 1,331
	tmp(400,j) = 100 ! temperature at row 400 and column 1 to 331
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! write the output in matrix form !!!!!!!!!!
!open(unit=3, file = "steadyState_1000.dat")
!write(3,*)
!do i = 1,n
!	do j = 1,m
!		write(3,'(f10.3$)') tmp(j,i)
!	end do
!	write(3,*)
!end do
!write(3,*)
!close(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tol = 1000 ! intial tolerance set very high
do while (tol > 0.01) !check the tolerance before going to next iteration (for updating temperature)
	iteration = iteration+1 !Counting the number of iterations
	!$OMP PARALLEL DO
	do j = 2,m-1
		do i = 2,n-1
			new_tmp(i,j) = (tmp(i+1,j)+tmp(i-1,j)+tmp(i,j+1)+tmp(i,j-1)+4*tmp(i,j))/8.0
		end do
	end do
!Setting the initial fixed temperatures
	!$OMP END PARALLEL DO
	new_tmp(200,500) = 100
	!$OMP PARALLEL DO
	do i = 1,n
		new_tmp(i,1) = 0
		new_tmp(i,m) = 0
	end do
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO
	do j = 1,m
		new_tmp(n,j) = 0
		new_tmp(1,j) = 100
	end do
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO
	do j = 1,331
		new_tmp(400,j) = 100
	end do
	!$OMP END PARALLEL DO
!UNCOMMENT FOR SAVING OUTPUT TO FILES
!	if (mod(iteration,3)==0) then
!		write(fileNum,'(i4)') iteration
!		fileName = "steadyState_"//fileNum//".dat"
!		open(unit=3, file = fileName)
!		write(3,*)
!		do i = 1,n
!			do j = 1,m
!				write(3,'(f10.3$)') new_tmp(j,i) 
!			end do
!			write(3,*)
!		end do
!		write(3,*)
!		close(3)
!	end if

	tol = maxval(abs(new_tmp-tmp)) !! Calcualte tolerance 
	tmp = new_tmp !set previous temperature to new temperature for updation 
	
end do
!Counting the cells having temperature above 50 degrees

!$OMP PARALLEL DO REDUCTION(+:tempabove50)  !tempabove50 is shared between the parallelized threads for updation
do j = 1,m
	do i = 1,n
		if (tmp(i,j) > 50) then
			tempabove50 = tempabove50+1
		end if
	end do
end do
!$OMP END PARALLEL DO
write(*,*), "Number of iterations taken to come to a steady state: ", iteration
write(*,*), "Number of cells with temperature above 50 degrees at the steady state: ",tempabove50


end program steadyStateParallel
