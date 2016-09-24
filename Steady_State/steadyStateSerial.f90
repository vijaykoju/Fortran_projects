!Author: *Filmon H. Gebreyesus
!	 *Vijay Koju
!COMS6100-Assignment #12 (Steady State Temperature Distribution)
!December 06, 2012


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program steadyStateSerial
implicit none
integer, parameter:: sgle = selected_real_kind(6,37)
integer::t,n=768,m=768,i,j,iteration=0, tempabove50=0
real (kind = sgle), dimension(768,768):: tmp,new_tmp
real (kind = sgle):: tol
integer (kind = sgle):: count, count_rate, count_max
!character (len = 4):: fileNum
!character (len = 25):: fileName

call system_clock(count,count_rate,count_max)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                    initial and fixed conditions                        !!!!!   
do j = 1,m
	do i = 1,n
		tmp(i,j) = 50 ! initial all the grid to 50 degrees
	end do
end do
new_tmp = tmp

tmp(200,500) = 100 ! temperature at grid (200,500)

do i = 1,n
	tmp(i,1) = 0 ! temperature at the bottom boundary
	tmp(i,m) = 0 ! temperature at the top boundary
end do

do j = 1,m
	tmp(n,j) = 0 ! temperature at the right boundary
	tmp(1,j) = 100 ! temperature at the left boundary
end do

do j = 1,331
	tmp(400,j) = 100 ! temperature at row 400 and column 1 to 331
end do
!!!!!!!!! write the output in matrix form !!!!!!!!!!
!open(unit=3, file = "steadyState_1000.txt")
!write(3,*)
!do i = 1,n
!	do j = 1,m
!		write(3,'(f10.3$)') tmp(j,i) !tmp(101,501), tmp(7,1), tmp(324, 523)
!	end do
!	write(3,*)
!end do
!write(3,*)
!close(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tol = 1000 ! intial tolerance set very high
do while (tol > 0.01) !check the tolerance before going to next iteration (for updating temperature)
	do j = 2,m-1
		do i = 2,n-1
			new_tmp(i,j) = (tmp(i+1,j)+tmp(i-1,j)+tmp(i,j+1)+tmp(i,j-1)+4*tmp(i,j))/8.0 
		end do
	end do
!Setting the initial fixed temperatures
	new_tmp(200,500) = 100
	do i = 1,n
		new_tmp(i,1) = 0
		new_tmp(i,m) = 0
	end do
	do j = 1,m
		new_tmp(n,j) = 0
		new_tmp(1,j) = 100
	end do
	do j = 1,331
		new_tmp(400,j) = 100
	end do
!UNCOMMENT FOR SAVING OUTPUT TO FILES
!	write(fileNum,'(i4)') iteration
!	fileName = "steadyState_"//fileNum//".txt"
!	open(unit=3, file = fileName)
!	write(3,*)
!	do i = 1,n
!		do j = 1,m
!			write(3,'(f10.3$)') new_tmp(j,i) !tmp(101,501), tmp(7,1), tmp(324, 523)
!		end do
!		write(3,*)
!	end do
!	write(3,*)
!	close(3)

	tol = maxval(abs(new_tmp-tmp))
	tmp = new_tmp
	iteration = iteration+1 !Counting the number of iterations
	
end do
!Counting the cells having temperature above 50 degrees
do j = 1,m
	do i = 1,n
		if (tmp(i,j) > 50) then
			tempabove50 = tempabove50+1
		end if
	end do
end do
write(*,*), "Number of iterations taken to come to a steady state: ", iteration
write(*,*), "Number of cells with temperature above 50 degrees at the steady state: ",tempabove50

end program steadyStateSerial
