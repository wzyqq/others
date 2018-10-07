program read
    implicit none
    integer*8 :: i, j
    integer*8 :: np, ips
    integer*8,parameter:: ndm = 28991130 
    real*4,allocatable :: x(:,:)
    real*8 :: rand
    integer*8 :: select
    real*4 :: ztp, omegat, lambdat, rLbox, xscale, vscale
    
    integer,dimension(8) :: values1,values2,values
	character(8)  :: dates1,dates2,dates
    character(10) :: times1,times2,times
    character(5)  :: zones1,zones2,zones

	call date_and_time(dates1,times1,zones1,values1)

    open(unit=10,file="/data/s1/simu/Jing6620/pos6620.5000.01",form="unformatted",status="old")
    open(unit=20,file="/data/s1/simu/Jing6620/pos6620.5000.02",form="unformatted",status="old")
    open(unit=30,file="/data/s1/simu/Jing6620/pos6620.5000.03",form="unformatted",status="old")
    open(unit=40,file="/data/s1/simu/Jing6620/pos6620.5000.04",form="unformatted",status="old")
    open(unit=50,file="select_random",status="replace")
    open(unit=60,file="position",status="replace")
    read(10)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
    write(*,*)np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
    !np =20000
    allocate(x(3,np))
    read(10)((x(j,i),j=1,3),i=1,np/4)
    read(20)((x(j,i),j=1,3),i=1+np/4,np/2)
    read(30)((x(j,i),j=1,3),i=1+np/2,np/4*3)
    read(40)((x(j,i),j=1,3),i=1+np/4*3,np)
    write(*,*)"done"
!    do i = 1,np
!        write(*,"(I5F15.6F15.6F15.6)")i,rLbox*x(1,i),rLbox*x(2,i),rLbox*x(3,i)
!    end do

!    call random_seed()
!    do i = 1,ndm
!        call random_number(rand)
!        select = rand * np + 0.5!In fortran, the array index cannot be 0,maybe some mistake
!        write(50,"(I16)")select
!        write(60,"(F15.6F15.6F15.6)")x(1,select)*rLbox,x(2,select)*rLbox,x(3,select)*rLbox
!    end do
    deallocate(x)
    
    close(10)
    close(20)
    close(30)
    close(40)
    close(50)
    close(60)

    call date_and_time(dates2,times2,zones2,values2)
	write(*,"(2(A12))")times1 , times2
    
    stop
end
