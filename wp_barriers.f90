program wp_prop

implicit none
logical                 ::      pbc
integer                 ::      i, step, step_print
integer, parameter      ::      un_input=101, un_out=102
integer                 ::      Nv, xPoints, Nsteps,Nprint, count_v=0
real(8)                 ::      a, k , x0, dt, t_init, dx, xInit, xFin, m, hBarra, w
real(8)                 ::      pi, dt_crit, vmax, norma
real(8)                 ::      start, finish
real(8)                 ::      init, fin, barrera
real(8),allocatable,dimension(:)        ::      psiR, psiI, x, v
logical,allocatable,dimension(:)        ::      vTypes

call cpu_time(start)
pi = dacos(-1.d0)
call open_input()
call read_parameters()

!! Check that input file is valid
do i=1,size(vTypes)
        if (vTypes(i).eqv..TRUE.) then
                count_v = count_v +1
        endif
enddo

if (count_v.gt.1) then
        print*, 'More than one potential selected. Exitting program...'
        call exit()
endif

if (count_v .eq. 0) then
        print*, 'No potential selected. Exitting program...'
        call exit()
endif

!!Allocate arrays
allocate(psiR(xPoints),psiI(xPoints), v(xPoints), x(xPoints))

!!Build potential
v(:) = 0.d0
! Check the potential
if (vtypes(1).eqv..FALSE.) then
        if (vtypes(2).eqv..TRUE.) then
                v = v_barrera()
        endif
endif
vmax = maxval(v)
dx = (xFin - xInit)/(xPoints-1)
dt_crit=(2*m*hBarra*dx**2)/(hBarra**2 + 2.d0*m*dx**2*vmax)

! Check valid dt
if (dt .ge. dt_crit ) then
        print*,'dt=',dt,'will result in an inestable solution'
        print*, 'dt must be smaller than',dt_crit
        print*,'Exitting program'
        call exit()
endif

! Build x array
do i=1, xPoints
        x(i) = xInit + (i-1)*dx
enddo

! Build initial WP
call build_wp()

! Open output file
open(unit=un_out,file='wp_trajectory.log',status='new')

! Write initial WP
call print_squared()

!Starting loop propagation
do step=1, int(NSteps/Nprint)
        do step_print=1,Nprint
                call euler()
        enddo
        call print_squared()
enddo

! Check final Norm
print*, 'Final norm:',norm()
close(un_out)

call cpu_time(finish)

print*,'CPU time:', finish - start

contains

subroutine open_input()
     
character(24)           ::      fName
integer                 ::      fStat

call get_command_argument(1,fName, status=fStat)

if (fStat /= 0) then
        print*,'Failed at reading input file. Exitting program...'
        call exit()
endif

open(unit=un_input,file=trim(fName), status='old')

end subroutine open_input


subroutine read_parameters()

integer ::      i

read(un_input,*) xInit
read(un_input,*) xFin
read(un_input,*) xPoints
read(un_input,*) x0
read(un_input,*) dt
read(un_input,*) t_init
read(un_input,*) Nsteps
read(un_input,*) Nprint
read(un_input,*) a
read(un_input,*) k
read(un_input,*) w
read(un_input,*) m
read(un_input,*) hBarra
read(un_input,*) Nv
allocate(vTypes(Nv))

do i=1,Nv
        read(un_input,*) vTypes(i)
enddo

read(un_input,*) pbc
read(un_input,*) init
read(un_input,*) fin
read(un_input,*) barrera
close(un_input)

end subroutine read_parameters


subroutine build_wp()

psiR= (a**2/(2*pi))**0.25d0*dexp(-a**2*(x(:)-x0)**2/4.d0)*dcos(k*x(:) - w*t_init)
psiI= (a**2/(2*pi))**0.25d0*dexp(-a**2*(x(:)-x0)**2/4.d0)*dsin(k*x(:) - w*t_init)

end subroutine build_wp


function v_barrera()

integer         ::      init_point, fin_point
real(8),dimension(xpoints)       ::      v_barrera

init_point = int((init-xinit)/dx)
fin_point = int((fin-xinit)/dx)
v_barrera(init_point:fin_point) = barrera

end function

subroutine print_descomp()

do i=1,xPoints
        write(un_out,*) x(i), psiR(i), psiI(i), psiR(i)**2+psiI(i)**2, v(i)
enddo
write(un_out,*)
write(un_out,*)

end subroutine print_descomp


subroutine print_squared()

do i=1,xPoints
        write(un_out,*) x(i), psiR(i)**2+psiI(i)**2, v(i)
enddo
write(un_out,*)
write(un_out,*)

end subroutine print_squared


function norm()

real(8)         ::      norm

norm = dx*sum(psir**2+psiI**2)

end function norm


subroutine euler()

real(8),dimension(xPoints)        ::      pasR, pasI

pasR(2:xPoints-1) = dt*(v(2:xPoints-1)*psiI(2:xPoints-1)/hBarra - &
(hBarra/(2.d0*m*dx**2))*(psiI(3:xPoints)-2*psiI(2:xPoints-1)+psiI(1:xPoints-2)))

pasI(2:xPoints-1) = -dt*(v(2:xPoints-1)*psiR(2:xPoints-1)/hBarra - &
(hBarra/(2.d0*m*dx**2))*(psiR(3:xPoints)-2*psiR(2:xPoints-1)+psiR(1:xPoints-2)))

!! Apply or not periodic boundary conditions
if (pbc.eqv..FALSE.) then

        pasR(1)=0.d0; pasR(xPoints)=0.d0
        pasI(1)=0.d0; pasI(xPoints)=0.d0

else

        pasR(1) = dt*(v(1)*psiI(1)/hBarra - &
        (hBarra/(2.d0*m*dx**2))*(psiI(2)-2*psiI(1)+psiI(xPoints)))

        pasI(1) = -dt*(v(1)*psiR(1)/hBarra - &
        (hBarra/(2.d0*m*dx**2))*(psiR(2)-2*psiR(1)+psiR(xPoints)))
        
        pasR(xPoints) = dt*(v(xPoints)*psiI(xPoints)/hBarra - &
        (hBarra/(2.d0*m*dx**2))*(psiI(1)-2*psiI(xPoints)+psiI(xPoints-1)))

        pasI(xPoints) = -dt*(v(xPoints)*psiR(xPoints)/hBarra - &
        (hBarra/(2.d0*m*dx**2))*(psiR(1)-2*psiR(xPoints)+psiR(xPoints-1)))

endif

psiR = psiR + pasR
psiI = psiI + pasI

end subroutine euler



end program wp_prop
