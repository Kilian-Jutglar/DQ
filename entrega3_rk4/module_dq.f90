subroutine build_wp(a, x0, k, w, xpoints, x, psi_r, psi_i)
implicit none
integer                       ::      i
integer                       ::      xpoints
real(8)                       ::      pi
real(8)                       ::      a, x0, k, w
real(8),dimension(xpoints)    ::      x
real(8),dimension(xpoints)    ::      psi_r, psi_i

pi = dacos(-1.d0)
do i = 1, xpoints, 1
        psi_r(i) = (a**2/(2.d0*pi))**0.25d0*dexp(-a**2*(x(i)-x0)**2/4.d0)*dcos(k*x(i))
        psi_i(i) = (a**2/(2.d0*pi))**0.25d0*dexp(-a**2*(x(i)-x0)**2/4.d0)*dsin(k*x(i))
enddo
end subroutine build_wp

subroutine euler(dt,dx,m,h_barra,xpoints,v,psi_r,psi_i,steps_print)
implicit none
integer                           ::      i, steps_print
integer                           ::      xpoints
real(8)                           ::      dt, dx, m, h_barra
real(8),dimension(xpoints)        ::      v
real(8),dimension(xpoints)        ::      psi_r, psi_i
real(8),dimension(xpoints)        ::      pas_r, pas_i

do i=1,steps_print
        pas_r(2:xpoints-1) = dt*(v(2:xpoints-1)*psi_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(3:xpoints)-2*psi_i(2:xpoints-1)+psi_i(1:xpoints-2)))
        
        pas_i(2:xpoints-1) = -dt*(v(2:xpoints-1)*psi_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(3:xpoints)-2*psi_r(2:xpoints-1)+psi_r(1:xpoints-2)))
        
        pas_r(1)=0.d0; pas_r(xpoints)=0.d0
        pas_i(1)=0.d0; pas_i(xpoints)=0.d0
        
        psi_r = psi_r + pas_r
        psi_i = psi_i + pas_i
enddo
end subroutine euler

subroutine euler_pbc(dt,dx,m,h_barra,xpoints,v,psi_r,psi_i,steps_print)
implicit none
integer                           ::      xpoints, i, steps_print
real(8)                           ::      dt, dx, m, h_barra
real(8),dimension(xpoints)        ::      v
real(8),dimension(xpoints)        ::      psi_r, psi_i
real(8),dimension(xpoints)        ::      pas_r, pas_i

do i=1,steps_print

        pas_r(2:xpoints-1) = dt*(v(2:xpoints-1)*psi_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(3:xpoints)-2*psi_i(2:xpoints-1)+psi_i(1:xpoints-2)))
        
        pas_i(2:xpoints-1) = -dt*(v(2:xpoints-1)*psi_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(3:xpoints)-2*psi_r(2:xpoints-1)+psi_r(1:xpoints-2)))
        
        
        pas_r(1) = dt*(v(1)*psi_i(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(2)-2*psi_i(1)+psi_i(xpoints)))
        
        pas_i(1) = -dt*(v(1)*psi_r(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(2)-2*psi_r(1)+psi_r(xpoints)))
        
        pas_r(xpoints) = dt*(v(xpoints)*psi_i(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(1)-2*psi_i(xpoints)+psi_i(xpoints-1)))
        
        pas_i(xpoints) = -dt*(v(xpoints)*psi_r(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(1)-2*psi_r(xpoints)+psi_r(xpoints-1)))
        
        
        psi_r = psi_r + pas_r
        psi_i = psi_i + pas_i
enddo
end subroutine euler_pbc

subroutine rk4(dt,dx,m,h_barra,xpoints,v,psi_r,psi_i,steps_print)
implicit none
integer                           ::      xpoints, i, steps_print
real(8)                           ::      dt, dx, m, h_barra
real(8),dimension(xpoints)        ::      v
real(8),dimension(xpoints)        ::      psi_r, psi_i
real(8),dimension(xpoints)        ::      k1_r, k1_i, k2_r, k2_i, k3_r, k3_i, k4_r, k4_i
real(8),dimension(xpoints)        ::      psi2_r, psi2_i, psi3_r, psi3_i, psi4_r, psi4_i

do i=1, steps_print
        k1_r(2:xpoints-1) = (v(2:xpoints-1)*psi_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(3:xpoints)-2*psi_i(2:xpoints-1)+psi_i(1:xpoints-2)))
        
        k1_i(2:xpoints-1) = -(v(2:xpoints-1)*psi_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(3:xpoints)-2*psi_r(2:xpoints-1)+psi_r(1:xpoints-2)))
        
        
        k1_r(1) = 0.d0
        k1_i(1) = 0.d0
        k1_r(xpoints) = 0.d0
        k1_i(xpoints) = 0.d0
        
        psi2_r=psi_r+dt*k1_r/2.d0
        psi2_i=psi_i+dt*k1_i/2.d0
        
        k2_r(2:xpoints-1) = (v(2:xpoints-1)*psi2_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_i(3:xpoints)-2*psi2_i(2:xpoints-1)+psi2_i(1:xpoints-2)))
        
        k2_i(2:xpoints-1) = -(v(2:xpoints-1)*psi2_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_r(3:xpoints)-2*psi2_r(2:xpoints-1)+psi2_r(1:xpoints-2)))
        
        k2_r(1) = 0.d0
        k2_i(1) = 0.d0
        k2_r(xpoints) = 0.d0
        k2_i(xpoints) = 0.d0
        
        psi3_r=psi_r+dt*k2_r/2.d0
        psi3_i=psi_i+dt*k2_i/2.d0
        
        k3_r(2:xpoints-1) = (v(2:xpoints-1)*psi3_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_i(3:xpoints)-2*psi3_i(2:xpoints-1)+psi3_i(1:xpoints-2)))
        
        k3_i(2:xpoints-1) = -(v(2:xpoints-1)*psi3_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_r(3:xpoints)-2*psi3_r(2:xpoints-1)+psi3_r(1:xpoints-2)))
        
        k3_r(1) = 0.d0
        k3_i(1) = 0.d0
        k3_r(xpoints) = 0.d0
        k3_i(xpoints) = 0.d0
        
        psi4_r=psi_r+dt*k3_r
        psi4_i=psi_i+dt*k3_i
        
        k4_r(2:xpoints-1) = (v(2:xpoints-1)*psi4_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_i(3:xpoints)-2*psi4_i(2:xpoints-1)+psi4_i(1:xpoints-2)))
        
        k4_i(2:xpoints-1) = -(v(2:xpoints-1)*psi4_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_r(3:xpoints)-2*psi4_r(2:xpoints-1)+psi4_r(1:xpoints-2)))
                
        k4_r(1) = 0.d0
        k4_i(1) = 0.d0
        k4_r(xpoints) = 0.d0
        k4_i(xpoints) = 0.d0

        psi_r = psi_r + dt/6.d0*(k1_r+2.d0*k2_r+2.d0*k3_r+k4_r)
        psi_i = psi_i + dt/6.d0*(k1_i+2.d0*k2_i+2.d0*k3_i+k4_i)
enddo
end subroutine rk4

subroutine rk4_pbc(dt,dx,m,h_barra,xpoints,v,psi_r,psi_i,steps_print)
implicit none
integer                           ::      xpoints, i, steps_print
real(8)                           ::      dt, dx, m, h_barra
real(8),dimension(xpoints)        ::      v
real(8),dimension(xpoints)        ::      psi_r, psi_i
real(8),dimension(xpoints)        ::      k1_r, k1_i, k2_r, k2_i, k3_r, k3_i, k4_r, k4_i
real(8),dimension(xpoints)        ::      psi2_r, psi2_i, psi3_r, psi3_i, psi4_r, psi4_i

do i=1, steps_print,1
        k1_r(2:xpoints-1) = (v(2:xpoints-1)*psi_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(3:xpoints)-2*psi_i(2:xpoints-1)+psi_i(1:xpoints-2)))
        
        k1_i(2:xpoints-1) = -(v(2:xpoints-1)*psi_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(3:xpoints)-2*psi_r(2:xpoints-1)+psi_r(1:xpoints-2)))
        
        k1_r(1) = (v(1)*psi_i(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(2)-2*psi_i(1)+psi_i(xpoints)))
        
        k1_i(1) = -(v(1)*psi_r(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(2)-2*psi_r(1)+psi_r(xpoints)))
        
        k1_r(xpoints) = (v(xpoints)*psi_i(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_i(1)-2*psi_i(xpoints)+psi_i(xpoints-1)))
        
        k1_i(xpoints) = -(v(xpoints)*psi_r(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi_r(1)-2*psi_r(xpoints)+psi_r(xpoints-1)))
        
        psi2_r=psi_r+dt*k1_r/2.d0
        psi2_i=psi_i+dt*k1_i/2.d0
        
        k2_r(2:xpoints-1) = (v(2:xpoints-1)*psi2_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_i(3:xpoints)-2*psi2_i(2:xpoints-1)+psi2_i(1:xpoints-2)))
        
        k2_i(2:xpoints-1) = -(v(2:xpoints-1)*psi2_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_r(3:xpoints)-2*psi2_r(2:xpoints-1)+psi2_r(1:xpoints-2)))
        
        k2_r(1) = (v(1)*psi2_i(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_i(2)-2*psi2_i(1)+psi2_i(xpoints)))
        
        k2_i(1) = -(v(1)*psi2_r(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_r(2)-2*psi2_r(1)+psi2_r(xpoints)))
        
        k2_r(xpoints) = (v(xpoints)*psi2_i(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_i(1)-2*psi2_i(xpoints)+psi2_i(xpoints-1)))
        
        k2_i(xpoints) = -(v(xpoints)*psi2_r(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi2_r(1)-2*psi2_r(xpoints)+psi2_r(xpoints-1)))
        
        psi3_r=psi_r+dt*k2_r/2.d0
        psi3_i=psi_i+dt*k2_i/2.d0
        
        k3_r(2:xpoints-1) = (v(2:xpoints-1)*psi3_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_i(3:xpoints)-2*psi3_i(2:xpoints-1)+psi3_i(1:xpoints-2)))
        
        k3_i(2:xpoints-1) = -(v(2:xpoints-1)*psi3_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_r(3:xpoints)-2*psi3_r(2:xpoints-1)+psi3_r(1:xpoints-2)))
        
        k3_r(1) = (v(1)*psi3_i(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_i(2)-2*psi3_i(1)+psi3_i(xpoints)))
        
        k3_i(1) = -(v(1)*psi3_r(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_r(2)-2*psi3_r(1)+psi3_r(xpoints)))
        
        k3_r(xpoints) = (v(xpoints)*psi3_i(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_i(1)-2*psi3_i(xpoints)+psi3_i(xpoints-1)))
        
        k3_i(xpoints) = -(v(xpoints)*psi3_r(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi3_r(1)-2*psi3_r(xpoints)+psi3_r(xpoints-1)))
        
        psi4_r=psi_r+dt*k3_r
        psi4_i=psi_i+dt*k3_i
        
        k4_r(2:xpoints-1) = (v(2:xpoints-1)*psi4_i(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_i(3:xpoints)-2*psi4_i(2:xpoints-1)+psi4_i(1:xpoints-2)))
        
        k4_i(2:xpoints-1) = -(v(2:xpoints-1)*psi4_r(2:xpoints-1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_r(3:xpoints)-2*psi4_r(2:xpoints-1)+psi4_r(1:xpoints-2)))
                
        k4_r(1) = (v(1)*psi4_i(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_i(2)-2*psi4_i(1)+psi4_i(xpoints)))
        
        k4_i(1) = -(v(1)*psi4_r(1)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_r(2)-2*psi4_r(1)+psi4_r(xpoints)))
        
        k4_r(xpoints) = (v(xpoints)*psi4_i(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_i(1)-2*psi4_i(xpoints)+psi4_i(xpoints-1)))
        
        k4_i(xpoints) = -(v(xpoints)*psi4_r(xpoints)/h_barra - &
        (h_barra/(2.d0*m*dx**2))*(psi4_r(1)-2*psi4_r(xpoints)+psi4_r(xpoints-1)))

        psi_r = psi_r + dt/6.d0*(k1_r+2.d0*k2_r+2.d0*k3_r+k4_r)
        psi_i = psi_i + dt/6.d0*(k1_i+2.d0*k2_i+2.d0*k3_i+k4_i)
enddo
end subroutine rk4_pbc
