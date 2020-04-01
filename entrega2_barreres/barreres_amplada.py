## Dinàmica Quàntica
from time import time
import f90
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def build_r(x,t=0.):
    return (a**2/(2.*np.pi))**0.25*np.exp(-a**2*(x-x0)**2/4.)*np.cos(k*x-w*t)

def build_i(x,t=0.):
    return (a**2/(2.*np.pi))**0.25*np.exp(-a**2*(x-x0)**2/4.)*np.sin(k*x-w*t)

def euler(psi_r,psi_i,steps=1):
    for i in range(steps):
        pas_r=np.zeros_like(psi_r); pas_i=np.zeros_like(psi_i); pas_r[0]=0.; pas_i[0]=0.
        
        pas_r[1:-1]= dt*(v[1:-1]*psi_i[1:-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_i[2:]-2*psi_i[1:-1] + psi_i[:-2]))
        pas_i[1:-1]=-dt*(v[1:-1]*psi_r[1:-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_r[2:]-2*psi_r[1:-1] + psi_r[:-2]))
    
        pas_r[-1]=0.; pas_i[-1]=0.
        psi_r += pas_r
        psi_i += pas_i
    return psi_r, psi_i

def euler_pbc(psi_r,psi_i,steps=1):
    for i in range(steps):
        pas_r=np.zeros_like(psi_r); pas_i=np.zeros_like(psi_i)
                
        pas_r[0]=dt*(v[0]*psi_i[0]/h_barra - (h_barra/(2.*m*dx**2))*(psi_i[1]-2*psi_i[0] + psi_i[-1]))
        pas_r[1:-1]= dt*(v[1:-1]*psi_i[1:-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_i[2:]-2*psi_i[1:-1] + psi_i[:-2]))
        pas_r[-1]=dt*(v[-1]*psi_i[-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_i[-2]-2*psi_i[-1] + psi_i[0]))

        pas_i[0]=-dt*(v[0]*psi_r[0]/h_barra - (h_barra/(2.*m*dx**2))*(psi_r[1]-2*psi_r[0] + psi_r[-1]))
        pas_i[1:-1]=-dt*(v[1:-1]*psi_r[1:-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_r[2:]-2*psi_r[1:-1] + psi_r[:-2]))
        pas_i[-1]=-dt*(v[-1]*psi_r[-1]/h_barra - (h_barra/(2.*m*dx**2))*(psi_r[-2]-2*psi_r[-1] + psi_r[0]))
                
        psi_r += pas_r
        psi_i += pas_i
    return psi_r, psi_i

def v_barrera(h,width):
    init = -width/2.
    fin = width/2.
    init_point = int((init-xinit)/dx)
    fin_point = int((fin-xinit)/dx)
    v[init_point:fin_point] = h
    v_mid = int((fin_point+init_point)/2)
    return v, v_mid

def v_caixa():
    h = 10
    v[:int(xpoints/4)]=h
    v[-int(xpoints/4):]=h
    return v

def v_dent_de_serra(d,h):
    indx_tram = int(d/dx)
    rect = np.linspace(0,d,indx_tram)
    count=0
    while count < xpoints:
        for j in range(indx_tram):
            v[count] = rect[j]
            count += 1
    return v

def v_harmonic():
    v = x**2
    return v


def n(psi2):
    n=dx*np.sum(psi2)
    return n

def r(psi2,v_mid):
    refl_val=dx*np.sum(psi2[:v_mid])
    return refl_val

def t(psi2,v_mid):
    trans_val=dx*np.sum(psi2[v_mid:])
    return trans_val

start=time()

global a, k, w, x0, dt, dx, xpoints, xinit, xfin, x, m, h_barra, v,v_mid, steps_print

steps_print=1000
steps=1.5e3
m=1 ; h_barra=1; a=5. ;k=10. ;w=1. ;x0=-7. # Paràmetres del paquet d'ones
dt = 1e-6 # ; dx=0.01
L=30
xinit=-L/2.; xfin=L/2.
xpoints = 1500
dx=(xfin-xinit)/xpoints
x = np.linspace(xinit,xfin,xpoints)
v_mid=0
v = np.zeros_like(x)
v_points=22
E=h_barra**2*k**2/(2*m)
min_v=0.
max_v=L/40.
v_values=np.linspace(min_v,max_v,v_points)
t_v= np.zeros(v_points)
r_v= np.zeros(v_points)
## Select the potential
h=9.*E/10.

# Manual without animation
for value in range(len(v_values)):
    v = np.zeros_like(x)
    v,v_mid = v_barrera(h,v_values[value])
    ## Build WP    
    dt_crit = (2.*m*h_barra*dx**2)/(h_barra**2 + 2.*m*dx**2*np.max(v))

    if dt >= dt_crit:
        print('dt should be bigger than',dt_crit)
        print('Exitting program')
        exit()
    #psi_r = build_r(x)
    #psi_i = build_i(x)
    psi_r=np.zeros_like(x)
    psi_i=np.zeros_like(x)
    f90.build_wp(a=a,x0=x0,k=k,w=w,xpoints=xpoints,x=x,psi_r=psi_r,psi_i=psi_i)
    psi2 = psi_r**2+psi_i**2

    norm = np.zeros(1+int(steps))
    refl = np.zeros(1+int(steps))
    trans = np.zeros(1+int(steps))
    t_pas = np.zeros(1+int(steps))

    norm[0]=n(psi2)
    refl[0]= r(psi2,v_mid)
    trans[0]= t(psi2,v_mid)
    t_pas[0]=0.

    for i in range(int(steps)):
        j = i + 1
        f90.euler(dt=dt,dx=dx,m=m,h_barra=h_barra,v=v,psi_r=psi_r,psi_i=psi_i,steps_print=steps_print)
        #psi_r, psi_i = euler_pbc(psi_r,psi_i,steps_print)
        psi2 = psi_r**2+psi_i**2
        t_pas_val = (j)*dt*steps_print
        norm_val = n(psi2)
        refl_val= r(psi2,v_mid)/norm_val
        trans_val= t(psi2,v_mid)/norm_val
    
        refl[j]=refl_val
        trans[j]=trans_val
        norm[j]=norm_val
        t_pas[j]= t_pas_val
    
       # if  trans[i]  > trans[j]:
    stop = i
    t_v[value]= trans[stop]
    r_v[value]= refl[stop]
        # break 
# Final figure Manual
'''norm_template = 'Norm = %.5fs'
norm=dx*np.sum(psi_r**2+psi_i**2)
norm_text = plt.text(0.05, 0.9, '', transform=plt.transAxes)
norm_text.set_text(norm_template % (norm))
plt.title("Paquet d'ones planes",fontsize=16)
plt.text(0.05,0.6,"Norma=")
plt.plot(x,psi_r**2 + psi_i**2,label='|ψ|²')
plt.plot(x,psi_r,label='Part Real')
plt.plot(x,psi_i,label='Part Imaginària')
plt.xlabel('x',fontsize=14)
plt.ylabel('|ψ|²',fontsize=14)
plt.legend(loc='best')
plt.savefig('ona_plana_desc_t0.eps', format='eps')'''

# Animation
def animate(i):
    global psi_r, psi_i, v, trans, refl, t_pas
    psi_r, psi_i = euler_pbc(psi_r,psi_i,steps_print)
    psi2 = psi_r**2+psi_i**2
    #real_line.set_data(x,psi_r)
    #imagine_line.set_data(x,psi_i)
    quadrat_line.set_data(x,psi2)
    t_pas_val = i*dt*steps_print
    t_pas_text.set_text(t_pas_template % (t_pas_val))
    norm_val = n(psi2)
    norm_text.set_text(norm_template % (norm_val))
    refl_val= r(psi2,v_mid)
    refl_text.set_text(refl_template % (refl_val))
    trans_val= t(psi2,v_mid)
    trans_text.set_text(trans_template % (trans_val))
    trans[i]=trans_val
    refl[i]=refl_val
    t_pas[i]= t_pas_val
    return quadrat_line,t_pas_text, norm_text, refl_text, trans_text, real_line, imagine_line

# Animation figure
'''
fig = plt.figure(figsize=(4,4),dpi=200)
ax=plt.axes(xlim=(xinit,xfin),ylim=(0.,4.0))
t_pas_template = 'Time  = %.5fs'
norm_template = 'Norm = %.5f'
refl_template = 'Refl    = %.5f'
trans_template = 'Trans  = %.5f'
t_pas_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
norm_text = ax.text(0.05, 0.85, '', transform=ax.transAxes)
refl_text = ax.text(0.05, 0.75, '', transform=ax.transAxes)
trans_text = ax.text(0.05, 0.7, '', transform=ax.transAxes)
plt.tight_layout()
#real_line, = ax.plot([],[],label='Real')
#imagine_line, = ax.plot([],[],label='Imagin')
quadrat_line, = ax.plot([],[],label='Modul al quadrat')
v_line, = ax.plot([],[],label='Potential')
v_line.set_data(x,v)
plt.legend(loc = 'best')
ani = FuncAnimation(fig,animate,frames=int(steps),interval=0,blit=True,repeat=False)
'''
# Figure R, T vs t
plt.ylim(0.,1.026)
plt.xlim(v_values[0],v_values[-1])
plt.title('Reflexion & Transmission of a WP vs Barrier width',fontsize=16)
plt.plot(v_values,t_v,label='WP Transmission')
plt.plot(v_values,r_v,label='WP Reflexion')
plt.xlabel('Barrier Width (a.u.)',fontsize=14)
plt.legend(loc = 'best')
plt.xticks([v_values[0],v_values[int(v_points/3)],v_values[int(2*v_points/3)],v_values[-1]],['0','L/60','L/30','L/20'])
plt.savefig('barrier_width.png',format='png',dpi=200)
plt.savefig('barrier_width.eps',format='eps')
print('CPU Time:',time()-start)
plt.show()
