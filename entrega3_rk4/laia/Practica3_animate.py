# PRACTICA 1 DINAMICA QUANTICA

# --------------------------------------------------------------------------------------------------------
# MODULS
# --------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --------------------------------------------------------------------------------------------------------
# FUNCIONS
# --------------------------------------------------------------------------------------------------------

# Definim la part real i la part imaginaria del paquet d'ones
def Psi_Real(x, t):
    return (a**2/(2.*np.pi))**(1/4) * np.exp(-a**2*(x-x0)**2/4.) * np.cos(k0*x - w0*t)    # Part real a temps = 0 (inicial)


def Psi_Imaginari(x, t):
    return (a**2/(2.*np.pi))**(1/4) * np.exp(-a**2*(x-x0)**2/4.) * np.sin(k0*x - w0*t)    # Part imaginaria a temps = 0 (inicial)


# Definim la solucio analitica
def Psi_sol_analitica(x, t):
    return (a**2/(2.*np.pi))**(1/4) * (sigma/(np.sqrt(sigma**2 + (1j*t*h_barra/m)))) * np.exp(1j*(k0*x-w0*t))*np.exp(-(x-x0-k0*h_barra*t/m)**2/(2*sigma**2+2*1j*t*h_barra/m))


# Definim els diferents tipus de potencial
def V_part_lliure(x):
    return np.zeros(len(x))         


# Metode Forward Euler
def Euler_fwd(PsiR, PsiI,V):
    PsiR_t = np.zeros(len(x))
    PsiI_t = np.zeros(len(x))
    
    # Condicions de contorn
    PsiR_t[0] = 0 ; PsiR_t[-1] = 0
    PsiR_t[0] = 0 ; PsiI_t[-1] = 0

    # Loop per calcular vectors Psi per tot el vector x a cada temps
    for i in range(1, len(x)-1):                                                # El primer punt i l'ultim no els tenim en compte, s'anul·la la funcio d'ona
        PsiR_t[i] = dt*(V[i]*PsiI[i]/h_barra - h_barra/(2.*m*dx**2)*(PsiI[i+1] - 2*PsiI[i] + PsiI[i-1]))
        PsiI_t[i] = -dt*(V[i]*PsiR[i]/h_barra - h_barra/(2.*m*dx**2)*(PsiR[i+1] - 2*PsiR[i] + PsiR[i-1]))
        
    # Actualitzem els vectors
    PsiR = PsiR + PsiR_t                                                    
    PsiI = PsiI + PsiI_t
        
    return PsiR, PsiI

# Metode Runge Kutta d'ordre 4
def RK4(PsiR, PsiI, V):  
    
    k1_I = np.zeros(len(x)) ; k1_R = np.zeros(len(x))
    k2_I = np.zeros(len(x)) ; k2_R = np.zeros(len(x))
    k3_I = np.zeros(len(x)) ; k3_R = np.zeros(len(x))
    k4_I = np.zeros(len(x)) ; k4_R = np.zeros(len(x))
    
    PsiR_2 = np.zeros(len(x)) ; PsiI_2 = np.zeros(len(x))
    PsiR_3 = np.zeros(len(x)) ; PsiI_3 = np.zeros(len(x))
    PsiR_4 = np.zeros(len(x)) ; PsiI_4 = np.zeros(len(x))
    
    # Condicions de contorn
    k1_I[0] = 0 ; k1_R[-1] = 0 
    k2_I[0] = 0 ; k2_R[-1] = 0
    k3_I[0] = 0 ; k3_R[-1] = 0 
    k4_I[0] = 0 ; k4_R[-1] = 0

    for i in range(1, len(x)-1):
        #print("PUTAAAAAA")
        k1_I[i] = (h_barra/(2*m*dx**2))*(PsiR[i+1] - 2*PsiR[i] + PsiR[i-1]) - PsiR[i]*V[i]/h_barra
        k1_R[i] = -(h_barra/(2*m*dx**2))*(PsiI[i+1] - 2*PsiI[i] + PsiI[i-1]) + PsiI[i]*V[i]/h_barra
        
    PsiR_2 = PsiR + k1_R*dt/2
    PsiI_2 = PsiI + k1_I*dt/2
    #print("k1_I", k1_I, "k1_R", k1_R)    
    #print("PsiR=", PsiR, "PsiR_2=", PsiR_2)
        
        #print(count)
        
    for i in range(1, len(x)-1):
        k2_I[i] = (h_barra/(2*m*dx**2))*(PsiR_2[i+1] - 2*PsiR_2[i] + PsiR_2[i-1]) - PsiR_2[i]*V[i]/h_barra
        k2_R[i] = -(h_barra/(2*m*dx**2))*(PsiI_2[i+1] - 2*PsiI_2[i] + PsiI_2[i-1]) + PsiI_2[i]*V[i]/h_barra
        
    PsiR_3 = PsiR + k2_R*dt/2
    PsiI_3 = PsiI + k2_I*dt/2
    
    for i in range(1, len(x)-1):
        k3_I[i] = (h_barra/(2*m*dx**2))*(PsiR_3[i+1] - 2*PsiR_3[i] + PsiR_3[i-1]) - PsiR_3[i]*V[i]/h_barra
        k3_R[i] = -(h_barra/(2*m*dx**2))*(PsiI_3[i+1] - 2*PsiI_3[i] + PsiI_3[i-1]) + PsiI_3[i]*V[i]/h_barra
        
    PsiR_4 = PsiR + k3_R*dt
    PsiI_4 = PsiI + k3_I*dt
        
    for i in range(1, len(x)-1):
        k4_I[i] = (h_barra/(2*m*dx**2))*(PsiR_4[i+1] - 2*PsiR_4[i] + PsiR_4[i-1]) - PsiR_4[i]*V[i]/h_barra
        k4_R[i] = -(h_barra/(2*m*dx**2))*(PsiI_4[i+1] - 2*PsiI_4[i] + PsiI_4[i-1]) + PsiI_4[i]*V[i]/h_barra

    PsiR = PsiR + (dt/6)*(k1_R + 2*k2_R + 2*k3_R + k4_R)
    PsiI = PsiI + (dt/6)*(k1_I + 2*k2_I + 2*k3_I + k4_I)
        
    #Psi2[i] = abs(PsiR[i])**2 + abs(1j*PsiI[i])**2

    #print("PsiR=", PsiR, "PsiI=", PsiI)
    return PsiR, PsiI


def animate(i):
    global PsiR, PsiI
    PsiR, PsiI = RK4(PsiR, PsiI, V_part_lliure(x))         # Posar el potencial que es vulgui
    #PsiR, PsiI = Euler_fwd(PsiR, PsiI, V_part_lliure(x))         # Posar el potencial que es vulgui
    Psi_analitica = Psi_sol_analitica(x, i*dt)
    
    norma_numerica = dx*np.sum(abs(PsiR)**2 + abs(1j*PsiI)**2)
    norma_analitica = dx*np.sum(abs(Psi_analitica)**2)

    real_line.set_data(x, PsiR)
    imagine_line.set_data(x, PsiI)
    quadrat_line.set_data(x, abs(PsiR)**2 + abs(1j*PsiI)**2)
    sol_analitica.set_data(x, abs(Psi_analitica)**2)
    
    time_text.set_text(time_template % (i*dt))
    dif_norma_text.set_text(dif_norma_template % abs(norma_numerica-norma_analitica))
    
    return real_line, imagine_line, quadrat_line, sol_analitica, time_text, dif_norma_text

global a, k0, w0, x0, dt, dx, npts, x_ini, x_fin, m, h_barra

# --------------------------------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------------------------------

# Dades inicials per començar el problema
a = 3                                                                           # 1/amplada
x_ini = -3
x_fin = 3
x0 = 0                                                                         # On esta centrat el paquet
npts = 180
x = np.linspace(x_ini, x_fin, npts)                                             # Vector posicio
dx = (x_fin - x_ini)/npts
w0 = 1
k0 = 10
m = 1
h_barra = 1
t0 = 0                                                                          # Temps inicial
dt = 5e-5
sigma = np.sqrt(2)/a

'''
# A TEMPS INICIAL (temps = 0) -------------------------------------------------
PsiR_0 = Psi_Real(x, t0)
PsiI_0 = Psi_Imaginari(x, t0)
Psi_analitica_0 = Psi_sol_analitica(x, t0)
# Part real
plt.xlabel('x')
plt.ylabel('$Ψ_{R}$')
plt.plot(x, PsiR_0, label = 'Part real')
plt.savefig('Practica3_PsiR0.png', dpi=200)
#plt.show()

# Part imaginaria
plt.xlabel('x')
plt.ylabel('$Ψ_{I}$')
plt.plot(x, PsiI_0, label = 'Part imaginaria')
plt.savefig('Practic31_PsiI0.png', dpi=200)
#plt.show()

# Comparacio solucio numerica i analitica
plt.xlabel('x')
plt.ylabel('$|Ψ|^{2}$')
plt.plot(x, abs(1j*PsiI_0)**2 + abs(PsiR_0)**2, label='Modul al quadrat', color='green')
plt.legend(loc='best')
plt.savefig('Practica3_modul.png', dpi=200)
#plt.show()

# Comparacio solucio numerica i analitica
plt.xlabel('x')
plt.ylabel('$|Ψ|^{2}$')
plt.plot(x, abs(1j*PsiI_0)**2 + abs(PsiR_0)**2, label='solucio numerica')
plt.plot(x, abs(Psi_analitica_0)**2, linestyle='dashed', color='r', label='solucio analitica')
plt.legend(loc='best')
plt.savefig('Practica3_numvsanalitic.png', dpi=200)
#plt.show()
'''
# EVOLUCIO DEL PAQUET D'ONES --------------------------------------------------
PsiR = Psi_Real(x, t0)
PsiI = Psi_Imaginari(x, t0)
Psi_analitica = Psi_sol_analitica(x, t0)

fig = plt.figure(figsize=(5,5), dpi=200)
ax = plt.axes(xlim=(x_ini, x_fin), ylim=(x_ini, x_fin))

time_template = 'temps = %.5fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

dif_norma_template = 'Δnorma = %.5f'
dif_norma_text = ax.text(0.05, 0.75, '', transform=ax.transAxes)

plt.tight_layout()
real_line, = ax.plot([],[],label='$Ψ_{R}$')
imagine_line, = ax.plot([],[],label='$Ψ_{I}$')
quadrat_line, = ax.plot([],[],label='$|Ψ|^{2}_{numerica}$')
sol_analitica, = ax.plot([],[],label='$|Ψ|^{2}_{analitica}$')
ani = FuncAnimation(fig, animate, frames =10000, interval=0, blit=True, repeat=True)
plt.legend(loc='best')
plt.show()


