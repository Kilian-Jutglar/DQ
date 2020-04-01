# PRACTICA 3 DINAMICA QUANTICA (RUNGE-KUTTA ORDRE 4)

# --------------------------------------------------------------------------------------------------------
# MODULS
# --------------------------------------------------------------------------------------------------------
import numpy as np

# --------------------------------------------------------------------------------------------------------
# FUNCIONS
# --------------------------------------------------------------------------------------------------------

# Definim la part real i la part imaginaria del paquet d'ones
def Psi_Real(x, t):
    return (a**2/(2.*np.pi))**(1/4) * np.exp(-a**2*(x-x0)**2/4.) * np.cos(k0*x - w0*t)    # Part real a temps = 0 (inicial)


def Psi_Imaginari(x, t):
    return (a**2/(2.*np.pi))**(1/4) * np.exp(-a**2*(x-x0)**2/4.) * np.sin(k0*x - w0*t)    # Part imaginaria a temps = 0 (inicial)


# Definim el tipus de potencial

# Potencial step function
def V_barrera(x0_b, xf_b):
    V = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] >= x0_b and x[i] <= xf_b:
            V[i] = 3*E/100
        else:
            V[i] = 0
    #print("V", V)
    #print("x0_b, xf_b:", x0_b,",", xf_b)
    all_x = np.argwhere(V == 3*E/100).T[0]                             # Retorna els indexs del vector V en que el potencial no es zero
    #print("all_x=", all_x)                                            # .T[0] es per transposar ja que argwhere retorna una matriu d'1 columna
    index_x0, index_xf = all_x[0], all_x[-1]                           # Necessitem els indexs per obtenir la reflexio i la transmissio
    #print("index_x0, index_xf:", index_x0,",", index_xf)
    return index_x0, index_xf, V
       

# Definim el metode amb el qual resolem les equacions diferencials
# Metode Runge Kutta d'ordre 4
def RK4(PsiR, PsiI, V):
    k1_I = np.zeros(len(x)) ; k1_R = np.zeros(len(x))
    k2_I = np.zeros(len(x)) ; k2_R = np.zeros(len(x))
    k3_I = np.zeros(len(x)) ; k3_R = np.zeros(len(x))
    k4_I = np.zeros(len(x)) ; k4_R = np.zeros(len(x))
    
    PsiR_2 = np.zeros(len(x)) ; PsiI_2 = np.zeros(len(x))
    PsiR_3 = np.zeros(len(x)) ; PsiI_3 = np.zeros(len(x))
    PsiR_4 = np.zeros(len(x)) ; PsiI_4 = np.zeros(len(x))

    for i in range(1, len(x)-1):
        #print("PUTAAAAAA")
        k1_I[i] = (h_barra/2*m*dx**2)*(PsiR[i+1] - 2*PsiR[i] + PsiR[i-1]) - PsiR[i]*V[i]/h_barra
        k1_R[i] = -(h_barra/2*m*dx**2)*(PsiI[i+1] - 2*PsiI[i] + PsiI[i-1]) + PsiI[i]*V[i]/h_barra
        
    PsiR_2 = PsiR + k1_R*dt/2
    PsiI_2 = PsiI + k1_I*dt/2
    #print("k1_I", k1_I, "k1_R", k1_R)    
    #print("PsiR=", PsiR, "PsiR_2=", PsiR_2)
        
        #print(count)
        
    for i in range(1, len(x)-1):
        k2_I[i] = (h_barra/2*m*dx**2)*(PsiR_2[i+1] - 2*PsiR_2[i] + PsiR_2[i-1]) - PsiR_2[i]*V[i]/h_barra
        k2_R[i] = -(h_barra/2*m*dx**2)*(PsiI_2[i+1] - 2*PsiI_2[i] + PsiI_2[i-1]) + PsiI_2[i]*V[i]/h_barra
        
    PsiR_3 = PsiR + k2_R*dt/2
    PsiI_3 = PsiI + k2_I*dt/2
    
    for i in range(1, len(x)-1):
        k3_I[i] = (h_barra/2*m*dx**2)*(PsiR_3[i+1] - 2*PsiR_3[i] + PsiR_3[i-1]) - PsiR_3[i]*V[i]/h_barra
        k3_R[i] = -(h_barra/2*m*dx**2)*(PsiI_3[i+1] - 2*PsiI_3[i] + PsiI_3[i-1]) + PsiI_3[i]*V[i]/h_barra
        
    PsiR_4 = PsiR + k3_R*dt/2
    PsiI_4 = PsiI + k3_I*dt/2
        
    for i in range(1, len(x)-1):
        k4_I[i] = (h_barra/2*m*dx**2)*(PsiR_4[i+1] - 2*PsiR_4[i] + PsiR_4[i-1]) - PsiR_4[i]*V[i]/h_barra
        k4_R[i] = -(h_barra/2*m*dx**2)*(PsiI_4[i+1] - 2*PsiI_4[i] + PsiI_4[i-1]) + PsiI_4[i]*V[i]/h_barra

    PsiR = PsiR + (dt/6)*(k1_R + 2*k2_R + 2*k3_R + k4_R)
    PsiI = PsiI + (dt/6)*(k1_I + 2*k2_I + 2*k3_I + k4_I)
        
    #Psi2[i] = abs(PsiR[i])**2 + abs(1j*PsiI[i])**2

    #print("PsiR=", PsiR, "PsiI=", PsiI)
    return PsiR, PsiI

# Metode forward Euler
def Euler_fwd(PsiR, PsiI, V):
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


# Definim la funcio per obtenir els coeficients de reflexio i transmissio
def get_refl_trans(func):
    reflexio = np.sum(func[0:index_xf])*dx            # Integrem des del primer valor de Psi2 fins el primer valor on el potencial no es zero
    transmissio = np.sum(func[index_xf:])*dx          # Integrem des del primer valor on el potencial no es zero fins l'ultim valor de Psi2
    return reflexio, transmissio


# --------------------------------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------------------------------

# Dades inicials per començar el problema
a = 3                                                                           # 1/amplada
x_ini = -5
x_fin = 5
x0 = -2.5                                                                         # On esta centrat el paquet
npts = 120
x = np.linspace(x_ini, x_fin, npts)                                             # Vector posicio
dx = (x_fin - x_ini)/npts
k0 = 10
m = 1
h_barra = 1
w0 = h_barra*k0**2/(2*m)
t0 = 0                                                                        # Temps inicial
sigma = np.sqrt(2)/a
E = h_barra*w0
x0_b, xf_b = x[int(2*len(x)/3)], x[int(3*len(x)/4)]
amplada = xf_b - x0_b
print("Amplada barrera:", amplada)

Emax = 0.03*E 


# SIMULACIO
index_x0, index_xf, V = V_barrera(x0_b, xf_b)
count=1

while V[index_x0] <= Emax and V[index_xf] <= Emax:
    
    time = 0
    dt = 1e-5
    Ntimesteps = 100000
    Nprint = 1000
    
    print("Alçada barrera:", V[index_x0])
    
    # VALORS INICIALS
    PsiR = Psi_Real(x, t0)
    PsiI = Psi_Imaginari(x, t0)
    
    #print("PsiRinicial=", PsiR, "PsIinicial=", PsiI)
    results_Psi = open('resultats_Psi_{}.dat'.format(count), 'w')
    results_rt = open('resultats_Reflexio_Transmisio_V_{}.dat'.format(count), 'w')
    
    results_Psi.write("# FUNCIO D'ONA\n\n")
    results_Psi.write("# x        PsiR        # PsiI        # Psi2        # Potencial"+'\n')
    
    results_rt.write("# E = %f ; Alçada = %.20f\n"%(E, V[index_x0]))
    results_rt.write("# Amplada = %f\n\n"%(amplada))
    results_rt.write("#Temps       #Reflexio       #Transmissio"+'\n')
    
    # EVOLUCIO DEL PAQUET D'ONES --------------------------------------------------
    for t in range(Ntimesteps):
        
        PsiR, PsiI = RK4(PsiR, PsiI, V)  
        #print("PsiR=", PsiR, "PsiI=", PsiI)
        
        Psi2 = abs(PsiR)**2 + abs(1j*PsiI)**2
        
        #print("temps", time)
        time = time + dt
        
        if t % Nprint == 0:                      # Calcul del residu. Per no printar cada calcul
            print("puta")
            reflexio, transmissio = get_refl_trans(Psi2)
            print(time, reflexio, transmissio)
            
            # RESULTATS
            for i in range(len(x)):
                results_Psi.write('%f %f %f %f %f\n'%(x[i], PsiR[i], PsiI[i], Psi2[i], V[i]))
                
            results_rt.write('%e %f %f\n'%(time, reflexio, transmissio))
    
    V[index_x0:index_xf+1] = V[index_x0] + 0.03*E
    count = count + 1

    results_Psi.close()
    results_rt.close()

print("Finished")
