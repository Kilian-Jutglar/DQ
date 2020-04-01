# Doble Gif de la funció d'Ona + els coeficients de Transmissió i Reflexió.
##############################################################################
# Fitxer de les Funcions d'Ona:
#	X, FR, FI, FNorma, Potencial
#
# Fitxer dels Coeficients:
#	Temps, T, R
##############################################################################
# Definim variables (tot aixó ho pots treure del input dels teus càlculs):
Ntimesteps=100  # El número de time steps que tinguis
Nspace=180     # Número de punts en l'espai en que has dividit la caixa
tmax=9.900100e-01      # Temps màxim al que arriba la teva simulació
###############################################################################
# SIMULACIÓ
###############################################################################
# Definim la terminal per fer un Gif i l'output
set terminal gif animate
set output 'gif_altura_1.gif'

# Definim el loop, el qual s'ha de fer tants cops com timesteps
do for [ii=1:Ntimesteps] {
ini = (ii-1)*Nspace + 4   	# En el fitxer de la funcio d'ona quin punt inicial ha d'agafar per aquest timestep
fin = ini + Nspace - 1          # En el fitxer de la funcio d'ona quin punt final ha d'agafar per aquest timestep 

# Plot de la Funcio d'Ona
set title "Funcio d'ona" font ',14'
set xlabel 'x'
set yrange [-1:3]  	# Aquest rang és per que la funció només pot agafar el valor 1 o -1 com a molt
set xrange [-3:6]       # Mida de la caixa
plot 'resultats_Psi_1.dat' every ::ini::fin u 1:2  w l t 'Ψ_{Real}', '' every ::ini::fin u 1:3 w l t 'Ψ_{Imaginaria}', '' every ::ini::fin u 1:4 w l t '|Ψ|^2', '' every ::ini::fin u ($1):($5) w l t 'V(x)'
# Fixa't que el 'every' serveix per especificar quin rang de punts en concret vols fer servir
}

unset terminal
