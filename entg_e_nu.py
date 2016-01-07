#simulation of entangled electrons with nuclear spins

import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#one electron was interacting with 3 electrons each
psi_nu1 = tensor(rand_ket(2).unit,rand_ket(2).unit,rand_ket(2).unit).unit()
psi_nu2 = tensor(rand_ket(2).unit,rand_ket(2).unit,rand_ket(2).unit).unit()


#initial electron spin state
# |0> state = basis(2,0)
# |1> state = basis(2,1)

psi_e = (tensor(basis(2,0),basis(2,0)) + tensor(basis(2,1),basis(2,1))).unit()

#exchange Hamiltonian
#qeye(2) -> 2 x 2 identity matrix
#sigmaz is Pauli z spin operator

#identity matrix for the three nuclear spins
nu_iden = tensor(qeye(2),qeye(2),qeye(2))
#identity matrix for a nuclear spin
e_iden = qeye(2)

I_z = tensor(sigmaz(),qeye(2),qeye(2)) + tensor(qeye(2),sigmaz(),qeye(2)) + tensor(qeye(2),qeye(2),sigmaz())
I_p = tensor(sigmap(),qeye(2),qeye(2)) + tensor(qeye(2),sigmap(),qeye(2)) + tensor(qeye(2),qeye(2),sigmap())
I_m = tensor(sigmam(),qeye(2),qeye(2)) + tensor(qeye(2),sigmam(),qeye(2)) + tensor(qeye(2),qeye(2),sigmam())


#Hamiltonian
# e x e x nu1 x nu2

#interaction of first set of nuclear spins with the first electron
H_nu1e = tensor(sigmaz(),e_iden,I_z,nu_iden) + 0.5 * (tensor(sigmap(),e_iden,I_m,nu_iden) + tensor(sigmam(),e_iden,I_p,nu_iden))
#interaction of second set of nuclear spins with the second electron
H_nu2e = tensor(e_iden,sigmaz(),nu_iden,I_z) + 0.5 * (tensor(e_iden,sigmap(),nu_iden,I_m) + tensor(e_iden,sigmam(),nu_iden,I_p))

#interaction of electron-electron
H_ee = tensor(sigmaz(),sigmaz(),nu_iden,nu_iden) + 0.5 * (tensor(sigmap(),sigman()) + tensor(sigman(),sigmap()))

#paramters with control the strength of interaction
#nu1 and e1
g1 = 1.0
#nu2 and e2
g2 = 1.0

#e and e
ge = 10.0
H = g1 * H_nu1e +  g2 * H_nu2e + ge * H_ee

#initial state 
psi0 = tensor(psi_e,psi_nu1,psi_nu2)


#number of simulation points
N = 100
tlist = np.linspace(0,100,N)
opt = solver.Options(store_states = True)
output = mesolve(H,psi0,tlist,[],[tensor(sigmaz(),qeye(128)),tensor(eye(2),sigmaz(),qeye(64))],options=opt)

#calculate the expectation values
#spin of electron 1
s1 = output.expect[0]
#spin of electron 2
s2 = output.expect[1]

#plot of electron spins vs time
plt.plot(s1)
plt.plot(s2)
plt.show()



