#simulation of exchange Hamiltonian for electron spin and ensemble of nuclear spins

import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#initial nuclear spin state
# 6 nuclear spins have been used

#each nuclear has a random spin
#psi_nu = tensor(rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit())

#the entire ensemble has a random state
psi_nu = rand_ket(64)

#all nuclei have the same spin |1> state
#psi_nu = tensor(basis(2,1),basis(2,1),basis(2,1),basis(2,1),basis(2,1),basis(2,1))

#initial electron spin state
# |0> state = basis(2,0)
# |1> state = basis(2,1)

psi_e = basis(2,0)


#creating the state vector of the entire system
psi0 = tensor(psi_e,psi_nu)
psi0 = psi0.unit()

#exchange Hamiltonian
#qeye(2) -> 2 x 2 identity matrix
#sigmaz is Pauli z spin operator

I_z = tensor(sigmaz(),qeye(2),qeye(2),qeye(2),qeye(2),qeye(2)) + tensor(qeye(2),sigmaz(),qeye(2),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),sigmaz(),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),sigmaz(),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmaz(),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),qeye(2),sigmaz())

I_p = tensor(sigmap(),qeye(2),qeye(2),qeye(2),qeye(2),qeye(2)) + tensor(qeye(2),sigmap(),qeye(2),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),sigmap(),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),sigmap(),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmap(),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),qeye(2),sigmap())

I_m = tensor(sigmam(),qeye(2),qeye(2),qeye(2),qeye(2),qeye(2)) + tensor(qeye(2),sigmam(),qeye(2),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),sigmam(),qeye(2),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),sigmam(),qeye(2),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),sigmam(),qeye(2)) +tensor(qeye(2),qeye(2),qeye(2),qeye(2),qeye(2),sigmam())

#overall strength of the hyperfine interaction
A = 1

#Magnetic Field
B = 0
H = A * (tensor(sigmaz(),I_z) + 0.5*tensor(sigmam(),I_p) + 0.5*tensor(sigmap(),I_m)) + B * tensor(sigmaz(),qeye(2),qeye(2),qeye(2),qeye(2),qeye(2),qeye(2)) 

#number of simulation points
N = 100
tlist = np.linspace(0,100,N)
opt = solver.Options(store_states = True)
output = mesolve(H,psi0,tlist,[],[tensor(sigmaz(),qeye(64)),tensor(qeye(2),I_z),tensor(sigmax(),qeye(64)),tensor(sigmay(),qeye(64))],options=opt)

#calculate the expectation values
S_z = output.expect[0]
I_exp = output.expect[1]

x,y,z = output.expect[2],output.expect[3],output.expect[0]

#plot of spins vs time

plt.figure()
plt.subplot(211)
plt.title('Collective Nuclear Spin vs Time')
plt.plot(I_exp)
plt.ylabel(r'$<A_z>$ (Nuclear Spin)',fontsize=16)
plt.xlabel(r'Time',fontsize=16)
plt.legend()
#plt.ylim([-1,1])

plt.subplot(212)
plt.title('Electron Spin vs Time')
plt.plot(S_z/S_z[0])
plt.ylabel(r'$<S_z>$ (Electron Spin)',fontsize=16)
plt.xlabel(r'Time',fontsize=16)
plt.legend()
#plt.ylim([-1,1])
plt.tight_layout()
plt.savefig('nuclear_e_spin.png',format='png')
plt.show()

states = output.states
ent = np.zeros(N)
#ent1 = np.zeros(N)
for i in range(N): 
    rho = ket2dm(states[i])
    rho_A = rho.ptrace([0])
    #rho_B = rho.ptrace([1])
    ent[i] = entropy_vn(rho_A)
    #ent1[i] = entropy_vn(rho_B)

#plot of entropy vs time
plt.plot(ent,label= r'Entanglement Entropy')
#plt.plot(ent1,label= r'Entanglement Entropy')
plt.ylabel(r'Tr($\rho_A \log \rho_A )$',fontsize=16)
plt.xlabel('Time',fontsize=16)
plt.ylim([0,1])
plt.legend()
plt.savefig('nuclear_e_entropy.png',format='png')
plt.show()

#fig = plt.figure()
#ax = Axes3D(fig,azim=-40,elev=30)
#sphere = Bloch(axes=ax)
#
#def animate(i):
#    sphere.clear()
#    sphere.vector_color = ['b']
#    sphere.add_vectors([x[i],y[i],z[i]])
#    sphere.make_sphere()
#    return ax
#
#def init():
#    #sphere.vector_color = ['g']
#    return ax
#
#ani = animation.FuncAnimation(fig, animate, np.arange(len(x)),init_func=init, blit=True, repeat=False)
#ani.save('nu_e_exp.mp4', fps=10)

