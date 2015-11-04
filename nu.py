import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#psi_nu = tensor(rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit(),rand_ket(2).unit())
psi_nu = rand_ket(64)
#psi_nu = tensor(basis(2,1),basis(2,1),basis(2,1),basis(2,1),basis(2,1),basis(2,1))

psi_e = basis(2,0)
psi0 = tensor(psi_e,psi_nu)
psi0 = psi0.unit()

#exchange Hamiltonian
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
    
S_z = output.expect[0]
I_exp = output.expect[1]

x,y,z = output.expect[2],output.expect[3],output.expect[0]

plt.figure()
plt.subplot(211)
plt.plot(I_exp,label='Collective Nuclear Spin')
plt.ylabel(r'$<I_z>$',fontsize=16)
plt.xlabel(r'Time')
plt.legend()
#plt.ylim([-1,1])

plt.subplot(212)
plt.plot(S_z/S_z[0],label='Electron Spin')
plt.ylabel(r'$<S_z>$',fontsize=16)
plt.xlabel(r'Time')
plt.legend()
#plt.ylim([-1,1])
plt.tight_layout()
plt.show()

states = output.states
ent = np.zeros(N)
for i in range(N): 
    rho = ket2dm(states[i])
    rho_A = rho.ptrace([1])
    print rho_A
    ent[i] = entropy_vn(rho_A)


plt.plot(ent,label= r'Entanglement Entropy')
plt.ylabel(r'Tr($\rho_A \log \rho_A )$')
plt.xlabel('Time')
plt.ylim([-1,1])
plt.legend()
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
