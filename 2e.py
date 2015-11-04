import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#1e
#a,b,c,d = 0,1,0,0
#a,b,c,d = np.random.random_sample() + np.random.random_sample()*1j,np.random.random_sample() + np.random.random_sample()*1j,np.random.random_sample() + np.random.random_sample()*1j,np.random.random_sample() + np.random.random_sample()*1j
#print a,b,c,d
#psi0 = a * tensor(basis(2,0),basis(2,0)) + b * tensor(basis(2,0),basis(2,1)) +c * tensor(basis(2,1),basis(2,0)) +d * tensor(basis(2,1),basis(2,1)) 
#psi0 = psi0.unit()


#maximally entangled state
#psi0 = tensor(basis(2,1),basis(2,0)) - tensor(basis(2,0),basis(2,1))
#psi0 = tensor(basis(2,1),basis(2,0)) + tensor(basis(2,0),basis(2,1))
#psi0 = psi0.unit()

#pure state z*z
psi0 = tensor(basis(2,1),basis(2,0))

#pure state x*x
#tmp1 = basis(2,0) + basis(2,1)
#tmp1 = tmp1.unit()
#tmp2 = basis(2,0) + basis(2,1)
#tmp2 = tmp1.unit()
#psi0 = tensor(tmp1,tmp2)



B = 0.00
#g = coupling constant
g = 1
#exchange Hamiltonian
H = g * 0.25 * (tensor(sigmaz(),sigmaz()) + 0.5*tensor(sigmap(),sigmam()) + 0.5*tensor(sigmam(),sigmap())) + B * tensor(sigmaz(),qeye(2)) 
print H
print H.eigenstates()
print "psi_0",psi0
#N = number of simulation poitns
N = 100
tlist = np.linspace(0,50,N)
opt = solver.Options(store_states = True)
output = mesolve(H,psi0,tlist,[],[tensor(sigmax(),qeye(2)),tensor(sigmay(),qeye(2)),tensor(sigmaz(),qeye(2)),tensor(qeye(2),sigmax()),tensor(qeye(2),sigmay()),tensor(qeye(2),sigmaz())],options=opt)
x1,y1,z1 = output.expect[0],output.expect[1],output.expect[2]
x2,y2,z2 = output.expect[3],output.expect[4],output.expect[5]

plt.plot(z1)
plt.plot(z2)
plt.show()

states = output.states
ent = np.zeros(N)
for i in range(N): 
    rho = ket2dm(states[i])
    rho_A = rho.ptrace(0)
    ent[i] = entropy_vn(rho_A)
plt.plot(ent,label= r'Entanglement Entropy')
plt.ylabel(r'Tr($\rho_A \log \rho_A )$')
plt.xlabel('Time')
plt.ylim([-1,1])
plt.legend()
plt.show()

#x,y,z = x1,y1,z1
#fig = plt.figure()
#ax = Axes3D(fig,azim=-40,elev=30)
#sphere = Bloch(axes=ax)
#
#def animate(i):
#    sphere.clear()
#    sphere.vector_color = ['b','r']
#    sphere.add_vectors([x1[i],y1[i],z1[i]])
#    sphere.add_vectors([x2[i],y2[i],z2[i]])
#    sphere.make_sphere()
#    return ax
#
#def init():
#    #sphere.vector_color = ['g']
#    return ax
#
#ani = animation.FuncAnimation(fig, animate, np.arange(len(x)),init_func=init, blit=True, repeat=False)
#ani.save('2e_exp.mp4', fps=20)
#
