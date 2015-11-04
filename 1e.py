import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#1 electron Hamiltonian
#examples of initial states
# 0.707 * (1,1)
psi0 = basis(2,0) + basis(2,1)
psi0 = psi0.unit()

# (0,1)
#psi0 = basis(2,1)

# random ket
#psi0 = rand_ket(2)

#non zero B produces Lamor precession
B = -4

H = sigmaz()*sigmaz() + sigmax()*sigmax() + sigmay()*sigmay() + B*sigmaz()

#N = number of simulation points
N = 25
tlist = np.linspace(0,20,N)
opt = solver.Options(store_states = True)
output = mesolve(H,psi0,tlist,[],[sigmax(),sigmay(),sigmaz()],options=opt)
x,y,z = output.expect[0],output.expect[1],output.expect[2]

states = output.states
ent = np.zeros(N)
for i in range(N): 
    rho = ket2dm(states[i])
    ent[i] = entropy_vn(rho)

plt.plot(ent)
plt.plot(ent,label= r'Von Neumann Entropy')
plt.ylabel(r'Tr($\rho \log \rho )$')
plt.xlabel('Time')
plt.ylim([-1,1])
plt.legend()
plt.show()

fig = plt.figure()
ax = Axes3D(fig,azim=-40,elev=30)
sphere = Bloch(axes=ax)

def animate_exp(i):
    sphere.clear()
    sphere.add_vectors([x[i],y[i],z[i]])
    sphere.make_sphere()
    return ax

def animate_states(i):
    sphere.clear()
    sphere.add_states(states[i])
    sphere.make_sphere()
    return ax
def init():
    sphere.vector_color = ['g']
    return ax

ani = animation.FuncAnimation(fig, animate_exp, np.arange(len(x)),init_func=init, blit=True, repeat=False)
ani.save('1e_exp.mp4', fps=10)
ani = animation.FuncAnimation(fig, animate_states, np.arange(len(x)),init_func=init, blit=True, repeat=False)
ani.save('1e_states.mp4', fps=10)

