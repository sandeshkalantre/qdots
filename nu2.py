import numpy as np 
from qutip import *
import matplotlib.pyplot as plt

psi_nu = tensor(basis(2,1),basis(2,1))
psi_e = basis(2,0)
psi0 = tensor(psi_e,psi_nu)

I_z = tensor(sigmaz(),qeye(2)) + tensor(qeye(2),sigmaz())
I_p = tensor(sigmap(),qeye(2)) + tensor(qeye(2),sigmap())
I_m = tensor(sigmam(),qeye(2)) + tensor(qeye(2),sigmam())

H = tensor(sigmaz(),I_z) + 0.5*tensor(sigmap(),I_m) + 0.5*tensor(sigmam(),I_p)

tlist = np.linspace(0,10,100)
output = mesolve(H,psi0,tlist,[],[tensor(sigmaz()),qeye(4),tensor(qeye(2),I_z)],options=opt)
