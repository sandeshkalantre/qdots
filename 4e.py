import numpy as np
from qutip import *

H1 = tensor(sigmaz(),sigmaz(),qeye(2),qeye(2)) + 0.5* ( tensor(sigmap(),sigmam(),qeye(2),qeye(2)) + tensor(sigmam(),sigmap(),qeye(2),qeye(2)))
H2 = tensor(qeye(2),sigmaz(),sigmaz(),qeye(2)) + 0.5 * ( tensor(qeye(2),sigmap(),sigmam(),qeye(2)) + tensor(qeye(2),sigmam(),sigmap(),qeye(2)))
H3 = tensor(qeye(2),qeye(2),sigmaz(),sigmaz()) + 0.5 * (tensor(qeye(2),qeye(2),sigmap(),sigmam()) + tensor(qeye(2),qeye(2),sigmam(),sigmap()))

H = H1 + H2 + H3

print H.eigenstates()

