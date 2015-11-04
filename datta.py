import numpy as np 
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

e1 = 0.1
e2 = 0.3
U = 0.0
t = 0.1
H = np.array([[2*e1 - U,0.0,t, t,0.0,0.0],[0,2*e2 + U,t,t,0,0],[t,t,e1+e2,0,0,0],[t,t,0,e1+e2,0,0],[0,0,0,0,e1+e2,0],[0,0,0,0,0,e1+e2]])
print H
H1 = Qobj(H)
print H1.eigenstates()



