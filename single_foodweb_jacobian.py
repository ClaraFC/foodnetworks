import numpy as np

food_web_filename = 'Aggregated_Baltic_Ecosystem_Wulff_1989'
crc = 'crc_f1a1c657'
n = 15
l = 12
b0 = np.array([1.066, 0.021418, 0.2102137, 0.06975778, 0.005467852, 0.202513, 3.448311, 0.2791077, 1.817268, 0.09461108, 0.585146, 0.5046653, 10.5455, 265.5, 608.5196])
p = np.array([208.24146, 6.507774, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.654, 2.04624, 0])
q = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.17041477698, 0.004313561616, 0, 0, 1.8116583156])
r = np.array([48.0501, 1.501794, 21.200135292, 29.987636238, 2.0271821976, 59.39471034, 11.252441844, 3.73705542, 9.171678852, 0.5167435644, 1.7694206334, 0.9619239042, 1.454085549, 0, 27.438868926])
F = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24.630357024, 28.06272, 0],
    [29.232, 0, 31.492952694, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 4.069905588, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [52.03296, 0, 0, 15.368661882, 0, 0, 0, 0, 0, 0, 0, 0, 53.10453204, 0, 0],
    [18.27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.049877786, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.24227162],
    [0, 0, 0, 0, 0, 0, 0, 1.1986702182, 2.7577406682, 0, 0, 0, 0, 0, 22.699915938],
    [0, 0, 0, 0, 0, 0, 0.22859920944, 0.18542793024, 0.366867081, 0.13100547348, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0.765882054, 2.0719510056, 0, 0, 1.0125533628, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0.12986696016, 0, 0.5738530266, 0, 0.20287753416, 0.0385738164, 0.44703036, 0, 0, 0, 0],
    [34.63992, 0, 0, 15.368661882, 1.1469730608, 54.96960672, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [26.01648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 5.00598, 0, 0, 0, 0, 13.265276976, 5.12115408, 13.144852098, 0.22562050518, 1.4635209078, 0.4259639538, 23.540277474, 0, 0]
])

#weights from models
sd=0.1
sr=0
sl = 0.1

#Consumption intensities
Cd = F/np.conjugate(b0).T
Cr = F/b0
Cl = F/(b0*b0)

#Setting biomass value
b=b0 #donor-controlled model

#Jacobian computation
J_offdiag= sd*Cd - sr*Cr + sl*(Cl*b) - sl*(Cl.T * b)

J_diag= sd * (np.diag(Cd)- np.diag(sd*np.sum(Cd, axis=0))) + sr * (np.diag(np.sum(Cr,axis=1))- np.diag(Cr))\
    +sl*np.diag(np.tile((np.conj(b).T *Cl),(n,1)))-sl* np.diag(np.tile((b*Cl),(n,1)).T)\
    -sl* np.diag(q/b + r/b)

J = J_offdiag - np.diag(J_offdiag) + J_diag

