### Calculates the rotational energy levels and 
### hyperfine splitting of HCN, assuming the rigid-rotor
### approximation

from astropy import units as u
from astropy.units import cds
import numpy as np
import matplotlib.pyplot as plt

# Import some constants
h = 6.62607e-34*u.J*u.s
k = 1.38065e-16*u.erg/u.K
#print k.value

# masses and bond distances
mh = 1.67e-24*u.g
mc = 1.99e-23*u.g
mn = 2.33e-23*u.g
rhc = 1.068e-10*u.m # From Herzberg 1966
rcn = 1.156e-10*u.m
rhn = rhc+rcn

# moment of inertia
I = (mh*mc*rhc**2 + mh*mn*rhn**2 + mn*mc*rcn**2)/(mn+mc+mh)
I = I.cgs
#print I

# Calculate first 15 rotational levels
J=np.arange(0,16)
E = (h**2 / ((8.*np.pi**2)*I)) * J * (J+1)
E = (E.cgs).value / k.value * u.K
print E

# Calculate the transition frequencies 
B = (h.cgs / ((8.*np.pi**2)*I))
v = 2*B*(J+1)
v = v.to(u.GHz)
print v

xx = []
for i in range(len(E)-1):
	x = E[i+1]-E[i]
	xx.append(x.value)
print xx

# Plot the transition frequencies versus 
# the change in energy between each pair of levels
fig = plt.figure(figsize=(5,7))
plt.scatter(xx, v[0:15])
plt.xlabel('$\Delta$E (J+1)->J [K]')
plt.ylabel('Transition Frequency [GHz]')
fig.savefig('freq.pdf', bbox_inches='tight')
plt.show()

# Make labels for the energy levels
labels = []
for i in J:
	l = 'J = '+str(i)
	labels.append(l)

# Plot the energy levels with labels
fig = plt.figure(figsize=(5,7))
plt.scatter(np.zeros(len(E)), E, marker='_', s=20000)
c = 0
for lab, x,y in zip(labels, np.zeros(len(E)), E):
	if c<3:
		s = 12
	else:
		s = 12		
	plt.annotate(lab, xy=(x, y.value), xytext=(80,-3), textcoords='offset points', size=s)
	c+=1
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.ylabel('Energy [K]')
plt.xlabel('Levels')
plt.ylim(-20,max(E.value)+20)
fig.savefig('levels.pdf', bbox_inches='tight')
plt.show()

# Calculate hyperfine splitting of J=1 level and plot
Q = -4.58*u.MHz
J = 1.
I = 1.
F = np.array([2.,1.,0.])
C = F*(F+1.) - I*(I+1.) - J*(J+1.)
de = -1.*(0.75*C*(C+1.) - I*(I+1.)*J*(J+1.)) / ( 2.*I*(2.*I - 1.)*(2.*J - 1.)*(2.*J + 3.) )
dE = Q*de
print dE.to(u.kHz)

labels = []
for i in F:
	l = 'F = '+str(int(i))
	labels.append(l)
fig = plt.figure()
plt.scatter(np.ones(len(dE)), dE.to(u.kHz), marker='_', s=20000, color='red')
plt.scatter([0.], [0.], marker='_', s=20000, color='blue')
plt.annotate('J = 1', xy=(0, 0), xytext=(-110,-3), textcoords='offset points', size=14, color='blue')
for lab, x,y in zip(labels, np.ones(len(dE)), dE.to(u.kHz)):
	plt.annotate(lab, xy=(x, y.value), xytext=(80,-3), textcoords='offset points', size=14, color='red')
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.xlim(-0.6,1.6)
plt.ylim(-2500, 2500)
plt.ylabel('$\Delta$E / h [kHz]')
fig.savefig('hyperfine1.pdf', bbox_inches='tight')
plt.show()

# Calculate hyperfine splitting of J=2 level and plot
Q = -4.58*u.MHz
J = 2.
I = 1.
F = np.array([3.,2.,1.])
C = F*(F+1.) - I*(I+1.) - J*(J+1.)
de = -1.*(0.75*C*(C+1.) - I*(I+1.)*J*(J+1.)) / ( 2.*I*(2.*I - 1.)*(2.*J - 1.)*(2.*J + 3.) )
dE = Q*de
print dE.to(u.kHz)

labels = []
for i in F:
	l = 'F = '+str(int(i))
	labels.append(l)
fig = plt.figure()
plt.scatter(np.ones(len(dE)), dE.to(u.kHz), marker='_', s=20000, color='red')
plt.scatter([0.], [0.], marker='_', s=20000, color='blue')
plt.annotate('J = 2', xy=(0, 0), xytext=(-110,-3), textcoords='offset points', size=14, color='blue')
for lab, x,y in zip(labels, np.ones(len(dE)), dE.to(u.kHz)):
	plt.annotate(lab, xy=(x, y.value), xytext=(80,-3), textcoords='offset points', size=14, color='red')
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.xlim(-0.6,1.6)
plt.ylim(-2500, 2500)
plt.ylabel('$\Delta$E / h [kHz]')
fig.savefig('hyperfine2.pdf', bbox_inches='tight')
plt.show()

# Calculate hyperfine splitting of J=3 level and plot
Q = -4.58*u.MHz
J = 3.
I = 1.
F = np.array([4.,3.,2.])
C = F*(F+1.) - I*(I+1.) - J*(J+1.)
de = -1.*(0.75*C*(C+1.) - I*(I+1.)*J*(J+1.)) / ( 2.*I*(2.*I - 1.)*(2.*J - 1.)*(2.*J + 3.) )
dE = Q*de
print dE.to(u.kHz)

labels = []
for i in F:
	l = 'F = '+str(int(i))
	labels.append(l)
fig = plt.figure()
plt.scatter(np.ones(len(dE)), dE.to(u.kHz), marker='_', s=20000, color='red')
plt.scatter([0.], [0.], marker='_', s=20000, color='blue')
plt.annotate('J = 3', xy=(0, 0), xytext=(-110,-3), textcoords='offset points', size=14, color='blue')
for lab, x,y in zip(labels, np.ones(len(dE)), dE.to(u.kHz)):
	plt.annotate(lab, xy=(x, y.value), xytext=(80,-3), textcoords='offset points', size=14, color='red')
plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.xlim(-0.6,1.6)
plt.ylim(-2500, 2500)
plt.ylabel('$\Delta$E / h [kHz]')
fig.savefig('hyperfine3.pdf', bbox_inches='tight')
plt.show()

