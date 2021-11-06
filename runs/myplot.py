import matplotlib.pyplot as plt

X, L, dL = [], [], []
for line in open('LEGENDRE.dat', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  L.append(values[1])
  dL.append(values[2])

plt.plot(X, L, X, dL)
plt.show()
