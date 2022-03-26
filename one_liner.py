# importing libraries
import pandas as pd
import chemparse
import matplotlib.pyplot as plt

# defining the periodic table
mass_number = {}
with open("neo p_table.txt") as f:
    for line in f:
        (key, val) = line.split(    )
        mass_number[key] = float(val)
periodic_table = list(mass_number)

# used for recognizing letters
alphabe = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ()'

# intro lines
print('Please enter your mixture like this:')
print("e.g. 40H2O+60CO2")
print('---------------------------------------------------')

# extracting mixture from input (and the % of compounds)
a = input('enter mixture below:\n').split('+')
mixture = []
prc = []
for i in a:
    j = i.lstrip('0123456789.')
    k = i[:3].rstrip(alphabe)
    #print(j)
    #print(k)
    mixture.append(j)
    prc.append(float(k))
print()
print('mixture: ', mixture)
print('%: ', prc)

def formula(n):
    return chemparse.parse_formula(n)

guns = []
bullets = []
for c in range(0, len(mixture)):
    formula_1 = formula(mixture[c])
    comp_const = list(formula_1)
    for i in comp_const:
        guns.append(i)
    for i in comp_const:
        bullets.append(formula_1[i] * prc[c])
    
# finding element ratios
elements = []
occ = []
j = 0
for i in guns:
    if i in elements:
        occ[elements.index(i)] = occ[elements.index(i)] + bullets[j]
        j += 1
        continue
    elements.append(i)
    occ.append(bullets[j])
    j += 1
print('elements: ', elements)
#print('numbers :', occ)

def f(n):
    return occ[n]/sum(occ)
print()
print('element ratios:')
for i in range(0, len(elements)):
    print('{:s} ---> {:.4f}'.format(elements[i], f(i)))

# findind atomic numbers
Z = []
for i in range(0, len(elements)):
    Z.append(periodic_table.index(elements[i]) + 1)
    i =+ 1

# finding energies
print()
b = input('choose energies in MeV below: (e.g. 0.02 0.03)\n').split()
energies = []
for i in b:
    energies.append(float(i))
#print()
#print('energies = ', energies) 

# output decoration
print('-----------------------------------\n')
print('\nEnergy (MeV) \t\t\t Z-eff\n')

zeff = []
# finding Z-eff for each given energy
for k in energies:
    u = []
# accessing the files
    for i in Z:
        websites = pd.read_csv("info/{:s}.txt".format(str(i)), 
                               delimiter = '\t', header = None)
        websites.columns = ['Photon-Energy (MeV)', 'Coh. Scatt.', 
                            'incoh. scatt.', 'photoelectric', 'Nuclear PP',
                            'Electron PP','Tot. w/', 'Tot. w/o', 'extra slot']
        websites.to_csv('test.csv', index = None)
        websites_output = pd.read_csv('test.csv')
        photon_str = list(websites['Photon-Energy (MeV)'][2:])
        photon = [float(x) for x in photon_str]
        MAC_str = list(websites['Tot. w/'][2:])
        MAC = [float(x) for x in MAC_str]   
# finding miu
        index = photon.index(k)
        miu = MAC[index]
        u.append(miu)
# calculation   
    num = []
    den = []
    
    for n in range(0, len(elements)):
        a = f(n) * mass_number[elements[n]] * u[n]
        num.append(a)
    for m in range(0, len(elements)):
        b = f(m) * (mass_number[elements[m]]/Z[m]) * u[m]
        den.append(b)
    Z_effective = sum(num)/sum(den)
    zeff.append(Z_effective)
    print(k, '\t\t\t', Z_effective)
    
#print(energies)
#print(zeff)    
plt.plot(energies, zeff)
plt.show()