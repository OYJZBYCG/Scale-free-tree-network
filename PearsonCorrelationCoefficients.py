import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams

config = {
    "font.family": 'serif',
    "font.size": 14,
    "mathtext.fontset": 'stix',
    "font.serif": ['SimSun'],
}
rcParams.update(config)


# Output directory
outdir = input(r'Please enter the image output path (e.g. D:\\Users\\test):')

def F(t):
    temp = 0
    for i in range(1, t-1):
        temp += 2**(2*i+1)*(5**(t-2-i)+5**(t-1-i))
    return 4**t+temp+2*5**(t-2)

def Gamma_1(t):
    return 6*F(t)+4*5**(t-1)

def Gamma_2(t):
    return 4*F(t)+2*5**t

def Gamma_3(t):
    temp = 0
    for i in range(1, t-1):
        temp += 2**(3*i+1)*(5**(t-2-i)+5**(t-1-i))
    return 8*(4*8**(t-1)+temp+2*5**(t-2))+18*5**(t-1)

def r(t):
    return (Gamma_1(t)/5**t-(Gamma_2(t)/(2*5**t))**2)/(Gamma_3(t)/(2*5**t)-(Gamma_2(t)/(2*5**t))**2)

def GammaStar_1(t):
    temp1 = 0
    for i in range(t):
        temp2 = 0
        for j in range(t - 1 - i):
            temp2 += 9 ** (t - 2 - j) / 5 ** (t - j) * 9 ** i
        temp1 += 2 * 9 ** (t - 1) / 5 ** t + 4 * temp2
    return 9 ** t / 5 ** t + 6 * temp1

def GammaStar_2(t):
    temp1 = 0
    for i in range(t):
        temp2 = 0
        for j in range(t - 1 - i):
            temp2 += 9 ** (t - 2 - j) / 5 ** (t - j)
        temp1 += 3 ** i * (2 * 9 ** (t - 1 - i) / 5 ** t + 4 * temp2)
    temp3 = 0
    for i in range(t):
        temp3 += 3 ** i / 5 ** (1 + i)
    return 3 ** t / 5 ** t + 3 * temp1 + 2 * temp3

def GammaStar_3(t):
    temp1 = 0
    for i in range(t):
        temp2 = 0
        for j in range(t - 1 - i):
            temp2 += 27 ** (t - 2 - j) / 5 ** (t - j)
        temp1 += 9 ** i * (2 * 27 ** (t - 1 - i) / 5 ** t + 4 * temp2)
    temp3 = 0
    for i in range(t):
        temp3 += 9 ** i / 5 ** (1 + i)
    return 9 ** t / 5 ** t + 9 * temp1 + 2 * temp3

def rStar(t):
    return (GammaStar_1(t) - GammaStar_2(t) ** 2) / (GammaStar_3(t) - GammaStar_2(t) ** 2)


t = np.arange(2, 11, dtype='int64')
r1 = []
r2 = []

for i in t:
    r1.append(r(i))
    r2.append(rStar(i))

# print(r1)
# print(r2)

fig = plt.figure(figsize=(6, 4), dpi=300)
plt.plot(t, r1, color='r', marker='v', label=r'$ r(t) $', linewidth=1.5)
plt.plot(t, r2, color='b', marker='o', label=r'$ r^\star(t) $', linewidth=1.5)
plt.xticks(np.arange(2, 11), fontname='Times New Roman')
plt.yticks(np.arange(-0.7, 0.2, 0.1), fontname='times new roman')
plt.grid(axis='y', linestyle='--')
plt.xlabel(r'$t$')
plt.ylabel(r'$r$')
plt.legend()
plt.savefig(outdir+'\\r&rStar.png')
plt.savefig(outdir+'\\r&rStar.jpg')

plt.show()