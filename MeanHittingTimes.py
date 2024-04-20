import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams

config = {
    "font.family": 'serif',
    "font.size": 16,
    "mathtext.fontset": 'stix',
    "font.serif": ['SimSun'],
    "axes.unicode_minus": False,
}
rcParams.update(config)


#################### Parameters ####################
outdir = input(r'Please enter the image output path (e.g. D:\\Users\\test):')      # Output directory
T0 = int(input('Please enter T0:'))
W0 = int(input('Please enter W0:'))


def H(t):      # H_t
    return 2*(15**t*W0-4*((15**t-5**t)*(T0-1)**2/10+(15**t-1)*(T0-1)/14))/(5**t*(T0-1)+1)*5**t

def HStar(t):      # H^\star_t
    return 2*(5**(t+1)*W0+10*t*5**(t-1)*(T0-1)**2-3*(5**t-1)*(T0-1)/2)/(5**t*(T0-1)+1)*5**(t-1)

def f2(m):
    return m * (2 * m + 1)

def f3(m):
    return m * (5 * m + 3)

def f4(m):
    return m * (3 * m + 2)

def HStarm(m, t):      # H^\star_{m;t}
    temp1 = 0
    for i in range(t):
        temp1 += (((T0 - 1) * (2 * m + 1) ** (t - 1 - i) + 1) / (2 * m + 1) ** (t - i)) ** 2
    temp2 = 0
    for i in range(t):
        temp2 += ((T0 - 1) * (2 * m + 1) ** (t - 1 - i) + 1) / (2 * m + 1) ** (t - i) / (2 * m + 1) ** (t - i)
    temp3 = 0
    for i in range(t):
        temp3 += ((2 * m + 1) ** i / (2 * m + 1) ** t) ** 2
    return 2 * (W0 + f2(m) * temp1 - f3(m) * temp2 + f4(m) * temp3) * (2 * m + 1) ** t / (
                (T0 - 1) * (2 * m + 1) ** t + 1) * (2 * m + 1) ** t

def f5(m):
    return 3*(m+1)**2
    
def g5(m):
    return (m-1)*(m+2)

def h5(m):
    return m+2

def TV(m, t):
    return T0*(m+1)**t

def HV(m, t):      # H^V_{m;t}
    return 2*(f5(m)**t/TV(m, t)*W0+g5(m)*sum([f5(m)**i*TV(m, t-1-i)**2/TV(m, t) for i in range(t)])+
                   h5(m)*sum([f5(m)**i*TV(m, t-1-i)/TV(m, t) for i in range(t)]))

def f_4(m):
    return 2*(m+2)**2

def g4(m):
    return m+2

def h4(m):
    return (m-1)*(m+2)

def l4(m):
    return m**2+2*m

def TT(m, t):
    return (T0-1)*(m+2)**t+1

def HT(m, t):
    return 2*(f_4(m)**t*W0/TT(m, t)-g4(m)*sum([f_4(m)**i*TT(m, t-1-i)**2/TT(m, t) for i in range(t)])-
                   h4(m)*sum([f_4(m)**i*TT(m, t-1-i)/TT(m, t) for i in range(t)])+l4(m)*sum([f_4(m)**i/TT(m, t) for i in range(t)]))

#################### H_t & H^\star_t ####################
t = np.arange(1, 11, dtype = 'int64')
fig = plt.figure(figsize = (15, 10), dpi = 300)
fig.subplots_adjust(wspace = 0.25, hspace = 0.4)

h = H(t)
hs = HStar(t)
# print(h1)
# print(hs1)
fig1 = plt.figure(figsize = (8, 6), dpi = 300)
plt.subplots_adjust(bottom = 0.2)
plt.semilogy(t, h, c = 'r', marker = 'v', label = r'${\langle \mathcal{H}_t \rangle}$')
plt.semilogy(t, hs, c = 'b', marker = 'o', label = r'${\langle \mathcal{H}_t^\star \rangle}$')
plt.xticks(np.arange(1, 11), fontname = 'Times New Roman')
plt.yticks(fontname = 'Times New Roman')
plt.title(r'${|\mathcal{T}_0|=%d, \mathcal{W}_0=%d}$'%(T0, W0), y = -0.22)
plt.xlabel(r'$t$')
plt.ylabel(r'${\langle \mathcal{H} \rangle}$')
plt.legend()
plt.grid(axis = 'y', linestyle = '--')
plt.savefig(outdir+'\\H_t&HStar_t(T0=%d,W0=%d).png'%(T0, W0))
plt.savefig(outdir+'\\H_t&HStar_t(T0=%d,W0=%d).jpg'%(T0, W0))

#################### H^\star_{m;t} ####################
t = np.arange(1, 11, dtype='int64')
m = np.arange(1, 6, dtype='int64')
tt, mm = np.meshgrid(t, m)
r, c = tt.shape

hsm = []
for i in range(r):
    tmp = []
    for j in range(c):
        tmp.append(HStarm(tt[i][j], mm[i][j]))
    hsm.append(tmp)
hsm = np.array(hsm, dtype='int64')
fig2 = plt.figure(figsize=(6, 6), dpi=300)
ax = plt.axes(projection='3d')
ax.plot_surface(tt, mm, np.log(hsm), rstride=1, cstride=1, cmap='rainbow')
ax.set_xticks(np.arange(1, 11))
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$m$')
ax.set_zlabel(r'${\ln \langle \mathcal{H}_{m;t}^\star \rangle}$')
ax.set_title(r'${|\mathcal{T}_0|=%d, \mathcal{W}_0=%d}$'%(T0, W0), y=-0.12)
plt.savefig(outdir+'\\HStar_{m;t}(T0=%d,W0=%d).png'%(T0, W0))
plt.savefig(outdir+'\\HStar_{m;t}(T0=%d,W0=%d).jpg'%(T0, W0))

#################### H^V_{m;t} & H^T_{m;t} ####################
t = np.arange(1, 10, dtype='int64')
m = np.arange(1, 6, dtype='int64')
tt, mm = np.meshgrid(t, m)
r, c = tt.shape

# H^V_{m;t}
hv = np.array([[HV(mm[i][j], tt[i][j]) for j in range(c)] for i in range(r)], dtype = 'float64')
# print(hv)
fig3 = plt.figure(figsize=(6, 6), dpi=300)
ax = plt.axes(projection='3d')
ax.plot_surface(tt, mm, np.log(hv), rstride=1, cstride=1, cmap='rainbow')
# ax.view_init(elev=30, azim=-70)
ax.set_xticks(np.arange(1, 10))
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$m$')
ax.set_zlabel(r'${\ln \langle \mathcal{H}_{m;t}^V \rangle}$')
ax.set_title(r'${|\mathcal{T}_0|=%d, \mathcal{W}_0=%d}$'%(T0, W0), y=-0.12)
plt.savefig(outdir+'\\HV(T0=%d,W0=%d).png'%(T0, W0))
plt.savefig(outdir+'\\HV(T0=%d,W0=%d).jpg'%(T0, W0))

# H^T_{m;t}
ht = np.array([[HT(mm[i][j], tt[i][j]) for j in range(c)] for i in range(r)], dtype = 'float64')
# print(ht)
fig4 = plt.figure(figsize=(6, 6), dpi=300)
ax = plt.axes(projection='3d')
ax.plot_surface(tt, mm, np.log(ht), rstride=1, cstride=1, cmap='rainbow')
# ax.view_init(elev=30, azim=-70)
ax.set_xticks(np.arange(1, 10))
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$m$')
ax.set_zlabel(r'${\ln \langle \mathcal{H}_{m;t}^T \rangle}$')
ax.set_title(r'${|\mathcal{T}_0|=%d, \mathcal{W}_0=%d}$'%(T0, W0), y=-0.12)
plt.savefig(outdir+'\\HT(T0=%d,W0=%d).png'%(T0, W0))
plt.savefig(outdir+'\\HT(T0=%d,W0=%d).jpg'%(T0, W0))

plt.show()
