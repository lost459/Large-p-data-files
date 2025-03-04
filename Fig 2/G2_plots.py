import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate as intg

plt.rcParams.update({"mathtext.fontset": 'cm',
                     'font.family':'times new roman'})

plt.rc('font', size=9) #controls default text size
plt.rc('axes', titlesize=9) #fontsize of the title
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rc('xtick', labelsize=9) #fontsize of the x tick labels
plt.rc('ytick', labelsize=9) #fontsize of the y tick labels
plt.rc('legend', fontsize=9) #fontsize of the legend
plt.rc('lines', linewidth=1) #linewidth of plot lines
plt.rc('lines', markersize=3.75) #size of datapoint
plt.rc('axes', labelpad=0.75)

###############################################################################

# load data, rows sorted according to G2, D, p, lam, q, t.
data = np.loadtxt(r'/Users/s5257291/Documents/Large_p_Paper/2024/Data Files/Fig 2/G2_large.txt', delimiter='\t', skiprows=1)

data = data[data[:,3].argsort()] # Sort according to lambda

q_data = data[0:202]
q_data = q_data[q_data[:,2].argsort()] # Sort q data according to p
q_lp = q_data[0:101]
q_hp = q_data[101:202]
q_lp = q_lp[q_lp[:,5].argsort()] # Sort q data according to taus
q_hp = q_hp[q_hp[:,5].argsort()] # Sort q data according to taus

lam_data = data[202:404]
lam_data = lam_data[lam_data[:,2].argsort()] # Sort lam data according to p
lam_lp = lam_data[0:101]
lam_hp = lam_data[101:202]
lam_lp = lam_lp[lam_lp[:,5].argsort()] # Sort lam data according to taus
lam_hp = lam_hp[lam_hp[:,5].argsort()] # Sort lam data according to taus

taus_lp = lam_lp[:,5]
taus_hp = lam_hp[:,5]

g2lam_lp = lam_lp[:,0]
g2lam_hp = lam_hp[:,0]
g2q_lp = q_lp[:,0]
g2q_hp = q_hp[:,0]

###############################################################################

def g2lam_an(mu,p,lam,ts):
    k = np.pi**2 * p / (4 * mu**2)
    g2 = 1 + lam * (lam - 1) * k * np.exp(- k * ts)
    return g2

def g2q_an(mu,p,q,ts):
    k = np.pi**2 * p / (4 * mu**2)
    g2 = 1 + q * (1/2 + q/4) * k * np.exp(- k * (1 + q/2) * ts)
    return g2

mu = 250
lp = 4.15
k_lp = (np.pi**2*lp)/(4*mu**2)
hp = 50
k_hp = (np.pi**2*hp)/(4*mu**2)
lam = 0.5
q = -1

ts_lp = np.linspace(0,taus_lp[-1],1000)
ts_hp = np.linspace(0,taus_hp[-1],1000)
g2la_lp = g2lam_an(mu,lp,0.5,ts_lp)
g2qa_lp = g2q_an(mu,lp,-1,ts_lp)
g2la_hp = g2lam_an(mu,hp,0.5,ts_hp)
g2qa_hp = g2q_an(mu,hp,-1,ts_hp)

###############################################################################

w_rnge_lp = 6e-4
w_rnge_hp = 6e-3

# Numerical Intensity Fluctuation Spectrum
w_lp = np.linspace(-w_rnge_lp, w_rnge_lp, len(g2lam_lp))
Swlam_lp = np.zeros(len(g2lam_lp))
Swq_lp = np.zeros(len(g2lam_lp))
w_hp = np.linspace(-w_rnge_hp, w_rnge_hp, len(g2lam_hp))
Swlam_hp = np.zeros(len(g2lam_hp))
Swq_hp = np.zeros(len(g2lam_hp))

for i in range(len(g2lam_lp)):
    Swlam_lp[i] = 1 + 2*np.real(intg.trapz(np.exp(w_lp[i]*taus_lp*1j)*(g2lam_lp-1),taus_lp))
    Swq_lp[i] = 1 + 2*np.real(intg.trapz(np.exp(w_lp[i]*taus_lp*1j)*(g2q_lp-1),taus_lp))
    Swlam_hp[i] = 1 + 2*np.real(intg.trapz(np.exp(w_hp[i]*taus_hp*1j)*(g2lam_hp-1),taus_hp))
    Swq_hp[i] = 1 + 2*np.real(intg.trapz(np.exp(w_hp[i]*taus_hp*1j)*(g2q_hp-1),taus_hp))
    
# Analytical Intensity Fluctuation Spectrum
w_an_lp = np.linspace(-w_rnge_lp, w_rnge_lp, 101)
Swlam_an_lp = 1 + lam*(lam - 1)*2*k_lp**2/(k_lp**2 + w_an_lp**2)
Swq_an_lp = 1 + q*(k_lp*(1+q/2))**2/((k_lp*(1+q/2))**2 + w_an_lp**2)
w_an_hp = np.linspace(-w_rnge_hp, w_rnge_hp, 101)
Swlam_an_hp = 1 + lam*(lam - 1)*2*k_hp**2/(k_hp**2 + w_an_hp**2)
Swq_an_hp = 1 + q*(k_hp*(1+q/2))**2/((k_hp*(1+q/2))**2 + w_an_hp**2)
    
###############################################################################

fig = plt.figure(figsize = (1*(3+3/8),2*(3+3/8)))

m = 0.1 # amount to take off axes width, in units of figure width
axs1 = fig.add_axes([0 + m, 3/4 + m, 1 - m, 1/4 - m]) # [left, bottom, width, height]

axs1.plot(taus_lp,g2q_lp-1, 'C3+', markersize = 5.0, label = r'$p,q$-Family Numerical')
axs1.plot(taus_lp,g2lam_lp-1, 'C2.', label = r'$p,\lambda$-Family Numerical')
axs1.plot(ts_lp,g2qa_lp-1, color = [0.25,0.25,0.25], label = r'$p,q$-Family Analytical')
axs1.plot(ts_lp,g2la_lp-1, '--', color = [0.25,0.25,0.25], label = r'$p,\lambda$-Family Analytical')
axs1.set_xlim([0 - 0.025*6e4,6e4])
axs1.grid('on')
axs1.set_ylim([-6e-5, 0 + 2.5e-6])
axs1.legend(labelspacing = 0.15)
axs1.set_xlabel(r'$\mathcal{N}t$')
axs1.set_ylabel(r'$g^{(2)}_{\rm ps}(t) - 1$')
axs1.ticklabel_format(axis = 'both', style='sci', scilimits=(0,0))
axs1.yaxis.get_offset_text().set_visible(False)
axs1.text(0.0, 1.025, r'$\times10^{-5}$', transform=axs1.transAxes)
axs1.xaxis.get_offset_text().set_visible(False)
axs1.text(0.9, -0.35, r'$\times10^{4}$', transform=axs1.transAxes)
# axs1.text(0.56, 0.45, r'$\mu = 250$' '\n' r'$p = 50$', transform=axs1.transAxes)

axs2 = fig.add_axes([0 + m, 1/4 + m, 1 - m, 1/4 - m]) # [left, bottom, width, height]

axs2.plot(taus_hp,g2q_hp-1, 'C3+', markersize = 5.0,)
axs2.plot(taus_hp,g2lam_hp-1, 'C2.')
axs2.plot(ts_hp,g2qa_hp-1, color = [0.25,0.25,0.25])
axs2.plot(ts_hp,g2la_hp-1, '--', color = [0.25,0.25,0.25])
axs2.set_xlim([0 - 0.025*6e3,6e3])
axs2.grid('on')
axs2.set_ylim([-6e-4, 0 + 2.5e-5])
axs2.set_xlabel(r'$\mathcal{N}t$')
axs2.set_ylabel(r'$g^{(2)}_{\rm ps}(t) - 1$')
axs2.ticklabel_format(axis = 'both', style='sci', scilimits=(0,0))
axs2.yaxis.get_offset_text().set_visible(False)
axs2.text(0.0, 1.025, r'$\times10^{-4}$', transform=axs2.transAxes)
axs2.xaxis.get_offset_text().set_visible(False)
axs2.text(0.9, -0.35, r'$\times10^{3}$', transform=axs2.transAxes)

axs3 = fig.add_axes([0 + m, 2/4 + m, 1 - m, 1/4 - m]) # [left, bottom, width, height]

axs3.plot(w_lp,Swq_lp, 'C3+', markersize = 5.0,)
axs3.plot(w_lp,Swlam_lp, 'C2.')
axs3.plot(w_an_lp,Swq_an_lp, color = [0.25,0.25,0.25])
axs3.plot(w_an_lp,Swlam_an_lp, '--', color = [0.25,0.25,0.25])
axs3.set_xlim([w_an_lp[0],w_an_lp[-1]])
axs3.set_ylim([-0.05,1.05])
axs3.grid('on')
axs3.set_xlabel(r'$\omega /\mathcal{N}$')
axs3.set_ylabel(r'$S(\omega)/\mathcal{N}$')
axs3.ticklabel_format(axis = 'x', style='sci', scilimits=(0,0))
axs3.xaxis.get_offset_text().set_visible(False)
axs3.text(0.9, -0.35, r'$\times10^{-4}}$', transform=axs3.transAxes)

axs4 = fig.add_axes([0 + m, 0 + m, 1 - m, 1/4 - m]) # [left, bottom, width, height]

axs4.plot(w_hp,Swq_hp, 'C3+', markersize = 5.0,)
axs4.plot(w_hp,Swlam_hp, 'C2.')
axs4.plot(w_an_hp,Swq_an_hp, color = [0.25,0.25,0.25])
axs4.plot(w_an_hp,Swlam_an_hp, '--', color = [0.25,0.25,0.25])
axs4.set_xlim([w_an_hp[0],w_an_hp[-1]])
axs4.set_ylim([-0.05,1.05])
axs4.grid('on')
axs4.set_xlabel(r'$\omega /\mathcal{N}$')
axs4.set_ylabel(r'$S(\omega)/\mathcal{N}$')
axs4.ticklabel_format(axis = 'x', style='sci', scilimits=(0,0))
axs4.xaxis.get_offset_text().set_visible(False)
axs4.text(0.9, -0.35, r'$\times10^{-3}}$', transform=axs4.transAxes)

axs1.text(0.025, 0.85, r'(a)', transform=axs1.transAxes)
axs1.text(0.2, 0.2, r'$\mu = 250$' '\n' r'$p = 4.15$', transform=axs1.transAxes)
axs2.text(0.025, 0.85, r'(c)', transform=axs2.transAxes)
axs2.text(0.2, 0.2, r'$\mu = 250$' '\n' r'$p = 50$', transform=axs2.transAxes)
axs3.text(0.025, 0.75, r'(b)', transform=axs3.transAxes)
axs3.text(0.2, 0.2, r'$\mu = 250$' '\n' r'$p = 4.15$', transform=axs3.transAxes)
axs4.text(0.025, 0.75, r'(d)', transform=axs4.transAxes)
axs4.text(0.2, 0.2, r'$\mu = 250$' '\n' r'$p = 50$', transform=axs4.transAxes)

# fig.savefig('G2_and_Spec.pdf', bbox_inches = 'tight')


