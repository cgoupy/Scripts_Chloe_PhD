import numpy as np
from matplotlib import pyplot as plt
from scipy.misc import derivative
import scipy.integrate as integrate
import scipy
from matplotlib.widgets import Button, Slider
import colorednoise as cn

tau1 = 1*1e-6 #s
tau2 = 2*1e-3 #s
A=20
t0=0.005

SamplingFreq=1*1e6 #Hz
Tmax=20*1e-3 # sec
x = np.linspace(0, 0.02, int(Tmax*SamplingFreq))
N = len(x)

noise_freqs = 2*np.pi*np.array([1.5*1e3, 100, 2e3, 2.3e3, 3e3, 2e4, 3e4]) #Hz
ampls = np.array([0.5, 0.9, 0.5, 0.5, 0.5,0.3, 0.1]) #Hz

def cdf(x, mu, sig):
    return 0.5*(1+scipy.special.erf((x-mu)/(np.sqrt(2)*sig)))

def pulse(x):
    res = cdf(x, t0, tau1)*A*(np.exp(-(x-t0)/tau2))
    return res


def pulse_pileup(x):
    delay=t0+2e-3
    res = pulse(x)+0.05*(cdf(x, delay, tau1)*A*(np.exp(-(x-delay)/tau2)))
    return res
    
#Deffing some colores
X_brown = cn.powerlaw_psd_gaussian(2, len(x))
X_white = cn.powerlaw_psd_gaussian(0, len(x),fmin=0.5)


def noise(x):
    res=0.6*X_brown+0.15*X_white
    for i in range(len(noise_freqs)):
        res += ampls[i]*np.sin(x*noise_freqs[i])
    return res
    
np.vectorize(noise);

pulse_f = scipy.fftpack.fft(pulse(x))
N = len(x)
# sample spacing
T = x[1]-x[0]
yf = scipy.fftpack.fft(pulse(x))
xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
yf_noise = scipy.fftpack.fft(noise(x))

w_base = np.linspace(10, 1000000, 100000)


## ------------------------ Interactive plot starts here ------------------------
fig_inter = plt.figure()
ax = fig_inter.add_subplot(121)
ax2 = fig_inter.add_subplot(122)

fig_inter.set_figheight(8)
fig_inter.set_figwidth(10)
fig_inter.subplots_adjust(top=0.95, left=0.05, right=0.95)

ax.loglog(xf, (2.0/N * np.abs(pulse_f[:N//2])), label="signal")
ax.loglog(xf, 2.0/N * np.abs(yf_noise[:N//2]), label="noise")
ax.set_xlabel('Frequency [Hz]')

ax2.set_xlim([4.8, 7])
ax2.set_xlabel('time [ms]')

ax2.plot(1e3*x, pulse(x), label='pulse')

# The parametrized function to be plotted

def ampliTF(G, R1, R2, C1, C4, R13, R14, R9, R10):

    w_base = np.linspace(10, 1000000, 100000)

    om0 = 1./(R1*C1+R13*C4)
    om1_pow2 = 1./(R1*C1*R13*C4)
    om2 = 1./(R1*C1*(1+R13/R14)+R13*C4*(1+R1/R2))
    om3_pow2 = 1./(R1*C1*R13*C4)

    totG = G*(R9+R10)/R9
    num=[totG/om1_pow2,totG/om0, totG]
    denom=[1./om3_pow2, 1./om2, (1+R1/R2)*(1+R13/R14)]
    
    TF_Amplifier = scipy.signal.TransferFunction(num, denom)
    
    w, mag_NuShaper, phase_NuShaper = scipy.signal.bode(TF_Amplifier, w_base)
    
    return TF_Amplifier, 10**(mag_NuShaper/20)
    
def unshapedTF(G, R1, R2, C1, C4, R13, R14, R9, R10):

    w_base = np.linspace(10, 1000000, 100000)

    om0 = 1./(R1*C1+R13*C4)
    om1_pow2 = 1./(R1*C1*R13*C4)
    om2 = 1./(R1*C1*(1+R13/R14)+R13*C4*(1+R1/R2))
    om3_pow2 = 1./(R1*C1*R13*C4)

    totG = G*(R9+R10)/R9
    num=[1/om1_pow2,1/om0, 1]
    denom=[1./om3_pow2, 1./om2, (1+R1/R2)*(1+R13/R14)]
    
    TF_Amplifier = scipy.signal.TransferFunction(num, denom)

    w, mag_NuShaper, phase_NuShaper = scipy.signal.bode(TF_Amplifier, w_base)
    
    return TF_Amplifier, 10**(mag_NuShaper/20)

# Define initial parameters
init_R1 = 1e6
init_R2 = 1e3
init_C1 = 10*1e-6

init_G = 11
init_R9= 100
init_R10= 359 #1e3
init_R13 = 770 #1e6
init_R14 = 461 #470
init_C4 = 100*1e-9

# Create the figure and the line that we will manipulate
line, = ax.loglog(w_base, ampliTF(init_G, init_R1, init_R2, init_C1, init_C4, init_R13, init_R14, init_R9, init_R10)[1],
                 label='amplified no shaper', color='tab:red')
line2, = ax2.plot(1e3*x, 0.01*scipy.signal.lsim(ampliTF(init_G, init_R1, init_R2, init_C1, init_C4, init_R13, init_R14,init_R9, init_R10)[0], pulse(x), x)[1],
                  label='0.01*amplified no shaper', color='tab:red')
                  
line3, = ax.loglog(w_base, unshapedTF(init_G, init_R1, init_R2, init_C1, init_C4, init_R13, init_R14, init_R9, init_R10)[1],
                 label='no shaper', color='tab:green')
line4, = ax2.plot(1e3*x, 1*scipy.signal.lsim(unshapedTF(init_G, init_R1, init_R2, init_C1, init_C4, init_R13, init_R14,init_R9, init_R10)[0], pulse(x), x)[1],
                 label='no shaper', color='tab:green')


# adjust the main plot to make room for the sliders
fig_inter.subplots_adjust(bottom=0.5)

# Make a horizontal slider to control the frequency.
axfreq = fig_inter.add_axes([0.25, 0.39, 0.5, 0.03])
R1_slider = Slider(
    ax=axfreq,
    label='log(R1 [Ohm])',
    valmin=1,
    valmax=10,
    valstep=np.linspace(1, 10, 19),
    valinit=np.log10(init_R1),
    color='tab:red'
)

axfreq = fig_inter.add_axes([0.25, 0.36, 0.5, 0.03])
R2_slider = Slider(
    ax=axfreq,
    label='R2 [Ohm]',
    valmin=10,
    valmax=1000,
    valstep=np.linspace(10, 1000, 100),
    valinit=init_R2,
    color='tab:red'
)

axfreq = fig_inter.add_axes([0.25, 0.33, 0.5, 0.03])
C1_slider = Slider(
    ax=axfreq,
    label='log(C1 [F])',
    valmin=-10,
    valmax=-5,
    valstep=np.linspace(-10, -1, 10),
    valinit=np.log10(init_C1),
    color='tab:red'
)


axfreq = fig_inter.add_axes([0.25, 0.14, 0.5, 0.03])
G_slider = Slider(
    ax=axfreq,
    label='Gain',
    valmin=-11,
    valmax=11,
    valstep=np.linspace(-11, 11, 10000),
    valinit=np.log10(init_G),
    color='tab:green'
)

axfreq = fig_inter.add_axes([0.25, 0.17, 0.5, 0.03])
R9_slider = Slider(
    ax=axfreq,
    label='R9 [Ohm]',
    valmin=10,
    valmax=1000,
    valstep=np.linspace(10, 1000, 100),
    valinit=init_R9,
    color='tab:green'
)

axfreq = fig_inter.add_axes([0.25, 0.20, 0.5, 0.03])
R10_slider = Slider(
    ax=axfreq,
    label='log(R10 [Ohm])',
    valmin=1,
    valmax=6,
    valstep=np.linspace(1, 6, 12),
    valinit=np.log10(init_R10),
     color='tab:green'
)

axfreq = fig_inter.add_axes([0.25, 0.11, 0.5, 0.03])
R13_slider = Slider(
    ax=axfreq,
    label='log(R13 [Ohm])',
    valmin=1,
    valmax=10,
    valstep=np.linspace(1, 10, 19),
    valinit=np.log10(init_R13),
    color='tab:green'
)

axfreq = fig_inter.add_axes([0.25, 0.08, 0.5, 0.03])
R14_slider = Slider(
    ax=axfreq,
    label='R14 [Ohm]',
    valmin=10,
    valmax=1000,
    valstep=np.linspace(10, 1000, 100),
    valinit=init_R14,
    color='tab:green'
)

axfreq = fig_inter.add_axes([0.25, 0.05, 0.5, 0.03])
C4_slider = Slider(
    ax=axfreq,
    label='log(C4 [F])',
    valmin=-10,
    valmax=-4,
    valstep=np.linspace(-10, -4, 7),
    valinit=np.log10(init_C4),
    color='tab:green'
)

# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(ampliTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(C4_slider.val), 10**(R13_slider.val), R14_slider.val, R9_slider.val, 10**(R10_slider.val))[1])
    line2.set_ydata(0.01*scipy.signal.lsim(ampliTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(C4_slider.val), 10**(R13_slider.val), R14_slider.val, R9_slider.val, 10**(R10_slider.val))[0], pulse(x), x)[1])
    line3.set_ydata(unshapedTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(C4_slider.val), 10**(R13_slider.val), R14_slider.val, R9_slider.val, 10**(R10_slider.val))[1])
    line4.set_ydata(scipy.signal.lsim(unshapedTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(C4_slider.val), 10**(R13_slider.val), R14_slider.val, R9_slider.val, 10**(R10_slider.val))[0], pulse(x), x)[1])
    
ax.legend()
ax2.legend()

# register the update function with each slider
R1_slider.on_changed(update)
R2_slider.on_changed(update)
C1_slider.on_changed(update)

G_slider.on_changed(update)
R9_slider.on_changed(update)
R10_slider.on_changed(update)
R13_slider.on_changed(update)
R14_slider.on_changed(update)
C4_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig_inter.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')

def reset(event):
    R1_slider.reset()
    R2_slider.reset()
    C1_slider.reset()
   
    G_slider.reset()
    R9_slider.reset()
    R10_slider.reset()
    R13_slider.reset()
    R14_slider.reset()
    C4_slider.reset()
button.on_clicked(reset)

plt.show()
