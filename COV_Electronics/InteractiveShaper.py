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
    
#Define some colors
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

w_base = np.linspace(1, 1000000, 1000000)


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

ax2.set_xlim([4.8, 10])
ax2.set_xlabel('time [ms]')

ax2.plot(1e3*x, pulse(x), label='pulse')

# The parametrized function to be plotted

def shaperTF(R1, R2, C1, R4, R6, C2):


    om1 = 1./(R1*C1*(R2+R4))
    om2 = 1./((R2+R1)*R6*C2+(R2+R1)*(1+R6/R2)*C2*R6+R1*C1*R2*(1+R6/R2))
    square_om3 = 1./((R2+R1)*R6*R6*C2*C2+R1*C1*R2*R6*C2+R1*C1*R2*C2*R6*(1+R6/R2))
    pow3_om4 = 1./(R1*C1*R2*R6*R6*C2*C2)

#    num=[1./om1, R2+R4]
#    denom=[1./pow3_om4, 1./square_om3, 1./om2, (R2+R1)*(1+R6/R2)]
#    TF_ItalianShaper = scipy.signal.TransferFunction(num, denom)

    num=[1.]
    denom=[R6*R6*C2*C2, R6*C2*(2+R6/R4), 1+R6/R4]
    TF_ItalianShaper = scipy.signal.TransferFunction(num, denom)
    
    w, mag_NuShaper, phase_NuShaper = scipy.signal.bode(TF_ItalianShaper, w_base*2*np.pi)
    
    return TF_ItalianShaper, 10**(mag_NuShaper/20)
        
def amplifiedshaperTF(G, R1, R2, C1, R4, R6, C2, R9, R10, R13, R14, C4):

    om0 = 1./(R1*C1*(R2+R4))
    om1 = 1./((R2+R1)*R6*C2+(R2+R1)*(1+R6/R2)*C2*R6+R1*C1*R2*(1+R6/R2))
    om2_pow2 = 1./((R2+R1)*R6*R6*C2*C2+R1*C1*R2*R6*C2+R1*C1*R2*C2*R6*(1+R6/R2))
    om3_pow3 = 1./(R1*C1*R2*R6*R6*C2*C2)

    om4 = 1./(1./om0 + (R4+R2)*R13*C4)
    om5_pow2 = 1./(R13*C4*1./om0)
    om6 = 1./((1+R13/R14)/om1+R13*C4*(R2*(1+R1/R2)*(1+R6/R2)))
    om7_pow2 = 1./((1+R13/R14)/om2_pow2 + R13*C4/om1)
    om8_pow3 = 1./((1+R13/R14)/om3_pow3 + R13*C4/om2_pow2)
    om9_pow4 = 1./(R13*C4/om3_pow3)

    totG = G*(R9+R10)/R9
    
    num=[totG/om5_pow2, totG/om4, totG*(R2+R4)]
    denom=[1./om9_pow4, 1./om8_pow3, 1./om7_pow2, 1./om6, (1+R13/R14)*(R2+R1)*(1+R6/R2)]
    TF_ItalianShaper = scipy.signal.TransferFunction(num, denom)

    w, mag_NuShaper, phase_NuShaper = scipy.signal.bode(TF_ItalianShaper, w_base*2*np.pi)
    
    return TF_ItalianShaper, 10**(mag_NuShaper/20)

# Define initial parameters
init_R1 = 1e6
init_R2 = 470
init_R4 = 1e3
init_C1 = 22*1e-9
init_R6 = 2.2e3
init_C2 = 100*1e-9

init_G = 11
init_R9= 470
init_R10= 1e3
init_R13 = 1e6
init_R14 = 470
init_C4 = 100*1e-9

# Create the figure and the line that we will manipulate
line3, = ax.loglog(w_base, shaperTF(init_R1, init_R2, init_C1, init_R4, init_R6, init_C2)[1],
                 label='shaper', color='tab:red')
line4, = ax2.plot(1e3*x, scipy.signal.lsim(shaperTF(init_R1, init_R2, init_C1, init_R4, init_R6, init_C2)[0], pulse(x), x)[1],
                  label='shaper', color='tab:red')

line, = ax.loglog(w_base, amplifiedshaperTF(init_G, init_R1, init_R2, init_C1, init_R4, init_R6, init_C2, init_R9, init_R10, init_R13, init_R14, init_C4)[1],
                 label='amplified shaper', color='tab:green')
line2, = ax2.plot(1e3*x, scipy.signal.lsim(amplifiedshaperTF(init_G, init_R1, init_R2, init_C1, init_R4, init_R6, init_C2, init_R9, init_R10, init_R13, init_R14, init_C4)[0], pulse(x), x)[1],
                 label='amplified shaper', color='tab:green')


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
    valstep=np.linspace(-10, -5, 10),
    valinit=np.log10(init_C1),
    color='tab:red'
)

axfreq = fig_inter.add_axes([0.25, 0.30, 0.5, 0.03])
R4_slider = Slider(
    ax=axfreq,
    label='log(R4 [Ohm])',
    valmin=1,
    valmax=6,
    valstep=np.linspace(1, 6, 12),
    valinit=np.log10(init_R4),
    color='tab:red'
)

axfreq = fig_inter.add_axes([0.25, 0.27, 0.5, 0.03])
R6_slider = Slider(
    ax=axfreq,
    label='log(R6 [Ohm])',
    valmin=1,
    valmax=6,
    valstep=np.linspace(1, 6, 12),
    valinit=np.log10(init_R6),
    color='tab:red'
)

axfreq = fig_inter.add_axes([0.25, 0.24, 0.5, 0.03])
C2_slider = Slider(
    ax=axfreq,
    label='log(C2 [F])',
    valmin=-10,
    valmax=-5,
    valstep=np.linspace(-10, -5, 10),
    valinit=np.log10(init_C2),
    color='tab:red'
)


axfreq = fig_inter.add_axes([0.25, 0.14, 0.5, 0.03])
G_slider = Slider(
    ax=axfreq,
    label='Gain',
    valmin=-100,
    valmax=100,
    valstep=np.linspace(-100, 100, 20),
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
    valstep=np.linspace(-10, -4, 10),
    valinit=np.log10(init_C4),
    color='tab:green'
)

# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(amplifiedshaperTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(R4_slider.val), 10**(R6_slider.val), 10**(C2_slider.val), R9_slider.val, 10**(R10_slider.val), 10**(R13_slider.val), R14_slider.val, 10**(C4_slider.val))[1])
    line2.set_ydata(10*scipy.signal.lsim(amplifiedshaperTF(G_slider.val, 10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(R4_slider.val), 10**(R6_slider.val), 10**(C2_slider.val), R9_slider.val, 10**(R10_slider.val), 10**(R13_slider.val), R14_slider.val, 10**(C4_slider.val))[0], pulse(x), x)[1])
    line3.set_ydata(shaperTF(10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(R4_slider.val), 10**(R6_slider.val), 10**(C2_slider.val))[1])
    line4.set_ydata(10*scipy.signal.lsim(shaperTF(10**(R1_slider.val), R2_slider.val, 10**(C1_slider.val), 10**(R4_slider.val), 10**(R6_slider.val), 10**(C2_slider.val))[0], pulse(x), x)[1])
    fig_inter.canvas.draw_idle()

ax.legend()
ax2.legend()

# register the update function with each slider
R1_slider.on_changed(update)
R2_slider.on_changed(update)
R4_slider.on_changed(update)
C1_slider.on_changed(update)
R6_slider.on_changed(update)
C2_slider.on_changed(update)

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
    R4_slider.reset()
    C1_slider.reset()
    R6_slider.reset()
    C2_slider.reset()
    
    G_slider.reset()
    R9_slider.reset()
    R10_slider.reset()
    R13_slider.reset()
    R14_slider.reset()
    C4_slider.reset()
button.on_clicked(reset)

plt.show()
