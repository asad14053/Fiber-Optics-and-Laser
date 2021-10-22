# The code running procedure; (Windows) Python
# 1. install pycharm:
# install in windows: https://www.jetbrains.com/pycharm/download/#section=windows
# install in mac: https://www.jetbrains.com/pycharm/download/#section=mac
# 2. Click community version downloads and install it.
#And, install python after downloading (https://www.python.org/downloads/)
# 3. open pycharm and create a project (any  name), paste this code into main.py.
# 4. To run code needs 2 packages: (a) matplotlib (b) Lightpipes
# 5. To install packages: go to view>tool windows>python packages, a popup window will appear and type those following 2 package names and click install.
#or, click terminal button on the bottom left and put command --> "pip install matplotlib" and "pip install Lightpipes" and enter
# if might see a bar for downloading like "==========================================================="
# 6. # After successful installing click run ( green  ">" button on upper right) button and see the results.
# The response time is too slow (Please wait atleast 10 seconds after each changes)

#https://opticspy.github.io/lightpipes/command-reference.html
#The Gain command introduces a simple single-layer model of a laser amplifier. The output field is
# F_out = F_in(x,y)e^((alpha)*L_gain), alpha = alpha0/(1+2*I(x,y)/I_sat), 2*alpha*L_gain = net round trip intensity gain, alpha0 = small signal intensity gain
# I_sat = is the saturation intensity of the gain medium with a length of resonator L_gain



###################################################################################### import necessary packages

import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use("TkAgg")
from scipy.stats import multivariate_normal
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import math
import time
import sys

if sys.version_info[0] < 3:
    from Tkinter import *
    import Tkinter as Tk
else:
    from tkinter import *
    import tkinter as Tk
from LightPipes import *

root = Tk.Tk()
root.wm_title("Laser with stable resonator")
root.wm_protocol("WM_DELETE_WINDOW", root.quit)




########################################################################### initialize necessary variables and taking input
power = []
roundtrip = 0


#OSA initialize
wavelength = 10600 * nm;  #Lasing wavelength
size = 20 * mm;    #field size- the square grid
N = 100; #field size- the grid dimension

#OSA gain
Isat = 131 * W / cm / cm; #Intensity saturation in laser
alpha = 0.0067 / cm;  #small gain threshold
Lgain = 30 * cm;  #gain medium length

# Lens
f1 = 2.0 * m  #mirror M1 radius
f2 = 5 * m    #mirror M2 radius
w0 = 2.4 * mm # lens aperture diameter

#gauss parameters
width1 = 0.0   # the initial width of the Gauss shape
mean1 = 0.0
var1 = 0.0

#cavity and coupling ratio
L = 30 * cm #resonator length
Reflect = 0.9; #coupling ratio


#tx = 0.00 * mrad;
#ty = 0.00 * mrad;
#xwire = 10.0 * mm
#ywire = 10.0 * mm


# Initialize Plotting
fig = plt.figure(figsize=(6, 6))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312, projection='3d')
ax3 = fig.add_subplot(313)
#ax2 = fig.add_subplot(212)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
v = StringVar()
v2 = StringVar()
dnt = IntVar()
LG = IntVar()

F = Begin(size, wavelength, N) #Initiates a field with a grid size, a wavelength and a grid dimension



################################################################################# The main code is here

###Taking input from GUI

#Block Diagram:   |Intensity -->M1-->M2|----------->SOA ---> Filter-->Coupler --->Atteneuation--| loop
                     #<<------------------------------------------------------------------------|


def TheExample():
    global F, f1, f2, L, w0, roundtrip, Reflect, width1, mean1, var1
    # start_time=time.time()
    w0 = float(scale_w0.get()) * mm / 2  # lens aperture diameter
    #xwire = float(scalexwire.get()) * mm
    #ywire = float(scaleywire.get()) * mm
    f1 = float(scale_f1.get() * cm) / 2   #Mirror M1 radius
    f2 = float(scale_f2.get() * cm) / 2   #Mirror M2 radius
    L = float(scale_L.get()) * cm        #Resonator length
    Reflect = float(scale_Reflect.get())  #Coupling ratio
    #tx = -float(scale_tx.get()) * mrad
    #ty = float(scale_ty.get()) * mrad
    alpha = float(scale_gain.get()) / cm  #Small gain threshold


### Executing the Methods and equipments
    F = RandomIntensity(time.time(), 1e-8, F)      # Adds random intensity to the field:  RandomIntensity(Fin, seed=123, noise=1.0)
    F = CircAperture(w0, 0, 0, F)                  # Inserts a circular aperture in the field:  CircAperture(Fin, R, x_shift=0.0, y_shift=0.0)
    #F = RectScreen(size, 0.2 * mm, 0.0, ywire, 0.0, F)
    #F = RectScreen(0.2 * mm, size, xwire, 0.0, 0.0, F)
    Iw = Intensity(0, F)                           # Calculates the intensity of the field, Intensity(Fin, flag=0)
    F = Lens(f2, 0, 0, F);                         # Propagates the field through an ideal, thin lens. If the input field is pure Gaussian, the ABCD matrix theory is used.
    F = Forvard(L, F);                             #Propagates the field using a FFT algorithm:  Forvard(Fin, z)
    F = Gain(Isat, alpha, Lgain, F);               # Propagates the field through a thin saturable gain sheet: Gain(Fin, Isat, alpha0, Lgain)
    F = Lens(f1, 0, 0, F);                         # Propagates the field through an ideal, thin lens. If the input field is pure Gaussian, the ABCD matrix theory is used.
    #F = Tilt(tx, ty, F)
    F = Forvard(L, F);                             #Propagates the field using a FFT algorithm:  Forvard(Fin, z)
    F = Gain(Isat, alpha, Lgain, F);               # Propagates the field through a thin saturable gain sheet: Gain(Fin, Isat, alpha0, Lgain)

    #Random filter
    F1 = PipFFT(1, F);
    F1 = GaussAperture(w0, 0, 0, F1)                  # Inserts a circular aperture in the field:  CircAperture(Fin, R, x_shift=0.0, y_shift=0.0)
    F1 = PipFFT(-1, F1);
    I1 = Intensity(1, F1);
    #F= I1


    F = IntAttenuator(Reflect, F)                  # Attenuates the intensity of the field, IntAttenuator(Fin, att=0.5)
    P = Power(F) * (1 - Reflect) * size / N * size / N   #Calculates the total power.
    power.append(P);                                     #append the power in the power array
    roundtrip = roundtrip + 1



#for realtime plot, needs to refresh the power arrays
    if (roundtrip > 500):
        power.pop(0)
    Iout = Isat * (alpha * Lgain - 0.5 * math.log(1 / Reflect)) * math.pi * w0 * w0  #intesity equation




##########################################################################     end of  main code

    ax1.clear()
    ax2.clear()
    ax3.clear()

    g1 = 1 - L / (2 * f1);   #gain from lens 1
    g2 = 1 - L / (2 * f2);   #gain from lens 2
    g = g1 * g2
    v.set("Power=%5.3f mW\n" % P +
          # "g1 = %5.3f\n" % g1 +
          # "g2 = %5.3f\n" % g2 +
          # "g  = %5.3f\n" % g +
          "Spatial width of Gauss  = %5.3f mm\n" %  width1 +
          "Mean of Gauss  = %5.3f\n" % mean1 +
          "Standard Deviation of Gauss  = %5.3f\n" % math.sqrt(var1) +
          "variance of Gauss  = %5.3f\n" % var1
          )

#################################################################  Plot 1
    ax1.imshow(Iw, cmap='rainbow');
    ax1.axis('off');
    ax1.axis('equal')
    ax1.set_title('laser mode')

    # ax2.imshow(I1, cmap='rainbow');
    # ax2.axis('off');
    # ax2.axis('equal')
    # ax2.set_title('Filter mode')

    # X = np.zeros(N)
    # Y = np.zeros(N)
    # I2 = np.array(Iw)
    # ax2.plot_surface(X, Y, I2, rstride=1, cstride=1, cmap='rainbow', linewidth=0.0)
    # #plt.axis('off');
    # #plt.title('Far-field intensity distribution')
    # #
    # # plt.show()
##################################################################  Plot 2

    #X = np.linspace(0, wavelength, 10)
    X = np.linspace(-5, 5, 50)
    Y = np.linspace(0, wavelength, 10)
    X, Y = np.meshgrid(X, Y)    #mesh 3D design
    X_mean = 0;
    Y_mean = 0
    X_var = 5;
    Y_var = 8
    pos = np.empty(X.shape + (2,))



    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    rv = multivariate_normal([X_mean, Y_mean], [[X_var, 0], [0, Y_var]])  #apply gauss equation
    gauss = rv.pdf(pos*P);                                                #apply laser power into gauss function
    mean, var = norm.stats(2*gauss)                                       #find normal distribution of the gauss
    #width = norm.interval(gauss)
    mean1 = mean[1][1]*P                                                  #find mean of the gauss
    var1 = var[1][1] * P                                                  #find variance of the gauss
    width1 = (2* math.sqrt(var[1][1]))*1/P                                #find width of gauss = 2*sigma of the gauss, sigma = sqrt (variance)

    #print()
    #print("-------------------")
    #print(var[0][0][0])
    ax2.plot_surface(X, Y, gauss, cmap="plasma")
    ax2.set_xlabel('Wavelength', fontsize=10)
    #ax2.set_ylabel('Spatial dimension', fontsize=10)
    ax2.set_zlabel('Power', fontsize=10)
    #ax2.axis('off');
    #ax2.axis('equal')
    ax2.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labelbottom='off'  # labels along the bottom edge are off)
    )

    #ax.plot(x, y, z)

   # plt.show()



    #plt.show()
    #################################################################  Plot 3

    ax3.plot(power);
    ax3.set_ylim(0, 10);
    ax3.set_xlim(0, 500)
    s = '%3.1f ns/div' % (2.0 * L / 2.988 * 1000.0)
    ax3.set_xlabel(s);
    ax3.set_ylabel('power [W]')
    ax3.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='on',  # ticks along the bottom edge are off
        top='on',  # ticks along the top edge are off
        labelbottom='off')
    ax3.grid()

    # ax3.imshow(Iw, cmap='rainbow');
    # ax3.axis('off');
    # ax3.axis('equal')
    # ax3.set_title('Filter mode')



    canvas.draw()

    # I2 = np.array(Iw)
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(X, Y, I2, rstride=1, cstride=1, cmap='rainbow', linewidth=0.0)
    # plt.axis('off');
    # plt.title('Far-field intensity distribution')
    #
    # plt.show()


#######################################root plot destroy process


def _quit():
    root.quit()
    root.destroy()

################################################     Eign modes TE, TM, TEM00, TE 01, .etc.
def _eigenmode():
    global F, f1, f2, L, w0
    g1 = 1 - L / (2 * f1);
    g2 = 1 - L / (2 * f2);
    g = g1 * g2
    z1 = L * g2 * (1 - g1) / (g1 + g2 - 2 * g1 * g2);
    z2 = L - z1;
    if (g > 0):
        w0 = math.sqrt(wavelength * L / math.pi) * (g1 * g2 * (1 - g1 * g2) / (g1 + g2 - 2 * g1 * g2) ** 2) ** 0.25;
    mode_m = int(order_m.get())   #get 1st input
    mode_n = int(order_n.get())   #get 2nd input

    if dnt.get():
        m = mode_m
        if m == 0:
            m = 1
        v2.set(r'Injected eigen mode: ' + f'dougnut LG{m},{mode_n}*')
        F = GaussBeam(F, w0, doughnut=True, m=m, n=mode_n)       #Gauss beam for laser
    else:
        if LG.get():
            F = GaussBeam(F, w0, LG=True, m=mode_m, n=mode_n)    #Gauss beam for laser
            v2.set(r'Injected eigen mode: ' + f'Laguerre-Gauss{mode_m},{mode_n}')
        else:
            F = GaussBeam(F, w0, LG=False, m=mode_m, n=mode_n)   #Gauss beam for laser
            v2.set(r'Injected eigen mode: ' + f'Hermite-Gauss{mode_m},{mode_n}')
    F = Forvard(z2, F);   #Propagates the field using a FFT algorithm:  Forvard(Fin, z)

##################################################### Framing for the hovering in each button
frame1 = Frame(root)
frame1.pack(side=Tk.BOTTOM)
frame2 = Frame(frame1)
frame2.pack(side=Tk.BOTTOM)
frame3 = Frame(frame2)
frame3.pack(side=Tk.BOTTOM)
frame4 = Frame(frame3)
frame4.pack(side=Tk.BOTTOM)
frame5 = Frame(frame4)
frame5.pack(side=Tk.BOTTOM)
frame6 = Frame(frame5)
frame6.pack(side=Tk.BOTTOM)
frame7 = Frame(frame6)
frame7.pack(side=Tk.BOTTOM)

Label(root, textvariable=v).pack(side=Tk.LEFT)
Label(root, textvariable=v2).pack(side=Tk.LEFT)

#scalexwire = Tk.Scale(frame1, orient='horizontal', label='x-wire position [mm]', length=200, from_=-size / 2 / mm, to=size / 2 / mm, resolution=0.001)
#scalexwire.pack(side=Tk.LEFT)
#scalexwire.set(xwire / mm)

#scaleywire = Tk.Scale(frame1, orient='horizontal', label='y-wire position [mm]', length=200, from_=-size / 2 / mm, to=size / 2 / mm, resolution=0.001)
#scaleywire.pack(side=Tk.LEFT)
#scaleywire.set(ywire / mm)

scale_w0 = Tk.Scale(frame2, orient='horizontal', label='aperture diameter [mm]', length=200, from_=0.0, to=size / mm,
                    resolution=0.01)
scale_w0.pack(side=Tk.LEFT)
scale_w0.set(2 * w0 / mm)

scale_Reflect = Tk.Scale(frame2, orient='horizontal', label='Reflectivity', length=200, from_=0.0, to=1.0,
                         resolution=0.01)
scale_Reflect.pack(side=Tk.LEFT)
scale_Reflect.set(Reflect)

scale_f1 = Tk.Scale(frame3, orient='horizontal', label='mirror M1 radius [cm]', length=200, from_=10.0, to=1000.0,
                    resolution=0.1)
scale_f1.pack(side=Tk.LEFT)
scale_f1.set(f1 / cm)

scale_f2 = Tk.Scale(frame3, orient='horizontal', label='mirror M2 radius [cm]', length=200, from_=10.0, to=1000.0,
                    resolution=0.1)
scale_f2.pack(side=Tk.LEFT)
scale_f2.set(f2 / cm)

scale_L = Tk.Scale(frame4, orient='horizontal', label='resonator length [cm]', length=200, from_=10.0, to=100.0,
                   resolution=0.01)
scale_L.pack(side=Tk.LEFT)
scale_L.set(L / cm)

scale_gain = Tk.Scale(frame4, orient='horizontal', label='gain [cm^-1]', length=200, from_=0.0, to=0.01,
                      resolution=0.0001)
scale_gain.pack(side=Tk.LEFT)
scale_gain.set(alpha * cm)

#scale_tx = Tk.Scale(frame5, orient='horizontal', label='mirror M2 x-tilt [mrad]', length=200, from_=-10.0, to=10.0, resolution=0.1)
#scale_tx.pack(side=Tk.LEFT)
#scale_tx.set(tx / mrad)

#scale_ty = Tk.Scale(frame5, orient='horizontal', label='mirror M2 y-tilt [mrad]', length=200, from_=-10.0, to=10.0, resolution=0.1)
#scale_ty.pack(side=Tk.LEFT)
#scale_ty.set(ty / mrad)

button_eigenmode = Tk.Button(frame6, width=18, text='eigen mode', command=_eigenmode)
button_eigenmode.pack(side=Tk.LEFT, pady=10)

order_m = Tk.Spinbox(frame6, width=1, from_=0, to=5)
order_m.pack(side=Tk.LEFT)

order_n = Tk.Spinbox(frame6, width=1, from_=0, to=5)
order_n.pack(side=Tk.LEFT, pady=10)

#doughnut = Tk.Checkbutton(frame6, text='doughnut', variable=dnt)
#doughnut.pack(side=Tk.LEFT, pady=10)

#Laguerre = Tk.Checkbutton(frame6, text='Laguerre Gauss', variable=LG)
#Laguerre.pack(side=Tk.LEFT, pady=10)

button_quit = Tk.Button(frame7, width=24, text='Quit', command=_quit)
button_quit.pack(side=Tk.LEFT, pady=10)


def task():
    TheExample()
    root.after(1, task)


root.after(1, task)
root.mainloop()
