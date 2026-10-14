#Intended to compare the experimentally determined relaxed helical repeat h_0 according to split band data of Soumya to the Marko Siggia prediction using the data available in 1994 (G/C ~ 2.4).
#This is done by performing the same curve-fitting (at least I suspect) done by MS to the experimental split band locations (via computing the helical repeat N_bp/Tw, where Tw is approximately a half-integer at the splits).
#Produces a plot with a green curve fit to the experimental data by Soumya, and comparison red curve predicted by Lipfert/Skoruppa et al. from magnetic torque tweezer experiments and orange curve I got from studying semicircularly bent DNA in oxDNA2.

#Initially written by Tommy to do the above^, I (Emma) added in the experimental data with the addition of EtBr (two different concentrations) and HMfB
#and also two different temps (8 deg C and 55 deg C), all fit with exponential decay

#For usage, run "python3 plot_exp_vs_TBC_h_0.py --help".

import matplotlib as mpl
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter, StrMethodFormatter
from scipy.optimize import curve_fit
from typing import List, Tuple
from fig_style import apply_style, fig_size_inches, MS, LW

#Parse CLI input parameters.
prog = "plot_exp_vs_TBC_h_0.py"
parser = argparse.ArgumentParser(prog = prog, description='''Intended to plot the experimentally measured helical repeat h_0 (N_bp/Tw) from the split band data (at time of writing, at 73, 83, 93, etc. base pairs with proposed Tw of 6.5, 7.5, 8.5, etc. turns).
                                 Also included are datapoints that MS borrowed from experimental data in the 90s (black datapoints).
                                 The red curve represents the approximate expected twist-bend coupling MS formula as recorded from magnetic torque tweezers experiments in Skoruppa et al. 2017 (G/C = 40/110, omega_0 = 1.75 rad/nm).
                                 The green fit is to the experimental split band data in current state and the MS borrowed datapoints, which implicates a much higher ratio G/C.
                                 Because of the nature of this plot, the program parameters is heavily hard-coded near the header unfortunately in current state.''')
parser.add_argument("plot_dir", type = str, nargs = 1, help = "Output plot directory (.png/.svg).")
parser.add_argument("-x", "--xlim", type = str, nargs = 1, dest = "xlim", help = "optional: Comma-separated x-limits for the plot (inclusive on both ends). Minimum and maximum datapoints currently are at 73 and 366 bp. Defaults to (60,370).")
parser.add_argument("-y", "--ylim", type = str, nargs = 1, dest = "ylim", help = "optional: Comma-separated y-limits for the plot (inclusive on both ends). Defaults to (10.2,11.6).")
args = parser.parse_args()
plot_dir: str = args.plot_dir[0]
xlim: Tuple[float, float]; ylim: Tuple[float, float]
if args.xlim:
    xlim = (float(args.xlim[0].split(",")[0]), float(args.xlim[0].split(",")[1]))
else:
    xlim = (50, 110)
if args.ylim:
    ylim = (float(args.ylim[0].split(",")[0]), float(args.ylim[0].split(",")[1]))
else:
    ylim = (10.4, 12.8)

#Program inputs.
#OxDNA2 (red) curve related parameters.
oxDNA2_ratio: float = 40./110. #Skoruppa et al. suggests torsional stiffness C ~ 40 nm and twist-bend coupling coefficient G ~ 110 I think by fitting to magnetic tweezers data of Lipfert 
oxDNA2_omega_0: float = 1.75 #approximate no curvature twist density in oxDNA2 model (approximately 10.55 bp/turn if convert)

#Split band data (N_bp and proposed Tw)
#no additives
turn_points: List[int] = [63, 73, 83, 93, 103, 123, 143, 164] #the turning points by which Lk_0 shifts from n to n + 1 in experiments (split bands); newest 123, 143, and 164 are separated by 2 (instead of 1)
min_Lk: int = 5 #determinant of the value of the twist number of the first split band location (i.e., turn_Lks[0] = min_Lk + 0.5)
turn_Lks: List[float] = [min_Lk + (i + 1./2.) for i in range(len(turn_points[0:5]))]; turn_Lks.append(11.5); turn_Lks.append(13.5); turn_Lks.append(15.5) #current split band twist numbers suspected are 5.5, 6.5, 7.5, 8.5, 9.5, 11.5, 13.5, and 15.5
turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(turn_points, turn_Lks)] #corresponding suspected helical repeats of form N_bp/Lk for split band data

#EtBr 5 ug/mL
EtBr5_turn_points: List[float] = [54, 64, 74, 84, 94, 105]
EtBr5_turn_Lks: List[float] = [4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
EtBr5_turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(EtBr5_turn_points, EtBr5_turn_Lks)]

#EtBr 12.5 ug/mL
EtBr12_turn_points: List[float] = [55, 65, 75, 86, 98]
EtBr12_turn_Lks: List[float] = [4.5, 5.5, 6.5, 7.5, 8.5]
EtBr12_turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(EtBr12_turn_points, EtBr12_turn_Lks)]

#HMfB
HMfB_turn_points: List[float] = [62, 71, 81, 91, 101]
HMfB_turn_Lks: List[float] = [5.5, 6.5, 7.5, 8.5, 9.5]
HMfB_turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(HMfB_turn_points, HMfB_turn_Lks)]

#8 deg C
temp8_turn_points: List[float] = [72.5, 82.5, 92.0, 102.5]
temp8_min_Lk: int = 6
temp8_turn_Lks: List[float] = [temp8_min_Lk + (i + 1./2.) for i in range(len(temp8_turn_points))] #6.5, 7.5, 8.5, 9.5
temp8_turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(temp8_turn_points, temp8_turn_Lks)]

#55 deg C
temp55_turn_points: List[float] = [74.0, 84.5, 94.0, 104.0]
temp55_min_Lk: int = 6
temp55_turn_Lks: List[float] = [temp55_min_Lk + (i + 1./2.) for i in range(len(temp55_turn_points))] #6.5, 7.5, 8.5, 9.5
temp55_turn_h_0s: List[float] = [N_bp/Lk for N_bp, Lk in zip(temp55_turn_points, temp55_turn_Lks)]

#Marko Siggia values used originally to extrapolate G/C ~ 2.4 through a similar curve-fitting routine (black datapoints). Copied verbatim from MS paper.
MS_N_bps: List[int] = [200, 250, 366]
MS_h_0s: List[float] = [10.54, 10.45, 10.44]

#Finally, parameters found from a separate curve-fitting of MS formula to runs I did of semicircularly bent DNA. Studied the helical repeat as a function of the radii of circles that spatially best fit the configuration.
semicirc_ratio: float = .665
semicirc_asymp_h_0: float = 10.543

#Converts the twist density omega (in units of rad/nm) to the helical repeat (in units of bp/turn).
#It's worth noting that this function is also its own inverse function (i.e., omega2h(h) yields omega). It's most evident writing it out that h = 1/omega * 2Pi/a implies omega = 1/h * 2Pi/a.
def omega2h(omega: float, helical_rise: float = .34) -> float:
	h: float = 1/omega * 2*np.pi/helical_rise #helical repeat (in bp/turn)
	
	return h

#Converts the number of base pairs in a dsDNA circles (i.e., the circumference in units of bp) to the corresponding approximate radius (in units of nm).
def N_bp2r(N_bp: int, helical_rise: float = .34) -> float:
    L: float = N_bp * helical_rise #circumference (in nm)
    r: float = L/(2*np.pi) #radius (in nm)

    return r

#Rewritten version of the equation on page 986 of the original Marko Siggia paper for omega + <Omega_3> 
#that gives the expected helical repeat of a duplex bent to the curvature of an N_bp minicircle.
#This loosely involves some conversion of units from r -> N_bp and omega -> h. 
def calc_TBC_h_0(N_bp: int, ratio: float, asymp_h_0: float) -> float:
   TBC_h_0: float = asymp_h_0/(1-1/2*(ratio*asymp_h_0/N_bp)**2)
    
   return TBC_h_0

#exponential decay func to fit data
def calc_exp_decay_h_0(N_bp: int, amplitude: float, decay: float, offset: float) -> float:
   exp_decay_h_0: float = offset + amplitude*np.exp(-decay*np.array(N_bp))

   return exp_decay_h_0

#Helper to run curve_fit and print results for each dataset
def fit_and_print(label: str, N_bps: List[float], h_0s: List[float]) -> Tuple[float, float]:
    out = curve_fit(calc_TBC_h_0, N_bps, h_0s, p0 = [2.4, 10.5])
    ratio: float = out[0][0]
    asymp_h_0: float = out[0][1]
    print(f"  Fit ratio G/C (unitless): {ratio}")
    print(f"  Fit asymptotic h_0 (bp/turn): {asymp_h_0} (in rad/nm: {omega2h(asymp_h_0)})")
    return ratio, asymp_h_0

#Helper to run an exponential-decay fit and print results.
def fit_exp_decay_and_print(label: str, N_bps: List[float], h_0s: List[float]) -> Tuple[float, float, float]:
    out = curve_fit(calc_exp_decay_h_0, N_bps, h_0s, p0 = [2.0, 0.05, 10.5], maxfev = 10000)
    amplitude: float = out[0][0]
    decay: float = out[0][1]
    offset: float = out[0][2]
    print(f"  Exponential amplitude: {amplitude}")
    print(f"  Exponential decay constant: {decay}")
    print(f"  Exponential offset h_0 (bp/turn): {offset} (in rad/nm: {omega2h(offset)})")
    return amplitude, decay, offset

#Attempt to curve-fit functional form of Marko Siggia paper to Soumya's experimental split band data and MS borrowed experimental datapoints.
print("No additives:")
concat_N_bps: List[int] = turn_points + MS_N_bps
concat_h_0s: List[float] = turn_h_0s + MS_h_0s
out = curve_fit(calc_TBC_h_0, concat_N_bps, concat_h_0s, p0 = [2.4, 10.5])
ratio_fit: float = out[0][0]
asymp_h_0_fit: float = out[0][1]

#Finally, print some minor diagnostics before plotting.
print(f"  Fit ratio G/C (unitless): {ratio_fit}")
print(f"  Fit asymptotic h_0 (bp/turn): {asymp_h_0_fit} (in rad/nm: {omega2h(asymp_h_0_fit)})")
print()
print("The below columns represent numerical calculations of the helical repeat at various N_bp. In order, the columns represent:")
print("1) DNA length (in bp)")
print("2) resulting helical repeat from fit of MS equation to Soumya's experimental split band data (blue 'x's) and the borrowed MS experiment (black 'o's) datapoints") 
print("3) resulting helical repeat from MS equation assuming G/C = 40/110 ~ .36 and omega_0 = 1.75 rad/nm as predicted by Lipfert/Skoruppa et al.") 
print("4) resulting helical repeat using parameters from semicircularly bent DNA in oxDNA2 by me (G/C ~ .665, asymp_h_0 ~ 10.543 bp/turn)")
print()

#EtBr 5 ug/mL
print("EtBr 5 ug/mL:")
#EtBr5_ratio_fit, EtBr5_asymp_h_0_fit = fit_and_print("EtBr 5 ug/mL", EtBr5_turn_points, EtBr5_turn_h_0s)
EtBr5_amp_fit, EtBr5_decay_fit, EtBr5_offset_fit = fit_exp_decay_and_print("EtBr 5 ug/mL", EtBr5_turn_points, EtBr5_turn_h_0s)
print()

#EtBr 12.5 ug/mL
print("EtBr 12.5 ug/mL:")
#EtBr12_ratio_fit, EtBr12_asymp_h_0_fit = fit_and_print("EtBr 12.5 ug/mL", EtBr12_turn_points, EtBr12_turn_h_0s)
EtBr12_amp_fit, EtBr12_decay_fit, EtBr12_offset_fit = fit_exp_decay_and_print("EtBr 12.5 ug/mL", EtBr12_turn_points, EtBr12_turn_h_0s)
print()

#HMfB
print("HMfB:")
#HMfB_ratio_fit, HMfB_asymp_h_0_fit = fit_and_print("HMfB", HMfB_turn_points, HMfB_turn_h_0s)
HMfB_amp_fit, HMfB_decay_fit, HMfB_offset_fit = fit_exp_decay_and_print("HMfB", HMfB_turn_points, HMfB_turn_h_0s)
print()

#8 deg C
print("8 deg C:")
temp8_ratio_fit, temp8_asymp_h_0_fit = fit_and_print("8 $\degree$ C", temp8_turn_points, temp8_turn_h_0s)
print()

#55 deg C
print("55 deg C:")
temp55_ratio_fit, temp55_asymp_h_0_fit = fit_and_print("55 $\degree$ C", temp55_turn_points, temp55_turn_h_0s)
print()

#Unit conversion prep before printing, and creating a list of all integer N_bp in provided limits of plot.
oxDNA2_asymp_h_0: float = omega2h(oxDNA2_omega_0)
N_bps: List[int] = list(range(int(xlim[0]), int(xlim[1]) + 1))
for N_bp in N_bps:
    print(N_bp, calc_TBC_h_0(N_bp, ratio_fit, asymp_h_0_fit), calc_TBC_h_0(N_bp, oxDNA2_ratio, oxDNA2_asymp_h_0), calc_TBC_h_0(N_bp, semicirc_ratio, semicirc_asymp_h_0))
print()

#Now that fit has been completed, and main diagnostics are finished, create the plot.

#Set plot parameters:
apply_style()

#If want control over size of image/aspect ratio (here, golden ratio):
#golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width, fig_height = fig_size_inches()

#Create figure/axes:
fig, ax = plt.subplots(figsize=(fig_width,fig_height))
ax2 = ax.secondary_yaxis('right', functions = (lambda x: 360/x, lambda x: 360/x)) 
ax3 = ax.twiny() #don't need a secondary axis for top x-axis as r scales linearly with N_bp it turns out

#Control gridlines:
ax.grid(which='major',axis='y',color='#404040',linestyle = '--',linewidth = LW.grid,alpha = 0.3, zorder = -20)
ax.grid(which='major',axis='x',color='#404040',linestyle = '--',linewidth = LW.grid,alpha = 0.3, zorder = -20)

#Edit the major and minor tick locations
for axis in [ax.xaxis, ax.yaxis]: #for x and y axis
    axis.set_tick_params(which='major', direction='out', top=False)
    axis.set_tick_params(which='minor', direction='out', top=False)
    axis.set_minor_formatter(NullFormatter()) #believe this makes it so that minor ticks are unlabeled 
#x-axis formatting and locators for major and minor axes
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
#y-axis formatting and locators
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(.5))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(.1))
ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))

#Second y-axis formatting and locators
ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(.25))
ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax2.tick_params(axis = 'y', which = 'major')
ax2.tick_params(axis = 'y', which = 'minor')

#Second x-axis formatting and locators
ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
ax3.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax3.tick_params(axis = 'x', which = 'major')
ax3.tick_params(axis = 'x', which = 'minor')

#Set bounds to plot from (if want full control of bounds)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ylim_phi: Tuple[float, float] = (360/ylim[0],360/ylim[1]) #right y-axis ymin and ymax (twist angle bounds)
ax2.set_ylim(*ylim_phi)
xlim_r: Tuple[float, float] = (N_bp2r(xlim[0]), N_bp2r(xlim[1])) #top x-axis xmin and xmax (radius of curvature bounds)
ax3.set_xlim(*xlim_r)

#Set axis labels:
ax.set_xlabel(r'$N_{\mathregular{bp}}$')
ax.set_ylabel(r'$h_0$ (bp/turn)', rotation = 90)
ax2.set_ylabel(r'$\phi_0$ (°)', rotation = 270, labelpad = 8)
ax3.set_xlabel(r'$r$ (nm)')

#Fit of MS equation (green line) to Soumya split band data (green 'x's) and MS borrowed experiment datapoints (black 'o's)
pts = np.linspace(xlim[0], xlim[1], 400)
ax.plot(turn_points, turn_h_0s, 'x', color = 'green',  markersize = MS.cross, zorder = 3, label = 'No additives') #green h_0 'x' markers (actual Soumya split band data)
#ax.plot(MS_N_bps, MS_h_0s, 'o', color = 'black', markersize = MS.large, zorder = 3) #black h_0 'o' markers (MS experiment cited data)
ax.plot(pts, calc_TBC_h_0(pts, ratio_fit, asymp_h_0_fit), '-', color = 'green', linewidth = LW.fit, zorder = 2, alpha = .75) #fit to above two sets of markers

#oxDNA2 parameter MS curve using G/C = 40/110 and omega_0 = 1.75 cited by Belgian group (red curve)
#ax.plot(pts, calc_TBC_h_0(pts, oxDNA2_ratio, oxDNA2_asymp_h_0), '-', color = 'red', linewidth = LW.fit, zorder = 1, alpha = .75)

#Semicircularly bent DNA oxDNA2 parameter MS curve with G/C ~ .665 and asymp_h_0 ~ 10.543 (orange curve)
#ax.plot(pts, calc_TBC_h_0(pts, semicirc_ratio, semicirc_asymp_h_0), '-', color = 'orange', linewidth = LW.fit, zorder = 1, alpha = .75)

#EtBr 5 ug/mL split band datapoints and fit
ax.plot(EtBr5_turn_points, EtBr5_turn_h_0s, 'x', color = 'cornflowerblue', markersize = MS.cross, zorder = 4, label = 'EtBr 5 ug/mL')
#ax.plot(pts, calc_TBC_h_0(pts, EtBr5_ratio_fit, EtBr5_asymp_h_0_fit), '-', color = 'cyan', linewidth = LW.fit, zorder = 2, alpha = .75)
ax.plot(pts, calc_exp_decay_h_0(pts, EtBr5_amp_fit, EtBr5_decay_fit, EtBr5_offset_fit), '-', color = 'cornflowerblue', linewidth = LW.fit, zorder = 2, alpha = .75)

#EtBr 12.5 ug/mL split band datapoints and fit
ax.plot(EtBr12_turn_points, EtBr12_turn_h_0s, 'x', color = 'blue', markersize = MS.cross, zorder = 3, label = 'EtBr 12.5 ug/mL')
#ax.plot(pts, calc_TBC_h_0(pts, EtBr12_ratio_fit, EtBr12_asymp_h_0_fit), '-', color = 'blue', linewidth = LW.fit, zorder = 2, alpha = .75)
ax.plot(pts, calc_exp_decay_h_0(pts, EtBr12_amp_fit, EtBr12_decay_fit, EtBr12_offset_fit), '-', color = 'blue', linewidth = LW.fit, zorder = 2, alpha = .75)

#HMfB split band datapoints and fit
ax.plot(HMfB_turn_points, HMfB_turn_h_0s, 'x', color = 'red', markersize = MS.cross, zorder = 3, label = 'HMfB')
#ax.plot(pts, calc_TBC_h_0(pts, HMfB_ratio_fit, HMfB_asymp_h_0_fit), '-', color = 'red', linewidth = LW.fit, zorder = 2, alpha = .75)
ax.plot(pts, calc_exp_decay_h_0(pts, HMfB_amp_fit, HMfB_decay_fit, HMfB_offset_fit), '-', color = 'red', linewidth = LW.fit, zorder = 2, alpha = .75)

#8 deg C split band datapoints and fit
ax.plot(temp8_turn_points, temp8_turn_h_0s, 'x', color = 'plum', markersize = MS.cross, zorder = 3, label = '8°C')
ax.plot(pts, calc_TBC_h_0(pts, temp8_ratio_fit, temp8_asymp_h_0_fit), '-', color = 'plum', linewidth = LW.fit, zorder = 2, alpha = .75)

#55 deg C split band datapoints and fit
ax.plot(temp55_turn_points, temp55_turn_h_0s, 'x', color = 'purple', markersize = MS.cross, zorder = 3, label = '55°C')
ax.plot(pts, calc_TBC_h_0(pts, temp55_ratio_fit, temp55_asymp_h_0_fit), '-', color = 'purple', linewidth = LW.fit, zorder = 2, alpha = .75)

legend = ax.legend(loc = 'upper right', ncol = 2, frameon = True, facecolor = 'white', edgecolor = 'black', framealpha = 1.0, handlelength = 1.2, columnspacing = 0.6, handletextpad = 0.4)
legend.get_frame().set_linewidth(0.5)

#HYPOTHETICAL DATAPOINTS
#Plotting the helical repeat if the split band had instead been observed closer to the expected value if h = 10.5 bp/turn.
#For example, the first set corresponds to turn_points = [72, 82, 92, ...] with turn_Lks = [6.5, 7.5, 8.5, ...].
#This is an attempt to get a better intuition of how big a deviance this really is relative to the predicted twist-bend coupling term.
#The datapoints get smaller and darker in shade as the deviance from what was actually measured in experiments is increased.
#print("Hypothetical datapoints shifts (in bp) and fit parameters (G/C, asymp_h_0):") 
#print(f"0 {ratio_fit} {asymp_h_0_fit}")
#for i in range(1,5):
    #turn_points_new = [turn_point - i for turn_point in turn_points] #shifted split band lengths to see what the implicated helical repeat would be if had observed them nearer to expectations
    #turn_h_0s_new = [N_bp/Lk for N_bp, Lk in zip(turn_points_new, turn_Lks)] #corresponding new values of the helical repeat
    #ax.plot(turn_points_new, turn_h_0s_new, 'x', color = (0, 0, 1 - i*.2),  markersize = MS.cross - i*.2, zorder = 3) #additional markers

    #Attempt to curve-fit functional form of Marko Siggia paper to Soumya's experimental split band data and MS borrowed experimental datapoints.
    #out = curve_fit(calc_TBC_h_0, turn_points_new, turn_h_0s_new, p0 = [1, 10.5])
    #hypo_ratio_fit: float = out[0][0]
    #hypo_asymp_h_0_fit: float = out[0][1]

    #ax.plot(pts, calc_TBC_h_0(pts, hypo_ratio_fit, hypo_asymp_h_0_fit), '-', color = (0, 0, 1 - i*.2), linewidth = LW.fit, zorder = 2, alpha = .75)
    #print(f"{-i} {hypo_ratio_fit} {hypo_asymp_h_0_fit}")
#print()    

#Tight layout to resize plot to fit plot + labels within canvas (bbox_inches = tight flag in savefig makes boundaries "stretch" to fit in the labels; this preserves canvas size as specified by fig_width and fig_height)
plt.tight_layout(pad = .2)

plt.savefig(plot_dir)

print("Saved plot to", plot_dir + ".")

plt.close()
