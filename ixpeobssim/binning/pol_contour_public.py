import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from ixpeobssim.binning.polarization import *

plt.rcParams["font.family"] = "Gill Sans"

#### THIS IS A POLARIZATION CONTOUR TOOL VER 3.0 ####
#### RELEASED JAN. 24. 2023. - DAWOON E. KIM ####

#### DEFAULT SETTING ####
#### C.L. = {"50.0%", "90.0%", "99.0%", "99.9%"} ####
#### CALCULATION BASED ON NORMALIZED Q & U ####

#### THIS IS A UNIT SCALING FACTOR FOR THE POLARIZATION DEGREE ####
UNIT = 100

#### PARAMETERS SETTING ####
def CHI_LIST(levels=None, sigma_unit=None):
    # sigma for 99% = 2.57583 in 1 D.O.F.
    # sigma for 99% = 3.03485492 in 2 D.O.F.
    from scipy.stats import chi2, norm

    if not levels:
        result = np.array([0.50, 0.90, 0.99, 0.999])
        labels = ["50.0 %", "90.0 %", "99.0 %", "99.9 %"]

    else:
        result=[]
        labels=[]

        for i in range(len(levels)):
            if not sigma_unit:
                globals()['std'+str(levels[i])] = levels[i]
                print("#### VALUES MUST BE SET WITHIN 0. - 1. ####")
                print("CONTOUR LEVELS :", globals()['std'+str(levels[i])])
            else :
                mean = 0
                SD = 1
                globals()['std'+str(levels[i])] = norm.cdf((levels[i])*SD, mean, SD) - norm.cdf(-(levels[i])*SD, mean, SD)
                print("CONTOUR LEVELS :", globals()['std'+str(levels[i])])

            result.append(globals()['std'+str(levels[i])])
            CL = format(globals()['std'+str(levels[i])]*100, ".4f")
            labels.append(f'{CL}$\%$')

    degree_freedom = 2

    ### CALCULATING CHI_2 FROM PERCENTAGE ###
    list_chi_2 = np.array(np.sqrt(chi2.ppf(result, degree_freedom)))

    print("CONTOUR LEVELS IN CHI_2, SQRT(CHI**2) :", list_chi_2)
    print("CONFIDENCE LEVELS :", labels)

    return list_chi_2, labels


################### POL PROPERTIES CALCULATION ####################
def POL_CALCULATION(QN, UN):
    PD = np.sqrt(QN**2 + UN**2)
    PA = 0.5 * np.arctan2(UN,QN)
    return PD, PA

def TRANS_AND_RADIUS(QN, UN, EPSILON, aspect=True):
    q0 = QN
    u0 = UN
    zeta = np.linspace(0, 2*np.pi, 5000)
    Q_Ec= q0 + EPSILON * np.cos(zeta)
    U_Ec= u0 + EPSILON * np.sin(zeta)
    PI, PSI = POL_CALCULATION(Q_Ec, U_Ec)
    PI_0, PSI_0 = POL_CALCULATION(q0, u0)

    if not aspect:
        for i in range(PSI.size):
            if PSI[i] < 0:
                PSI[i] = PSI[i] + np.pi
        if PSI_0 < 0:
            PSI_0 = PSI_0 +np.pi
    else:
        PI = PI
        PI_0 = PI_0

    return PI_0, PSI_0, PI, PSI

################### PLOT ###################
def POL_CONTOUR_PLOT(QN, UN, QN_ERR, UN_ERR, aspect=True, ax=None, levels=None, sigma_unit=None, legend=True, text=True, rmax=None,global_color=None):
    from scipy.stats import chi2, norm

    print("QN(%): ", np.round(QN*UNIT ,2),", UN(%): ", np.round(UN*UNIT ,2))
    print("QN_ERR(%): ", np.round(QN_ERR*UNIT ,2),", UN_ERR(%): ", np.round(UN_ERR*UNIT ,2))

    list_chi_2, labels = CHI_LIST(levels=levels, sigma_unit=None)
    sigma= np.abs(QN_ERR)
    epsilon = list_chi_2 * sigma

    if not ax:
        cm = 1/2.54
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5, 5), tight_layout=True)

    # MDP #
    MDP_99 = np.sqrt(chi2.ppf(0.99,2)) * np.abs(QN_ERR) * UNIT

    # MAIN PLOT #
    for i in range(len(epsilon)):
        pi_0, psi_0, pi, psi = TRANS_AND_RADIUS(QN, UN, epsilon[i], aspect=aspect)
        pi_0 = pi_0 * UNIT
        pi = pi * UNIT

        if not global_color:
            color_input = ['red','orange','limegreen','tab:blue']
            ax.scatter(psi, pi, marker='.', label=labels[i], color=color_input[i], s=0.02, alpha=0.2)
            ax.scatter(psi_0, pi_0, marker='o', s=20, color='k', edgecolors='k', alpha=0.2)
        else:
            ax.scatter(psi, pi, marker='.', label=labels[i], s=0.02, color=global_color, alpha=0.2)
            ax.scatter(psi_0, pi_0, marker='o', s=20, color=global_color, edgecolors=global_color, alpha=0.2)

    POL_DEGREE= np.round(pi_0, 2)
    POL_ANGLE= np.round(np.rad2deg(psi_0), 2)

    stat_signif = np.round(pi_0/(np.abs(sigma)*UNIT), 3)
    print("==================== SUMMARY ====================")
    print("Statistical Significance", np.round(stat_signif,2))
    print("MDP_99, 2 D.O.F. (%): ", np.round(MDP_99,2))

    Detection_CL= np.round(chi2.cdf(stat_signif**2, 2)*UNIT, 2)
    print("Detection Significance (%)", Detection_CL)

    print("PD(%):", POL_DEGREE, "PD_ERR_1D(±)", np.round(np.abs(QN_ERR)*UNIT,2), 'PA(°):', POL_ANGLE, "PA_ERR_1D(±)", np.round(np.rad2deg(np.abs(QN_ERR)*UNIT)/(2*pi_0),2))

    ####### MAXIMA PD 100 % ######
    if np.max(pi[-1]) >= 100:
        ax.set_rmax(100)

    elif rmax:
        ax.set_rmax(rmax)

    ax.set_theta_zero_location("N")

    ####### ASPECT CHANGE ######
    if not aspect:
        ax.set_thetamax(180)
        tex_info = [[0.69, 1.05, 'N'], [0.69, -0.05,  'S'], [0.14, 0.5, 'E']]
        # title_pos =-0.2
        legend_pos = [0.3,0.97]
        grid_lines, grid_labels = plt.thetagrids(range(0, 181, 30))
        ax.text(0.85,0.75,'$\Pi$\n(%)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, fontsize='large')
    elif aspect=='half left':
        ax.set_thetamax(90)
        tex_info = [[1., 1.1, 'N'], [-0.13, 0., 'E']]
        # title_pos =-0.1
        legend_pos = [0.3,.95]
        grid_lines, grid_labels = plt.thetagrids(range(0, 91, 30))
    elif aspect=='half right':
        ax.set_thetamin(-90)
        ax.set_thetamax(0)
        tex_info = [[0., 1.2, 'N'], [1.25, 0, 'W']]
        # title_pos =-0.27
        legend_pos = [.98,0.92]
        grid_lines, grid_labels = plt.thetagrids(range(-90, 1, 30))
    else:
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        tex_info = [[0.45, 0.8, 'N'], [1.14, 0.25, 'W'], [-0.12, 0.25, 'E']]
        # title_pos =0.1
        legend_pos = [0.3,.95]
        grid_lines, grid_labels = plt.thetagrids(range(-90, 91, 30))

    ### LEGEND SETTING ###
    if legend:
        ax.legend(markerscale=30, fontsize='small', bbox_to_anchor=legend_pos, bbox_transform=plt.gcf().transFigure)

    ### DIRECTION TEXT ###
    if text:
        for i in range(len(tex_info)):
            ax.text(tex_info[i][0],tex_info[i][1],tex_info[i][2], horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, fontsize='large')

    ax.set_rlim(0)
    ax.grid(True, lw=0.4, zorder=0, alpha=0.5, ls=':')

    return ax

def polarization_contour(input, bkg_input, output=None, aspect=True, ax=None, levels=None, sigma_unit=None,legend=False, text=False, rmax=None, global_color=None):
    src_pcube = xBinnedPolarizationCube.from_file_list(input)
    bkg_pcube = xBinnedPolarizationCube.from_file_list(bkg_input)
    print("SRC_COUNTS", float(src_pcube.COUNTS), "BKG_COUNTS", float(bkg_pcube.COUNTS))

    bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()
    src_pcube -= bkg_pcube
    print('BACKSCALE_SRC', src_pcube.backscal(), 'BACKSCALE_BKG', bkg_pcube.backscal(), 'BACKSCALE_RATIO', np.round(src_pcube.backscal() / bkg_pcube.backscal(), 2))

    EBINS_NUM=  src_pcube.hdu_list[1].header['NAXIS2']
    for num in range(EBINS_NUM):
        print("EBINS_NUM", f'{num+1} / {EBINS_NUM}')
        ax_out = POL_CONTOUR_PLOT(float(src_pcube.QN[num]), float(src_pcube.UN[num]), float(src_pcube.QN_ERR[num]), float(src_pcube.UN_ERR[num]), aspect=aspect, ax=ax, levels=levels, sigma_unit=sigma_unit, legend=legend, text=text, rmax=rmax, global_color=global_color)

    if output:
        plt.savefig(output, dpi=300, transparent=True, bbox_inches='tight')

    return ax_out



######## EXAMPLE ########
# src_list =  [PCUBE_DU1.FITS, PCUBE_DU2.FITS, PCUBE_DU3.FITS]
# bkg_list =  [PCUBE_DU1_BKG.FITS, PCUBE_DU2_BKG.FITS, PCUBE_DU3_BKG.FITS]
# polarization_contour(src_list, bkg_list, aspect=False, rmax=20, global_color='k')
# plt.show()

