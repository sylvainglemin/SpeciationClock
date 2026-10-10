# -*- coding: utf-8 -*-
"""

# Script to infer simple 2D demographic models from GBS data
# Speciation clock project
# Sylvain Glémin - December 2025
# Same scenario but adding inbreeding with F fixed from the data
"""

import os
import sys

# NLopt is the optimization libary dadi uses
import nlopt
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
#import math
import numpy
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

# To extract list of files
import glob

import random

# To read and manipulate csv files
import csv





################
# SCENARIOS ####
################



def split(params, ns, pts):
    """
    params = (nu1,nu2,T)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

	F1: Inbreeding coefficient of population 1
	F2: Inbreeding coefficient of population 2
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    F1,F2,nu1,nu2,T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))

    return fs


def split_mig(params, ns, pts):
    """
    params = (nu1,nu2,T,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.
    
	F1: Inbreeding coefficient of population 1
	F2: Inbreeding coefficient of population 2
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    F1,F2,nu1,nu2,T,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))

    return fs


def split_ancmig(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

	F1: Inbreeding coefficient of population 1
	F2: Inbreeding coefficient of population 2    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    F1,F2,nu1,nu2,T0,T1,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=m, m21=m)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))

    return fs

def sec_contact(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

	F1: Inbreeding coefficient of population 1
	F2: Inbreeding coefficient of population 2   nu1: Size of population 1 after split.
    nu2: Size of population 2 after split."
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    F1,F2,nu1,nu2,T0,T1,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))

    return fs

def mig_twoepoch(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m0,m1)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    F1: Inbreeding coefficient of population 1
    F2: Inbreeding coefficient of population 2  nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    F1,F2,nu1,nu2,T0,T1,m0,m1 = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=m0, m21=m0)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m1, m21=m1)
    
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    
    return fs


def perturbation(demo_model_ex,params,lower_bound,upper_bound,ll_max,popt_max,rep1=30,rep2=20,rep3=10,rep4=10,MAXEVAL=100000,FTOL=1e-5,XTOL=1e-5):
    # First round of perturbation - strong perturbation
    random.seed(1975)
    for i in range(rep1):
        p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)
        p0[0] = F1 # No perurbation F
        p0[1] = F2
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
        									lower_bound=lower_bound, upper_bound=upper_bound,
        									ftol_abs=FTOL, xtol_abs=XTOL,
        									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
        if ll_model > ll_max:
            ll_max = ll_model
            popt_max = popt
        print(ll_max)
    # Second round of perturbation - strong perturbation starting from the best estimate of the previous round
    params = popt_max
    for i in range(rep2):
        p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)
        p0[0] = F1
        p0[1] = F2
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
        									lower_bound=lower_bound, upper_bound=upper_bound,
        									ftol_abs=FTOL, xtol_abs=XTOL,
        									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
        if ll_model > ll_max:
            ll_max = ll_model
            popt_max = popt
        print(ll_max)
    # # Third round of perturbation - intermediate perturbation, of the previous best result
    params = popt_max
    for i in range(rep3):
        p0 = dadi.Misc.perturb_params(popt_max, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
        p0[0] = F1
        p0[1] = F2        
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
        									lower_bound=lower_bound, upper_bound=upper_bound,
        									ftol_abs=FTOL, xtol_abs=XTOL,
        									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
        if ll_model > ll_max:
            ll_max = ll_model
            popt_max = popt
        print(ll_max)
    # # Fourth round of perturbation - weak perturbation of the previous best result
    params = popt_max
    for i in range(rep4):
        p0 = dadi.Misc.perturb_params(popt_max, fold=0.5, upper_bound=upper_bound,lower_bound=lower_bound)
        p0[0] = F1
        p0[1] = F2
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l,
        									lower_bound=lower_bound, upper_bound=upper_bound,
        									ftol_abs=FTOL, xtol_abs=XTOL,
        									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
        if ll_model > ll_max:
            ll_max = ll_model
            popt_max = popt
        print(ll_max)
    # Calculating theta
    model_fs = demo_model_ex(popt_max, ns, pts_l)
    theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
    # Bootstrap and standard error
    # all_boot =[fs.sample() for i in range(100)]
    # uncert = dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, all_boot, popt_max, fs)
    # Use of Fisher information matrix instead
    uncert = dadi.Godambe.FIM_uncert(demo_model_ex, pts_l, popt_max, fs)
    return [ll_max,theta,popt_max,uncert,model_fs]




############################
# MAIN - RUNING THE MODELS #
############################


#os.chdir('/Users/sglemin/Documents/Boulot/Recherche/Thematiques/SystRepro/SPECIATIONCLOCK/ArticleData/SpeciationClock/data_analyses/data/sfs')
os.chdir('/scratch/sglemin/speciation_clock/sfs/free_inbreeding')

MIN = 0.01 #Minimum boundary for parameters
MAX = 50 #Maximum boundary for parameters
MAXEVAL = 100000
FTOL = 1e-5
XTOL = 1e-5

# To send file in argument and to the loop outside the python program
list_file = sys.argv[1:]

# To do the loop within the Python program
# list_file = glob.glob('*.sfs')
# list_file.sort()

#list_file = list_file[5:len(list_file)]
#list_file = list_file[26:len(list_file)]

# read the csv file with inbreeding coefficients

# Create a dictionary to store the data from the CSV file for quick lookup
file_data = {}
# Open the CSV file and read all rows
with open('F_per_sfs_file.csv', mode='r', encoding='utf-8') as file:
    csv_reader = csv.DictReader(file, delimiter=';')
    # Store each row's data in the dictionary with File as the key
    for row in csv_reader:
        file_name = row['File']
        f1_value = float(row['F1'])
        f2_value = float(row['F2'])
        file_data[file_name] = {'F1': f1_value, 'F2': f2_value}

for f in list_file:

    FILENAME = f.split('.')[0]
    print(FILENAME)
    # Loading the dataset
    fs = dadi.Spectrum.from_file(FILENAME+'.sfs') + 0.01 # to avoid problem with numbers < 1
    fs = numpy.ceil(fs) # to avoid problem with numbers < 1
    ns = fs.sample_sizes
    print(ns)
    # Loading the Fis values
    if ns[0] % 2 > 0:
        n0 = ns[0] - 1
        if ns[1] % 2 > 0:
            n1 = ns[1] - 1
            fs = fs.project([n0, n1])
        else:
            n1 = ns[1]
            fs = fs.project([n0, n1])
    else:
        if ns[1] % 2 > 0:
            n0 = ns[0]
            n1 = ns[1] - 1
            fs = fs.project([n0,n1])

    ns = fs.sample_sizes

    F1 = min(file_data[f]['F1'],0.9) # To avoid to high F values that may leads to numerical issues
    F2 = min(file_data[f]['F2'],0.9)
    
    # Define the grid points based on the sample size.
    # For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
    # For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]
    
    # Statistics on the dataset: help for defining starting points
    fs0 = fs.marginalize([1])
    fs1 = fs.marginalize([0])
    pi0 = fs0.pi()
    pi1 = fs1.pi()
    D0 = fs0.Tajima_D()
    D1 = fs1.Tajima_D()
    fst = fs.Fst()

    #Creating the output file"
    output = open(FILENAME+'_free_inbreeding.out','w')
    output.write('dataset pi0 pi1 D0 D1 fst model lnL theta F1 F2 N1 N2 Tsplit Tinter Manc Mrec stheta sF1 sF2 sN1 sN2 sTsplit sTinter sManc sMrec\n')
    output.close()

    #################################
    #### DIVERGENCE IN ALLOPATRY ####
    #################################

    demo_model = split
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = [F1,F2,1, 1, 1]
    lower_bound = [MIN,MIN,MIN, MIN, MIN]
    upper_bound = [0.95,0.95,MAX, MAX, MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l,
    									lower_bound=lower_bound, upper_bound=upper_bound,
    									ftol_abs=FTOL, xtol_abs=XTOL,
    									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
    ll_max = ll_model
    popt_max = popt
    
    # Optimization with iterative perturbation of the starting point
    optim = perturbation(demo_model_ex, params, lower_bound, upper_bound, ll_max, popt_max)

    # Data to be exported
    lnL_split = optim[0]
    theta = optim[1]
    popt_split = optim[2].tolist()
    se_split = optim[3].tolist()

    # Writing results in the output file
    a = f'{FILENAME} {pi0} {pi1} {D0} {D1} {fst} split {lnL_split} {theta} '
    b = f'{popt_split[0]} {popt_split[1]} {popt_split[2]} {popt_split[3]} {popt_split[4]} NA 0 0 '
    c = f'{se_split[5]} {se_split[0]} {se_split[1]} {se_split[2]} {se_split[3]} {se_split[4]} NA 0 0\n'
    output = open(FILENAME+'_free_inbreeding.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = True)
    plt.savefig(FILENAME+'_free_inbreeding_allo.png',dpi=200)
    
    del(optim,theta,a,b,c)


    ###############################
    #### SYMMETRICAL MIGRATION ####
    ###############################

    demo_model = split_mig
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = [F1,F2] + popt_split[0:2] + [min(2*popt_split[2],50)] + [1] # As a starting point: double divergence time and add migration
    #params = [1,1,1,1]
    lower_bound = [MIN,MIN,MIN,MIN,MIN,MIN]
    upper_bound = [0.95,0.95,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        ftol_abs=FTOL, xtol_abs=XTOL,
    									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
    ll_max = ll_model
    popt_max = popt

    # Optimization with iterative perturbation of the starting point
    optim = perturbation(demo_model_ex, params, lower_bound, upper_bound, ll_max,popt_max)

    # Data to be exported
    lnL_splitmig = optim[0]
    theta = optim[1]
    popt_splitmig = optim[2].tolist()
    se_splitmig = optim[3].tolist()

    # Writing results in output file
    a = f'{FILENAME} {pi0} {pi1} {D0} {D1} {fst} split_mig {lnL_splitmig} {theta} '
    b = f'{popt_splitmig[0]} {popt_splitmig[1]} {popt_splitmig[2]} {popt_splitmig[3]} {popt_splitmig[4]} NA {popt_splitmig[5]} {popt_splitmig[5]}  '
    c = f'{se_splitmig[6]} {se_splitmig[0]} {se_splitmig[1]} {se_splitmig[2]} {se_splitmig[3]} {se_splitmig[4]} NA {se_splitmig[5]} {se_splitmig[5]} \n'
    output = open(FILENAME+'_free_inbreeding.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_free_inbreeding_splitmig.png',dpi=200)
    
    del(optim,theta,a,b,c)
    
    #############################
    #### ANCESTRAL MIGRATION ####
    #############################

    demo_model = split_ancmig
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = [F1,F2] + popt_splitmig[0:2] + [max(popt_splitmig[2]/2,2e-2),max(popt_splitmig[2]/2,2e-2)] + [popt_splitmig[3]] # Keep the same total divergence time but split into two epochs with and without migration
    #params = [1,1,0.5,0.5,1]
    lower_bound = [MIN,MIN,MIN,MIN,MIN,MIN,MIN]
    upper_bound = [0.95,0.95,MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        ftol_abs=FTOL, xtol_abs=XTOL,
    									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
    popt_max = popt

    # Optimization with iterative perturbation of the starting point
    optim = perturbation(demo_model_ex, params, lower_bound, upper_bound, ll_max, popt_max)
        
    # Data to be exported
    lnL_splitancmig = optim[0]
    theta = optim[1]
    popt_splitancmig = optim[2].tolist()
    se_splitancmig = optim[3].tolist()

    # Writing results in output file
    a = f'{FILENAME} {pi0} {pi1} {D0} {D1} {fst} split_ancmig {lnL_splitancmig} {theta} '
    b = f'{popt_splitancmig[0]} {popt_splitancmig[1]} {popt_splitancmig[2]} {popt_splitancmig[3]} {popt_splitancmig[4]} {popt_splitancmig[5]} {popt_splitancmig[6]} 0 '
    c = f'{se_splitancmig[7]} {se_splitancmig[0]} {se_splitancmig[1]} {se_splitancmig[2]} {se_splitancmig[3]} {se_splitancmig[4]} {se_splitancmig[5]} {se_splitancmig[6]} 0\n'
    output = open(FILENAME+'_free_inbreeding.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_free_inbreeding_splitancmig.png',dpi=200)
    
    del(optim,theta,a,b,c)

    ###########################
    #### SECONDARY CONTACT ####
    ###########################

    demo_model = sec_contact
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = popt_splitancmig # Keep the same estimates as a starting point
    #params = [1,1,0.5,0.5,1]
    lower_bound = [MIN,MIN,MIN,MIN,MIN,MIN,MIN]
    upper_bound = [0.95,0.95,MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        ftol_abs=FTOL, xtol_abs=XTOL,
    									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
    ll_max = ll_model
    popt_max = popt

    # Optimization with iterative perturbation of the starting point
    optim = perturbation(demo_model_ex, params, lower_bound, upper_bound, ll_max, popt_max)

    # Data to be exported
    lnL_seccontact = optim[0]
    theta = optim[1]
    popt_seccontact = optim[2].tolist()
    se_seccontact = optim[3].tolist()

    # Writing results in output files
    a = f'{FILENAME} {pi0} {pi1} {D0} {D1} {fst} sec_contact {lnL_seccontact} {theta} '
    b = f'{popt_seccontact[0]} {popt_seccontact[1]} {popt_seccontact[2]} {popt_seccontact[3]} {popt_seccontact[4]} {popt_seccontact[5]} 0 {popt_seccontact[6]} '
    c = f'{se_seccontact[7]} {se_seccontact[0]} {se_seccontact[1]} {se_seccontact[2]} {se_seccontact[3]} {se_seccontact[4]} {se_seccontact[5]} 0 {se_seccontact[6]}\n'
    output = open(FILENAME+'_free_inbreeding.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_free_inbreeding_seccontact.png',dpi=200)
    
    del(optim,theta,a,b,c)

    ####################
    #### TWO EPOCHs ####
    ####################

    demo_model = mig_twoepoch
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = popt_seccontact + [popt_seccontact[4]]
    #params = [1,1,0.5,0.5,1,1]
    lower_bound = [MIN,MIN,MIN,MIN,MIN,MIN,MIN,MIN]
    upper_bound = [0.95,0.95,MAX,MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        ftol_abs=FTOL, xtol_abs=XTOL,
    									algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=1)
    ll_max = ll_model
    popt_max = popt

    # Optimization with iterative perturbation of the starting point
    optim = perturbation(demo_model_ex, params, lower_bound, upper_bound, ll_max, popt_max)

    # Data to be exported
    lnL_twoepoch = optim[0]
    theta = optim[1]
    popt_twoepoch = optim[2].tolist()
    se_twoepoch = optim[3].tolist()

    # Writing results in output file
    a = f'{FILENAME} {pi0} {pi1} {D0} {D1} {fst} two_epoch {lnL_twoepoch} {theta} '
    b = f'{popt_twoepoch[0]} {popt_twoepoch[1]} {popt_twoepoch[2]} {popt_twoepoch[3]} {popt_twoepoch[4]} {popt_twoepoch[5]} {popt_twoepoch[6]} {popt_twoepoch[7]} '
    c = f'{se_twoepoch[8]} {se_twoepoch[0]} {se_twoepoch[1]} {se_twoepoch[2]} {se_twoepoch[3]} {se_twoepoch[4]} {se_twoepoch[5]} {se_twoepoch[6]} {se_twoepoch[7]}\n'
    output = open(FILENAME+'_free_inbreeding.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_free_inbreeding_twoepoch.png',dpi=200)
    
    del(optim,theta,a,b,c)

exit(0)

