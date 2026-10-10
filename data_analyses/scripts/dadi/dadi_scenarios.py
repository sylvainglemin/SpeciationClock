# -*- coding: utf-8 -*-
"""

# Script to infer simple 2D demographic models from GBS data
# Speciation clock project
# Sylvain Glémin - September 2022
"""

import os
# Pickle is used to save variables as files for future use --> Maybe note usefull for me
# import pickle
# NLopt is the optimization libary dadi uses
import nlopt
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
#import math
import dadi
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

# To extract list of files
import glob

import random
random.seed(1975)




################
# SCENARIOS ####
################



def split(params, ns, pts):
    """
    params = (nu1,nu2,T)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_mig(params, ns, pts):
    """
    params = (nu1,nu2,T,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def split_ancmig(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T0,T1,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=m, m21=m)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def sec_contact(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T0,T1,m = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def mig_twoepoch(params, ns, pts):
    """
    params = (nu1,nu2,T0,T1,m0,m1)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T0: Time between the split and the arrest of migration (in units of 2*Na generations) 
    T1: Time between the arrest of migration and present
    The Time of the split before present is T = T0+T1
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T0,T1,m0,m1 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T0, nu1, nu2, m12=m0, m21=m0)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m1, m21=m1)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def perturbation(demo_model_ex,params,lower_bound,upper_bound,ll_max,popt_max,rep1=20,rep2=20,rep3=10,rep4=10,MAXEVAL=20):
    # First round of perturbation - strong perturbation
    # for i in range(rep1):
    #     p0 = dadi.Misc.perturb_params(params, fold=3, upper_bound=upper_bound,lower_bound=lower_bound)
    #     popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=100)
    #     if ll_model > ll_max:
    #         ll_max = ll_model
    #         popt_max = popt
    #     print(ll_max)
    # Second round of perturbation - intermediate perturbation, starting from the best estimate of the previous round
    # params = popt_max
    for i in range(rep2):
        p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=100)
        if ll_model > ll_max:
            ll_max = ll_model
            popt_max = popt
        print(ll_max)
    # # Third round of perturbation - intermediate perturbation, of the previous best result
    # for i in range(rep3):
    #     p0 = dadi.Misc.perturb_params(popt_max, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)
    #     popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=100)
    #     if ll_model > ll_max:
    #         ll_max = ll_model
    #         popt_max = popt
    #     print(ll_max)
    # # Fourth round of perturbation - weak perturbation, of the previous best result
    # for i in range(rep4):
    #     p0 = dadi.Misc.perturb_params(popt_max, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
    #     popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=MAXEVAL, verbose=100)
    #     if ll_model > ll_max:
    #         ll_max = ll_model
    #         popt_max = popt
    #     print(ll_max)
    # Calculating theta
    model_fs = demo_model_ex(popt_max, ns, pts_l)
    theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
    # Bootstrap and standard error
    all_boot =[fs.sample() for i in range(100)]
    uncert = dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, all_boot, popt_max, fs)
    return [ll_max,theta,popt_max,uncert,model_fs]




#####################
# RUNING THE MODELS #
#####################


# os.chdir('/Users/sglemin/Documents/Boulot/Recherche/Thematiques/SystRepro/SPECIATIONCLOCK/Abel/Data/2022_09_30_SFS_SpCl_GBS')
os.chdir('/scratch/sglemin/speciation_clock')


MIN = 1e-3 #Minimum boundary for parameters
MAX = 100 #Maximum boundary for parameters

list_file = glob.glob('*.sfs')
list_file.sort(reverse=True)
#list_file = list_file[5:len(list_file)]

list_file = list_file[26:len(list_file)]


for f in list_file:

    FILENAME = f.split('.')[0]
    print(FILENAME)
    # Loading the dataset
    fs = dadi.Spectrum.from_file(FILENAME+'.sfs')
    ns = fs.sample_sizes
    
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
    output = open(FILENAME+'.out','w')
    output.write('dataset pi0 pi1 D0 D1 fst model lnL theta N0 N1 Tsplit Tinter Manc Mrec stheat sN0 sN1 sTsplit sTinter sManc sMrec\n')
    output.close()

    #################################
    #### DIVERGENCE IN ALLOPATRY ####
    #################################

    demo_model = split
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = [1, 1, 1]
    lower_bound = [MIN, MIN, MIN]
    upper_bound = [MAX, MAX, MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=20, verbose=1)
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
    b = f'{popt_split[0]} {popt_split[1]} {popt_split[2]} NA 0 0 '
    c = f'{se_split[3]} {se_split[0]} {se_split[1]} {se_split[2]} NA 0 0\n'
    output = open(FILENAME+'.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = True)
    plt.savefig(FILENAME+'_allo.png',dpi=200)
    
    del(optim,theta,a,b,c)

    ###############################
    #### SYMMETRICAL MIGRATION ####
    ###############################

    demo_model = split_mig
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = popt_split[0:2] + [2*popt_split[2]] + [1] # As a starting point: double divergence time and add migration
    #params = [1,1,1,1]
    lower_bound = [MIN,MIN,MIN,MIN]
    upper_bound = [MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=20, verbose=1)
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
    b = f'{popt_splitmig[0]} {popt_splitmig[1]} {popt_splitmig[2]} NA {popt_splitmig[3]} {popt_splitmig[3]} '
    c = f'{se_splitmig[4]} {se_splitmig[0]} {se_splitmig[1]} {se_splitmig[2]} NA {se_splitmig[3]} {se_splitmig[3]}\n'
    output = open(FILENAME+'.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_splitmig.png',dpi=200)
    
    del(optim,theta,a,b,c)
    
    #############################
    #### ANCESTRAL MIGRATION ####
    #############################

    demo_model = split_ancmig
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # Starting parameters and boundaries
    params = popt_splitmig[0:2] + [max(popt_splitmig[2]/2,2e-2),max(popt_splitmig[2]/2,2e-2)] + [popt_splitmig[3]] # Keep the same total divergence time but split into two epochs with and without migration
    #params = [1,1,0.5,0.5,1]
    lower_bound = [MIN,MIN,MIN,MIN,MIN]
    upper_bound = [MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=20, verbose=1)
    ll_max = ll_model
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
    b = f'{popt_splitancmig[0]} {popt_splitancmig[1]} {popt_splitancmig[2]} {popt_splitancmig[3]} {popt_splitancmig[4]} 0 '
    c = f'{se_splitancmig[5]} {se_splitancmig[0]} {se_splitancmig[1]} {se_splitancmig[2]} {se_splitancmig[3]} {se_splitancmig[4]} 0\n'
    output = open(FILENAME+'.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_splitancmig.png',dpi=200)
    
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
    lower_bound = [MIN,MIN,MIN,MIN,MIN]
    upper_bound = [MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=20, verbose=1)
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
    b = f'{popt_seccontact[0]} {popt_seccontact[1]} {popt_seccontact[2]} {popt_seccontact[3]} 0 {popt_seccontact[4]} '
    c = f'{se_seccontact[5]} {se_seccontact[0]} {se_seccontact[1]} {se_seccontact[2]} {se_seccontact[3]} 0 {se_seccontact[4]}\n'
    output = open(FILENAME+'.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_seccontact.png',dpi=200)
    
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
    lower_bound = [MIN,MIN,MIN,MIN,MIN,MIN]
    upper_bound = [MAX,MAX,MAX,MAX,MAX,MAX]

    # First estimation with the initial starting point 
    popt, ll_model = dadi.Inference.opt(params, fs, demo_model_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, algorithm=nlopt.LN_BOBYQA, maxeval=20, verbose=1)
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
    b = f'{popt_twoepoch[0]} {popt_twoepoch[1]} {popt_twoepoch[2]} {popt_twoepoch[3]} {popt_twoepoch[4]} {popt_twoepoch[5]} '
    c = f'{se_twoepoch[6]} {se_twoepoch[0]} {se_twoepoch[1]} {se_twoepoch[2]} {se_twoepoch[3]} {se_twoepoch[4]} {se_twoepoch[5]}\n'
    output = open(FILENAME+'.out','a')
    output.write(a+b+c)
    output.close()

    # Saving plot
    fs_opt = optim[4]
    dadi.Plotting.plot_2d_comp_multinom(optim[4], fs,vmin=fs_opt.min()+ 0.01,vmax=fs_opt.max(),show = False)
    plt.savefig(FILENAME+'_twoepoch.png',dpi=200)
    
    del(optim,theta,a,b,c)




exit(0)

