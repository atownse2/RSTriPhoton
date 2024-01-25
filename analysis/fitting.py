import os
import sys

import json

import multiprocessing

import numpy as np
import awkward as ak
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

from .utils import sample_info as si
from .utils import selections as sel

top = si.top_dir

# Save signal events into arrays for unbinned fits
array_dir = top+'/outputs/arrays/preselection/'
def get_arrays(signal_point: tuple,
                tuple_version=si.get_tuple_version(),
                remake=False):
    
    from coffea.nanoevents import NanoEventsFactory, BaseSchema

    # signal_point = si.get_mass_point(signal_filetag)
    signal_filetag = si.get_signal_filetag(signal_point)

    out_file = array_dir + f'{signal_filetag}_{tuple_version}.parquet'
    if os.path.exists(out_file) and not remake:
        return ak.from_parquet(out_file)
    
    # print(f'Processing {signal_filetag}')

    ## Load events
    files = si.get_filelist('signal', subset=signal_filetag, tuple_version=tuple_version)

    fdict = {fn : "/flattener/tree" for fn in files}
    # fdict = {filepath : "/flattener/tree"}

    events = NanoEventsFactory.from_root(
    fdict,
    schemaclass=BaseSchema,
    metadata={'dataset': 'signal'},
    permit_dask=True,
    ).events()

    candidate = sel.get_candidates(events)

    out_arr = ak.zip({'triphoton_mass': candidate.triphoton.mass,
                      'diphoton_mass': candidate.diphoton.mass,
                      'alpha': candidate.alpha,})

    # Save arrays
    ak.to_parquet(out_arr, out_file)
    return out_arr

def get_all_arrays(tuple_version=si.get_tuple_version,
                   mass_grid=si.get_mass_grid(),
                   remake=False):
    from coffea.nanoevents import NanoEventsFactory, BaseSchema

    arguments = [(signal_point, tuple_version, remake) for signal_point in mass_grid]
    with multiprocessing.Pool(9) as pool:
        results = pool.starmap(get_arrays, arguments)
    return results

# Fitting
fit_dir = top + '/plots/signal/fits'
fit_info_dir = top + '/analysis/metadata/'
def arr_to_tree(arr):
    """ Converts a flat awkward array to a TTree """
    import ROOT
    from array import array

    tree = ROOT.TTree('tree', 'tree')
    x = array('d', [0.0])
    tree.Branch('x', x, 'x/D')
    for val in arr:
        x[0] = val
        tree.Fill()
    return tree

def crystal_ball_fit(awk_arr, name=None, save_fit=False):
    """
    Do a crystal ball fit on the given array.
    """
    import ROOT

    min_val = 1e-1
    max_val = 1e7

    arr = awk_arr.to_numpy()
    if len(arr) == 0:
        print(f'No data for fit: {name}')
        return
    
    mean = arr.mean()
    std = arr.std()

    # Estimate mode
    bin_values, bin_edges = np.histogram(arr, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    mode = bin_centers[np.argmax(bin_values)]

    # Set mean to mode and set std using values around the mode
    mean = mode

    arr = arr[(arr > mode - 2*std) & (arr < mode + 2*std)]
    std = np.std(arr[(arr > mode - 1*std) & (arr < mode + 1*std)])
    
    tree = arr_to_tree(arr)
    xrange = (min(arr), max(arr))

    x = ROOT.RooRealVar('x', 'x', *xrange)
    mean = ROOT.RooRealVar('mean', 'mean', mean, 0.5*mean, 1.5*mean)
    sigma = ROOT.RooRealVar('sigma', 'sigma', 0.8*std, 0.2*std, 1.3*std)
    alphaL = ROOT.RooRealVar('alphaL', 'alphaL', 1, min_val, 50)
    nL = ROOT.RooRealVar('nL', 'nL', 2, min_val, 100)
    alphaR = ROOT.RooRealVar('alphaR', 'alphaR', 1, min_val, 50)
    nR = ROOT.RooRealVar('nR', 'nR', 2, min_val, 100)
    # Create RooFit crystal ball
    cb = ROOT.RooCrystalBall('cb', 'cb', x, mean, sigma, alphaL, nL, alphaR, nR)
    # Create RooFit dataset
    data = ROOT.RooDataSet('data', 'data', tree, ROOT.RooArgSet(x))
    # Fit
    # Documentation https://root.cern.ch/doc/master/classRooAbsPdf.html#ab0721374836c343a710f5ff92a326ff5
    result = cb.fitTo(data, ROOT.RooFit.Minimizer("Minuit2", "minimize"), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1));

    # Get CI
    cdf = cb.createCdf(ROOT.RooArgSet(x))

    cdf_min = 0.025
    cdf_max = 0.975

    lower = cdf.findRoot(x, *xrange, cdf_min)
    upper = cdf.findRoot(x, *xrange, cdf_max)

    # Plot
    frame = x.frame()

    # Set uniform binning
    data.plotOn(frame, ROOT.RooFit.Binning(64))
    cb.plotOn(frame)

    # Vertical lines for lower, upper
    frame.addObject(ROOT.TLine(lower, 0, lower, 1e7))
    frame.addObject(ROOT.TLine(upper, 0, upper, 1e7))

    if save_fit:
        c = ROOT.TCanvas()
        frame.Draw()

        # Put fit parameters in legend
        legend = ROOT.TLegend(0.1, 0.4, 0.4, 0.9)
        legend.AddEntry(0, f'mean: {mean.getVal():.2f}', '')
        legend.AddEntry(0, f'sigma_range: {sigma.getMin():.2f} - {sigma.getMax():.2f}', '')
        legend.AddEntry(0, f'sigma: {sigma.getVal():.2f}', '')
        legend.AddEntry(0, f'alphaL: {alphaL.getVal():.2f}', '')
        legend.AddEntry(0, f'nL: {nL.getVal():.2f}', '')
        legend.AddEntry(0, f'alphaR: {alphaR.getVal():.2f}', '')
        legend.AddEntry(0, f'nR: {nR.getVal():.2f}', '')
        legend.AddEntry(0, f'nll: {result.minNll():.2f}', '')


        legend.Draw()

        c.SaveAs(f'{fit_dir}/{name}_crystal_ball.png')
    
    info = {'lower': lower, 'upper': upper, 'mean': mean.getVal(), 'sigma': sigma.getVal(), 'n': len(arr)}

    return info

def do_signal_fit(signal_point: tuple,
                  tuple_version=si.get_tuple_version()):
    
    signal_filetag = si.get_signal_filetag(signal_point)

    arr = ak.from_parquet(array_dir + f'{signal_filetag}_{tuple_version}.parquet')

    # Do fits
    triphoton_fit = crystal_ball_fit(arr.triphoton_mass, name=f'{signal_filetag}_triphoton_mass', save_fit=True);
    diphoton_fit = crystal_ball_fit(arr.diphoton_mass, name=f'{signal_filetag}_diphoton_mass', save_fit=True)

    fits = {'BKK Mass': signal_point[0], 'Radion Mass': signal_point[1]}
    fits.update({f"triphoton_mass {key}": val for key, val in triphoton_fit.items()})
    fits.update({f"diphoton_mass {key}": val for key, val in diphoton_fit.items()})

    return fits

def get_signal_fit_info(tuple_version=si.get_tuple_version(),
                        mass_grid=si.get_mass_grid(),
                        refit=False):
    outfile = fit_info_dir + f'signal_fits.csv'
    if os.path.exists(outfile) and not refit:
        return pd.read_csv(outfile)

    import ROOT

    arguments = [(signal_point, tuple_version) for signal_point in mass_grid]
    with multiprocessing.Pool(6) as pool:
        results = pool.starmap(do_signal_fit, arguments)
    
    # Save fit results to csv
    df = pd.DataFrame(results)
    df.to_csv(outfile, index=False)

    return df

# Regress upper and lower CL on masses
def regress_signal_CL(tuple_version=si.get_tuple_version(), refit=False):
    outfile = fit_info_dir + f'signal_regression_coeffs.json'
    if os.path.exists(outfile) and not refit:
        # Load json
        with open(outfile, 'r') as f:
            regression_coeffs = json.load(f)
        return regression_coeffs

    from sklearn.linear_model import LinearRegression

    # Load fit info
    fit_info = get_signal_fit_info(tuple_version=tuple_version, refit=refit)

    regression_coeffs = {'triphoton_mass': {}, 'diphoton_mass': {}}

    for var_name, coef_dict in regression_coeffs.items():
        if var_name == 'triphoton_mass':
            v = 'BKK Mass'
        elif var_name == 'diphoton_mass':
            v = 'Radion Mass'

        X = fit_info[v]

        for bound in ['lower', 'upper']:
            y = fit_info[f'{var_name} {bound}']

            # Fit
            reg = LinearRegression().fit(X.values.reshape(-1, 1), y)

            # Save coefficients
            coef_dict[bound] = (reg.intercept_, reg.coef_[0])
    
    # Save to json
    with open(outfile, 'w') as f:
        json.dump(regression_coeffs, f, indent=4)
    
    return regression_coeffs

def plot_signal_CL_regression(tuple_version=si.get_tuple_version(),
                              refit=False, save=True):
    # Load fit info
    fit_info = get_signal_fit_info(tuple_version=tuple_version, refit=refit)

    # Load regression coefficients
    regression_coeffs = regress_signal_CL(tuple_version=tuple_version, refit=refit)

    # Plot
    for var_name, coef_dict in regression_coeffs.items():
        fig, ax = plt.subplots(2, 2, figsize=(12, 12),
                               gridspec_kw={'height_ratios': [3, 1]},
                               sharex=True)
    
        if var_name == 'triphoton_mass':
            v = 'BKK Mass'
        elif var_name == 'diphoton_mass':
            v = 'Radion Mass'

        X = fit_info[v]

        for j, bound in enumerate(['lower', 'upper']):
            y = fit_info[f'{var_name} {bound}']

            # Plot simulation
            ax[0, j].scatter(X, y, label='Simulation', color='black')

            # Plot regression
            x_range = (min(X), max(X))
            x = np.linspace(*x_range, 100)

            fn = lambda x: coef_dict[bound][0] + coef_dict[bound][1] * x
            ax[0, j].plot(x, fn(x), label='Regression', color='blue')
            ax[0, j].set_title(f'{var_name.capitalize()} {bound.capitalize()} 95% Crystal Ball CL',
                                 fontsize=20)
            
            # Plot residuals
            residuals = y - fn(X)

            ax[1, j].scatter(X, residuals, label='Residuals', color='red')
            ax[1, j].set_xlabel(f"{var_name.capitalize()} Mass", fontsize=16)
            ax[1, j].set_ylabel("Residuals", fontsize=16)

        fig.tight_layout()
        fig.savefig(f'{fit_dir}/signal_CL_fits.png')

    return

def get_signal_CL(signal_point,
                  tuple_version=si.get_tuple_version(),
                  regression=False,
                  refit=False):
    
    # Get mass values
    BKK_mass = signal_point[0]
    Radion_mass = signal_point[1]

    CL = {}
    if regression:
        regression_coeffs = regress_signal_CL(tuple_version=tuple_version, refit=refit)
    else:
        fit_info = get_signal_fit_info(tuple_version=tuple_version, refit=refit)
        signal_point_info = fit_info[(fit_info['BKK Mass']==BKK_mass) & (fit_info['Radion Mass']==Radion_mass)]

    for i, var_name in enumerate(['triphoton_mass', 'diphoton_mass']):
        CL[var_name] = {}
        for bound in ['lower', 'upper']:
            if regression:
                c = regression_coeffs[var_name][bound]
                CL[var_name][bound] = c[0] + c[1]*signal_point[i]

            else:
                CL[var_name][bound] = signal_point_info[f'{var_name} {bound}'].values[0]
    
    return CL

def plot_signal_CL(signal_points,
                   tuple_version=si.get_tuple_version(),
                   regression=False,
                   refit=False,
                   save=True):
    if isinstance(signal_points, tuple):
        signal_points = [signal_points]

    for var_name in ['triphoton_mass', 'diphoton_mass']:
        fig, ax = plt.subplots(figsize=(8, 8))

        x_min = 9000
        x_max = 0
        for i, signal_point in enumerate(signal_points):
            CL = get_signal_CL(signal_point,
                    tuple_version=tuple_version,
                    regression=regression,
                    refit=refit)
            
            bounds = CL[var_name]
            
            vals = get_arrays(signal_point, tuple_version=tuple_version, remake=False)[var_name].to_numpy()

            x_min = min(x_min, min(vals))
            x_max = max(x_max, max(vals))

            ax.hist(vals, bins=100, label=f'{signal_point}', histtype='step', color=f'C{i}')
            ax.axvline(bounds['lower'], color=f'C{i}', linestyle='--')
            ax.axvline(bounds['upper'], color=f'C{i}', linestyle='--')
        
        ax.set_title(f'{var_name.replace("_", " ")} 95% CL', fontsize=20)
        ax.set_xlabel(f'{var_name.replace("_", " ")}', fontsize=16)
        ax.legend()
        ax.set_xlim(x_min, x_max)

    
        if save:
            fig.savefig(f'{fit_dir}/signal_CL_{var_name}.png')

# Calculate signal ellipse for a given signal point using fits by default (not regression)
def calculate_signal_ellipse(signal_point,
                             tuple_version=si.get_tuple_version(),
                             regression=False,
                             refit=False):

    CL = get_signal_CL(signal_point,
                       tuple_version=tuple_version,
                       regression=regression,
                       refit=refit)

    h = 0.5*(CL['triphoton_mass']['upper'] + CL['triphoton_mass']['lower']) # triphoton center
    k = 0.5*(CL['diphoton_mass']['upper'] + CL['diphoton_mass']['lower']) # diphoton center
    a = 0.5*(CL['triphoton_mass']['upper'] - CL['triphoton_mass']['lower']) # triphoton radius
    b = 0.5*(CL['diphoton_mass']['upper'] - CL['diphoton_mass']['lower']) # diphoton radius

    return (h,k,a,b)

def in_ellipse(signal_point, triphoton_mass, diphoton_mass):
    h, k, a, b = calculate_signal_ellipse(signal_point)
    in_it = ((triphoton_mass - h)/a)**2 + ((diphoton_mass - k)/b)**2 <= 1
    return in_it
    

def plot_signal_CL_ellipse(signal_points,
                           tuple_version=si.get_tuple_version(),
                           regression=False,
                           refit=False,):
    if isinstance(signal_points, tuple):
        signal_points = [signal_points]
    
    fig, ax = plt.subplots()

    xrange = (9000, 0)
    yrange = (9000, 0)

    for i, signal_point in enumerate(signal_points):

        # Plot the data first
        arrs = get_arrays(signal_point, tuple_version=tuple_version, remake=False)

        triphoton_mass = arrs.triphoton_mass.to_numpy()
        diphoton_mass = arrs.diphoton_mass.to_numpy()

        ax.scatter(triphoton_mass, diphoton_mass, label=f'{signal_point}', color=f'C{i}', s=1, alpha=0.2)

        # Plot the ellipse
        h, k, a, b = calculate_signal_ellipse(signal_point,
                                              tuple_version=tuple_version,
                                              regression=regression,
                                              refit=refit)

        ellipse = mpl.patches.Ellipse((h, k), 2*a, 2*b, color=f'C{i}', fill=False)
        ax.add_patch(ellipse)
    
        # Update ranges
        xrange = (min(xrange[0], min(triphoton_mass)), max(xrange[1], max(triphoton_mass)))
        yrange = (min(yrange[0], min(diphoton_mass)), max(yrange[1], max(diphoton_mass)))
    
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)

    ax.set_xlabel('Triphoton Mass')
    ax.set_ylabel('Diphoton Mass')

    ax.legend(scatterpoints=10000)

    fig.savefig(f'{fit_dir}/signal_CL_ellipse.png')

