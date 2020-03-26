"""
  Date:     200325
  History:  
    - Originally written as a python notebook by Lukas Heinrich, and adapted for the 3L-RJ likelihood validation by Giordon Stark (https://github.com/kratsg/3L-RJ-mimic-likelihood-validation/blob/master/OutlierPlot.ipynb)
    - Adapted and generalized to script by Danika MacDonell [March 25, 2020]
    
  Details:  Script to visualize the fractional size of the systematics for pyhf likelihoods at each mass point, where the systematics are added together in quadrature. Should be run in a directory containing the background-only json likelihood file, along with a patch json likelihood file for each signal point. Assumes that the model has two variable mass points, and the signals in the json patch files are named as such.
  
  Usage:
  >> python OutlierPlot.py --signal_template <signal_template_{}_{}_for_masses> --v_max <max colourbar amplitude> --x_label <x axis label> --y_label <y axis label>
  
  Example for 1Lbb Wh analysis (https://glance.cern.ch/atlas/analysis/analyses/details.php?id=2969):
  
  >> python OutlierPlot.py --signal_template C1N2_Wh_hbb_{}_{} --v_max 10 --x_label '$m(\tilde{\chi}_{1}^{\pm}/\tilde{\chi}_{2}^{0})$ [GeV]' --y_label '$m(\tilde{\chi}_{1}^{0})$ [GeV]'
  
  The signal template is the name of an arbitrary signal in the json patch files, with the signal masses left as {}. Given a background-only file and a patch file, the signal name can be found under "samples" in the output of:
  >> jsonpatch BkgOnly.json patch_XXX.json | pyhf inspect
  
  If, for example, one of the signals is called C1N2_Wh_hbb_550_200, where 550 and 200 are the variable model masses, the signal template would be C1N2_Wh_hbb_{}_{}.
  """
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from parse import *
import click
import pyhf
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

def handle_deltas(delta_up, delta_dn):
    nom_is_center = np.bitwise_or(
        np.bitwise_and(delta_up > 0, delta_dn > 0),
        np.bitwise_and(delta_up <= 0, delta_dn <= 0)
    )
    span    = (delta_dn + delta_up)
    maxdel  = np.maximum(np.abs(delta_dn),np.abs(delta_dn))
    abs_unc = np.where(nom_is_center, span, maxdel)
    return abs_unc

def process_patch(p):
    nom = np.asarray(p['value']['data'])


    #histosys
    hid = np.asarray([m['data']['hi_data'] for m in p['value']['modifiers'] if m['type'] == 'histosys'])
    lod = np.asarray([m['data']['lo_data'] for m in p['value']['modifiers'] if m['type'] == 'histosys'])
    delta_up = hid-nom
    delta_dn = nom-lod
    A = handle_deltas(delta_up,delta_dn)

    hi = np.asarray([m['data']['hi'] for m in p['value']['modifiers'] if m['type'] == 'normsys'])
    lo = np.asarray([m['data']['lo'] for m in p['value']['modifiers'] if m['type'] == 'normsys'])
    delta_up = np.asarray([h*nom-nom for h in hi])
    delta_dn = np.asarray([l*nom-nom for l in lo])
    B = handle_deltas(delta_up,delta_dn)


    delta = np.asarray([m['data'] for m in p['value']['modifiers'] if m['type'] == 'staterror'])
    delta_up = (delta/2.)
    delta_dn = (delta/2.)
    C = handle_deltas(delta_up,delta_dn)
    C = np.zeros_like(A)

    systs = np.concatenate([A,B,C])
    inquad = np.sqrt(np.sum(np.square(systs), axis=0))
    rel =  inquad/nom
    rel = np.where(nom == 0, np.ones_like(nom), rel)
    return rel,nom

@click.command()
@click.option(
              "--signal_template",
              help="Signal name template, with 2 signal masses as variables (eg. signal_{}_{}_more_info)",
              )
@click.option(
              "--v_max",
              help="Maximum amplitude of the colourbar for plotting interpolated fractional syst sizes",
              default=10,
              required=False
              )
@click.option(
              "--x_label",
              help="x label for interpolated plot of fractional systematics",
              default=None,
              required=False
              )

@click.option(
              "--y_label",
              help="y label for interpolated plot of fractional systematics",
              default=None,
              required=False
              )
def outlier_plot(signal_template, v_max, x_label, y_label):

    patches = [json.load(open(x)) for x in glob.glob('patch*.json')]
    grouped_patches = {x[0]['value']['name']:{p['path']:p for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}


    data = {x[0]['value']['name']:{p['path']:process_patch(p) for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}
    #print(data)

    # Make the mapping of json channel names to analysis region names
    listOfPatches = glob.glob("patch_*.json")
    spec_sig = json.load(open(listOfPatches[0]))
    spec_bkg = json.load(open("BkgOnly.json"))

    names_json=[]
    for ichan, channel in enumerate(spec_bkg['channels']): names_json.append(channel['name'])

    channels_json=[]
    for ichan, channel in enumerate(spec_sig): channels_json.append(channel['path'])

    channel_names = dict(zip(channels_json, names_json))

    # Make the mapping of number of json channel name to number of bins
    workspace_bkg = pyhf.Workspace(json.load(open("BkgOnly.json")))
    channel_bins = {}
    for key, value in channel_names.items():
        channel_bins[key] = workspace_bkg.channel_nbins[value]
        
    outliers = []
    for k,v in data.items():
        for kk,vv in v.items():
            for b,(r,n) in enumerate(zip(*vv)):
                if r > 1.0:
                    outliers.append((k,kk,b,r,n))

    #print("Outliers (> 1.0):")
    #for o in list(reversed(sorted(outliers, key = lambda x: x[-1]))):
    #    print('\t',o[-1],o[-2],o[0],channel_names[o[1]],o[2])
    """
    missing_signal = []
    # missing signal in signal region
    print("Missing signal in signal region:")
    for k, v in data.items():
        if not '/channels/2/samples/5' in v or not '/channels/2/samples/5' in v:
            missing_signal.append(k)
            print('\t',k)
    """


    data = {x[0]['value']['name']:{p['path']:process_patch(p)[0] for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}

    sig_name_template = signal_template

    for ichan, channel in enumerate(channels_json):
        frac_systs = np.asarray(
            [[float(m) for m in [parse(sig_name_template, k)[0], parse(sig_name_template, k)[1]]] + v.get(channel, np.array([0.0]*channel_bins[channel])).tolist() for k,v in data.items()
            ]
        )
        x_min, x_max = min(frac_systs[:,0]), max(frac_systs[:,0])
        y_min, y_max = min(frac_systs[:,1]), max(frac_systs[:,1])
        x,y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        for jbin in range(2,2+channel_bins[channel]):
            f = plt.figure(figsize =(12, 6)); ax=f.add_subplot(1,1,1)
            if x_label != None: ax.set_xlabel(x_label, fontsize=20)
            if y_label != None: ax.set_ylabel(y_label, fontsize=20)
            ax.set_xlim(x_min-25, x_max+25)
            ax.set_ylim(y_min-25, y_max+25)
            bin_number = jbin - 1
            frac_systs = frac_systs[frac_systs[:, jbin] != 0]   # Remove any points with zero amplitude
            
            z = scipy.interpolate.griddata(frac_systs[:,:2],frac_systs[:,jbin],(x,y))
            #print(frac_systs)
            
            vmin,vmax = 0, v_max
            #z = np.log(z)
            ax.scatter(frac_systs[:,0],frac_systs[:,1], c = frac_systs[:,jbin], edgecolors='w', vmin = vmin, vmax = vmax)
            im = ax.contourf(x,y,z, levels = np.linspace(vmin,vmax,100))
            cb = plt.colorbar(im, ax=ax)
            cb.set_label(label='$\oplus$ (histosys, normsys, staterr)', fontsize=18)
            if channel_bins[channel] < 2: ax.set_title(channel_names[channel], fontsize=20)
            else: ax.set_title(f'{channel_names[channel]} (Bin {bin_number})', fontsize=20)
            
            outliers_chan = np.asarray(
                [[float(parse(sig_name_template, o[0])[0]), float(parse(sig_name_template, o[0])[1])] + [o[-2]] for o in outliers
                 if o[1] == channel and o[2] == jbin-2]
            )
            
            #print(outliers_chan)
            
            if outliers_chan.shape[0]:
                ax.scatter(outliers_chan[:,0],outliers_chan[:,1], c = outliers_chan[:,2], vmin = 0, vmax = 20, cmap = 'cool')
            for o in outliers_chan:
                ax.text(o[0]+5,o[1]+5,'{:.2f}'.format(o[2]), c='r')
        
            plt.tight_layout()
            plt.savefig(f'{channel_names[channel]}_bin{bin_number}.png')
            plt.close()

if __name__ == "__main__":
    outlier_plot()
