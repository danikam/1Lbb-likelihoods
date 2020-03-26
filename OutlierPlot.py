import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from parse import *
import click
import pyhf

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

def outlier_plot(signal_template, v_max):

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
        x,y = np.mgrid[min(frac_systs[:,0]):max(frac_systs[:,0]):100j,min(frac_systs[:,1]):max(frac_systs[:,1]):100j]
        for jbin in range(2,2+channel_bins[channel]):
            f = plt.figure(figsize =(12, 6)); ax=f.add_subplot(1,1,1)
            bin_number = jbin - 1
            frac_systs = frac_systs[frac_systs[:, jbin] != 0]   # Remove any points with zero amplitude
            
            z = scipy.interpolate.griddata(frac_systs[:,:2],frac_systs[:,jbin],(x,y))
            #print(frac_systs)
            
            vmin,vmax = 0, v_max
            #z = np.log(z)
            ax.scatter(frac_systs[:,0],frac_systs[:,1], c = frac_systs[:,jbin], edgecolors='w', vmin = vmin, vmax = vmax)
            im = ax.contourf(x,y,z, levels = np.linspace(vmin,vmax,100))
            #f.colorbar(im, cax=ax[1])
            f.colorbar(im, ax=ax, label='$\oplus$ (histosys, normsys, staterr)')
            if channel_bins[channel] < 2: ax.set_title(channel_names[channel])
            else: ax.set_title(f'{channel_names[channel]} (Bin {bin_number})')
            
            outliers = np.asarray(
                [list(map(lambda x: float(x.replace('p','.')),o[0].split('hbb')[-1].split('_')[1:3])) + [o[-2]] for o in outliers
                 if o[1] == channel and o[2] == jbin-2]
            )
            #print(outliers)
            
            if outliers.shape[0]:
                ax.scatter(outliers[:,0],outliers[:,1], c = outliers[:,2], vmin = 0, vmax = 20, cmap = 'cool')
            for o in outliers:
                ax.text(o[0]+5,o[1]+5,'{:.2f}'.format(o[2]), c='r')
        
            plt.savefig(f'{channel_names[channel]}_bin{bin_number}.png')
            plt.close()

if __name__ == "__main__":
    outlier_plot()
