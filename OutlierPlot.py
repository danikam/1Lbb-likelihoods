import json
import glob
import numpy as np
import matplotlib
import copy
import matplotlib.pyplot as plt

patches = [json.load(open(x)) for x in glob.glob('patch*.json')]
grouped_patches = {x[0]['value']['name']:{p['path']:p for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}

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

process_patch(grouped_patches['C1N2_Wh_hbb_1000_0']['/channels/2/samples/9'])

data = {x[0]['value']['name']:{p['path']:process_patch(p) for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}
#print(data)

channel_names = {
    '/channels/0/samples/9': 'Region_0',
    '/channels/1/samples/9': 'Region_1',
    '/channels/2/samples/9': 'Region_2',
    '/channels/3/samples/9': 'Region_3',
    '/channels/4/samples/9': 'Region_4',
    '/channels/5/samples/7': 'Region_5',
    '/channels/6/samples/8': 'Region_6',
    '/channels/7/samples/9': 'Region_7', 
}

outliers = []
for k,v in data.items():
    for kk,vv in v.items():
        for b,(r,n) in enumerate(zip(*vv)):
            if r > 1.0:
                outliers.append((k,kk,b,r,n))

#print("Outliers (> 1.0):")
#for o in list(reversed(sorted(outliers, key = lambda x: x[-1]))):
#  print('\t',o[-1],o[-2],o[0],channel_names[o[1]],o[2])
"""
missing_signal = []
# missing signal in signal region
print("Missing signal in signal region:")
for k, v in data.items():
    if not '/channels/2/samples/5' in v or not '/channels/2/samples/5' in v:
        missing_signal.append(k)
        print('\t',k)
"""

import scipy.interpolate
data = {x[0]['value']['name']:{p['path']:process_patch(p)[0] for p in x if p['op'] == 'add'} for x in patches if 'value' in x[0]}

f,axarr = plt.subplots(8,1)
for ichan, channel in enumerate(['/channels/0/samples/9']):#,'/channels/1/samples/9', '/channels/2/samples/9', '/channels/3/samples/9', '/channels/4/samples/9', '/channels/5/samples/7', '/channels/6/samples/8', '/channels/7/samples/9']):
    print(ichan)
    what = np.asarray(
        [[float(x.replace('p','.')) for x in k.split('hbb')[-1].split('_')[1:3]] + v.get(channel, np.array([0.0])).tolist() for k,v in data.items()
        ]
    )
    x,y = np.mgrid[100:500:100j,0:300:100j]
    for jbin in range(2,3):
        print(jbin)
        print(what[:,:2])
        print(what[:,jbin])
        z = scipy.interpolate.griddata(what[:,:2],what[:,jbin],(x,y))
        """
        vmin,vmax = 0, 1#np.log(0.001),np.log(20)
        #z = np.log(z)
        axarr[ichan].scatter(what[:,0],what[:,1], c = what[:,jbin], edgecolors='w', vmin = vmin, vmax = vmax)
        im = axarr[ichan].contourf(x,y,z, levels = np.linspace(vmin,vmax,100))
        #f.colorbar(im, cax=axarr[ichan][1])
        f.colorbar(im, ax=axarr[ichan], label='$\oplus$ (histosys, normsys, staterr)')
        axarr[ichan].set_title(channel_names[channel])
        
        out_here = np.asarray(
            [list(map(lambda x: float(x.replace('p','.')),o[0].split('WZ')[-1].split('_')[1:3])) + [o[-2]] for o in outliers
             if o[1] == channel and o[2] == jbin-2]
        )
#         print(out_here)
        
#         print(out_here.shape)
        if out_here.shape[0]:
            axarr[ichan].scatter(out_here[:,0],out_here[:,1], c = out_here[:,2], vmin = 0, vmax = 20, cmap = 'cool')
        for o in out_here:
            axarr[ichan].text(o[0]+5,o[1]+5,'{:.2f}'.format(o[2]), c='r')
        """

f.set_size_inches(15,15)
f.savefig('plot.png')
plt.close()
