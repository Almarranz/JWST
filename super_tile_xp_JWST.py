# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:46:07 2026

@author: amartinez
"""

import sys
sys.path.append("/Users/amartinez/Desktop/pythons_imports/")
from sys import exit as stop
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.colors import LogNorm
from scipy import stats
import compare_lists as complists
from astropy.table import vstack
from astropy.table import unique
import time
import IPython
from ds9_region import region
from itertools import combinations
from sklearn.neighbors import KDTree
# %%
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')    

def tile_filename(base_path, prefix, fields, ext=".txt"):
    fields_str = "_".join(f"f{f}" for f in fields)
    return f"{base_path}{prefix}_{fields_str}{ext}"

check_plots = 'yes'
# optimze_tiles = 'yes'
optimze_tiles = 'no'
# survey = 1
# survey = 2

  
    
# fields= [1,2,3,4]
# fields= ['a1','a2','a3','a4', 'b1','b2','b3','b4']
detectors = ['a1','a2','a3','a4', 'b1','b2','b3','b4']
# detectors = ['a1','a2','a3','a4']
pointings = np.arange(1,7)



band = 'H'
pixel_scale = 0.031 # hijos de puta!!
uncer_in_mas = 1 # express the positions uncertainties in mas
# uncer_in_mas = 0 #Otherhise it woulb be in deg


poserr = 1 #mas position uncertaintiy in arcsec
# poserr = 10 #mas position uncertaintiy in arcsec
if (band == 'J'):
    mag_bright = 13
    mag_faint = 19
if (band == 'H'):
    mag_bright = 10
    mag_faint = 22
if (band == 'Ks'):
    mag_bright = 13
    mag_faint = 16
tolerance = 23 * u.mas  # 
# tolerance = 1e6 * u.mas  #

tg_point = SkyCoord(l = 359.9443, b = -0.0462, unit = 'degree', frame = 'galactic').fk5


fig, ax = plt.subplots(1,1)
ax.invert_xaxis()
ax.axis('equal')
ax.set_xlabel('xp')
ax.set_ylabel('yp')



# subfix = tile_filename('', '', fields)
# x_off = np.loadtxt(inpath  + f'x_offsets{subfix}').reshape(len(fields),4)
# y_off = np.loadtxt(inpath  + f'y_offsets{subfix}').reshape(len(fields),4)
# m_off = np.loadtxt(inpath  + f'm_offsets{subfix}').reshape(len(fields),4)



d_chip = {}
ids = 0
total = 0

tile_a = Table(
    names=('x','y','f','sx','sy','sf','RA','Dec','m','sm','id','field'),
    dtype=(float, float, float, float, float, float,
           float, float, float, float, str, str)
)

folder = f'/Users/amartinez/Desktop/Projects/JWST/data/rainer_ls/'

for cf, de in enumerate(detectors):
   
   
   
    # for po in range(1,3,1):
    for po in pointings:
        print(de,po)
        
        
        
        cat = Table.read(folder + f'jw10678098001_02101_0000{po}_nrc{de}_destreak_o098_crf_stars_cal.txt', format = 'ascii')
        # cat = Table.read(folder + f'0{c}_nrc{f}_stars_calib.txt', format = 'ascii')
        cat['id'] = [f"{i}d{de}p{po}" for i in range(len(cat))]
        cat['field'] = [f"d{de}p{po}" for i in range(len(cat))]
        
        
        if uncer_in_mas:
            cat['sRA'] = cat['sx']*0.031*1000
            cat['sDec'] =  cat['sy']*0.031*1000
        else:
            cat['sRA'] = cat['sx']*0.031/3600
            cat['sDec'] =  cat['sy']*0.031/3600
        # #select stars in calibration magnitude range
        
        
# =============================================================================
#         Quality cuts
# =============================================================================
# =============================================================================
#         mag_m = np.nonzero((cat['m'] > mag_bright) & (cat['m'] < mag_faint) 
#                            & (cat['sx'] < poserr)& (cat['sy'] < poserr))[0]
#         cat = cat[mag_m]
#         total += len(cat)
#         
# =============================================================================
        coord1 = SkyCoord(ra = cat['RA'], dec = cat['Dec'], unit='deg')
        
        
        offset_frame = tg_point.skyoffset_frame()
        
        xyp = coord1.transform_to(offset_frame)

        
        cat['xp'] = xyp.lon 
        cat['yp'] = xyp.lat 
        if uncer_in_mas:
            cat['sxp'] = cat['sx']*0.031*1000
            cat['syp'] = cat['sx']*0.031*1000
        else:
            cat['sxp'] = cat['sx']*0.031/3600
            cat['syp'] = cat['sx']*0.031/3600
        
        
        # if optimze_tiles == 'yes':
        #     print(30*'<' + '\nApplied optimization\n' + 30*'>')
        #     cat['xp'] =  cat['xp'] +  x_off[cf, c-1]
        #     cat['yp'] =  cat['yp'] + y_off[cf, c-1]
        #     cat['m'] =  cat['m'] + m_off[cf, c-1]
        
        xypC = SkyCoord(lon = cat['xp'], lat = cat['yp'], frame = offset_frame)
        cat['xypC'] = xypC
        
        # print(30*'*'+ f'\nField{cf}, Chip{c}: x_off = {x_off[cf, c-1]}\n' + 30*'_')
        # stop(161)
        
        # ax.scatter(lb.l.wrap_at('180d'), cat['b'], s = 0.1,alpha = 0.3, label = f'field {f}, c {c}' )
        ax.scatter(cat['xp'], cat['yp'], s = 0.1,alpha = 0.3, label = f'dctr {de}, pntg {po}' )
        
        
        lg = ax.legend() 
        for handle in lg.legend_handles:
            handle.set_sizes([10])
            handle.set_alpha(1)
            
        # d_chip[f'f{f}c{c}'] = cat['xp','yp','sx','sy','fm']
        d_chip[f'd{de}p{po}'] = cat
        
        tile_a = vstack([tile_a, cat])
        
        # break
plt.show()





# stop(190)
# %%


chip_keys = list(d_chip.keys())
# Example: ['f2c1', 'f2c2', 'f2c3', 'f2c4', 'f19c1', ...]
key_to_index = {k: i for i, k in enumerate(chip_keys)}
# %%

from collections import defaultdict

# Dictionary: {star_id: list of (subfield_key, index_in_subfield)}
star_to_subfields = defaultdict(list)

for k in chip_keys:
    for i, row in enumerate(d_chip[k]):
        star_to_subfields[row['id']].append((k, i))
# %%
multi_subfield_stars = {star_id: subfields for star_id, subfields in star_to_subfields.items() if len(subfields) > 1}
# %%


tile_m1 = Table(
    names=('x','y','xp','yp','sxp','syp','f','sx','sy','sf','RA','Dec','sRA','sDec','m','sm','id','field'),
   
dtype=(float, float, float, float, float, float,
           float, float, float, float,   # up to 'sf'
           float, float, float, float,   # ra, dec, sra, sdec
           float, float,                 # m, smag
           str, str)                     # id, field
)



tile_m2 = Table(
    names=('x','y','xp','yp','sxp','syp','f','sx','sy','sf','RA','Dec','sRA','sDec','m','sm','id','field'),
   
dtype=(float, float, float, float, float, float,
           float, float, float, float,   # up to 'sf'
           float, float, float, float,   # ra, dec, sra, sdec
           float, float,                 # m, smag
           str, str)                     # id, field
)


tile_m12 = Table(
    names=('x','y','xp','yp','sxp','syp','f','sx','sy','sf','RA','Dec','sRA','sDec','m','sm','id','field'),
   
dtype=(float, float, float, float, float, float,
           float, float, float, float,   # up to 'sf'
           float, float, float, float,   # ra, dec, sra, sdec
           float, float,                 # m, smag
           str, str)                     # id, field
)


dic_m = {}

for k1, k2 in combinations(chip_keys, 2):    # all unordered chip pairs
    
    
    # print(f'{k1}, {len(d_chip[k1])}')
    # print(f'{k2}, {len(d_chip[k2])}')
    i = key_to_index[k1]
    j = key_to_index[k2]
    
    # Extract xy lists for each chip
    l1 = np.vstack([d_chip[k1]['xp'], d_chip[k1]['yp']]).T
    l2 = np.vstack([d_chip[k2]['xp'], d_chip[k2]['yp']]).T

    # comp = complists.compare_lists(l1, l2, tolerance.to(u.arcsec).value)
    comp = complists.compare_lists(l1, l2, tolerance.to(u.degree).value)
   
    if len(comp['ind_1']) < 1:
        print('No overlapping between',k1, k2, len(comp['ind_1']))
        continue
    print(k1,k2, len(comp['ind_1']))
    m1 = d_chip[k1][comp['ind_1']]
    m2 = d_chip[k2][comp['ind_2']]

   
    
    
    match_s = np.logical_not(np.isin(tile_a['id'], m1['id']) | (np.isin(tile_a['id'], m2['id'])))
    tile_a = tile_a[match_s]
    
   
    tile_m1 = vstack([tile_m1, m1])
    tile_m2 = vstack([tile_m2, m2])
    
    tile_um = vstack([tile_m1, tile_m2])




# %%
ta = time.time()

fig, ax = plt.subplots(1,1)
ax.set_title(f'Non matched {len(tile_a)} F{pointings}\nMetodo group-crossmatching {len(tile_a)}')
ax.scatter(tile_a['xp'], tile_a['yp'], color = 'lime',s = 1, alpha = 0.03)
# h2 = ax.hexbin(tile_a['xp'], tile_a['yp'], norm = LogNorm())
# plt.colorbar(h2[3], ax = ax)
ax.axis('equal')   

fig, (ax,ax2) = plt.subplots(1,2)

ax.hexbin(tile_a['m'], np.sqrt(tile_a['sxp']**2 + tile_a['syp']**2), norm = LogNorm())
ax2.hexbin(tile_a['m'], tile_a['sm'], norm = LogNorm())
ax.set_xlabel('[m]')
ax2.set_xlabel('[m]')
ax2.set_ylabel('$\sigma_{m}$[mag]')

if uncer_in_mas:
    ax.set_ylabel('$\sigma_{Pos}$[mas]]')
else:
    ax.set_ylabel('$\sigma_{Pos}$[deg]')
fig.tight_layout()
    



tile_um = unique(tile_um, keys = 'id')
tg = tile_um.group_by('field')
grupos = tg.groups.keys

# cols = ['x', 'y', 'f', 'sx', 'sy', 'sf', 'ra', 'dec', 'm', 'smag']
# cols = ['x', 'y', 'f', 'sx', 'sy', 'sf', 'ra', 'dec', 'm',]
cols = ['xp', 'yp', 'f', 'RA', 'Dec', 'm']
rows_to_add = []
for g in grupos:
    tile_gi = tile_um[tile_um['field'] == g[0]]
    tile_gj = tile_um[np.logical_not(tile_um['field'] == g[0])]
    print(f'Group {g[0]} has {len(tile_gi)} stars ')
    print(f'Rest of groups {len(tile_gj)} stars ')
    total = len(tile_gi) + len(tile_gj)
    print(f'Total: {total}')
    
    if len(tile_gi) > 0 and len(tile_gj):
    # if len(tile_gi) > 0 and len(tile_gj) >0:
    
        # for s in range(len(tile_gi)):
        count = 0
        

        # Precompute once per group
        p_gi = np.column_stack([tile_gi['xp'], tile_gi['yp']])
        p_gj = np.column_stack([tile_gj['xp'], tile_gj['yp']])
        
        tree = KDTree(p_gj)
        
        # Query all stars at once
        ind, dist = tree.query_radius(p_gi, r = tolerance.to(u.degree).value, return_distance=True)
        
        to_remove_ids = set()
        # to_remove_ids = []
        for s, matches in enumerate(ind):
            if len(matches) == 0:
                continue
            
            
            
            vals_RA = np.concatenate(([tile_gi['RA'][s]], tile_gj['RA'][matches]))
            tile_gi['sRA'][s] = np.std(vals_RA) / np.sqrt(len(vals_RA)) #* 3600
            
            vals_Dec = np.concatenate(([tile_gi['Dec'][s]], tile_gj['Dec'][matches]))
            tile_gi['sDec'][s] = np.std(vals_Dec) / np.sqrt(len(vals_Dec))# * 3600
            
            vals_xp = np.concatenate(([tile_gi['xp'][s]], tile_gj['xp'][matches]))
            tile_gi['sxp'][s] = np.std(vals_xp) / np.sqrt(len(vals_xp)) #* 3600
            
            vals_yp = np.concatenate(([tile_gi['yp'][s]], tile_gj['yp'][matches]))
            tile_gi['syp'][s] = np.std(vals_yp) / np.sqrt(len(vals_yp)) #* 3600
            
            vals_m = np.concatenate(([tile_gi['m'][s]], tile_gj['m'][matches]))
            tile_gi['sm'][s] = np.std(vals_m) / np.sqrt(len(vals_m))
            
            vals_f = np.concatenate(([tile_gi['f'][s]], tile_gj['f'][matches]))
            tile_gi['sf'][s] = np.std(vals_f) / np.sqrt(len(vals_f))
            
            
            # tile_gi['sRA'][s] = np.std(np.concatenate(([tile_gi['RA'][s]],  tile_gj['RA'][matches])))
            # tile_gi['sDec'][s] = np.std(np.concatenate(([tile_gi['Dec'][s]],  tile_gj['Dec'][matches])))
            # tile_gi['sxp'][s] = np.std(np.concatenate(([tile_gi['xp'][s]],  tile_gj['xp'][matches])))
            # tile_gi['syp'][s] = np.std(np.concatenate(([tile_gi['yp'][s]],  tile_gj['yp'][matches])))
            # tile_gi['sm'][s] = np.std(np.concatenate(([tile_gi['m'][s]],  tile_gj['m'][matches])))
            # tile_gi['sf'][s] = np.std(np.concatenate(([tile_gi['f'][s]],  tile_gj['f'][matches])))
            for col in cols:
                tile_gi[col][s] = np.mean(np.concatenate(([tile_gi[col][s]], tile_gj[col][matches])))
        
            rows_to_add.append(tile_gi[s])
            
            
            # to_remove_ids.append(tile_gi['id'][s])
            # to_remove_ids.append(tile_gj['id'][matches])
            to_remove_ids.add(tile_gi['id'][s])
            to_remove_ids.update(tile_gj['id'][matches])

            
            # rem_m = tile_um['id'] == tile_gi[s]['id']
            # tile_um.remove_rows(np.where(rem_m))    
            # rem_m2 = np.isin(tile_um['id'], tile_gj[matches]['id'])
            # tile_um.remove_rows(np.where(rem_m2))
        tile_um = tile_um[~np.isin(tile_um['id'], list(to_remove_ids))]
 
tile_a = vstack([tile_a, Table(rows_to_add)])  
tb = time.time()

print(f'loop took {(tb-ta)/60: .2f} min')

# filename = tile_filename(tiles_f, f"tile_gns{survey}_pro", fields)
# tile_a.write(filename, format="ascii", overwrite=True) 


# tile_a = unique(tile_a, keys = 'id')
%%


if check_plots == 'yes':
    fig, ax = plt.subplots(1,1)
    ax.set_title(f'Pntgn{pointings} dect{detectors}\nMetodo group-crossmatching {len(tile_a)}')
    # ax.scatter(tile['xp'], tile['yp'], color = 'lime',s = 0.5, alpha = 0.01)
    h2 = ax.hexbin(tile_a['xp'], tile_a['yp'], norm = LogNorm() )
    # h2 = ax.hexbin(tile_a['xp'], tile_a['yp'] )
    # h2 = ax.hist2d(tile['xp'], tile['yp'],bins = (60,60))
    # plt.colorbar(h2, ax = ax)
    ax.axis('equal')   
# idxc, group_md, d2d, d3d = ap_coor.search_around_sky(
#  
   

# %%


    # IPython.get_ipython().run_line_magic('matplotlib', 'auto')
    # IPython.get_ipython().run_line_magic('matplotlib', 'inline')    
    # tile_a = unique(tile_a, keys = 'id')
    fig, ax = plt.subplots(1,1)
    # ax.set_title(f'F{fields}')
    ax.set_title(f'Pntgn{pointings} dect{detectors}\nMetodo group-crossmatching {len(tile_a)}')
    ax.scatter(tile_a['xp'], tile_a['yp'], color = 'k',s = 0.001, alpha = 1)
    # h2 = ax.hexbin(tile_a['sl'], tile_a['sb'], norm = LogNorm() )
    # h2 = ax.hist2d(tile['xp'], tile['yp'],bins = (60,60))
    # plt.colorbar(h2, ax = ax)
    ax.axis('equal')   
    
    # ax.scatter(tile_gi['xp'], tile_gi['yp'], edgecolor = 'red', lw = 1,s = 10, alpha = 1, facecolor = 'none' ) 
    # ax.scatter(tile_gj['xp'], tile_gj['yp'], edgecolor = 'blue', lw = 1,s = 10, alpha = 1, facecolor = 'none' ) 
    
    # region(tile_gj,'RA', 'Dec', name = 'extra_detection',
    #        save_in = '/Users/amartinez/Desktop/for_people/for_Lucia/pruebas/',
    #        )
# %%

# %%


    fig, ax = plt.subplots(1,1)
    ax.set_title(f'$\sigma$ position.Crs-Math dis = {tolerance.to(u.arcsec):.2f}')
    ax.axis('equal')
    ax.set_xlabel('xp')
    ax.set_ylabel('yp')
    num_bins = 60
    s_xy = np.sqrt(tile_a['sxp']**2 + tile_a['syp']**2)
    # s_xy = np.sqrt(tile_a['sRA']**2 + tile_a['sDec']**2)
    
    # s_xy = np.sqrt(tile_a['sm']**2 )
    statistic, x_edges, y_edges, binnumber = stats.binned_statistic_2d(tile_a['xp'], tile_a['yp'], s_xy, statistic='mean', bins=(num_bins,int(num_bins/2)))
    # statistic, x_edges, y_edges, binnumber = stats.binned_statistic_2d(tile_a['RA'], tile_a['Dec'], s_xy, statistic='mean', bins=(num_bins,int(num_bins/2)))
    X, Y = np.meshgrid(x_edges, y_edges)
    
    # c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r',norm = LogNorm() )
    c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r')
    # c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r')
    if uncer_in_mas:
        cb1 = fig.colorbar(c, ax=ax, label=r'$\sqrt{sxp^2 +syp^2}$ [mas]', shrink = 0.8, )
    else:
        cb1 = fig.colorbar(c, ax=ax, label=r'$\sqrt{sxp^2 +syp^2}$ [deg]', shrink = 0.8, )
    # ax.colorbar(c, ax=ax, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.5)
# %%
    

    fig, ax = plt.subplots(1,1)
    ax.set_title(f'Pntgn{pointings} dect{detectors}')
    ax.axis('equal')
    ax.set_xlabel('xp')
    ax.set_ylabel('yp')
    num_bins = 60
    s_m = tile_a['sm']
    
    # s_xy = np.sqrt(tile_a['sm']**2 )
    statistic, x_edges, y_edges, binnumber = stats.binned_statistic_2d(tile_a['xp'], tile_a['yp'], s_m, statistic='mean', bins=(num_bins,int(num_bins/2)))
    X, Y = np.meshgrid(x_edges, y_edges)
    
    # c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r',norm = LogNorm() )
    c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r' )
    # c = ax.pcolormesh(X, Y, statistic.T, cmap='Spectral_r')
    cb1 = fig.colorbar(c, ax=ax, label=r'$\sigma$[m] mag', shrink = 0.8,)
    # ax.colorbar(c, ax=ax, label=r'$\sqrt{sl^2 +sb^2}$ [arcsec]', shrink = 0.5)

# %%
    

    fig, (ax,ax2) = plt.subplots(1,2)
    
    ax.hexbin(tile_a['m'], s_xy, norm = LogNorm())
    ax2.hexbin(tile_a['m'], s_m, norm = LogNorm())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    