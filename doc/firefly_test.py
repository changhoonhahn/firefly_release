import sys
import os.path
import numpy as np
import astropy.cosmology as co
import GalaxySpectrumFIREFLY as gs
import StellarPopulationModel as spm
import astropy.units as u 
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
# --- 
from fomospec import util as UT 
from fomospec import spectra as Spec

import matplotlib.pyplot as plt 

def fit_test(galid): 
    cosmo = co.Planck13
    #spec = gs.GalaxySpectrumFIREFLY(None, milky_way_reddening=True)

    # read in source spectra
    f_s = ''.join([UT.dat_dir(), 'Lgal/templates/', 'gal_spectrum_', str(galid), '_BGS_template_BC03_Stelib.fits'])
    f_inspec = fits.open(f_s)
    hdr = f_inspec[0].header
    specin = f_inspec[1].data
    meta = {}
    for k in hdr.keys(): 
        meta[k] = hdr[k]
    zred = meta['REDSHIFT']

    f_input = ''.join([UT.dat_dir(), 'Lgal/gal_inputs/', 
        'gal_input_', str(galid), '_BGS_template_FSPS_uvmiles.csv']) 
    gal_input = Table.read(f_input, delimiter=' ')

    spec_in = {}
    spec_in['wave'] = specin['wave']
    #flux = specin['flux_dust_nonoise'] * u.Watt / u.AA / u.m**2 #from W/A/m2 to 10e-17 erg/s/cm2/A
    #flux = flux.to(1e-17*u.erg/u.s/u.cm**2/u.AA) 
    spec_in['flux_dust_nonoise'] = specin['flux_dust_nonoise'] * 1e-4 * 1e7 *1e17 
    spec_in['flux_nodust_nonoise'] = specin['flux_nodust_nonoise'] * 1e-4 * 1e7 *1e17

    gspec = Spec.GSfirefly()
    gspec.generic(spec_in['wave'], spec_in['flux_nodust_nonoise'],  redshift=zred)
    gspec.path_to_spectrum = UT.dat_dir()
    
    outputfile = 'test.fits'
    bestfit = spm.StellarPopulationModel(gspec, outputfile, cosmo, 
            models = 'm11', 
            model_libs = ['MILES'], 
            imfs = ['cha'], 
            hpf_mode = 'on', 
            age_limits = [0., 15.], 
            downgrade_models = False,
            data_wave_medium = 'vacuum',
            Z_limits = [0.001, 4.], 
            wave_limits = [3350., 9000.],
            use_downgraded_models=False, 
            write_results=True)
    bestfit.fit_models_to_data()
    keys = ['age_lightW','age_massW','metallicity_lightW','metallicity_massW','stellar_mass']#, 'HIERARCH stellar_mass']
    units = ['Age_unit','Age_unit','Metallicity_unit','Metallicity_unit','Mass_unit']#, 'Mass_unit']
    
    #for k, un in zip(keys, units): 
    #    print('%s, %s = %f' % (k, bestfit.thdulist[1].header[un], bestfit.thdulist[1].header[k]))
    
    print('redshift = %f' % zred)
    print('input log M*= %f' % np.log10(np.sum(gal_input['sfh_disk']) + np.sum(gal_input['sfh_bulge'])))
    print(np.log10(bestfit.averages['stellar_mass']))
    print(bestfit.thdulist[1].header['stellar_mass'])
    print(bestfit.thdulist[1].header['total_mass'])

    plt.plot(bestfit.thdulist[1].data['wavelength'], bestfit.thdulist[1].data['original_data']) 
    plt.plot(bestfit.thdulist[1].data['wavelength'], bestfit.thdulist[1].data['firefly_model'], ls='--') 
    plt.ylim([0., 2.]) 
    plt.savefig('test'+str(galid)+'.png') 
    plt.close()
    return None 


if __name__=='__main__': 
    for gid in [2322, 3086, 2675]: 
        fit_test(gid)
