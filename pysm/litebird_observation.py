import numpy as np
import healpy as hp
import pysm
import os
from pysm.litebird_models import models

def make_litebird_instrument(nside=512, duration=1., seed=None, out_dir='./', prefix='test'):
    freqs = np.array([40, 50, 60, 68, 78, 89, 100, 119, 140, 166, 195, 235, 280, 337, 402])
    beams = np.array([ 60, 56, 48, 43, 39, 35, 29, 25, 23, 21, 20, 19, 24, 20, 17])
    sens_P = np.array([39.76, 25.76, 20.69, 12.72, 10.39, 8.95, 6.43, 4.3, 4.43, 4.86, 5.44, 9.72, 12.91, 19.07, 43.53])
    if seed is None:
        seed = np.random.randint(1)
    litebird_inst = {
    'nside' = nside,
    'frequencies' : freqs,
    'use_smoothing' : True,
    'beams' : beams,
    'add_noise' : True,
    'sens_I' : sens_P/np.sqrt(2.)*np.sqrt(duration),
    'sens_P' : sens_P*np.sqrt(duration),
    'noise_seed' : seed,
    'use_bandpass' : False,
    'output_units' : 'uK_CMB',
    'output_directory' : out_dir,
    'output_prefix' : prefix
    }

def make_litebird_instrument_noisless(nside=512, out_dir='./', prefix='test_noiseless'):
    freqs = np.array([40, 50, 60, 68, 78, 89, 100, 119, 140, 166, 195, 235, 280, 337, 402])
    beams = np.array([ 60, 56, 48, 43, 39, 35, 29, 25, 23, 21, 20, 19, 24, 20, 17])
    litebird_inst = {
    'nside' = nside,
    'frequencies' : freqs,
    'use_smoothing' : True,
    'beams' : beams,
    'add_noise' : False,
    'use_bandpass' : False,
    'output_units' : 'uK_CMB',
    'output_directory' : out_dir,
    'output_prefix' : prefix
    }

def make_cmb(cl, out_dir='./', prefix='cmb_input'):
    fname = prefix+'.fits'
    path = os.path.join(out_dir, fname)
    cmb = hp.synfast(cl, new=True)
    hp.write_map(path, cmb, overwrite=True)


def make_skys(nside=512, dust_model='LBd0', synch_model='LBs0', cmb_model='LBc0', cmb_in=None):
    sky_fg = pysm.Sky({"dust":models(dust_model, nside),
        "synchrotron":models(synch_model, nside)})
    sky_cmb = pysm.Sky({"cmb":models(cmb_model, nside, cmb_in)})
    return sky_fg, sky_cmb

def observe_sky(sky, instrument):
    instr = pysm.Instrument(instrument)
    sky_obs = instr.observe(sky_cmb)
    return sky_obs
