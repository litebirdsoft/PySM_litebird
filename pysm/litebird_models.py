from __future__ import absolute_import
from .common import read_map, loadtxt
import numpy as np
from healpy import nside2npix
import os

LB_data_dir = os.path.join(os.path.dirname(__file__), 'lb_template')
LB_template = lambda x: os.path.join(LB_data_dir, x)

def models(key, nside, pixel_indices=None, mpi_comm=None, cmb_in=None):
    if key=='LBc0':
        model = eval(key)(nside, pixel_indices=pixel_indices, mpi_comm=mpi_comm, cmb_in=cmb_in)
    else:
        model = eval(key)(nside, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    for m in model:
        m['pixel_indices'] = pixel_indices # include pixel indices in the model dictionary
        m['nside'] = nside
    return model


def LBd0(nside, pixel_indices=None, mpi_comm=None):
    A_I = read_map(LB_template('dust_T_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'modified_black_body',
        'nu_0_I': 545.,
        'nu_0_P': 353.,
        'A_I': A_I,
        'A_Q': read_map(LB_template('dust_Q_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(LB_template('dust_U_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': np.ones(len(A_I)) * 1.53,
        'temp': np.ones(len(A_I)) * 19.6,
        'add_decorrelation': False,
    }]

def LBd1(nside, pixel_indices=None, mpi_comm=None):
    A_I = read_map(LB_template('dust_T_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'modified_black_body',
        'nu_0_I': 545.,
        'nu_0_P': 353.,
        'A_I': A_I,
        'A_Q': read_map(LB_template('dust_Q_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(LB_template('dust_U_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': read_map(LB_template('beta_dust_ns512_1deg.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'temp': read_map(LB_template('temperature_dust_ns512_1deg.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'add_decorrelation': False,
    }]

def LBs0(nside, pixel_indices=None, mpi_comm=None):
    A_I = read_map(LB_template('synch_T_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'power_law',
        'nu_0_I': 0.408,
        'nu_0_P': 23.,
        'A_I': A_I,
        'A_Q': read_map(LB_template('synch_Q_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(LB_template('synch_U_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': np.ones(len(A_I)) * -3.1,
    }]

def LBs1(nside, pixel_indices=None, mpi_comm=None):
    A_I = read_map(LB_template('synch_T_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm)
    return [{
        'model': 'power_law',
        'nu_0_I': 0.408,
        'nu_0_P': 23.,
        'A_I': A_I,
        'A_Q': read_map(LB_template('synch_Q_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(LB_template('synch_U_ns512.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': read_map(LB_template('beta_synch_ns512_1deg.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
    }]



def LBf0(nside, pixel_indices=None, mpi_comm=None):
    return [{
        'model': 'power_law',
        'nu_0_I': 30.,
        'A_I': read_map(LB_template('freefree_T.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'spectral_index': -2.14,
    }]

def LBa0(nside, pixel_indices=None, mpi_comm=None):
    return [{
        'model': 'spdust',
        'nu_0_I': 22.8,
        'nu_0_P': 22.8,
        'A_I': read_map(LB_template('ame_T.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'nu_peak_0': 30.,
        'emissivity': loadtxt(LB_template('ame_emissivity.txt'), mpi_comm=mpi_comm, unpack=True),
        'nu_peak': 18.95,
    }, {
        'model': 'spdust',
        'nu_0_I': 41.0,
        'nu_0_P': 41.0,
        'A_I': read_map(LB_template('ame2_T.fits'), nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'nu_peak_0': 30.,
        'emissivity': loadtxt(LB_template('ame_emissivity.txt'), mpi_comm=mpi_comm, unpack=True),
        'nu_peak': 33.35,
    }]


def LBc0(nside, pixel_indices=None, mpi_comm=None, cmb_in=None):
    if cmb_in:
        cmb_map = cmb_in
    else:
        cmb_map = LB_template('cmb_lens_Planck2018_r0.fits')
    return [{
        'model': 'pre_computed',
        'A_I': read_map(cmb_map, nside, field=0, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_Q': read_map(cmb_map, nside, field=1, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'A_U': read_map(cmb_map, nside, field=2, pixel_indices=pixel_indices, mpi_comm=mpi_comm),
        'nside': nside
    }]
