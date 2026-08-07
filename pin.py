"""
Suzaku XIS extractor
"""
from datetime import datetime
import glob
import subprocess
from .xselect import Xselect
from .spec_util import *
import os
import re
from astropy.io import fits
import urllib.request

xis_names = {'xi0': 'XIS0', 'xi1': 'XIS1', 'xi2': 'XIS2', 'xi3': 'XIS3'}

class PinExtractor(object):
    #
    # class to extract data products from Suzaku PIN observations
    #

    def __init__(self, obsdir, evl_dir='reproc', run_reduction=False, suffix=None):
        self.obsdir = obsdir
        self.hxddir = obsdir + '/hxd'
        self.pindir = self.hxddir + '/pin'

        self.stem = 'ae%s' % obsdir

        self.suffix = suffix

        self.evlsdir = self.pindir + '/reproc'
        self.specdir = self.pindir + '/spectra'
        self.lcdir = self.pindir + '/lightcurves'

        self.evls = self.find_evls()
        if self.evls['evl'] is None or self.evls['pse'] is None:
            if run_reduction:
                self.reprocess()
                self.evls = self.find_evls()
            else:
                raise ValueError('Processed event lists not found in selected directory. Do you need to run aepipline?')

        with fits.open(self.evls['evl']) as f:
            self.date_obs = datetime.strptime(f[0].header['DATE-OBS'].split('.')[0], "%Y-%m-%dT%H:%M:%S")

    def reprocess(self):
        tmp_reproc_dir = 'tmp_reproc_%s_pin' % self.obsdir
        if not os.path.exists(tmp_reproc_dir):
            os.mkdir(tmp_reproc_dir)

        args = ['aepipeline',
                'indir=%s' % self.obsdir,
                'outdir=%s' % tmp_reproc_dir,
                'steminput=%s' % self.stem,
                'entry_stage=1',
                'exit_stage=2',
                'clobber=yes',
                'instrume=PIN'
                ]

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'aepipeline']).wait()
        proc = subprocess.Popen(args).wait() 

        if not os.path.exists(self.pindir):
            os.mkdir(self.pindir)
        os.rename(tmp_reproc_dir, self.pindir+'/reproc')

    def find_evls(self):
        evls = {}
        evl = glob.glob(self.evlsdir + '/%shxd_*_pinno_cl.evt' % self.stem)
        evls['evl'] = evl[0] if len(evl) > 0 else None
        pse = glob.glob(self.evlsdir + '/%shxd_*_pse_cl.evt' % self.stem)
        evls['pse'] = pse[0] if len(pse) > 0 else None
        bkg = glob.glob(self.evlsdir + '/%shxd_pinbgd.evt*' % self.stem)
        evls['bkg'] = bkg[0] if len(bkg) > 0 else None
        return evls

    def eV2pha(self, eV):
        return int(eV/375. - 1)

    def get_background(self, bkgver=None, bkgfile=None):
        if bkgver is None:
            bkgver = '2.2' if self.date_obs >= datetime(2012, 8, 1) else '2.0'
            print('Using background version %s based on observation date %s' % (bkgver, self.date_obs.strftime('%Y-%m-%d')))
        if bkgfile is None:
            bkgfile = self.evlsdir + '/%shxd_pinbgd.evt.gz' % self.stem

        bkg_url = 'https://heasarc.gsfc.nasa.gov/FTP/suzaku/data/background/pinnxb_ver%s_tuned/%04d_%02d/%s_hxd_pinbgd.evt.gz' % (bkgver, self.date_obs.year, self.date_obs.month, self.stem)
        print('Downloading background event list from %s' % bkg_url)
        urllib.request.urlretrieve(bkg_url, bkgfile)
        if not os.path.exists(bkgfile):
            raise ValueError('Failed to download background file: %s')
        self.evls['bkg'] = bkgfile

    def get_spectrum(self, groupmin=20, cxb_fname='CALC'):
        if not os.path.exists(self.specdir):
            os.mkdir(self.specdir)

        if self.evls['bkg'] is None or not os.path.exists(self.evls['bkg']):
            print('Background file not found. Downloading background file...')
            self.get_background()

        args = ['hxdpinxbpi',
                os.path.relpath(self.evls['evl'], start=self.specdir),
                os.path.relpath(self.evls['pse'], start=self.specdir),
                os.path.relpath(self.evls['bkg'], start=self.specdir),
                self.stem,
                'cxb_fname=%s' % cxb_fname
        ]
        if groupmin is not None:
            args += ['groupmin=%d' % groupmin]

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'hxdpinxbpi']).wait()
        proc = subprocess.Popen(args, cwd=self.specdir).wait()

    