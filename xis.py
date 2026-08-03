"""
Suzaku XIS extractor
"""
import glob
import subprocess
from .xselect import Xselect
from .spec_util import *
import os
import re
from astropy.io import fits

xis_names = {'xi0': 'XIS0', 'xi1': 'XIS1', 'xi2': 'XIS2', 'xi3': 'XIS3'}

class XisExtractor(object):
    #
    # class to extract data products from Suzaku XIS observations
    #

    def __init__(self, obsdir, filt_level=1, ccd=None, mode=None, evl_dir='reproc', run_reduction=False, suffix=None):
        self.obsdir = obsdir
        self.xisdir = obsdir + '/xis'

        self.stem = 'ae%s' % obsdir

        self.suffix = suffix

        self.evlsdir = self.xisdir + '/%s' % evl_dir + ('_%s' % suffix if suffix is not None else '')
        self.specdir = self.xisdir + '/spectra' + ('_%s' % suffix if suffix is not None else '')
        self.lcdir = self.xisdir + '/lightcurves' + ('_%s' % suffix if suffix is not None else '')
        self.regiondir = self.xisdir + '/regions'
        self.gtidir = self.xisdir + '/gti'
        self.ufdir = self.xisdir + '/event_uf'
        self.auxdir = self.obsdir + '/auxil'

        self.evls = self.find_evls()
        # if len(self.evls) == 0:
        #     print('Filter level %d event list not available, falling back to level 1' % filt_level)
        #     self.evls = self.find_evls(filt_level='', ccd=ccd, mode=mode)
        #     self.filt_level = 1

        if not os.path.exists(self.regiondir):
            os.mkdir(self.regiondir)

        self.regions = self.populate_regions()

        with fits.open(self.evls['xi0'][0]) as f:
            self.ra_nom = f[0].header['RA_NOM']
            self.dec_nom = f[0].header['DEC_NOM']

    def reprocess(self, filt_level=2, ccd=None, mode=None):
        tmp_reproc_dir = 'tmp_reproc_%s' % self.obsdir
        if not os.path.exists(tmp_reproc_dir):
            os.mkdir(tmp_reproc_dir)

        args = ['aepipeline',
                'indir=%s' % self.obsdir,
                'outdir=%s' % tmp_reproc_dir,
                'steminput=%s' % self.stem,
                'entry_stage=1',
                'exit_stage=2',
                'clobber=yes',
                'instrume=XIS'
                ]

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'aepipeline']).wait()
        proc = subprocess.Popen(args).wait() 

        os.rename(tmp_reproc_dir, self.xisdir+'/reproc')

    def find_evls(self):
        evls = {}
        for inst in ['xi0', 'xi1', 'xi2', 'xi3']:
            evls[inst] = sorted(glob.glob(self.evlsdir + '/%s%s_*_cl.evt*' % (self.stem, inst)))
        return evls


    def populate_regions(self):
        regions = {}
        for inst in ['xi0', 'xi1', 'xi2', 'xi3']:
            regions[inst] = {}
            regions[inst]['src'] = sorted(glob.glob(self.regiondir + '/src_%s.reg' % inst))
            regions[inst]['bkg'] = sorted(glob.glob(self.regiondir + '/bkg_%s.reg' % inst))

            regions[inst]['src'] = regions[inst]['src'][0] if len(regions[inst]['src']) > 0 else None
            regions[inst]['bkg'] = regions[inst]['bkg'][0] if len(regions[inst]['bkg']) > 0 else None
        return regions

    def eV2pha(self, eV):
        return int(eV/3.65)

    def extract_spectrum(self, evl, spec_file=None, src_region=None, bkg_file=None, bkg_region=None):
        if os.path.exists(spec_file):
            os.remove(spec_file)

        with Xselect(mission='SUZAKU') as xsl:
            if isinstance(evl, list):
                xsl.command('read event')
                xsl.command(os.path.dirname(evl[0]))
                xsl.command(','.join([os.path.basename(e) for e in evl]))
            else:
                xsl.read_event(evl)
            xsl.command('filter region %s' % src_region)
            xsl.command('extract spectrum')
            xsl.command('save spectrum %s resp=no group=no' % spec_file)
            if bkg_file is not None:
                xsl.command('clear region')
                xsl.command('filter region %s resp=no group=no' % bkg_region)
                xsl.command('extract spectrum')
                xsl.command('save spectrum %s resp=no group=no' % bkg_file)

    def make_rmf(self, outfile=None, spec_file=None):
        args = ['xisrmfgen',
                'phafile=%s' % spec_file,
                'outfile=%s' % outfile,
                'clobber=yes',
                'mode=h'
        ]

        proc = subprocess.Popen(['punlearn', 'xisrmfgen']).wait()
        proc = subprocess.Popen(args).wait()

    def make_arf_pointsource(self, outfile=None, inst='xi0', ra=None, dec=None, spec_file=None, rmf_file=None, gti_file=None, region_file=None, numphoton=400000):
        if ra is None:
            ra = self.ra_nom
        if dec is None:
            dec = self.dec_nom

        args = ['xissimarfgen',
                'instrume=%s' % xis_names[inst],
                'pointing=AUTO',
                'source_mode=J2000',
                'source_ra=%0.5f' % ra,
                'source_dec=%0.5f' % dec,
                'num_region=1',
                'region_mode=SKYREG',
                'regfile1=%s' % region_file,
                'arffile1=%s' % outfile,
                'limit_mode=NUM_PHOTON',
                'num_photon=%d' % numphoton,
                'phafile=%s' % spec_file,
                'detmask=none',
                'gtifile=%s' % (gti_file if gti_file is not None else spec_file),
                'attitude=%s' % (self.auxdir + '/%s.att' % self.stem),
                'rmffile=%s' % rmf_file,
                'estepfile=default',
                'clobber=yes',
                'mode=h'
        ]

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'xaarfgen']).wait()
        proc = subprocess.Popen(args).wait()
        # workaround for a bug in xaarfgen command line processing
        #print(' '.join(args))
        #proc = subprocess.Popen(' '.join(args), shell=True).wait()

    def get_spectrum(self, instruments=['xi0', 'xi1', 'xi2', 'xi3'], src_region=None, bkg_region=None, ra=None, dec=None, suffix=None, extract_spectrum=True, make_rmf=True, make_arf=True, link_resp=True, opt_bin=True):
        if not os.path.exists(self.specdir):
            os.mkdir(self.specdir)

        for inst in instruments:
            name_arr = [self.stem, inst]
            if suffix is not None:
                name_arr.append(suffix)
            
            spec_filename = '_'.join(name_arr) + '_src.pha'
            spec_file = self.specdir + '/' + spec_filename

            bkg_filename = '_'.join(name_arr) + '_bkg.pha'
            bkg_file = self.specdir + '/' + bkg_filename

            rmf_filename = '_'.join(name_arr) + '_src.rmf'
            rmf_file = self.specdir + '/' + rmf_filename

            arf_filename = '_'.join(name_arr) + '_src.arf'
            arf_file = self.specdir + '/' + arf_filename

            grp_file = spec_file.replace('.pha', '_opt.grp')

            if len(self.evls[inst]) == 0:
                print('No event list found for %s, skipping' % inst)
                continue

            if src_region is None:
                src_region = self.regions[inst]['src']
                print('Using default source region for %s: %s' % (inst, src_region))
            if bkg_region is None:
                bkg_region = self.regions[inst]['bkg']
                print('Using default background region for %s: %s' % (inst, bkg_region))

            if extract_spectrum:
                self.extract_spectrum(self.evls[inst], spec_file=spec_file, src_region=src_region, bkg_file=bkg_file, bkg_region=bkg_region)
            if make_rmf:
                self.make_rmf(rmf_file, spec_file=spec_file)
            if make_arf:
                self.make_arf_pointsource(outfile=arf_file, inst=inst, ra=None, dec=None, spec_file=spec_file, rmf_file=rmf_file, gti_file=None, region_file=src_region, numphoton=400000)
            if link_resp:
                link_spectra(spec_file, bkg=bkg_file, rmf=rmf_file, arf=arf_file)
            if opt_bin:
                group_spec(grp_file, spec_file, rmffile=rmf_file, grptype='opt')

    def extract_lightcurve(self, evl, lc_file=None, tbin=128.0, exposure=0.0, energy=(300, 12000), bkg_file=None, src_region=None, bkg_region=None, suffix=None):
        if evl is None:
            evl = self.evls[0]
        if lc_file is None:
            name_arr = ['%sxtd' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            if suffix is not None:
                name_arr.append(suffix)
            name_arr.append('tbin%g' % tbin)
            if energy is not None:
                name_arr.append('en%g-%g' % energy)

            lc_filename = '_'.join(name_arr) + '_src.lc'
            lc_file = self.lcdir + '/' + lc_filename
            bkg_filename = '_'.join(name_arr) + '_bkg.lc'
            bkg_file = self.lcdir + '/' + bkg_filename

        if os.path.exists(lc_file):
            os.remove(lc_file)

        with Xselect(mission='SUZAKU') as xsl:
            xsl.read_event(evl)
            xsl.command('set binsize %g' % tbin)
            if energy is not None:
                xsl.command('filter pha_cutoff %d %d' % (self.eV2pha(energy[0]), self.eV2pha(energy[1]) - 1))
            xsl.command('filter region %s' % src_region)
            xsl.command('extract curve exposure=%g' % exposure)
            xsl.command('save curve %s' % lc_file)
            if bkg_file is not None:
                xsl.command('clear region')
                xsl.command('filter region %s' % bkg_region)
                xsl.command('extract curve exposure=%g' % exposure)
                xsl.command('save curve %s' % bkg_file)

    def get_lightcurve(self, tbin=128.0, exposure=0.0, energy=(300, 12000), suffix=None, extract_dir=None, src_region=None, bkg_region=None):
        for evl in self.evls:
            name_arr = ['%sxtd' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.path.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            if suffix is not None:
                name_arr.append(suffix)
            name_arr.append('tbin%g' % tbin)
            if energy is not None:
                name_arr.append('en%g-%g' % energy)

            if src_region is None:
                if len(self.evls) > 1 and os.path.exists(self.regiondir + '/src_%s.reg' % evl_name):
                    src_region = self.regiondir + '/src_%s.reg' % evl_name
                else:
                    src_region = self.regiondir + '/src.reg'
            if bkg_region is None:
                if len(self.evls) > 1 and os.path.exists(self.regiondir + '/bkg_%s.reg' % evl_name):
                    bkg_region = self.regiondir + '/bkg_%s.reg' % evl_name
                else:
                    bkg_region = self.regiondir + '/bkg.reg'

            if extract_dir is None:
                extract_dir = self.lcdir

            if not os.path.exists(extract_dir):
                os.mkdir(extract_dir)

            lc_filename = '_'.join(name_arr) + '_src.lc'
            lc_file = self.lcdir + '/' + lc_filename
            bkg_filename = '_'.join(name_arr) + '_bkg.lc'
            bkg_file = self.lcdir + '/' + bkg_filename

            self.extract_lightcurve(evl, lc_file, tbin, exposure, energy, bkg_file, src_region, bkg_region)

    
