"""
XRISM Resolve extractor
"""
import glob
import subprocess
from .xselect import Xselect
from .spec_util import *
import os
import re
from astropy.io import fits

class XtendExtractor(object):
    #
    # class to extract data products from XRISM Xtend observations
    #

    def __init__(self, obsdir, filt_level=1, ccd=None, mode=None, evl_dir='event_cl', run_reduction=False, suffix=None):
        self.obsdir = obsdir
        self.xtddir = obsdir + '/xtend'

        self.stem = 'xa%s' % obsdir

        self.suffix = suffix

        self.evlsdir = self.xtddir + '/%s' % evl_dir + ('_%s' % suffix if suffix is not None else '')
        self.specdir = self.xtddir + '/spectra' + ('_%s' % suffix if suffix is not None else '')
        self.lcdir = self.xtddir + '/lightcurves' + ('_%s' % suffix if suffix is not None else '')
        self.regiondir = self.xtddir + '/regions'
        self.gtidir = self.xtddir + '/gti'
        self.ufdir = self.xtddir + '/event_uf'
        self.scratchdir = self.xtddir + '/scratch'

        if not os.path.exists(self.scratchdir):
            os.mkdir(self.scratchdir)

        self.ehk = glob.glob(self.obsdir + '/auxil/%s.ehk*' % self.stem)[0]

        self.evls = self.find_evls(filt_level=filt_level, ccd=ccd, mode=mode)
        if len(self.evls) == 0:
            print('Filter level %d event list not available, falling back to level 1' % filt_level)
            self.evls = self.find_evls(filt_level='', ccd=ccd, mode=mode)
            self.filt_level = 1

        with fits.open(self.evls[0]) as f:
            self.ra_nom = f[0].header['RA_NOM']
            self.dec_nom = f[0].header['DEC_NOM']

    def find_evls(self, filt_level='', ccd=None, mode=None):
        if filt_level == 1:
            filt_level = ''
        evls = [evl for evl in sorted(glob.glob(self.evlsdir + '/%sxtd_p*_cl%s.evt*' % (self.stem, filt_level)))]
        if ccd == 1 or ccd == 2:
            if mode == 'window':
                evls = [e for e in evls if 'p0311' in e]
            elif mode == 'full_burst':
                evls = [e for e in evls if 'p0312' in e]
            elif mode == 'window_burst':
                evls = [e for e in evls if 'p0313' in e]
            else:
                evls = [e for e in evls if 'p031' in e]
        if ccd == 3 or ccd == 4:
            evls = [e for e in evls if 'p0320' in e]
        return evls

    def filter(self, rise_time=True, anomalous_Ls=True, frame_events=True, adj_channels=False, filt_level=2):
        pass

    def eV2pha(self, eV):
        return int(0.1667 * eV)

    def extract_spectrum(self, evl=None, spec_file=None, bkg_file=None, src_region=None, bkg_region=None, suffix=None):
        if evl is None:
            evl = self.evls[0]
        if spec_file is None:
            name_arr = ['%sxtd' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            name_arr.append(''.join(grade))
            if suffix is not None:
                name_arr.append(suffix)

            spec_filename = '_'.join(name_arr) + '_src.pha'
            spec_file = self.specdir + '/' + spec_filename
            bkg_filename = '_'.join(name_arr) + '_bkg.pha'
            bkg_file = self.specdir + '/' + bkg_filename

        if os.path.exists(spec_file):
            os.remove(spec_file)

        with Xselect(mission='XRISM') as xsl:
            xsl.read_event(evl)
            xsl.command('filter region %s' % src_region)
            xsl.command('extract spectrum')
            xsl.command('save spectrum %s' % spec_file)
            if bkg_file is not None:
                xsl.command('clear region')
                xsl.command('filter region %s' % bkg_region)
                xsl.command('extract spectrum')
                xsl.command('save spectrum %s' % bkg_file)

    def make_rmf(self, outfile=None, spec_file=None):
        if spec_file is None:
            evl - self.evls[0]
            name_arr = ['%sxtd' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            name_arr.append(''.join(grade))
            if suffix is not None:
                name_arr.append(suffix)

            src_spec_filename = '_'.join(name_arr) + '_src.pha'
            spec_file = self.specdir + '/' + src_spec_filename

        if outfile is None:
            outfile = spec_file.replace('.pha', '.rmf')

        args = ['xtdrmf',
                spec_file,
                outfile]

        proc = subprocess.Popen(['punlearn', 'xtdrmf']).wait()
        proc = subprocess.Popen(args).wait()

    def exposure_map(self, outfile=None, evl=None):
        if evl is None:
            evl = self.evls[0]

        name_arr = ['%sxtd' % self.stem]
        if len(self.evls) > 1:
            evl_name = os.path.basename(evl).split('_')[1]
            name_arr.append(evl_name)
        if outfile is None:
            outfile = self.specdir + '/' + '_'.join(name_arr) + '.expo'

        modestr = re.search('(p[0-9]+)', os.path.basename(evl)).group(1)
        bagimgfile = glob.glob(self.ufdir + '/%sxtd_%s.bimg*' % (self.stem, modestr))[0]

        args = ['xaexpmap',
                'ehkfile=%s' % self.ehk,
                'gtifile=%s' % evl,
                'instrume=XTEND',
                'badimgfile=%s' % bagimgfile,
                'pixgtifile=NONE',
                'outfile=%s' % outfile,
                'outmaptype=EXPOSURE',
                'delta=20.0',
                'numphi=1']

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'xaexpmap']).wait()
        proc = subprocess.Popen(args).wait()

    def make_arf_pointsource(self, outfile=None, xrtevtfile=None, ra=None, dec=None, expomap=None, regionfile=None, rmffile=None, erange='0.3 15.0 0 0', numphoton=600000, minphoton=100, seed=7):
        if ra is None:
            ra = self.ra_nom
        if dec is None:
            dec = self.dec_nom

        if xrtevtfile is None:
            xrtevtfile = self.scratchdir + '/' + os.path.basename(outfile).replace('.arf', '_raytrace_pt.evt')

        if expomap is None or not os.path.exists(expomap):
            raise AssertionError('Exposure map does not exist')
        if regionfile is None or not os.path.exists(regionfile):
            raise AssertionError('Region file does not exist')
        if rmffile is None or not os.path.exists(rmffile):
            raise AssertionError('RMF does not exist')

        args = ['xaarfgen',
                'xrtevtfile=%s' % xrtevtfile,
                'source_ra=%0.5f' % ra,
                'source_dec=%0.5f' % dec,
                'telescop=XRISM',
                'instrume=XTEND',
                'emapfile=%s' % expomap,
                'regmode=RADEC',
                'regionfile=%s' % regionfile,
                'sourcetype=POINT',
                'rmffile=%s' % rmffile,
                'erange="%s"' % erange,
                'outfile=%s' % outfile,
                'numphoton=%d' % numphoton,
                'minphoton=%d' % minphoton,
                'teldeffile=CALDB',
                'qefile=CALDB',
                'contamifile=CALDB',
                'obffile=CALDB',
                'fwfile=CALDB',
                'gatevalvefile=CALDB',
                'onaxisffile=CALDB',
                'onaxiscfile=CALDB',
                'mirrorfile=CALDB',
                'obstructfile=CALDB',
                'frontreffile=CALDB',
                'backreffile=CALDB',
                'pcolreffile=CALDB',
                'scatterfile=CALDB',
                'imgfile=NONE',
                'seed=%d' % seed,
                'clobber=yes',
                'mode=h']

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'xaarfgen']).wait()
        # proc = subprocess.Popen(args).wait()
        # workaround for a bug in xaarfgen command line processing
        print(' '.join(args))
        proc = subprocess.Popen(' '.join(args), shell=True).wait()

    def get_spectrum(self, src_region=None, bkg_region=None, ra=None, dec=None, suffix=None, extract_spectrum=True, make_rmf=True, make_arf=True, link_resp=True, opt_bin=True):
        for evl in self.evls:
            name_arr = ['%sxtd' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.path.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            if suffix is not None:
                name_arr.append(suffix)

            if not os.path.exists(self.specdir):
                os.mkdir(self.specdir)

            spec_filename = '_'.join(name_arr) + '_src.pha'
            spec_file = self.specdir + '/' + spec_filename

            bkg_filename = '_'.join(name_arr) + '_bkg.pha'
            bkg_file = self.specdir + '/' + bkg_filename

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

            if not os.path.exists(src_region):
                raise AssertionError('Source region file does not exist')
            if not os.path.exists(bkg_region):
                bkg_region = None

            rmf_file = spec_file.replace('.pha', '.rmf')

            expomap_file = self.scratchdir + '/' + '_'.join(name_arr) + '.expo'
            region_file = self.scratchdir + '/' + self.stem + '_DET.reg'
            arf_file = self.specdir + '/' + '_'.join(name_arr) + '_src.arf'

            grp_file = spec_file.replace('.pha', '_opt.grp')

            if extract_spectrum:
                self.extract_spectrum(evl, spec_file, bkg_file, src_region, bkg_region, suffix)
            if make_rmf:
                self.make_rmf(rmf_file, spec_file=spec_file)
            if make_arf:
                self.exposure_map(expomap_file, evl)
                self.make_arf_pointsource(arf_file, ra=ra, dec=dec, expomap=expomap_file, regionfile=src_region, rmffile=rmf_file)
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

        with Xselect(mission='XRISM') as xsl:
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

    
