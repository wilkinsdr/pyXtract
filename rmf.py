import numpy as np
import itertools
import operator
import os
from astropy.io import fits


def unpack_rmf(rmf_file):
    with fits.open(rmf_file) as f:
        elow = f['MATRIX'].data['ENERG_LO']
        ehigh = f['MATRIX'].data['ENERG_HI']
        fchan = f['MATRIX'].data['F_CHAN']
        nchan = f['MATRIX'].data['N_CHAN']
        mat = f['MATRIX'].data['MATRIX']
        chan = f['EBOUNDS'].data['CHANNEL']

    en = 0.5 * (elow + ehigh)

    rmf = np.zeros((len(en), len(chan)))

    for i in range(len(en)):
        en_chan = np.concatenate([f + np.arange(n) for f, n in zip(fchan[i], nchan[i])])
        en_mat = mat[i]
        rmf[i][en_chan] = en_mat

    return rmf, en, chan


def pack_rmf(rmf):
    ngrp = []
    fchan = []
    nchan = []
    matrix = []

    for i in range(rmf.shape[0]):
        chan_groups = [[i for i, value in it] for key, it in
                       itertools.groupby(enumerate(rmf[i] > 0), key=operator.itemgetter(1)) if key != 0]

        ngrp.append(len(chan_groups))
        fchan.append([c[0] for c in chan_groups])
        nchan.append([len(c) for c in chan_groups])
        matrix.append(rmf[i][rmf[i] > 0])

    ngrp = np.array(ngrp)
    fchan = np.array(fchan, dtype=object)
    nchan = np.array(nchan, dtype=object)
    matrix = np.array(matrix, dtype=object)

    return ngrp, fchan, nchan, matrix


def read_sparse_matrix(rmf_file, mat_ext_idx):
    """Read the sparse matrix columns from one FITS extension of one RMF file."""
    with fits.open(rmf_file) as f:
        return {
            'fchan':  f[mat_ext_idx].data['F_CHAN'].copy(),
            'nchan':  f[mat_ext_idx].data['N_CHAN'].copy(),
            'matrix': f[mat_ext_idx].data['MATRIX'].copy(),
        }


def sum_sparse_matrices(sparse_data, weights, n_chan):
    """
    Accumulate a weighted sum of sparse RMF matrices row by row.

    sparse_data : list of dicts, each with 'fchan', 'nchan', 'matrix' keys
                  as returned by read_sparse_matrix
    weights     : sequence of floats, one per entry in sparse_data
    n_chan      : number of detector channels (width of the dense working row)

    Returns packed sparse arrays (ngrp, fchan, nchan, matrix) ready to be
    written as FITS columns.
    """
    n_energy = len(sparse_data[0]['fchan'])

    ngrp_out   = np.zeros(n_energy, dtype=np.int16)
    fchan_out  = np.empty(n_energy, dtype=object)
    nchan_out  = np.empty(n_energy, dtype=object)
    matrix_out = np.empty(n_energy, dtype=object)
    row = np.zeros(n_chan, dtype=np.float64)

    for i in range(n_energy):
        if i % 5000 == 0:
            print(f'  Processing row {i}/{n_energy}...')
        row[:] = 0.0
        for sd, w in zip(sparse_data, weights):
            groups = [fc + np.arange(nc) for fc, nc in zip(sd['fchan'][i], sd['nchan'][i])]
            if groups:
                row[np.concatenate(groups)] += w * sd['matrix'][i]

        chan_groups = [
            [j for j, _ in it]
            for key, it in itertools.groupby(enumerate(row > 0), key=operator.itemgetter(1))
            if key != 0
        ]
        ngrp_out[i]   = len(chan_groups)
        fchan_out[i]  = np.array([c[0] for c in chan_groups], dtype=np.int32)
        nchan_out[i]  = np.array([len(c) for c in chan_groups], dtype=np.int32)
        matrix_out[i] = row[row > 0].astype(np.float32)

    return ngrp_out, fchan_out, nchan_out, matrix_out


def add_combined_rmf(rmf_files, weights=None, spec_files=None, outfile='src_comb.rmf', weight_mode='exposure'):
    """
    Add multiple XRISM combined RMF files, each containing multiple MATRIX/EBOUNDS
    extension pairs (PRIMARY, MATRIX, EBOUNDS, MATRIX, EBOUNDS, ...).

    Processes each MATRIX component row by row in the sparse domain to avoid
    building a dense n_energy x n_chan matrix in memory (which would be ~14 GB
    for the 60000-channel XRISM Resolve fine-grid component).

    EBOUNDS extensions are copied directly from the first input file, since they
    are identical across observations from the same instrument configuration.
    """
    if weights is None:
        if spec_files is not None:
            weights = []
            for spec_file in spec_files:
                with fits.open(spec_file) as f:
                    if weight_mode == 'exposure':
                        weights.append(f['SPECTRUM'].header['EXPOSURE'])
                    elif weight_mode == 'counts':
                        weights.append(np.sum(f['SPECTRUM'].data['COUNTS']))
            weights = np.array(weights, dtype=float)
            weights /= np.sum(weights)
        else:
            weights = np.ones(len(rmf_files)) / len(rmf_files)

    with fits.open(rmf_files[0]) as f:
        matrix_indices = [i for i, hdu in enumerate(f) if hdu.name == 'MATRIX']
        ebounds_indices = [i for i, hdu in enumerate(f) if hdu.name == 'EBOUNDS']

    output_hdus = [fits.PrimaryHDU()]

    for mat_ext_idx, eb_ext_idx in zip(matrix_indices, ebounds_indices):
        with fits.open(rmf_files[0]) as f:
            elow = f[mat_ext_idx].data['ENERG_LO'].copy()
            ehigh = f[mat_ext_idx].data['ENERG_HI'].copy()
            n_chan = len(f[eb_ext_idx].data['CHANNEL'])
            src_matrix_hdr = f[mat_ext_idx].header.copy()
            ebounds_hdu = f[eb_ext_idx].copy()

        sparse_data = [read_sparse_matrix(rmf_file, mat_ext_idx) for rmf_file in rmf_files]
        ngrp_out, fchan_out, nchan_out, matrix_out = sum_sparse_matrices(sparse_data, weights, n_chan)

        elow_col   = fits.Column(name='ENERG_LO', format='E', array=elow)
        ehigh_col  = fits.Column(name='ENERG_HI', format='E', array=ehigh)
        ngrp_col   = fits.Column(name='N_GRP',    format='I', array=ngrp_out)
        fchan_col  = fits.Column(name='F_CHAN',   format='PJ', array=fchan_out)
        nchan_col  = fits.Column(name='N_CHAN',   format='PJ', array=nchan_out)
        matrix_col = fits.Column(name='MATRIX',  format='PE', array=matrix_out)

        matrix_hdu = fits.BinTableHDU.from_columns(
            [elow_col, ehigh_col, ngrp_col, fchan_col, nchan_col, matrix_col]
        )

        skip_keys = {'XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'PCOUNT',
                     'GCOUNT', 'TFIELDS', 'CHECKSUM', 'DATASUM', ''}
        for key, value, comment in src_matrix_hdr.cards:
            if (key in skip_keys
                    or key.startswith(('TTYPE', 'TFORM', 'TUNIT', 'TLMIN', 'TLMAX'))
                    or key in ('HISTORY', 'COMMENT')):
                continue
            out_key = ('HIERARCH ' + key) if len(key) > 8 else key
            matrix_hdu.header[out_key] = (value, comment)

        matrix_hdu.header['DETCHANS'] = n_chan
        matrix_hdu.header['NUMGRP']   = int(np.sum(ngrp_out))
        matrix_hdu.header['NUMELT']   = int(sum(len(c) for c in matrix_out))
        matrix_hdu.header['TLMIN4']   = 0
        matrix_hdu.header['TLMAX4']   = n_chan - 1

        output_hdus.append(matrix_hdu)
        output_hdus.append(ebounds_hdu)

    hdulist = fits.HDUList(output_hdus)

    if os.path.exists(outfile):
        os.remove(outfile)
    hdulist.writeto(outfile)
    print('Done')


def add_rmf(rmf_files, weights=None, spec_files=None, outfile='src_comb.rmf', weight_mode='exposure'):
    if weights is None:
        if spec_files is not None:
            weights = []
            for spec_file in spec_files:
                with fits.open(spec_file) as f:
                    if weight_mode == 'exposure':
                        weights.append(f['SPECTRUM'].header['EXPOSURE'])
                    elif weight_mode == 'counts':
                        weights.append(np.sum(f['SPECTRUM'].data['COUNTS']))
            weights = np.array(weights, dtype=float)
            weights /= np.sum(weights)
        else:
            weights = np.ones(len(rmf_files)) / len(rmf_files)

    with fits.open(rmf_files[0]) as f:
        chan = f['EBOUNDS'].data['CHANNEL']
        emin = f['EBOUNDS'].data['E_MIN']
        emax = f['EBOUNDS'].data['E_MAX']
        elow = f['MATRIX'].data['ENERG_LO']
        ehigh = f['MATRIX'].data['ENERG_HI']

        ebounds_hdr = f['EBOUNDS'].header
        matrix_hdr = f['MATRIX'].header

    sparse_data = [read_sparse_matrix(rmf_file, 'MATRIX') for rmf_file in rmf_files]
    ngrp, fchan, nchan, matrix = sum_sparse_matrices(sparse_data, weights, len(chan))

    chan_col = fits.Column(name='CHANNEL', format='J', array=chan)
    emin_col = fits.Column(name='E_MIN', format='E', array=emin)
    emax_col = fits.Column(name='E_MAX', format='E', array=emax)
    elow_col = fits.Column(name='ENERG_LO', format='E', array=elow)
    ehigh_col = fits.Column(name='ENERG_HI', format='E', array=ehigh)
    ngrp_col = fits.Column(name='N_GRP', format='I', array=ngrp)
    fchan_col = fits.Column(name='F_CHAN', format='PI', array=fchan)
    nchan_col = fits.Column(name='N_CHAN', format='PI', array=nchan)
    matrix_col = fits.Column(name='MATRIX', format='PE', array=matrix)

    pri_hdu = fits.PrimaryHDU()

    ebounds_hdu = fits.BinTableHDU.from_columns([chan_col, emin_col, emax_col])
    ebounds_hdu.header['EXTNAME'] = 'EBOUNDS'
    ebounds_hdu.header['TELESCOP'] = ebounds_hdr['TELESCOP']
    ebounds_hdu.header['INSTRUME'] = ebounds_hdr['INSTRUME']
    try:
        ebounds_hdu.header['FILTER'] = ebounds_hdr['FILTER']
    except:
        pass
    ebounds_hdu.header['CHANTYPE'] = ebounds_hdr['CHANTYPE']
    ebounds_hdu.header['DETCHANS'] = len(chan)
    ebounds_hdu.header['HDUCLASS'] = 'OGIP'
    ebounds_hdu.header['HDUCLAS1'] = 'RESPONSE'
    ebounds_hdu.header['HDUCLAS2'] = 'EBOUNDS'
    ebounds_hdu.header['HDUVERS'] = '1.2.0'

    matrix_hdu = fits.BinTableHDU.from_columns([elow_col, ehigh_col, ngrp_col, fchan_col, nchan_col, matrix_col])
    matrix_hdu.header['EXTNAME'] = 'MATRIX'
    matrix_hdu.header['TELESCOP'] = matrix_hdr['TELESCOP']
    matrix_hdu.header['INSTRUME'] = matrix_hdr['INSTRUME']
    try:
        matrix_hdu.header['FILTER'] = matrix_hdr['FILTER']
    except:
        pass
    matrix_hdu.header['CHANTYPE'] = matrix_hdr['CHANTYPE']
    matrix_hdu.header['DETCHANS'] = len(chan)
    matrix_hdu.header['HDUCLASS'] = 'OGIP'
    matrix_hdu.header['HDUCLAS1'] = 'RESPONSE'
    matrix_hdu.header['HDUCLAS2'] = 'RSP_MATRIX'
    matrix_hdu.header['HDUVERS'] = '1.3.0'
    matrix_hdu.header['TLMIN4'] = np.min(chan)
    matrix_hdu.header['TLMAX4'] = np.max(chan)
    matrix_hdu.header['NUMGRP'] = np.sum(ngrp)
    matrix_hdu.header['NUMELT'] = np.sum([np.sum(c) for c in nchan])
    try:
        matrix_hdu.header['HDUCLAS3'] = matrix_hdr['HDUCLAS3']
    except:
        pass

    hdulist = fits.HDUList([pri_hdu, ebounds_hdu, matrix_hdu])

    if (os.path.exists(outfile)):
        os.remove(outfile)

    hdulist.writeto(outfile)
    print("Done")
