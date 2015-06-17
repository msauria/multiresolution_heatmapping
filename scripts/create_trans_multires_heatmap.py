#!/usr/bin/env python

import sys
import argparse as ap
import struct
import binascii

import numpy
from scipy import weave

import hifive


def main():
    parser = generate_parser()
    args = parser.parse_args()
    n_res = (numpy.log(args.lres) - numpy.log(args.hres)) / numpy.log(args.zoom)
    if not numpy.isclose(n_res, numpy.round(n_res)):
        print >> sys.stderr, ("The high to low resolution ratio must an integer power of the zoom factor (H / L = Z**N).")
        sys.exit(1)
    n_res = int(numpy.round(n_res)) + 1
    zoom2 = args.zoom * args.zoom
    hic = hifive.HiC(args.project)
    # Fetch highest level resolution data
    data, mapping1, mapping2 = hic.trans_heatmap(chrom1=args.chrom, chrom2=args.chrom2, binsize=args.hres,
                                    datatype=args.dtype, arraytype='full', returnmapping=True)
    span1 = mapping1[-1, 1] - mapping1[0, 0]
    bins1 = (span1 - 1) / args.lres + 1
    mapping_start1 = max(0, int(round((mapping1[-1, 1] + mapping1[0, 0]) / 2 - (bins1 * args.lres) / 2.0)))
    mapping1[:, :2] -= mapping_start1
    mapping_stop1 = bins1 * args.lres + mapping_start1
    mids1 = (mapping1[:, 0] + mapping1[:, 1]) / 2
    span2 = mapping2[-1, 1] - mapping2[0, 0]
    bins2 = (span2 - 1) / args.lres + 1
    mapping_start2 = max(0, int(round((mapping2[-1, 1] + mapping2[0, 0]) / 2 - (bins2 * args.lres) / 2.0)))
    mapping2[:, :2] -= mapping_start2
    mapping_stop2 = bins2 * args.lres + mapping_start2
    mids2 = (mapping2[:, 0] + mapping2[:, 1]) / 2
    all_data = []
    all_indices = []
    # Create top-level resolution map as rectangle
    res = args.lres
    n_bins = int(mapping1[-1, 1] / res + 1)
    m_bins = int(mapping2[-1, 1] / res + 1)
    n = int(data.shape[0])
    m = int(data.shape[1])
    stop1 = n_bins * res
    stop2 = m_bins * res
    minobs = args.minobs
    new_data = numpy.zeros(n_bins * m_bins, dtype=numpy.float32)
    new_indices = numpy.zeros(new_data.shape[0], dtype=numpy.int32) - 1
    temp_indices = numpy.zeros(new_data.shape[0], dtype=numpy.int32)
    indices1 = numpy.searchsorted(mids1, numpy.linspace(0, stop1, n_bins + 1)).astype(numpy.int64)
    indices2 = numpy.searchsorted(mids2, numpy.linspace(0, stop2, m_bins + 1)).astype(numpy.int64)
    nan = numpy.nan
    code = """
        long long int index0, index1, index2, index3, i, j, k, l;
        long double expected;
        long int observed;
        for(i = 0; i < n_bins; i++){
            index0 = i * m_bins;
            for(j = 0; j < m_bins; j++){
                index1 = index0 + j;
                temp_indices[index1] = index1;
                expected = 0.0;
                observed = 0;
                for(k = indices1[i]; k < indices1[i + 1]; k++){
                    index2 = k * m * 2;
                    for(l = indices2[j]; l < indices2[j + 1]; l++){
                        index3 = index2 + l * 2;
                        observed += data[index3];
                        expected += data[index3 + 1];
                    }
                }
                if(observed < minobs){
                    new_data[index1] = nan;
                } else {
                    new_data[index1] = log2(observed / expected);
                }
            }
        }
    """
    weave.inline(code, ['indices1', 'indices2', 'data', 'new_data', 'n_bins', 'm_bins',
                        'temp_indices', 'nan', 'm', 'minobs'])
    all_data = [new_data]
    all_indices = [new_indices]
    pos = n_bins * m_bins
    zoom = args.zoom
    # Create all lower-level maps
    code = """
        long long int pos2 = 0;
        long int zoom2 = zoom * zoom;
        long long int observed, index0, index1, index2, index3, start1, start2, stop1, stop2, valid, i, j, k, l, p;
        long double expected;
        for(i = 0; i < t_bins; i++){
            if(prev_data[i] != nan){
                valid = 0;
                // find positions of lower-resolution bin
                start1 = (temp_indices[i] / old_m_bins) * zoom;
                stop1 = start1 + zoom;
                start2 = (temp_indices[i] % old_m_bins) * zoom;
                stop2 = start2 + zoom;
                // zero out temp_data array
                for(j = 0; j < zoom2; j++){
                    temp_data[j] = 0.0;
                }
                // find square grid of values
                for(j = start1; j < stop1; j++){
                    index2 = (j - start1) * zoom;
                    for(k = start2; k < stop2; k++){
                        index3 = index2 + k - start2;
                        observed = 0;
                        expected = 0.0;
                        for(l = indices1[j]; l < indices1[j + 1]; l++){
                            index0 = l * m * 2;
                            for(p = indices2[k]; p < indices2[k + 1]; p++){
                                index1 = index0 + p * 2;
                                observed += data[index1];
                                expected += data[index1 + 1];
                            }
                        }
                        if(observed < minobs){
                            temp_data[index3] = nan;
                        }else{
                            temp_data[index3] = log2(observed / expected);
                            valid += 1;
                        }
                    }
                }
                // if there are valid values at this resolution, add the values to the current resolution data
                if(valid > 0){
                    for(j = 0; j < zoom2; j++){
                        new_data[pos2] = temp_data[j];
                        new_temp_indices[pos2] = (start1 + j / zoom) * m_bins + start2 + j % zoom;
                        pos2++;
                    }
                    // fill in the higher resolution index at the lower resolution bin
                    prev_indices[i] = pos;
                    pos += zoom2;
                }else{
                    prev_indices[i] = -1;
                }
            }
        }
    """
    for h in range(1, n_res):
        new_temp_indices = numpy.zeros(all_data[-1].shape[0] * args.zoom ** 2, dtype=numpy.int32)
        new_temp_indices.fill(-1)
        new_data = numpy.zeros(all_data[-1].shape[0] * args.zoom ** 2, dtype=numpy.float32)
        res /= args.zoom
        old_m_bins = m_bins
        n_bins = stop1 / res
        m_bins = stop2 / res
        t_bins = all_data[-1].shape[0]
        indices1 = numpy.searchsorted(mids1, numpy.linspace(0, stop1, n_bins + 1)).astype(numpy.int64)
        indices2 = numpy.searchsorted(mids2, numpy.linspace(0, stop2, m_bins + 1)).astype(numpy.int64)
        prev_data = all_data[-1]
        prev_indices = all_indices[-1]
        temp_data = numpy.zeros(zoom ** 2, dtype=numpy.float32)
        weave.inline(code, ['data', 'new_data', 'new_temp_indices', 'prev_data', 'prev_indices', 'indices1',
                            'indices2', 'm', 'temp_data', 'temp_indices', 'pos', 'zoom', 't_bins', 'm_bins',
                            'old_m_bins', 'nan', 'minobs'])
        where = numpy.where(new_temp_indices >= 0)[0]
        pos += where.shape[0]
        temp_indices = new_temp_indices[where]
        all_data.append(new_data[where])
        if h < n_res - 1:
            all_indices.append(numpy.zeros(all_data[-1].shape[0], dtype=numpy.int32))
            all_indices[-1].fill(-1)
    data_bins = 0
    minscore = numpy.inf
    maxscore = -numpy.inf
    for i in range(len(all_data)):
        data_bins += all_data[i].shape[0]
        minscore = min(numpy.amin(all_data[i][numpy.where(numpy.logical_not(numpy.isnan(all_data[i])))]), minscore)
        maxscore = max(numpy.amax(all_data[i][numpy.where(numpy.logical_not(numpy.isnan(all_data[i])))]), maxscore)
    total_bins = data_bins
    for i in range(len(all_indices)):
        total_bins += all_indices[i].shape[0]
    # Header contains:
    # magic number - 42054205 (8 btyes)
    # Header size - int32 (either 56 or 68 bytes for cis or trans, respectively)
    # Top resolution - int32
    # Bottom resolution - int32
    # Zoom factor - int32
    # Minimum observation cutoff - int32
    # Minimum score - int32
    # Maximum score - int32
    # Trans flag - int32
    # Minimum coordinate - int32
    # if trans: Second chromosome minimum coordinate - int32
    # Maximum coordinate - int32
    # if trans: Second chromosome maximum coordinate - int32
    # Top level number of divisions - int32
    # if trans: Second chromosome top level number of bins - int32
    # Total data bins (indices offset) - int32
    # Total bins (data + indices) - int32
    output = open(args.output, 'wb')
    output.write(binascii.a2b_hex('42054205'))
    output.write(struct.pack('iiiiiffiiiiiiiii', 68, args.lres, args.hres, args.zoom, args.minobs,
                             minscore, maxscore, 1, mapping_start1, mapping_start2, mapping_stop1,
                             mapping_stop2, bins1, bins2, data_bins, total_bins))
    for i in range(len(all_data)):
        output.write(all_data[i].tostring())
    for i in range(len(all_indices)):
        output.write(all_indices[i].tostring())
    output.close()

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Produce a multi-resolution heatmap file for a chromosome from HiC data using HiFive"
    parser = ap.ArgumentParser(description=description)
    parser.add_argument("-L", "--lowest-resolution", dest="lres", required=False, type=int, default=1000000,
        action='store', help="The lowest resolution bin size to create. [default: %(default)s]")
    parser.add_argument("-H", "--highest-resolution", dest="hres", required=False, type=int, default=1600,
        action='store', help="The highest resolution bin size to create. [default: %(default)s]")
    parser.add_argument("-z", "--zoom", dest="zoom", required=False, type=int, default=5,
        action='store', help="The factor by which each higher resolution heatmap is zoomed in. [default: %(default)s]")
    parser.add_argument("-m", "--min-observations", dest="minobs", required=False, type=int, default=10,
        action='store', help="The minimum number of observations needed for a bin to be considered valid. [default: %(default)s]")
    parser.add_argument("-d", "--datatype", dest="dtype", required=False, default='fend',
        choices=['fend', 'enrichment'], help="The type of data to return. [default: %(default)s]")
    parser.add_argument("-c", "--chromosome", dest="chrom", required=True, type=str,
        action='store', help="The chromosome to construct a heatmap for.")
    parser.add_argument("-C", "--chromosome2", dest="chrom2", required=True, type=str,
        action='store', help="The second chromosome to construct a heatmap for.")
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive HiC project to construct the heatmap from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write the multi-resolution HiFive heatmap to.")
    return parser


if __name__ == "__main__":
    main()

