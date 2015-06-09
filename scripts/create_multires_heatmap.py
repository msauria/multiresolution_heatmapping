#!/usr/bin/env python

import sys
import argparse as ap
import struct

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
    data, mapping = hic.cis_heatmap(chrom=args.chrom, binsize=args.hres,
                                    datatype=args.dtype, arraytype='full', returnmapping=True)
    span = mapping[-1, 1] - mapping[0, 0]
    bins = (span - 1) / args.lres + 1
    mapping_start = int(round((mapping[-1, 1] + mapping[0, 0]) / 2 - (bins * args.lres) / 2.0))
    mapping[:, :2] -= mapping_start
    mids = (mapping[:, 0] + mapping[:, 1]) / 2
    all_data = []
    all_indices = []
    # Create top-level resolution map as upper-triangle including diagonal
    res = args.lres
    n_bins = int(mapping[-1, 1] / res + 1)
    n = int(data.shape[0])
    stop = n_bins * res
    minobs = args.minobs
    new_data = numpy.zeros((n_bins * (n_bins + 1)) / 2, dtype=numpy.float32)
    new_indices = numpy.zeros(new_data.shape[0], dtype=numpy.int32) - 1
    temp_indices = numpy.zeros(new_data.shape[0], dtype=numpy.int32)
    indices = numpy.searchsorted(mids, numpy.linspace(0, stop, n_bins + 1)).astype(numpy.int64)
    nan = numpy.nan
    code = """
        long long int index0, index1, index2, index3, i, j, k, l;
        long double expected;
        long int observed;
        for(i = 0; i < n_bins; i++){
            index0 = i * (n_bins - 1) - (i * (i - 1)) / 2;
            for(j = i; j < n_bins; j++){
                index1 = index0 + j;
                temp_indices[index1] = i * n_bins + j;
                expected = 0.0;
                observed = 0;
                for(k = indices[i]; k < indices[i + 1]; k++){
                    index2 = k * n * 2;
                    for(l = indices[j]; l < indices[j + 1]; l++){
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
    weave.inline(code, ['indices', 'data', 'new_data', 'n_bins', 'temp_indices', 'nan', 'n', 'minobs'])
    all_data = [new_data]
    all_indices = [new_indices]
    pos = (n_bins * (n_bins + 1)) / 2
    zoom = args.zoom
    t_bins = data.shape[0]
    # Create all lower-level maps
    code = """
        long long int pos2 = 0;
        long int zoom2 = zoom * zoom;
        long long int observed, index0, index1, index2, index3, start1, start2, stop1, stop2, valid, i, j, k, l, m;
        long double expected;
        for(i = 0; i < n; i++){
            if(prev_data[i] != nan){
                valid = 0;
                // find positions of lower-resolution bin
                start1 = (temp_indices[i] / old_n_bins) * zoom;
                stop1 = start1 + zoom;
                start2 = (temp_indices[i] % old_n_bins) * zoom;
                stop2 = start2 + zoom;
                // if on the diagonal, don't calculate values below the diagonal
                if(start1 == start2){
                    // zero out temp_diag array
                    for(j = 0; j < (zoom * (zoom + 1)) / 2; j++){
                        temp_data[j] = 0.0;
                    }
                    for(j = start1; j < stop1; j++){
                        index2 = (j - start1) * (zoom - 1) - ((j - start1) * (j - start1 - 1)) / 2;
                        for(k = j; k < stop2; k++){
                            index3 = index2 + k - start2;
                            observed = 0;
                            expected = 0.0;
                            for(l = indices[j]; l < indices[j + 1]; l++){
                                index0 = l * t_bins * 2;
                                for(m = indices[k]; m < indices[k + 1]; m++){
                                    index1 = index0 + m * 2;
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
                        l = 0;
                        for(j = 0; j < zoom; j++){
                            for(k = j; k < zoom; k++){
                                new_data[pos2] = temp_data[l];
                                new_temp_indices[pos2] = (start1 + j) * n_bins + start2 + k;
                                pos2++;
                                l++;
                            }
                        }
                        // fill in the higher resolution index at the lower resolution bin
                        prev_indices[i] = pos;
                        pos += (zoom * (zoom + 1)) / 2;
                    }else{
                        prev_indices[i] = -1;
                    }
                }else{
                    // zero out temp_data array
                    for(j = 0; j < zoom2; j++){
                        temp_data[j] = 0.0;
                    }
                    // find square grid of values off the diagonal
                    for(j = start1; j < stop1; j++){
                        index2 = (j - start1) * zoom;
                        for(k = start2; k < stop2; k++){
                            index3 = index2 + k - start2;
                            observed = 0;
                            expected = 0.0;
                            for(l = indices[j]; l < indices[j + 1]; l++){
                                index0 = l * t_bins * 2;
                                for(m = indices[k]; m < indices[k + 1]; m++){
                                    index1 = index0 + m * 2;
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
                            new_temp_indices[pos2] = (start1 + j / zoom) * n_bins + start2 + j % zoom;
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
        }
    """
    for h in range(1, n_res):
        new_temp_indices = numpy.zeros(all_data[-1].shape[0] * args.zoom ** 2, dtype=numpy.int32)
        new_temp_indices.fill(-1)
        new_data = numpy.zeros(all_data[-1].shape[0] * args.zoom ** 2, dtype=numpy.float32)
        res /= args.zoom
        old_n_bins = n_bins
        n_bins = stop / res
        n = all_data[-1].shape[0]
        indices = numpy.searchsorted(mids, numpy.linspace(0, stop, n_bins + 1)).astype(numpy.int64)
        prev_data = all_data[-1]
        prev_indices = all_indices[-1]
        temp_data = numpy.zeros(zoom ** 2, dtype=numpy.float32)
        temp_diag = numpy.zeros((zoom * (zoom + 1)) / 2, dtype=numpy.float32)
        weave.inline(code, ['data', 'new_data', 'new_temp_indices', 'prev_data', 'prev_indices', 'indices', 'n',
                            'temp_data', 'temp_diag', 'temp_indices', 'pos', 'zoom', 't_bins', 'n_bins', 'old_n_bins',
                            'nan', 'minobs'])

        where = numpy.where(new_temp_indices >= 0)[0]
        pos += where.shape[0]
        temp_indices = new_temp_indices[where]
        all_data.append(new_data[where])
        old_n_bins = n_bins
        if h < n_res - 1:
            all_indices.append(numpy.zeros(all_data[-1].shape[0], dtype=numpy.int32))
            all_indices[-1].fill(-1)
    data_bins = 0
    for i in range(len(all_data)):
        data_bins += all_data[i].shape[0]
    total_bins = data_bins
    for i in range(len(all_indices)):
        total_bins += all_indices[i].shape[0]
    # Header contains:
    # magic number - 42054205 (8 btyes)
    # Top resolution - int32
    # Bottom resolution - int32
    # Zoom factor - int32
    # Minimum observation cutoff - int32
    # First bin coordinate - int32
    # Top level number of bins - int32
    # Total data bins (indices offset) - int32
    # Total bins (data + indices) - int32
    output = open(args.output, 'wb')
    output.write(struct.pack('iiiiiiiii', int('42054205', 16), args.lres, args.hres, args.zoom, args.minobs,
                             mapping_start, all_data[0].shape[0], data_bins, total_bins))
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
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive HiC project to construct the heatmap from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write the multi-resolution HiFive heatmap to.")
    return parser


if __name__ == "__main__":
    main()

