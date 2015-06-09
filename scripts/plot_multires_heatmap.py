#!/usr/bin/env python

import sys
import argparse as ap
import struct

import numpy
from PIL import Image

pixels = 400

def main():
    parser = generate_parser()
    args = parser.parse_args()
    temp = load_file(args.heatmap)
    if temp is None:
        return None
    header, infile = temp
    if args.start is None:
        args.start = header['start']
    if args.start2 is None:
        args.start2 = args.start
    if args.width is None:
        args.width = header['start'] + header['n'] * header['lres'] - min(args.start2, args.start)
    if args.maxres is None:
        args.maxres = args.width / float(pixels * 0.5)
    canvas = numpy.zeros((pixels, pixels, 2), dtype=numpy.float32)
    args.bpp = args.width / float(pixels)
    if args.start <= args.start2:
        paint_canvas(canvas, header, args, infile)
    else:
        args.start, args.start2 = args.start2, args.start
        print args.start, args.start2
        paint_canvas(canvas, header, args, infile)
        canvas = canvas.transpose(1, 0, 2)
    infile.close()
    plot_canvas(canvas, args.output)

def plot_canvas(data, outfname):
    valid = numpy.where(data[:, :, 1] > 0.0)
    data[valid[0], valid[1], 0] -= numpy.amin(data[valid[0], valid[1], 0])
    data[valid[0], valid[1], 0] /= numpy.amax(data[valid[0], valid[1], 0]) * 0.5
    data[valid[0], valid[1], 0] -= 1.0
    img = numpy.zeros((data.shape[0], data.shape[0]), dtype=numpy.uint32)
    img.fill(int('ff999999', 16))
    where = numpy.where((data[:, :, 1] > 0.0) * (data[:, :, 0] >= 0.0))
    img[where] = ((256 ** 3 + 1) * 255 + (256 ** 2 + 256) *
                  (255 - numpy.round(255.0 * data[where[0], where[1], 0]).astype(numpy.int32)))
    where = numpy.where((data[:, :, 1] > 0.0) * (data[:, :, 0] < 0.0))
    img[where] = ((256 ** 3 + 256 ** 2) * 255 + (257) *
                  (255 + numpy.round(255.0 * data[where[0], where[1], 0]).astype(numpy.int32)))
    img = Image.frombuffer('RGBA', img.shape, img, 'raw', 'RGBA', 0, 1)
    img.save(outfname)

def paint_canvas(canvas, header, args, infile):
    start_pos1 = max(0, (args.start - header['start']) / header['lres'])
    end_pos1 = min(header['n'], (args.start + args.width - header['start']) / header['lres'] + 1)
    start_pos2 = max(0, (args.start2 - header['start']) / header['lres'])
    end_pos2 = min(header['n'], (args.start2 + args.width - header['start']) / header['lres'] + 1)
    resolution = header['lres']
    for i in range(start_pos1, end_pos1):
        # Find position in file for data with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * (header['n'] - 1) - (i * (i - 1)) / 2 + max(i, start_pos2)) * 4)
        data = numpy.fromstring(infile.read((end_pos2 - max(start_pos2, i)) * 4), dtype=numpy.float32)
        # Find position in file for indices with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * (header['n'] - 1) - (i * (i - 1)) / 2 + header['d_bins'] +
                    max(i, start_pos2)) * 4)
        indices = numpy.fromstring(infile.read((end_pos2 - max(start_pos2, i)) * 4), dtype=numpy.int32)
        for j in range(max(start_pos2, i), end_pos2):
            k = j - max(start_pos2, i)
            if not numpy.isnan(data[k]):
                start1 = i * resolution + header['start']
                start2 = j * resolution + header['start']
                paint_square(canvas, header, args, data[k], resolution, start1, start2)
                if indices[k] != -1:
                    paint_lower_level(canvas, header, args, infile, indices[k],
                                      resolution / header['zoom'], start1, start2)
    if args.start2 < args.start + args.width:
        span = int(round((args.start + args.width - args.start2) / args.bpp))
        for i in range(span):
            canvas[(canvas.shape[0] - span + i):, i, :] = canvas[canvas.shape[0] - span + i, i:span, :]
    return None

def paint_lower_level(canvas, header, args, infile, index, resolution, start1, start2):
    if resolution < args.maxres:
        return None
    infile.seek(header['offset'] + index * 4)
    if start1 == start2:
        data = numpy.fromstring(infile.read(2 * header['zoom'] * (header['zoom'] + 1)), dtype=numpy.float32)
        if index < header['i_bins']:
            infile.seek(header['offset'] + (index + header['d_bins']) * 4)
            indices = numpy.fromstring(infile.read(2 * header['zoom'] * (header['zoom'] + 1)), dtype=numpy.int32)
        else:
            indices = None
        for i in range(header['zoom']):
            for j in range(i, header['zoom']):
                k = i * (header['zoom'] - 1) - (i * (i - 1)) / 2 + j
                if not numpy.isnan(data[k]):
                    start1b = start1 + i * resolution
                    start2b = start2 + j * resolution
                    paint_square(canvas, header, args, data[k], resolution, start1b, start2b)
                    if not indices is None and indices[k] != -1:
                        paint_lower_level(canvas, header, args, infile, indices[k], resolution / header['zoom'],
                                          start1b, start2b)
    else:
        data = numpy.fromstring(infile.read(4 * header['zoom'] ** 2), dtype=numpy.float32)
        if index < header['i_bins']:
            infile.seek(header['offset'] + (index + header['d_bins']) * 4)
            indices = numpy.fromstring(infile.read(4 * header['zoom'] ** 2), dtype=numpy.int32)
        else:
            indices = None
        for i in range(header['zoom'] ** 2):
            if not numpy.isnan(data[i]):
                start1b = start1 + (i / header['zoom']) * resolution
                start2b = start2 + (i % header['zoom']) * resolution
                paint_square(canvas, header, args, data[i], resolution, start1b, start2b)
                if not indices is None and indices[i] != -1:
                    paint_lower_level(canvas, header, args, infile, indices[i], resolution / header['zoom'],
                                      start1b, start2b)
    return None

def paint_square(canvas, header, args, data, resolution, start1, start2):
    n = canvas.shape[0]
    index1a = max(0, int(round((start1 - args.start) / args.bpp)))
    index2a = max(0, int(round((start2 - args.start2) / args.bpp)))
    index1b = min(n, int(round((start1 + resolution - args.start) / args.bpp)))
    index2b = min(n, int(round((start2 + resolution - args.start2) / args.bpp)))
    if index1b - index1a > 0 and index2b - index2a > 0:
        canvas[index1a:index1b, index2a:index2b, 0] = data
        canvas[index1a:index1b, index2a:index2b, 1] = 1
    return None

def load_file(fname):
    infile = open(fname, 'rb')
    magic_number = hex(struct.unpack('i', infile.read(4))[0])
    if magic_number != '0x42054205':
        print >> sys.stderr, ('File does not appear to be a multi-resolution heatmap file.\n'),
        return None
    lres, hres, zoom, minobs, start, n_bins, d_bins, t_bins = struct.unpack('iiiiiiii', infile.read(32))
    header = {'lres':lres, 'hres':hres, 'zoom':zoom, 'start':start, 'n_bins':n_bins, 'd_bins':d_bins,
              't_bins':t_bins, 'offset':36}
    header['n_levels'] = int(numpy.round(numpy.log(lres / float(hres)) / numpy.log(zoom))) + 1
    header['n'] = int((0.25 + 2 * n_bins) ** 0.5 - 0.5)
    header['i_bins'] = header['t_bins'] - header['d_bins']
    return header, infile

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Produce a multi-resolution heatmap image from a HiC multi-resolution heatmap file."
    parser = ap.ArgumentParser(description=description)
    parser.add_argument("-s", "--start", dest="start", required=False, type=int, default=None,
        action='store', help="The first region start coordinate to plot. If no value is passed, this will be set to the first bin position in the heatmap. [default: %(default)s]")
    parser.add_argument("-S", "--start2", dest="start2", required=False, type=int, default=None,
        action='store', help="The second region start coordinate to plot. If no value is passed, this will be set to the same value as start. [default: %(default)s]")
    parser.add_argument("-w", "--width", dest="width", required=False, type=int, default=None,
        action='store', help="The width of the regions to plot. If no value is set, this will be set to the distance between min(start, start2) and the end of the last bin in the heatmap. [default: %(default)s]")
    parser.add_argument("-m", "--max-resolution", dest="maxres", required=False, type=int, default=None,
        action='store', help="A maximum resolution bound for plotting. [default: %(default)s]")
    parser.add_argument(dest="heatmap", type=str,
        help="The name of a HiFive multi-resolution heatmap file to construct the image from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write the multi-resolution HiFive heatmap image to.")
    return parser

if __name__ == "__main__":
    main()
