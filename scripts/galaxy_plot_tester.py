import argparse as ap
from math import ceil, floor

import urllib2
from PIL import Image
import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    if args.trans:
        trans = '&trans=True'
    else:
        trans = ''
    temp = urllib2.urlopen("http://localhost:8080/api/datasets/%s?data_type=raw_data&provider=json&start1=%i&start2=%i&window_size=%i&min_resolution=%i&max_resolution=%s&header=true%s" %
        (args.dataset, args.start1, args.start2, args.window, args.lres, args.hres, trans)).read()
    header = eval(temp)['data'][0]
    header['score_span'] = header['maxscore'] - header['minscore']
    temp = urllib2.urlopen("http://localhost:8080/api/datasets/%s?data_type=raw_data&provider=json&start1=%i&start2=%i&window_size=%i&min_resolution=%i&max_resolution=%i%s" %
        (args.dataset, args.start1, args.start2, args.window, args.lres, args.hres, trans)).read()
    data = eval(temp)['data']
    plot_data(data, args, header)
    return None

def plot_data(data, args, header):
    canvas = numpy.zeros((args.canvas, args.canvas), dtype=numpy.float32)
    for square in data:
        start1 = int(round((square['x1'] - args.start1) / float(args.window) * args.canvas))
        #cstart1 = int(ceil(start1))
        #fstart1 = int(floor(start1))
        stop1 = int(round((square['x2'] - args.start1) / float(args.window) * args.canvas))
        #fstop1 = int(floor(stop1))
        start2 = int(round((square['y1'] - args.start2) / float(args.window) * args.canvas))
        #cstart2 = int(ceil(start2))
        #fstart2 = int(floor(start2))
        stop2 = int(round((square['y2'] - args.start2) / float(args.window) * args.canvas))
        #fstop2 = int(floor(stop2))
        value = (square['value'] - header['minscore']) / header['score_span'] * 2.0 - 1.0
        canvas[start1:stop1, start2:stop2] = value
        """
        if start2 > 0.0:
            frac = cstart2 - start2
            canvas[cstart1:fstop1, fstart2] += value * frac
            if start1 > 0.0:
                canvas[fstart1, fstart2] += value * frac * (cstart1 - start1)
            if stop1 < args.canvas:
                canvas[fstop1, fstart2] += value * frac * (stop1 - fstop1)
        if stop2 < args.canvas:
            frac = stop2 - fstop2
            canvas[cstart1:fstop1, fstop2] += value * frac
            if start1 > 0.0:
                canvas[fstart1, fstop2] += value * frac * (cstart1 - start1)
            if stop1 < args.canvas:
                canvas[fstop1, fstop2] += value * frac * (stop1 - fstop1)
        if start1 > 0.0:
            frac = cstart1 - start1
            canvas[fstart1, cstart2:fstop2] += value * frac
        if stop1 < args.canvas:
            frac = stop1 - fstop1
            canvas[fstop1, cstart2:fstop2] += value * frac
        """
    canvas = (canvas * 255)[::-1, :].astype(numpy.int32)
    img = numpy.zeros(canvas.shape, dtype=numpy.uint32)
    where = numpy.where(canvas >= 0)
    img[where[1], where[0]] = (256 ** 3 + 1) * 255 + (255 - canvas[where]) * (256 ** 2 + 256)
    where = numpy.where(canvas < 0)
    img[where[1], where[0]] = (256 ** 3 + 256 ** 2) * 255 + (255 + canvas[where]) * (256 + 1)
    pilImage = Image.frombuffer('RGBA', img.shape, img, 'raw', 'RGBA', 0, 1)
    pilImage.save(args.outfile)
    return None

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Produce a multi-resolution heatmap file for a chromosome from HiC data using HiFive"
    parser = ap.ArgumentParser(description=description)
    parser.add_argument("-L", "--lowest-resolution", dest="lres", required=False, type=int, default=10000000,
        action='store', help="The lowest resolution bin size to create. [default: %(default)s]")
    parser.add_argument("-H", "--highest-resolution", dest="hres", required=False, type=int, default=1000000,
        action='store', help="The highest resolution bin size to create. [default: %(default)s]")
    parser.add_argument("-s", "--start1", dest="start1", required=False, type=int, default=20000000,
        action='store', help="First start position. [default: %(default)s]")
    parser.add_argument("-S", "--start2", dest="start2", required=False, type=int, default=20000000,
        action='store', help="Second start position. [default: %(default)s]")
    parser.add_argument("-w", "--window", dest="window", required=False, type=int, default=1000000,
        action='store', help="Window size of region to get. [default: %(default)s]")
    parser.add_argument("-c", "--canvas", dest="canvas", required=False, type=int, default=500,
        action='store', help="Canvas size in pixels. [default: %(default)s]")
    parser.add_argument("-d", "--dataset", dest="dataset", required=False, type=str, default='df7a1f0c02a5b08e',
        action='store', help="Galaxy dataset to use. [default: %(default)s]")
    parser.add_argument("-t", "--trans", dest="trans", required=False,
        action='store_true', help="Dataset is trans. [default: %(default)s]")
    parser.add_argument("-o", "--outfile", dest="outfile", required=False, type=str, default='test.png',
        action='store', help="File name to write plot to [default: %(default)s].")
    return parser

if __name__ == "__main__":
    main()