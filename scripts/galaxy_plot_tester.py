import argparse as ap

import urllib2
from PIL import Image
import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    args.span1 = args.end1 - args.start1
    args.span2 = args.end2 - args.start2
    temp = urllib2.urlopen("http://localhost:8080/api/datasets/%s?data_type=raw_data&provider=json&chrom1=%s&chrom2=%s&start1=%i&start2=%i&stop1=%i&stop2=%i&min_resolution=%i&max_resolution=%s&chromosomes=true" %
        (args.dataset, args.chrom1, args.chrom2, args.start1, args.start2, args.end1, args.end2, args.lres, args.hres)).read()
    chroms = eval(temp.replace('false', 'False').replace('true', 'True'))['data'][0]
    print chroms
    temp = urllib2.urlopen("http://localhost:8080/api/datasets/%s?data_type=raw_data&provider=json&chrom1=%s&chrom2=%s&start1=%i&start2=%i&stop1=%i&stop2=%i&min_resolution=%i&max_resolution=%s&header=true" %
        (args.dataset, args.chrom1, args.chrom2, args.start1, args.start2, args.end1, args.end2, args.lres, args.hres)).read()
    header = eval(temp.replace('false', 'False').replace('true', 'True'))['data'][0]
    print header
    header['score_span'] = header['maxscore'] - header['minscore']
    temp = urllib2.urlopen("http://localhost:8080/api/datasets/%s?data_type=raw_data&provider=json&chrom1=%s&chrom2=%s&start1=%i&start2=%i&stop1=%i&stop2=%i&min_resolution=%i&max_resolution=%i" %
        (args.dataset, args.chrom1, args.chrom2, args.start1, args.start2, args.end1, args.end2, args.lres, args.hres)).read()
    data = eval(temp)['data']
    plot_data(data, args, header)
    return None

def plot_data(data, args, header):
    canvas = numpy.zeros((args.canvas, args.canvas * args.span2 / args.span1), dtype=numpy.float32)
    canvas.fill(-10)
    for square in data:
        start1 = int(round((square['x1'] - args.start1) / float(args.span1) * args.canvas))
        stop1 = int(round((square['x2'] - args.start1) / float(args.span1) * args.canvas))
        start2 = int(round((square['y1'] - args.start2) / float(args.span1) * args.canvas))
        stop2 = int(round((square['y2'] - args.start2) / float(args.span1) * args.canvas))
        value = (square['value'] - header['minscore']) / header['score_span'] * 2.0 - 1.0
        canvas[start1:stop1, start2:stop2] = value
    canvas = (canvas * 255)[::-1, :].astype(numpy.int32)
    img = numpy.zeros(canvas.shape, dtype=numpy.uint32)
    img.shape = (img.shape[1], img.shape[0])
    img.fill(int("ff555555", 16))
    where = numpy.where(canvas >= 0)
    img[where[1], where[0]] = (256 ** 3 + 1) * 255 + (255 - canvas[where]) * (256 ** 2 + 256)
    where = numpy.where((canvas < 0) * (canvas > -10))
    img[where[1], where[0]] = (256 ** 3 + 256 ** 2) * 255 + (255 + canvas[where]) * (256 + 1)
    pilImage = Image.frombuffer('RGBA', canvas.shape, img, 'raw', 'RGBA', 0, 1)
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
    parser.add_argument("-c", "--chrom1", dest="chrom1", required=False, type=str, default='19',
        action='store', help="First chromosome. [default: %(default)s]")
    parser.add_argument("-C", "--chrom2", dest="chrom2", required=False, type=str, default='19',
        action='store', help="Second chromosome. [default: %(default)s]")
    parser.add_argument("-s", "--start1", dest="start1", required=False, type=int, default=20000000,
        action='store', help="First start position. [default: %(default)s]")
    parser.add_argument("-S", "--start2", dest="start2", required=False, type=int, default=20000000,
        action='store', help="Second start position. [default: %(default)s]")
    parser.add_argument("-e", "--end1", dest="end1", required=False, type=int, default=30000000,
        action='store', help="First stop position. [default: %(default)s]")
    parser.add_argument("-E", "--end2", dest="end2", required=False, type=int, default=30000000,
        action='store', help="Second stop position. [default: %(default)s]")
    parser.add_argument("-w", "--canvas", dest="canvas", required=False, type=int, default=500,
        action='store', help="Canvas width in pixels. [default: %(default)s]")
    parser.add_argument("-d", "--dataset", dest="dataset", required=False, type=str, default='c9468fdb6dc5c5f1',
        action='store', help="Galaxy dataset to use. [default: %(default)s]")
    parser.add_argument("-t", "--trans", dest="trans", required=False,
        action='store_true', help="Dataset is trans. [default: %(default)s]")
    parser.add_argument("-o", "--outfile", dest="outfile", required=False, type=str, default='test.png',
        action='store', help="File name to write plot to [default: %(default)s].")
    return parser

if __name__ == "__main__":
    main()