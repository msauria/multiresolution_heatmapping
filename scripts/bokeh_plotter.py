#!/usr/bin/env python

from math import ceil
import argparse

from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, output_file, show
import numpy

from complete_multires_heatmap_datahandler import MrhSquareDataProvider

def main():
    parser = generate_parser()
    args = parser.parse_args()
    output_file( args.outfile )
    dataprovider = MrhSquareDataProvider( source=args.datafile )
    header = dataprovider.get_data( header=True, chrom1=args.chrom1, chrom2=args.chrom2 )
    if args.minscore is None:
        args.minscore = header['minscore']
    if args.maxscore is None:
        args.maxscore = header['maxscore']

    temp_data = dataprovider.get_data(
            chrom1=args.chrom1, start1=args.start1, stop1=args.end1,
            chrom2=args.chrom2, start2=args.start2, stop2=args.end2,
            min_resolution=args.minres, max_resolution=args.maxres,
            minscore=args.minscore, maxscore=args.maxscore
        )
    TOOLS = "pan, wheel_zoom, reset"
    p1 = figure( title="Interaction Data", background_fill="#555555", tools=TOOLS )
    canvas = plot_data( temp_data, args )
    p1.image_rgba( image=[canvas], x=[args.start1], y=[args.start2],
              dw=[args.end1 - args.start1],
              dh=[args.end2 - args.start2] )

    show(p1)

def plot_data( data, args ):
    width1 = float( args.end1 - args.start1 )
    width2 = float( args.end2 - args.start2 )
    width = width1 / args.maxres
    height = int( ceil( width2 / width1 * width ) )
    width1 /= width
    width2 /= height
    canvas = numpy.empty( ( width, height ), dtype=numpy.uint32 )
    canvas.fill( int( "ff555555", 16 ) )
    for i in range( len( data[ 'x1' ] ) ):
        x1 = int( round( ( data[ 'x1' ][ i ] - args.start1 ) / width1 ) )
        x2 = int( round( ( data[ 'x2' ][ i ] - args.start1 ) / width1 ) )
        y1 = int( round( ( data[ 'y1' ][ i ] - args.start2 ) / width2 ) )
        y2 = int( round( ( data[ 'y2' ][ i ] - args.start2 ) / width2 ) )
        canvas[ x1:x2, y1:y2 ] = data[ 'color' ][i]
    return canvas

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Produce a multi-resolution heatmap file for a chromosome from HiC data using HiFive"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-L", "--lowest-resolution", dest="minres", required=False, type=int, default=10000000,
        action='store', help="The lowest resolution bin size to create. [default: %(default)s]")
    parser.add_argument("-H", "--highest-resolution", dest="maxres", required=False, type=int, default=1000000,
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
    parser.add_argument("-m", "--minscore", dest="minscore", required=False, type=int, default=None,
        action='store', help="Minimum score cutoff for plotting. [default: %(default)s]")
    parser.add_argument("-M", "--maxscore", dest="maxscore", required=False, type=int, default=None,
        action='store', help="Maximum score cutoff for plotting. [default: %(default)s]")
    parser.add_argument(dest="datafile", type=str,
        action='store', help="Multi-resolution heatmap file to use.")
    parser.add_argument(dest="outfile", type=str, default='temp.html',
        action='store', help="File name to write plot html to [default: %(default)s].")
    return parser

if __name__ == "__main__":
        main()
