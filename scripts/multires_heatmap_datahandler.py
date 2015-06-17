#!/usr/bin/env python

import sys
import struct
import math
import binascii

from pyx import *

def main():
    start1, start2, window, maxres, minres, fname = sys.argv[1:7]
    start1, start2, window, maxres, minres = int(start1), int(start2), int(window), int(maxres), int(minres)
    temp = load_file(fname)
    if temp is None:
        return None
    header, infile = temp
    print header
    if header['trans']:
        if sys.argv[7] == '0':
            transpose = False
        else:
            transpose = True
        args = {'start1':start1, 'start2':start2, 'stop1':start1 + window, 'stop2':start2 + window,
                'window':window, 'minres':minres, 'maxres':maxres, 'transpose':transpose}
        outdata = paint_trans_canvas(header, args, infile)
    else:
        if start1 <= start2:
            transpose = False
        else:
            transpose = True
            start1, start2 = start2, start1
        args = {'start1':start1, 'start2':start2, 'stop1':start1 + window, 'stop2':start2 + window,
                'window':window, 'minres':minres, 'maxres':maxres, 'transpose':transpose}
        if args['stop1'] > args['start2']:
            args['overlap'] = True
        else:
            args['overlap'] = False
        outdata = paint_cis_canvas(header, args, infile)
    #print_squares(outdata, args)
    plot_squares(outdata, args, header, sys.argv[8])
    infile.close()

def paint_cis_canvas(header, args, infile):
    outdata = []
    start_pos1 = max(0, (args['start1'] - header['start']) / header['lres'])
    end_pos1 = min(header['n'], (args['stop1'] - header['start']) / header['lres'] + 1)
    start_pos2 = max(0, (args['start2'] - header['start']) / header['lres'])
    end_pos2 = min(header['n'], (args['stop2'] - header['start']) / header['lres'] + 1)
    resolution = header['lres']
    for i in range(start_pos1, end_pos1):
        # Find position in file for data with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * (header['n'] - 1) - (i * (i - 1)) / 2 + max(i, start_pos2)) * 4)
        data = struct.unpack('f' * (end_pos2 - max(start_pos2, i)), infile.read((end_pos2 - max(start_pos2, i)) * 4))
        # Find position in file for indices with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * (header['n'] - 1) - (i * (i - 1)) / 2 + header['d_bins'] +
                    max(i, start_pos2)) * 4)
        indices = struct.unpack('i' * (end_pos2 - max(start_pos2, i)),
                                infile.read((end_pos2 - max(start_pos2, i)) * 4))
        for j in range(max(start_pos2, i), end_pos2):
            k = j - max(start_pos2, i)
            if not math.isnan(data[k]):
                start1 = i * resolution + header['start']
                start2 = j * resolution + header['start']
                if (start1 >= args['stop1'] or start1 + resolution <= args['start1'] or
                    start2 >= args['stop2'] or start2 + resolution <= args['start2']):
                    continue
                else:
                    if indices[k] != -1:
                        new_outdata, valid = paint_cis_lower_level(header, args, infile, indices[k],
                                                                   resolution / header['zoom'], start1, start2)
                    else:
                        valid = 0
                    if start1 == start2:
                        if resolution <= args['minres'] and valid < header['zoom_d']:
                            outdata.append([start1, start2, resolution, data[k]])
                            if (args['overlap'] and start1 != start2 and start2 + resolution > args['start1'] and
                                start1 < args['stop2']):
                                outdata.append([start2, start1, resolution, data[k]])
                    else:
                        if resolution <= args['minres'] and valid < header['zoom2']:
                            outdata.append([start1, start2, resolution, data[k]])
                            if (args['overlap'] and start1 != start2 and start2 + resolution > args['start1'] and
                                start1 < args['stop2']):
                                outdata.append([start2, start1, resolution, data[k]])
                    if valid > 0:
                        outdata += new_outdata
    return outdata

def paint_trans_canvas(header, args, infile):
    outdata = []
    start_pos1 = max(0, (args['start1'] - header['start']) / header['lres'])
    end_pos1 = min(header['n'], (args['stop1'] - header['start']) / header['lres'] + 1)
    start_pos2 = max(0, (args['start2'] - header['start2']) / header['lres'])
    end_pos2 = min(header['m'], (args['stop2'] - header['start2']) / header['lres'] + 1)
    resolution = header['lres']
    for i in range(start_pos1, end_pos1):
        # Find position in file for data with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * header['m'] + start_pos2) * 4)
        data = struct.unpack('f' * (end_pos2 - start_pos2), infile.read((end_pos2 - start_pos2) * 4))
        # Find position in file for indices with 'i' as upstream interaction
        infile.seek(header['offset'] + (i * header['m'] + start_pos2 + header['d_bins']) * 4)
        indices = struct.unpack('i' * (end_pos2 - start_pos2), infile.read((end_pos2 - start_pos2) * 4))
        for j in range(start_pos2, end_pos2):
            k = j - start_pos2
            if not math.isnan(data[k]):
                start1 = i * resolution + header['start']
                start2 = j * resolution + header['start2']
                if (start1 >= args['stop1'] or start1 + resolution <= args['start1'] or
                    start2 >= args['stop2'] or start2 + resolution <= args['start2']):
                    continue
                else:
                    if indices[k] != -1:
                        new_outdata, valid = paint_trans_lower_level(header, args, infile, indices[k],
                                                                     resolution / header['zoom'], start1, start2)
                    else:
                        valid = 0
                if resolution <= args['minres'] and valid < header['zoom2']:
                    outdata.append([start1, start2, resolution, data[k]])
                if valid > 0:
                    outdata += new_outdata
    return outdata

def paint_cis_lower_level(header, args, infile, index, resolution, start1, start2):
    if resolution < args['maxres']:
        return [], 0
    valid = 0
    outdata = []
    infile.seek(header['offset'] + index * 4)
    if start1 == start2:
        data = struct.unpack('f' * (header['zoom'] * (header['zoom'] + 1) / 2),
                             infile.read(2 * header['zoom'] * (header['zoom'] + 1)))
        if index < header['i_bins']:
            infile.seek(header['offset'] + (index + header['d_bins']) * 4)
            indices = struct.unpack('i' * (header['zoom'] * (header['zoom'] + 1) / 2),
                                    infile.read(header['zoom'] * (header['zoom'] + 1) * 2))
        else:
            indices = None
        for i in range(header['zoom']):
            for j in range(i, header['zoom']):
                k = i * (header['zoom'] - 1) - (i * (i - 1)) / 2 + j
                if not math.isnan(data[k]):
                    start1b = start1 + i * resolution
                    start2b = start2 + j * resolution
                    if (start1b >= args['stop1'] or start1b + resolution <= args['start1'] or
                        start2b >= args['stop2'] or start2b + resolution <= args['start2']):
                        valid += 1
                        continue
                    else:
                        if not indices is None and indices[k] != -1:
                            new_outdata, new_valid = paint_cis_lower_level(header, args, infile, indices[k],
                                                                           resolution / header['zoom'], start1b,
                                                                           start2b)
                        else:
                            new_valid = 0
                        if resolution <= args['minres']:
                            if start1b == start2b:
                                if new_valid < header['zoom_d']:
                                    outdata.append([start1b, start2b, resolution, data[k]])
                            else:
                                if new_valid < header['zoom2']:
                                    outdata.append([start1b, start2b, resolution, data[k]])
                                    if (start1b != start2b and
                                        args['overlap'] and start2b + resolution > args['start1'] and
                                        start1b < args['stop2']):
                                        outdata.append([start2b, start1b, resolution, data[k]])
                        if new_valid > 0:
                            outdata += new_outdata
                        valid += 1
    else:
        data = struct.unpack('f' * header['zoom2'], infile.read(4 * header['zoom2']))
        if index < header['i_bins']:
            infile.seek(header['offset'] + (index + header['d_bins']) * 4)
            indices = struct.unpack('i' * header['zoom2'], infile.read(header['zoom2'] * 4))
        else:
            indices = None
        for i in range(header['zoom2']):
            if not math.isnan(data[i]):
                start1b = start1 + (i / header['zoom']) * resolution
                start2b = start2 + (i % header['zoom']) * resolution
                if (start1b >= args['stop1'] or start1b + resolution <= args['start1'] or
                    start2b >= args['stop2'] or start2b + resolution <= args['start2']):
                    valid += 1
                    continue
                else:
                    if not indices is None and indices[i] != -1:
                        new_outdata, new_valid = paint_cis_lower_level(header, args, infile, indices[i],
                                                                  resolution / header['zoom'], start1b, start2b)
                    else:
                        new_valid = 0
                    if resolution <= args['minres'] and new_valid < header['zoom2']:
                        outdata.append([start1b, start2b, resolution, data[i]])
                        if (args['overlap'] and start2b + resolution > args['start1'] and
                            start1b < args['stop2']):
                            outdata.append([start2b, start1b, resolution, data[i]])
                    if new_valid > 0:
                        outdata += new_outdata
                    valid += 1
    return outdata, valid

def paint_trans_lower_level(header, args, infile, index, resolution, start1, start2):
    if resolution < args['maxres']:
        return [], 0
    valid = 0
    outdata = []
    infile.seek(header['offset'] + index * 4)
    data = struct.unpack('f' * header['zoom2'], infile.read(4 * header['zoom2']))
    if index < header['i_bins']:
        infile.seek(header['offset'] + (index + header['d_bins']) * 4)
        indices = struct.unpack('i' * header['zoom2'], infile.read(header['zoom2'] * 4))
    else:
        indices = None
    for i in range(header['zoom2']):
        if not math.isnan(data[i]):
            start1b = start1 + (i / header['zoom']) * resolution
            start2b = start2 + (i % header['zoom']) * resolution
            if (start1b > args['stop1'] or start1b + resolution < args['start1'] or
                start2b > args['stop2'] or start2b + resolution < args['start2']):
                valid += 1
                continue
            else:
                if not indices is None and indices[i] != -1:
                    new_outdata, new_valid = paint_trans_lower_level(header, args, infile, indices[i],
                                                                     resolution / header['zoom'], start1b, start2b)
                else:
                    new_valid = 0
                if resolution <= args['minres'] and new_valid < header['zoom2']:
                    outdata.append([start1b, start2b, resolution, data[i]])
                if new_valid > 0:
                    outdata += new_outdata
                valid += 1
    return outdata, valid

def print_squares(data, args):
    print """{\n\t"squares":["""
    for square in data:
        if args['transpose']:
            print """\t\t{"x1":%f,""" % min(max(square[1], args['start2']), args['stop2'])
            print """\t\t"x2":%f,""" % min(max(square[1] + square[2], args['start2']), args['stop2'])
            print """\t\t"y1":%f,""" % min(max(square[0], args['start1']), args['stop1'])
            print """\t\t"y2":%f,""" % min(max(square[0] + square[2], args['start1']), args['stop1'])
        else:
            print """\t\t{"y1":%f,""" % min(max(square[1], args['start2']), args['stop2'])
            print """\t\t"y2":%f,""" % min(max(square[1] + square[2], args['start2']), args['stop2'])
            print """\t\t"x1":%f,""" % min(max(square[0], args['start1']), args['stop1'])
            print """\t\t"x2":%f,""" % min(max(square[0] + square[2], args['start1']), args['stop1'])
        print """\t\t"value":%f},""" % square[3]
    print """\t]\n}"""
    return None

def plot_squares(data, args, header, fname):
    window = float(args['window'])
    c = canvas.canvas()
    g1 = color.gradient.WhiteRed
    g2 = color.gradient.WhiteBlue
    for square in data:
        if args['transpose']:
            y = min(10, max(0, (square[0] - args['start1']) / window * 10))
            x = min(10, max(0, (square[1] - args['start2']) / window * 10))
            y_width = min(10, max(0, (square[0] + square[2] - args['start1']) / window * 10)) - y
            x_width = min(10, max(0, (square[1] + square[2] - args['start2']) / window * 10)) - x
        else:
            x = min(10, max(0, (square[0] - args['start1']) / window * 10))
            y = min(10, max(0, (square[1] - args['start2']) / window * 10))
            x_width = min(10, max(0, (square[0] + square[2] - args['start1']) / window * 10)) - x
            y_width = min(10, max(0, (square[1] + square[2] - args['start2']) / window * 10)) - y
        score = (square[3] - header['minscore']) / (header['maxscore'] - header['minscore']) * 2.0 - 1.0
        if score >= 0.0:
            c.fill(path.rect(x, y, x_width, y_width), [g1.getcolor(score)])
        else:
            c.fill(path.rect(x, y, x_width, y_width), [g2.getcolor(-score)])
    c.writePDFfile(fname)

def load_file(fname):
    try:
        infile = open(fname, 'rb')
    except:
        return None
    magic_number = infile.read(4)
    if magic_number != binascii.a2b_hex('42054205'):
        print >> sys.stderr, ('File does not appear to be a multi-resolution heatmap file.\n'),
        return None
    offset, lres, hres, zoom, minobs, minscore, maxscore, trans = struct.unpack('iiiiiffi', infile.read(32))
    header = {'offset':offset, 'lres':lres, 'hres':hres, 'zoom':zoom,
              'minscore':minscore, 'maxscore':maxscore, 'trans':trans}
    if trans:
        start, start2, stop, stop2, n, m, d_bins, t_bins = struct.unpack('iiiiiiii', infile.read(32))
        header.update({'start':start, 'start2':start2, 'stop':stop, 'stop2':stop2, 'n':n,
                       'm':m, 'd_bins':d_bins, 't_bins':t_bins, 'n_bins':n * m})
    else:
        start, stop, n, d_bins, t_bins = struct.unpack('iiiii', infile.read(20))
        header.update({'start':start, 'stop':stop, 'n':n, 'd_bins':d_bins,
                       't_bins':t_bins, 'n_bins':(n * (n - 1)) / 2})
    header['n_levels'] = int(round(math.log(lres / float(hres)) / math.log(zoom))) + 1
    header['i_bins'] = header['t_bins'] - header['d_bins']
    header['zoom2'] = header['zoom'] ** 2
    header['zoom_d'] = (header['zoom'] * (header['zoom'] + 1)) / 2
    return header, infile


if __name__ == "__main__":
    main()
