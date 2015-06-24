#!/usr/bin/env python

import sys
import struct
import math
import binascii

import numpy
from PIL import Image
from pyx import *

def main():
    chrom1, start1, stop1, chrom2, start2, stop2, maxres, minres, fname = sys.argv[1:10]
    start1, stop1, start2, stop2, maxres, minres = (
        int(start1), int(stop1), int(start2), int(stop2), int(maxres), int(minres))
    try:
        infile = open(fname, 'rb')
    except:
        print >> sys.stderr, ("Couldn't open file.\n"),
        return None
    header = load_header(infile, chrom1, chrom2)
    if header is None:
        return None
    args = {'start1':start1, 'start2':start2, 'stop1':stop1, 'stop2':stop2,
            'window1':stop1 - start1, 'window2':stop2 - start2, 'minres':minres,
            'maxres':maxres, 'transpose':header['transpose']}
    if header['trans']:
        outdata = paint_trans_canvas(header, args, infile)
    else:
        if stop1 <= start2 or stop2 <= start1:
            args['overlap'] = False
        else:
            args['overlap'] = True
        outdata = paint_cis_canvas(header, args, infile)
    #print_squares(outdata, args)
    #plot_squares(outdata, args, header, sys.argv[10])
    plot_data(outdata, args, header, sys.argv[10])
    infile.close()

def paint_cis_canvas(header, args, infile):
    start_pos1 = max(0, (args['start1'] - header['start']) / header['lres'])
    end_pos1 = min(header['n'], (args['stop1'] - header['start']) / header['lres'] + 1)
    start_pos2 = max(0, (args['start2'] - header['start']) / header['lres'])
    end_pos2 = min(header['n'], (args['stop2'] - header['start']) / header['lres'] + 1)
    print start_pos1, end_pos1, start_pos2, end_pos2, args['overlap']
    rev_args = {'start1':args['start2'], 'stop1':args['stop2'], 'start2':args['start1'], 'stop2':args['stop1'],
                'minres':args['minres'], 'maxres':args['maxres'], 'overlap':args['overlap']}
    if args['overlap']:
        outdata = []
        mid_coord = min(args['stop1'], args['stop2'])
        mid_start = min(header['n'], (mid_coord - header['start']) / header['lres'])
        mid_stop = min(header['n'], (mid_coord - header['start']) / header['lres'] + 1)
        if start_pos1 <= start_pos2:
            outdata += paint_cis_toplevel(header, args, infile, start_pos1, mid_stop,
                                          start_pos2, mid_stop, False)
        else:
            outdata += paint_cis_toplevel(header, rev_args, infile, start_pos2, mid_stop,
                                          start_pos1, mid_stop, True)
        if args['stop1'] > mid_coord:
            outdata += paint_cis_toplevel(header, rev_args, infile, start_pos2, end_pos2,
                                          mid_start, end_pos1, True)
        elif args['stop2'] > mid_coord:
            outdata += paint_cis_toplevel(header, args, infile, start_pos1, end_pos1,
                                          mid_start, end_pos2, False)
    elif args['start2'] < args['start1']:
        outdata = paint_cis_toplevel(header, rev_args, infile, start_pos2, end_pos2,
                                      start_pos1, end_pos1, True)
    else:
        outdata = paint_cis_toplevel(header, args, infile, start_pos1, end_pos1,
                                      start_pos2, end_pos2, False)
    return outdata

def paint_cis_toplevel(header, args, infile, start_pos1, end_pos1, start_pos2, end_pos2, transpose):
    resolution = header['lres']
    outdata = []
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
                    else:
                        if resolution <= args['minres'] and valid < header['zoom2']:
                            outdata.append([start1, start2, resolution, data[k]])
                            if (args['overlap'] and start2 + resolution > args['start1'] and
                                start1 < args['stop2']):
                                outdata.append([start2, start1, resolution, data[k]])
                    if valid > 0:
                        outdata += new_outdata
    if transpose:
        for i in range(len(outdata)):
            outdata[i] = [outdata[i][1], outdata[i][0]] + outdata[i][2:]
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
        data = struct.unpack('f' * header['zoom_d'], infile.read(4 * header['zoom_d']))
        if index < header['i_bins']:
            infile.seek(header['offset'] + (index + header['d_bins']) * 4)
            indices = struct.unpack('i' * header['zoom_d'], infile.read(4 * header['zoom_d']))
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
                                    if (args['overlap'] and start2b + resolution > args['start1'] and
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
    width = 10.0
    window = float(args['window1']) / width
    height = width * float(args['window2']) / args['window1']
    c = canvas.canvas()
    g1 = color.gradient.WhiteRed
    g2 = color.gradient.WhiteBlue
    for square in data:
        if args['transpose']:
            y = min(width, max(0, (square[0] - args['start1']) / window))
            x = min(height, max(0, (square[1] - args['start2']) / window))
            y_width = min(width, max(0, (square[0] + square[2] - args['start1']) / window)) - y
            x_width = min(height, max(0, (square[1] + square[2] - args['start2']) / window)) - x
        else:
            x = min(width, max(0, (square[0] - args['start1']) / window))
            y = min(height, max(0, (square[1] - args['start2']) / window))
            x_width = min(width, max(0, (square[0] + square[2] - args['start1']) / window)) - x
            y_width = min(height, max(0, (square[1] + square[2] - args['start2']) / window)) - y
        score = (square[3] - header['minscore']) / (header['maxscore'] - header['minscore']) * 2.0 - 1.0
        if score >= 0.0:
            c.fill(path.rect(x, y, x_width, y_width), [g1.getcolor(score)])
        else:
            c.fill(path.rect(x, y, x_width, y_width), [g2.getcolor(-score)])
    d = canvas.canvas()
    d.insert(c, [trafo.scale(-1, -1)])
    d.writePDFfile(fname)

def plot_data(data, args, header, fname):
    args['canvas'] = 500
    header['score_window'] = header['maxscore'] - header['minscore']
    canvas = numpy.zeros((args['canvas'], args['canvas'] * args['window2'] / args['window1']), dtype=numpy.float32)
    for square in data:
        start1 = int(round((square[0] - args['start1']) / float(args['window1']) * args['canvas']))
        #cstart1 = int(ceil(start1))
        #fstart1 = int(floor(start1))
        stop1 = int(round((square[0] + square[2] - args['start1']) / float(args['window1']) * args['canvas']))
        #fstop1 = int(floor(stop1))
        start2 = int(round((square[1] - args['start2']) / float(args['window1']) * args['canvas']))
        #cstart2 = int(ceil(start2))
        #fstart2 = int(floor(start2))
        stop2 = int(round((square[1] + square[2] - args['start2']) / float(args['window1']) * args['canvas']))
        #fstop2 = int(floor(stop2))
        value = (square[3] - header['minscore']) / header['score_window'] * 2.0 - 1.0
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
    img.shape = (img.shape[1], img.shape[0])
    where = numpy.where(canvas >= 0)
    img[where[1], where[0]] = (256 ** 3 + 1) * 255 + (255 - canvas[where]) * (256 ** 2 + 256)
    where = numpy.where(canvas < 0)
    img[where[1], where[0]] = (256 ** 3 + 256 ** 2) * 255 + (255 + canvas[where]) * (256 + 1)
    pilImage = Image.frombuffer('RGBA', canvas.shape, img, 'raw', 'RGBA', 0, 1)
    pilImage.save(fname)
    return None

def load_header(infile, chrom1, chrom2):
    magic_number = infile.read(4)
    if magic_number != binascii.a2b_hex('42054205'):
        print >> sys.stderr, ('File does not appear to be a multi-resolution heatmap file.\n'),
        return None
    num_chroms, includetrans = struct.unpack('ii', infile.read(8))
    includetrans = bool(includetrans)
    if not includetrans and chrom1 != chrom2:
        return None
    if includetrans:
        num_chrom_pairings = (num_chroms * (num_chroms + 1)) / 2
        repeats = 2
    else:
        num_chrom_pairings = num_chroms
        repeats = 1
    chr2int = {}
    for n in range(num_chroms):
        chr2int[''.join(struct.unpack('c' * 10, infile.read(10))).strip(' ')] = n
    if chrom1 not in chr2int:
        return None
    chrint1 = chr2int[chrom1]
    if chrom1 == chrom2:
        transpose = False
        trans = False
        if includetrans:
            pairing_index = chrint1 * num_chroms - (chrint1 * (chrint1 - 1)) / 2
        else:
            pairing_index = chrint1
        infile.seek(pairing_index * 4, 1)
        index_start, index_stop = struct.unpack('ii', infile.read(8))
        print pairing_index, index_start, index_stop
        if index_start == index_stop:
            return None
        infile.seek((num_chrom_pairings - pairing_index - 1 + chrint1) * 4, 1)
        n = struct.unpack('i', infile.read(4))[0]
        infile.seek((repeats * num_chroms - chrint1 - 1 + pairing_index) * 4, 1)
        d_bins = struct.unpack('i', infile.read(4))[0]
        infile.seek((num_chrom_pairings - pairing_index - 1 + chrint1) * 4, 1)
        start1 = struct.unpack('i', infile.read(4))[0]
        infile.seek((repeats * num_chroms - 1) * 4, 1)
        stop1 = struct.unpack('i', infile.read(4))[0]
        infile.seek((repeats * num_chroms - chrint1 - 1 + pairing_index) * 4, 1)
        minscore = struct.unpack('f', infile.read(4))[0]
        infile.seek((num_chrom_pairings - 1) * 4, 1)
        maxscore = struct.unpack('f', infile.read(4))[0]
        infile.seek((num_chrom_pairings - pairing_index - 1) * 4, 1)
        lres = struct.unpack('i', infile.read(4))[0]
        if includetrans:
            infile.seek(4, 1)
        hres = struct.unpack('i', infile.read(4))[0]
        if includetrans:
            infile.seek(4, 1)
        zoom = struct.unpack('i', infile.read(4))[0]
        if includetrans:
            infile.seek(4, 1)
        n_bins = (n * (n + 1)) / 2
    else:
        if chrom2 not in chr2int:
            return None
        trans = True
        chrint2 = chr2int[chrom2]
        if chrint1 > chrint2:
            chrint2, chrint1 = chrint1, chrint2
            transpose = True
        else:
            transpose = False
        pairing_index = chrint1 * (num_chroms - 1) - (chrint1 * (chrint1 - 1)) / 2 + chrint2
        infile.seek(pairing_index * 4, 1)
        index_start, index_stop = struct.unpack('ii', infile.read(8))
        if index_start == index_stop:
            return None
        infile.seek((num_chrom_pairings - pairing_index - 1 + num_chroms + chrint1) * 4, 1)
        n = struct.unpack('i', infile.read(4))[0]
        infile.seek((chrint2 - chrint1 - 1) * 4, 1)
        m = struct.unpack('i', infile.read(4))[0]
        infile.seek((num_chroms - chrint2 - 1 + pairing_index) * 4, 1)
        d_bins = struct.unpack('i', infile.read(4))[0]
        infile.seek((num_chrom_pairings - pairing_index - 1 + num_chroms + chrint1) * 4, 1)
        start1 = struct.unpack('i', infile.read(4))[0]
        infile.seek((chrint2 - chrint1 - 1) * 4, 1)
        start2 = struct.unpack('i', infile.read(4))[0]
        infile.seek((num_chroms * 2 - chrint2 - 1 + chrint1) * 4, 1)
        stop1 = struct.unpack('i', infile.read(4))[0]
        infile.seek((chrint2 - chrint1 - 1) * 4, 1)
        stop2 = struct.unpack('i', infile.read(4))[0]
        infile.seek((num_chroms - chrint2 - 1 + pairing_index) * 4, 1)
        minscore = struct.unpack('f', infile.read(4))[0]
        infile.seek((num_chrom_pairings - 1) * 4, 1)
        maxscore = struct.unpack('f', infile.read(4))[0]
        infile.seek((num_chrom_pairings - pairing_index) * 4, 1)
        lres = struct.unpack('i', infile.read(4))[0]
        infile.seek(4, 1)
        hres = struct.unpack('i', infile.read(4))[0]
        infile.seek(4, 1)
        zoom = struct.unpack('i', infile.read(4))[0]
        n_bins = n * m
    t_bins = (index_stop - index_start) / 4
    i_bins = t_bins - d_bins
    zoom2 = zoom ** 2
    zoom_d = (zoom * (zoom + 1)) / 2

    print [index_start, lres, hres, zoom, minscore, maxscore, start1, stop1, n, n_bins, d_bins, t_bins,
                       i_bins, zoom2, zoom_d, trans, transpose]
    n_levels = int(round(math.log(lres / float(hres)) / math.log(zoom))) + 1
    header = dict(zip(['offset', 'lres', 'hres', 'zoom', 'minscore', 'maxscore', 'start', 'stop', 'n', 'n_bins',
                       'd_bins', 't_bins', 'i_bins', 'n_levels', 'zoom2', 'zoom_d', 'trans', 'transpose'],
                      [index_start, lres, hres, zoom, minscore, maxscore, start1, stop1, n, n_bins, d_bins, t_bins,
                       i_bins, n_levels, zoom2, zoom_d, trans, transpose]))
    if trans:
        header.update(zip(['m', 'start2', 'stop2'], [m, start2, stop2]))
    return header


if __name__ == "__main__":
    main()
