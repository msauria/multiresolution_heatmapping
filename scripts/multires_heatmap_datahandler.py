#!/usr/bin/env python

import sys
import struct
import math


def main():
    start1, start2, width, maxres, fname = sys.argv[1:6]
    start1, start2, width, maxres = int(start1), int(start2), int(width), int(maxres)
    if start1 <= start2:
        transpose = False
    else:
        transpose = True
        start1, start2 = start2, start1
    args = {'start1':start1, 'start2':start2, 'stop1':start1 + width, 'stop2':start2 + width,
            'width':width, 'maxres':maxres, 'transpose':transpose}
    if args['stop1'] > args['start2']:
        args['overlap'] = True
    else:
        args['overlap'] = False
    temp = load_file(fname)
    if temp is None:
        return None
    header, infile = temp
    outdata = paint_canvas(header, args, infile)
    print_squares(outdata, args)
    infile.close()

def paint_canvas(header, args, infile):
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
                if indices[k] != -1:
                    valid, new_outdata = paint_lower_level(header, args, infile, indices[k],
                                                           resolution / header['zoom'], start1, start2)
                    if valid < header['zoom2']:
                        outdata.append([start1, start2, resolution, data[k]])
                        if (args['overlap'] and start1 != start2 and start2 + resolution > args['start1'] and
                            start1 < args['stop2']):
                            outdata.append([start2, start1, resolution, data[k]])
                    outdata += new_outdata
                else:
                    outdata.append([start1, start2, resolution, data[k]])
    return outdata

def paint_lower_level(header, args, infile, index, resolution, start1, start2):
    if resolution < args['maxres']:
        return 0, []
    outdata = []
    valid = 0
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
                    if (start1b > args['stop1'] or start1b + resolution < args['start1'] or
                        start2b > args['stop2'] or start2b + resolution < args['start2']):
                        valid += 1
                    else:
                        if not indices is None and indices[k] != -1:
                            new_valid, new_outdata = paint_lower_level(header, args, infile, indices[k],
                                                     resolution / header['zoom'], start1b, start2b)
                            if start1b == start2b:
                                if new_valid < header['zoom'] * (header['zoom'] + 1) / 2:
                                    outdata.append([start1b, start2b, resolution, data[k]])
                            else:
                                if new_valid < header['zoom2']:
                                    outdata.append([start1b, start2b, resolution, data[k]])
                                    if (args['overlap'] and start2b + resolution > args['start1'] and
                                        start1b < args['stop2']):
                                        outdata.append([start2b, start1b, resolution, data[k]])
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
                if (start1b > args['stop1'] or start1b + resolution < args['start1'] or
                    start2b > args['stop2'] or start2b + resolution < args['start2']):
                    valid += 1
                else:
                    if not indices is None and indices[i] != -1:
                        new_valid, new_outdata = paint_lower_level(header, args, infile, indices[i],
                                                 resolution / header['zoom'], start1b, start2b)
                        if new_valid < header['zoom2']:
                            outdata.append([start1b, start2b, resolution, data[i]])
                            if args['overlap'] and start2b + resolution > args['start1'] and start1b < args['stop2']:
                                outdata.append([start2b, start1b, resolution, data[i]])
                        outdata += new_outdata
                        valid += 1
    return valid, outdata

def print_squares(data, args):
    width = float(args['width'])
    print """{\n\t"squares":["""
    for square in data:
        if args['transpose']:
            print """\t\t{"x1":%f,""" % min(max((square[1] - args['start2']) / width, 0.0), 1.0)
            print """\t\t"x2":%f,""" % min(max((square[1] - args['start2'] + square[2]) / width, 0.0), 1.0)
            print """\t\t"y1":%f,""" % min(max((square[0] - args['start1']) / width, 0.0), 1.0)
            print """\t\t"y2":%f,""" % min(max((square[0] - args['start1'] + square[2]) / width, 0.0), 1.0)
        else:
            print """\t\t{"x1":%f,""" % min(max((square[0] - args['start1']) / width, 0.0), 1.0)
            print """\t\t"x2":%f,""" % min(max((square[0] - args['start1'] + square[2]) / width, 0.0), 1.0)
            print """\t\t"y1":%f,""" % min(max((square[1] - args['start2']) / width, 0.0), 1.0)
            print """\t\t"y2":%f,""" % min(max((square[1] - args['start2'] + square[2]) / width, 0.0), 1.0)
        print """\t\t"value":%f},""" % square[3]
    print """\t]\n}"""
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
    header['n_levels'] = int(round(math.log(lres / float(hres)) / math.log(zoom))) + 1
    header['n'] = int((0.25 + 2 * n_bins) ** 0.5 - 0.5)
    header['i_bins'] = header['t_bins'] - header['d_bins']
    header['zoom2'] = header['zoom'] ** 2
    return header, infile


if __name__ == "__main__":
    main()
