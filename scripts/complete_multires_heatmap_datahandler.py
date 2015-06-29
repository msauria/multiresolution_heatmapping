#!/usr/bin/env python

import struct
import math
import binascii



# -----------------------------------------------------------------------------
class MrhSquareDataProvider( object ):
    """
    """

    def __init__( self, source ):
        """
        :param source: mrh file name
        :type source: str.
        """
        self.source = open( source, 'rb' )
        self.header = None


    def get_data( self, chrom1='', start1=0, stop1=0, chrom2='', start2=0, stop2=0, min_resolution=0,
                  max_resolution=0, header=False, chromosomes=False, minscore=None, maxscore=None ):
        """
        :param chrom1: chromosome for sequence1
        :type chrom1: str.

        :param start1: starting bp in sequence 1
        :type start1: long

        :param stop1: stopping bp in sequence 1
        :type stop1: long

        :param chrom2: chromosome for sequence2
        :type chrom2: str.
        
        :param start2: starting bp in sequence 2
        :type start2: long

        :param stop2: stopping bp in sequence 2
        :type stop2: long

        :param min_resolution: largest bin size to return
        :type min_resolution: long

        :param max_resolution: smallest bin size to return
        :type max_resolution: long

        :param header: T/F send chromosome-specific header data instead
        :type resolution: bool

        :param chromosomes: T/F send chromosome list instead
        :type resolution: bool
        """
        if chromosomes:
            return self.load_chrom_data( )
        elif header and chrom1 != '' and chrom2 != '':
            self.chrom1 = chrom1.strip( 'chr' )
            self.chrom2 = chrom2.strip( 'chr' )
            return self.load_header( )
        elif chrom1 != '' and chrom2 != '':
            self.chrom1 = chrom1.strip( 'chr' )
            self.chrom2 = chrom2.strip( 'chr' )
            self.start1 = int( start1 )
            self.start2 = int( start2 )
            self.stop1 = int( stop1 )
            self.stop2 = int( stop2 )
            self.minres = int( min_resolution )
            self.maxres = int( max_resolution )
            self.window1 = self.stop1 - self.start1
            self.window2 = self.stop2 - self.start2
            if chrom1 == chrom2 and (stop1 > start2 or stop2 > start1):
                self.overlap = True
            else:
                self.overlap = False
            if self.header is None:
                self.load_header( )
            if minscore is None:
                self.minscore = self.header[ 'minscore' ]
            else:
                self.minscore = minscore
            if maxscore is None:
                self.maxscore = self.header[ 'maxscore' ]
            else:
                self.maxscore = maxscore
            squares = {'x1':[], 'x2':[], 'y1':[], 'y2':[], 'value':[], 'color':[]}
            if self.chrom1 != self.chrom2:
                for square in self.paint_trans_canvas( ):
                    self.interpolate_square( square, squares )
                return squares
            else:
                for square in self.paint_cis_canvas( ):
                    self.interpolate_square( square, squares )
                return squares
        return None

    def load_chrom_data( self ):
        """
        """
        self.header = {}
        # get and validate magic number for mrh type
        mrh_magic_number = '42054205'
        mrh_magic_number_size = 4
        self.source.seek( 0 )
        magic_number = binascii.b2a_hex( self.source.read( mrh_magic_number_size ) )
        if magic_number != mrh_magic_number:
            raise TypeError( 'File does not appear to be a multi-resolution heatmap file' )

        # get the number of chromosomes and whether the file includes inter-chromosome maps
        int_float32_size = 4
        trans, num_chroms = struct.unpack( 'ii', self.source.read( int_float32_size * 2 ) )
        self.trans = trans
        self.num_chroms = num_chroms
        name_sizes = struct.unpack( 'i' * num_chroms, self.source.read( int_float32_size * num_chroms ) )

        # create dictionary of chromosome names and indices, retrieve chrom indices for requested data
        self.chr2int = {}
        for i in range( num_chroms ):
            self.chr2int[ ''.join( struct.unpack( 'c' * name_sizes[i], self.source.read( name_sizes[i] ) ) ).strip( ) ] = i
        return { 'chromosomes':self.chr2int.keys( ), 'includes_trans':bool( trans ) }

    def load_header( self ):
        int_float32_size = 4
        self.load_chrom_data( )
        trans = self.trans
        num_chroms = self.num_chroms
        if self.chrom1 != self.chrom2 and not trans:
            raise TypeError( 'File does not appear to contain inter-chromosome data' )
        if trans:
            num_pairings = ( num_chroms * ( num_chroms + 1 ) ) / 2
        else:
            num_pairings = num_chroms
        if self.chrom1 in self.chr2int and self.chrom2 in self.chr2int:
            chr_index1 = self.chr2int[ self.chrom1 ]
            chr_index2 = self.chr2int[ self.chrom2 ]
        else:
            raise KeyError( 'File does not appear to contain data for the requested chromosome(s)' )
        self.transpose = False
        if trans:
            if chr_index1 > chr_index2:
                self.transpose = True
                chr_index1, chr_index2 = chr_index2, chr_index1
            pairing_index = chr_index1 * ( num_chroms - 1 ) - ( chr_index1 * ( chr_index1 - 1 ) ) / 2 + chr_index2
        else:
            pairing_index = chr_index1

        # get index of chrom/chrom-pair data starting byte and ending byte
        self.source.seek( pairing_index * int_float32_size, 1 )
        start_index, stop_index = struct.unpack( 'i' * 2, self.source.read( int_float32_size * 2 ) )
        self.header[ 'offset' ] = start_index
        self.header[ 'total_bins' ] = ( stop_index - start_index ) / int_float32_size
        skip = num_pairings - pairing_index - 1

        # get number of top-layer partitions
        if chr_index1 != chr_index2:
            # skip intra-chromosomal partitions and find inter-chromosome partitions for requested chroms
            self.source.seek( ( skip + chr_index1 + num_chroms ) * int_float32_size, 1 )
            self.header[ 'n' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            self.source.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'm' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            # if chroms are reversed from data order, switch dimensions
            if self.transpose:
                self.header['m'], self.header['n'] = self.header['n'], self.header['m']
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome partitions for requested chrom
            self.source.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'n' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal paritions
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get total number of data bins for chrom/chrom-pair
        self.source.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'data_bins' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
        self.header[ 'index_bins' ] = self.header[ 'total_bins' ] - self.header[ 'data_bins' ]
        skip = num_pairings - pairing_index - 1

        # get starting coordinate(s)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal start coordinates and find inter-chromosome start coordinates for requested chroms
            self.source.seek( ( skip + num_chroms + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'start1' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            self.source.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'start2' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome start coordinate for requested chrom
            self.source.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'start' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal start coordinates
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get stopping coordinate(s)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal stop coordinates and find inter-chromosome stop coordinates for requested chroms
            self.source.seek( ( skip + num_chroms + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'stop1' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            self.source.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'stop2' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome stop coordinate for requested chrom
            self.source.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'stop' ] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal stop coordinates
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get minimum score
        self.source.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'minscore' ] = struct.unpack( 'f', self.source.read( int_float32_size ) )[0]
        skip = num_pairings - pairing_index - 1

        # get maximum score
        self.source.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'maxscore' ] = struct.unpack( 'f', self.source.read( int_float32_size ) )[0]
        skip = num_pairings - pairing_index - 1

        # get maximum bin size (lowest resolution)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal largest bin size
            skip += 1
        self.source.seek( skip * int_float32_size, 1 )
        self.header[ 'lres'] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]

        # get maximum bin size (highest resolution)
        if trans:
            # skip intra- or interchromosomal value (whichever is unneeded)
            self.source.seek( int_float32_size, 1 )
        self.header[ 'hres'] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]

        # get zooming factor
        if trans:
            # skip intra- or interchromosomal value (whichever is unneeded)
            self.source.seek( int_float32_size, 1 )
        self.header[ 'zoom'] = struct.unpack( 'i', self.source.read( int_float32_size ) )[0]
        self.header[ 'zoom2' ] = self.header[ 'zoom' ] ** 2
        self.header[ 'zoom_d' ] = ( self.header[ 'zoom' ] * ( self.header[ 'zoom' ] + 1 ) ) / 2

        return self.header

    def paint_cis_canvas( self ):
        """
        """
        if ( min( self.stop1, self.stop2 ) <= self.header[ 'start' ] or
             max( self.start1, self.start2 ) >= self.header[ 'stop' ] ):
            return []
        # ensure positions are within bounds of data present in file
        start_pos1 = max( 0, ( self.start1 - self.header[ 'start' ] ) / self.header[ 'lres' ] )
        end_pos1 = min( self.header[ 'n' ], ( self.stop1 - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
        start_pos2 = max( 0, ( self.start2 - self.header[ 'start' ] ) / self.header[ 'lres' ] )
        end_pos2 = min( self.header[ 'n' ], ( self.stop2 - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
        
        if self.overlap:
            # if data overlap, break query into two parts, since data only covers upper-triangle
            outdata = []
            mid_coord = min(self.stop1, self.stop2)
            mid_start = min( self.header[ 'n' ], ( mid_coord - self.header[ 'start' ] ) / self.header[ 'lres' ] )
            mid_stop = min( self.header[ 'n' ], ( mid_coord - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
            if self.start1 <= self.start2:
                self.reverse = False
                self.eff_start1, self.eff_start2 = self.start1, self.start2
                self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( start_pos1, mid_stop, start_pos2, mid_stop )
            else:
                self.reverse = True
                self.eff_start2, self.eff_start1 = self.start1, self.start2
                self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( start_pos2, mid_stop, start_pos1, mid_stop )
            if self.stop1 > mid_coord:
                self.revese = True
                self.eff_start2, self.eff_start1 = self.start1, self.start2
                self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( start_pos2, end_pos2, mid_start, end_pos1 )
            elif self.stop2 > mid_coord:
                self.reverse = False
                self.eff_start1, self.eff_start2 = self.start1, self.start2
                self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( start_pos1, end_pos1, mid_start, end_pos2 )

        # if data don't overlap, figure out which range is upstream
        elif self.start2 < self.start1:
            self.reverse = True
            self.eff_start2, self.eff_start1 = self.start1, self.start2
            self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
            outdata = self.paint_cis_upper_level( start_pos2, end_pos2, start_pos1, end_pos1 )
        else:
            self.reverse = False
            self.eff_start1, self.eff_start2 = self.start1, self.start2
            self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
            outdata = self.paint_cis_upper_level( start_pos1, end_pos1, start_pos2, end_pos2 )
        return outdata

    def paint_cis_upper_level( self, start_pos1, end_pos1, start_pos2, end_pos2 ):
        """
        """
        outdata = []
        resolution = self.header[ 'lres' ]
        for i in range( start_pos1, end_pos1 ):
            # Find position in file for data with 'i' as upstream interaction
            if end_pos2 - max( i, start_pos2 ) <= 0:
                continue
            self.source.seek( self.header[ 'offset' ] + ( i * ( self.header[ 'n' ] - 1 ) -
                        ( i * ( i - 1 ) ) / 2 + max( i, start_pos2 ) ) * 4 )
            data = struct.unpack( 'f' * ( end_pos2 - max( start_pos2, i ) ),
                   self.source.read( ( end_pos2 - max( start_pos2, i ) ) * 4 ) )
            # Find position in file for indices with 'i' as upstream interaction
            self.source.seek( self.header[ 'offset' ] + ( i * ( self.header[ 'n' ] - 1 ) -
                        ( i * ( i - 1 ) ) / 2 + self.header[ 'data_bins' ] + max( i, start_pos2 ) ) * 4 )
            indices = struct.unpack( 'i' * ( end_pos2 - max( start_pos2, i ) ),
                      self.source.read( ( end_pos2 - max( start_pos2, i ) ) * 4 ) )
            for j in range( max( start_pos2, i ), end_pos2 ):
                k = j - max( start_pos2, i )
                # if data is valid, seek down to lower levels for smaller bin sizes
                if not math.isnan( data[k] ):
                    start1 = i * resolution + self.header[ 'start' ]
                    start2 = j * resolution + self.header[ 'start' ]
                    # if valid lower-level data exists, try retrieving it
                    if indices[k] != -1:
                        valid, new_outdata = self.paint_cis_lower_level( indices[k],
                                             resolution / self.header[ 'zoom' ], start1, start2 )
                    else:
                        valid = 0
                    if start1 == start2:
                        # if data is on the diagonal and not completely covered by lower-level data,
                        # only add one square
                        if resolution <= self.minres and valid < self.header[ 'zoom_d' ]:
                            outdata.append( [start1, start2, resolution, data[k]] )
                    else:
                        # if data is off-diagonal and not completely covered by lower-level data, add it to dataset
                        if resolution <= self.minres and valid < self.header[ 'zoom2' ]:
                            outdata.append( [start1, start2, resolution, data[k]] )
                            # if data is off-diagonal, check if the mirror image is also needed
                            if self.overlap and start2 + resolution > self.eff_start1 and start1 < self.eff_stop2:
                                outdata.append( [start2, start1, resolution, data[k]] )
                    if valid > 0:
                        outdata += new_outdata
        # if sequence 2 was upstream of sequence 1, flip the x and y coordinates
        if self.reverse:
            for i in range( len( outdata ) ):
                outdata[i] = [ outdata[i][1], outdata[i][0] ] + outdata[i][2:]
        return outdata

    def paint_cis_lower_level( self, index, resolution, start1, start2 ):
        """
        """
        # don't return data the is higher resolution than requested
        if resolution < self.maxres:
            return 0, []
        outdata = []
        valid = 0
        self.source.seek( self.header[ 'offset' ] + index * 4 )
        # if data is on the diagonal, lower level bin arrays will be upper-triangular, not square
        if start1 == start2:
            data = struct.unpack( 'f' * self.header[ 'zoom_d' ], self.source.read( 4 * self.header[ 'zoom_d' ] ) )
            if index < self.header[ 'index_bins' ]:
                self.source.seek( self.header[ 'offset' ] + ( index + self.header[ 'data_bins' ] ) * 4 )
                indices = struct.unpack( 'i' * self.header[ 'zoom_d' ], self.source.read( self.header[ 'zoom_d' ] * 4 ) )
            else:
                indices = None
            # for each bin, find if it has valid data, and if lower levels need to be queried
            for i in range( self.header[ 'zoom' ] ):
                for j in range( i, self.header[ 'zoom' ] ):
                    k = i * ( self.header[ 'zoom' ] - 1 ) - ( i * ( i - 1 ) ) / 2 + j
                    if not math.isnan( data[k] ):
                        start1b = start1 + i * resolution
                        start2b = start2 + j * resolution
                        # if the bin is valid but out of the query range, do return data
                        if ( start1b >= self.eff_stop1 or start1b + resolution <= self.eff_start1 or
                             start2b >= self.eff_stop2 or start2b + resolution <= self.eff_start2 ):
                            valid += 1
                        else:
                            # if the bin is valid and has valid data at a lower level, go down a level
                            if indices is not None and indices[k] != -1:
                                new_valid, new_outdata = self.paint_cis_lower_level( indices[k],
                                                         resolution / self.header[ 'zoom' ], start1b, start2b )
                            else:
                                new_valid = 0
                            # if data is high enough resolution and isn't completely covered by
                            # lower-level bins, return bin data
                            if resolution  <= self.minres:
                                if start1b == start2b:
                                    if new_valid < self.header[ 'zoom_d' ]:
                                        outdata.append( [start1b, start2b, resolution, data[k]] )
                                else:
                                    if new_valid < self.header[ 'zoom2' ]:
                                        outdata.append( [start1b, start2b, resolution, data[k]] )
                                        if ( self.overlap and start2b + resolution > self.eff_start1 and
                                             start1b < self.eff_stop2 ):
                                            outdata.append( [start2b, start1b, resolution, data[k]] )
                            if new_valid > 0:
                                outdata += new_outdata
                            valid += 1
        else:
            data = struct.unpack( 'f' * self.header[ 'zoom2' ], self.source.read( 4 * self.header[ 'zoom2' ] ) )
            if index < self.header[ 'index_bins' ]:
                self.source.seek( self.header[ 'offset' ] + ( index + self.header[ 'data_bins' ] ) * 4 )
                indices = struct.unpack( 'i' * self.header[ 'zoom2' ], self.source.read( self.header[ 'zoom2' ] * 4 ) )
            else:
                indices = None
            # for each bin, find if it has valid data, and if lower levels need to be queried
            for i in range( self.header[ 'zoom2' ] ):
                if not math.isnan( data[i] ):
                    start1b = start1 + ( i / self.header[ 'zoom' ] ) * resolution
                    start2b = start2 + ( i % self.header[ 'zoom' ] ) * resolution
                    # if the bin is valid but out of the query range, do return data
                    if ( start1b > self.eff_stop1 or start1b + resolution < self.eff_start1 or
                         start2b > self.eff_stop2 or start2b + resolution < self.eff_start2 ):
                        valid += 1
                    else:
                        # if the bin is valid and has valid data at a lower level, go down a level
                        if indices is not None and indices[i] != -1:
                            new_valid, new_outdata = self.paint_cis_lower_level( indices[i],
                                                     resolution / self.header[ 'zoom' ], start1b, start2b )
                        else:
                            new_valid = 0
                        # if data is high enough resolution and isn't completely covered by
                        # lower-level bins, return bin data
                        if resolution <= self.minres and new_valid < self.header[ 'zoom2' ]:
                            outdata.append( [start1b, start2b, resolution, data[i]] )
                            if self.overlap and start2b + resolution > self.eff_start1 and start1b < self.eff_stop2:
                                outdata.append( [start2b, start1b, resolution, data[i]] )
                        if new_valid > 0:
                            outdata += new_outdata
                        valid += 1
        return valid, outdata

    def paint_trans_canvas( self ):
        """
        """
        outdata = []
        start_pos1 = max( 0, ( self.start1 - self.header[ 'start1' ] ) / self.header[ 'lres' ] )
        end_pos1 = min( self.header[ 'n' ], ( self.stop1 - self.header[ 'start1' ] ) / self.header[ 'lres' ] + 1 )
        start_pos2 = max( 0, ( self.start2 - self.header[ 'start2' ] ) / self.header[ 'lres' ] )
        end_pos2 = min( self.header[ 'm' ], ( self.stop2 - self.header[ 'start2' ] ) / self.header[ 'lres' ] + 1 )
        resolution = self.header[ 'lres' ]

        for i in range( start_pos1, end_pos1 ):
            # Find position in file for data with 'i' as upstream interaction
            self.source.seek( self.header[ 'offset' ] + i * self.header[ 'm' ] * 4 )
            data = struct.unpack( 'f' * ( end_pos2 - start_pos2 ), self.source.read( ( end_pos2 - start_pos2 ) * 4 ) )
            # Find position in file for indices with 'i' as upstream interaction
            self.source.seek( self.header[ 'offset' ] + ( i * self.header[ 'm' ] +
                        self.header[ 'data_bins' ] + start_pos2 ) * 4 )
            indices = struct.unpack( 'i' * ( end_pos2 - start_pos2 ), self.source.read( ( end_pos2 - start_pos2 ) * 4 ) )
            for j in range( start_pos2, end_pos2 ):
                k = j - start_pos2
                # if data is valid, seek down to lower levels for smaller bin sizes
                if not math.isnan( data[k] ):
                    start1 = i * resolution + self.header[ 'start1' ]
                    start2 = j * resolution + self.header[ 'start2' ]
                    # if valid lower-level data exists, try retrieving it
                    if indices[k] != -1:
                        new_valid, new_outdata = self.paint_trans_lower_level( indices[k],
                                                 resolution / self.header[ 'zoom' ], start1, start2 )
                    else:
                        new_valid = 0
                    # if top-level square not completely covered by lower-level squares, add it to dataset
                    if resolution <= self.minres and new_valid < self.header[ 'zoom2' ]:
                        outdata.append( [start1, start2, resolution, data[k]] )
                    if new_valid > 0:
                        outdata += new_outdata
        return outdata

    def paint_trans_lower_level( self, index, resolution, start1, start2 ):
        """
        """
        # don't return data the is higher resolution than requested
        if resolution < self.maxres:
            return 0, []
        outdata = []
        valid = 0
        self.source.seek( self.header[ 'offset' ] + index * 4 )
        data = struct.unpack( 'f' * self.header[ 'zoom2' ], self.source.read( 4 * self.header[ 'zoom2' ] ) )
        if index < self.header[ 'index_bins' ]:
            self.source.seek( self.header[ 'offset' ] + ( index + self.header[ 'data_bins' ] ) * 4 )
            indices = struct.unpack( 'i' * self.header[ 'zoom2' ], self.source.read( self.header[ 'zoom2' ] * 4 ) )
        else:
            indices = None
        for i in range( self.header[ 'zoom2' ] ):
            # for each bin, find if it has valid data, and if lower levels need to be queried
            if not math.isnan( data[i] ):
                start1b = start1 + ( i / self.header[ 'zoom' ] ) * resolution
                start2b = start2 + ( i % self.header[ 'zoom' ] ) * resolution
                # if the bin is valid but out of the query range, do return data
                if ( start1b > self.stop1 or start1b + resolution < self.start1 or
                     start2b > self.stop2 or start2b + resolution < self.start2 ):
                    valid += 1
                else:
                    if indices is not None and indices[i] != -1:
                        new_valid, new_outdata = self.paint_trans_lower_level( indices[i],
                                                 resolution / self.header[ 'zoom' ], start1b, start2b )
                    else:
                        new_valid = 0
                    # if data is high enough resolution and isn't completely covered by
                    # lower-level bins, return bin data
                    if resolution <= self.minres and new_valid < self.header[ 'zoom2' ]:
                        outdata.append( [start1b, start2b, resolution, data[i]] )
                    if new_valid > 0:
                        outdata += new_outdata
                    valid += 1
        return valid, outdata

    def interpolate_square( self, square, squares ):
        """
        """
        if self.start1 == self.stop1 or self.start2 == self.stop2:
            return
        if self.transpose:
            squares['y1'].append( min( max( square[0], self.start1 ), self.stop1 ) )
            squares['x1'].append( min( max( square[1], self.start2 ), self.stop2 ) )
            squares['y2'].append( min( max( square[0] + square[2], self.start1 ), self.stop1 ) )
            squares['x2'].append( min( max( square[1] + square[2], self.start2 ), self.stop2 ) )
        else:
            squares['x1'].append( min( max( square[0], self.start1 ), self.stop1 ) )
            squares['y1'].append( min( max( square[1], self.start2 ), self.stop2 ) )
            squares['x2'].append( min( max( square[0] + square[2], self.start1 ), self.stop1 ) )
            squares['y2'].append( min( max( square[1] + square[2], self.start2 ), self.stop2 ) )
        squares['value'].append( square[3] )

        color = min( 255, max(-255, int( round( 255 * ( ( square[3] - self.minscore ) /
                ( self.maxscore - self.minscore ) * 2.0 - 1.0 ) ) ) ) )
        if color < 0:
            squares['color'].append( int( "ffff%02x%02x" % (255 + color, 255 + color ), 16 ) )
            #squares['color'].append( "#%02x%02xff" % (255 + color, 255 + color ) )
        else:
            squares['color'].append( int( "ff%02x%02xff" % (255 - color, 255 - color ), 16 ) )
            #squares['color'].append( "#ff%02x%02x" % (255 - color, 255 - color ) )
       
        return None
