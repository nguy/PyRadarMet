#! /usr/bin/env python
"""
An exectutable program that outputs characteristics of raw Sigmet format
files of radar data.  An optional output file may be created.

EXAMPLE USAGE::
To create a list echoed to terminal:
listfile_raw.py "/Users/nickguy/data/iphex/noxp/140615/raw/NOX*.RAW*" dump

To create a list echoed to terminal and output file
listfile_raw.py -m "/Users/nickguy/data/iphex/noxp/140615/raw/NOX*.RAW*" info_140615.txt

HISTORY::
19 Jun 2014   Nick Guy (NRC; NOAA/NSSL)

"""
import pyart.io as fio
import argparse
import os
from glob import glob


if __name__ == '__main__':

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Show file characteristics of raw Sigmet files.')
    parser.add_argument('searchstring', type=str, help='radar file(s) to list, if more than one file quotations are needed')
    parser.add_argument('outFile', type=str, help='output file name to store information')

    igroup = parser.add_argument_group(
        title='output characteristic file, optional',
        description=('Decide whether to create output file:'))

    parser.add_argument('-m', '--make_file', action='store_true',
                        help='create an output text file with characteristics')

    args = parser.parse_args()

    # Search for the file(s)
    flist = glob(args.searchstring)

    # Create an output file if requested
    if args.make_file:
        OutTx = open(args.outFile, 'w')

    for file in flist:
        fLName, fileExt = os.path.splitext(file)
        fName = os.path.basename(file)

        # Get the file size
        fSize = os.path.getsize(file)

        # read in the file
        radar = fio.read_sigmet(file, file_field_names=True)

        # Create a latitude string
        if radar.latitude['units'] == 'degrees_north':
            rlat = str(radar.latitude['data'])[1:-1] + 'N'
        if radar.latitude['units'] == 'degrees_south':
            rlat = str(radar.latitude['data'])[1:-1] + 'S'

        # Create a longitude string
        if radar.longitude['units'] == 'degrees_east':
            rlon = str(radar.longitude['data'])[1:-1] + 'E'
        if radar.longitude['units'] == 'degrees_west':
            rlon = str(radar.longitude['data'])[1:-1] + 'W'

        # Create field names, check for extended header
        if radar.metadata['sigmet_extended_header'] == 'true':
            Vars = "'Xhdr', " + str(radar.fields.keys())[1:-1]
        else:
            Vars = str(radar.fields.keys())[1:-1]

        # Create Altitude variable
        Alt = str(radar.altitude['data'])[1:-1] + ' ' + radar.altitude['units']

        # Create PRF variable
        PRF = str(1./radar.instrument_parameters['prt']['data'][0]) + ' Hz'

        # Create Pulse width variable
        pWid = str(radar.instrument_parameters['pulse_width']['data'][0]) + ' us'

        txVar = (fName+" "+radar.metadata['instrument_name'] + " " +
                 radar.metadata['sigmet_task_name'] + " File Size: " + str(fSize) +
                 "  Alt: " + Alt + "  PRF: " + PRF +
                 " Pulse Width: " + pWid + " " + Vars + " " + rlat + " " + rlon)

        # Echo this data to the terminal window
        print(txVar)

        # Add to output file if desired
        if args.make_file:
            OutTx.write(txVar + "\n")

        del fName, fileExt, radar, rlat, rlon, Vars, Alt, PRF, pWid, txVar

    if args.make_file:
        OutTx.close()
