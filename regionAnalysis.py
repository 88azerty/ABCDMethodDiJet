#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Analyse data from RegionsHV.log files for complete analysis.')
parser.add_argument('-i', '--inputfile',   help='Input File (required).',  required = True, type = argparse.FileType('r'))
parser.add_argument('-o', '--outputfile',  help='Output File (required).', required = True, type = argparse.FileType('w'))
parser.add_argument('-H', '--horizontal',   help='Position of horizontal boundary.', required = True, type = float)
parser.add_argument('-V', '--vertical',   help='Position of vertical boundary.', required = True, type = float)
args=parser.parse_args()


with args.inputfile as f:
	InputContent = f.readlines()		#cargar en memoria archivo de entrada y cerrarlo
    f.close()

data = []

for i in InputContent[1:]:
    data.append([   str(i.split()[0]),      #DataSetName
                    int(i.split()[1]),      #RegionA
                    int(i.split()[2]),      #RegionB
                    int(i.split()[3]),      #RegionC
                    int(i.split()[4]),      #RegionD
                    float(i.split()[5]),    #RegionAW
                    float(i.split()[6]),    #RegionBW
                    float(i.split()[7]),    #RegionCW
                    float(i.split()[8]),    #RegionDW
                    float(i.split()[9]),    #RegionAE
                    float(i.split()[10]),   #RegionBE
                    float(i.split()[11]),   #RegionCE
                    float(i.split()[12])    #RegionDE
                ])

RegionA = sum(data[][1])
RegionB = sum(data[][2])
RegionC = sum(data[][3])
RegionD = sum(data[][4])
RegionAW = sum(data[][5])
RegionBW = sum(data[][6])
RegionCW = sum(data[][7])
RegionDW = sum(data[][8])

with args.outputfile as g:
#   g.write('Boundaries' + '\t' + 'RegionA' + '\t' + 'RegionB'+ '\t' + 'RegionC'+ '\t' + 'RegionD'+ '\t' + 'RegionAW' + '\t' + 'RegionBW' + '\t' + 'RegionCW' + '\t' + 'RegionDW')
    g.seek(0, 2)
    g.write(        'H' + str(args.horizontal) + 'V' + str(args.vertical) + '\t' +
                    str(RegionA) + '\t' +
                    str(RegionB) + '\t' +
                    str(RegionC) + '\t' +
                    str(RegionD) + '\t' +
                    str(RegionAW) + '\t' +
                    str(RegionBW) + '\t' +
                    str(RegionCW) + '\t' +
                    str(RegionDW) + '\t'
            )
