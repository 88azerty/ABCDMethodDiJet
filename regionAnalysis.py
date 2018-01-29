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
	data.append([   str(i.split('\t')[0]),      #DataSetName
	int(i.split('\t')[1]),
	int(i.split('\t')[2]),
	int(i.split('\t')[3]),
	int(i.split('\t')[4]),
	float(i.split('\t')[5]),
	float(i.split('\t')[6]),
	float(i.split('\t')[7]),
	float(i.split('\t')[8]),
	float(i.split('\t')[9]),
	float(i.split('\t')[10]),
	float(i.split('\t')[11]),
	float(i.split('\t')[12])
	])

columnA = [e[1] for e in data]
columnB = [e[2] for e in data]
columnC = [e[3] for e in data]
columnD = [e[4] for e in data]
columnAW = [e[5] for e in data]
columnBW = [e[6] for e in data]
columnCW = [e[7] for e in data]
columnDW = [e[8] for e in data]

RegionA = sum( columnA )
RegionB = sum( columnB )
RegionC = sum( columnC )
RegionD = sum( columnD )
RegionAW = sum( columnAW )
RegionBW = sum( columnBW )
RegionCW = sum( columnCW )
RegionDW = sum( columnDW )

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
