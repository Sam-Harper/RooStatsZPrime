#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Write all workspaces in specified files to a single file')
parser.add_argument("inputFiles",type=str,nargs="+",help="input files")
parser.add_argument("--outFile",type=str,required=True,help="output filename")
args = parser.parse_args()
#print args
from ROOT import RooWorkspace
from ROOT import TFile

outFile = TFile(args.outFile,"RECREATE")

for filename in args.inputFiles:
    rootFile = TFile(filename,"READ")
    for key in rootFile.GetListOfKeys():
 #       help(key)
   #     print key.Class(),key.Class_Name(),key.GetName()
        obj = rootFile.Get(key.GetName())
        if(obj.Class_Name()=="RooWorkspace"):
            print "writing ",obj.GetName()
            outFile.WriteTObject(obj)

outFile.Close()
