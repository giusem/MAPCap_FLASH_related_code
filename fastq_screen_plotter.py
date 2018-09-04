#!/usr/bin/env python

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

def main():
    parser =argparse.ArgumentParser(description='Plot the output of fastqc_screen')
    parser.add_argument('-1', '--inputfile', help='The txt file from fastq_screen', required=True)
    args = parser.parse_args()

    plotFastqScreen(args.inputfile)

def plotFastqScreen(fname):
    species=[]
    ohol=[]
    mhol=[]
    ohml=[]
    mhml=[]
    for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
        if(len(line) == 0) :
            break
        if(line[0].startswith("#")) :
            continue
        if(line[0].startswith("Library")) :
            continue
        species.append(line[0])
        ohol.append(float(line[5]))
        mhol.append(float(line[7]))
        ohml.append(float(line[9]))
        mhml.append(float(line[11]))

    ohol = np.array(ohol)
    mhol = np.array(mhol)
    ohml = np.array(ohml)
    mhml = np.array(mhml)

    ind = np.arange(len(species))
    p1 = plt.bar(ind, tuple(ohol), color="#0000FF")
    p2 = plt.bar(ind, tuple(mhol), color="#6699FF", bottom=tuple(ohol))
    p3 = plt.bar(ind, tuple(ohml), color="#FF0000", bottom=tuple(ohol+mhol))
    p4 = plt.bar(ind, tuple(mhml), color="#FF6699", bottom=tuple(ohol+mhol+ohml))

    plt.title("%s" % fname.replace("_R1_screen.txt","").split("/")[-1])
    plt.ylabel("%")
    plt.ylim((0,105))
    plt.xticks(ind+0.4, species, rotation="vertical")
    plt.yticks(np.arange(0,110,10))
    plt.legend((p4[0], p3[0], p2[0], p1[0]), ("repeat", "conserved", "multimap", "unique"))
    plt.tight_layout()
    plt.savefig("%s.png" % fname.replace("_screen.txt","_screen"))
    plt.close()


if __name__ == "__main__":
    main()
