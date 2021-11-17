import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC_skew
import datetime

parser = argparse.ArgumentParser(description="""GC_GC-Skew_FASTA  -  Meriam Guellil  -  July 2021 v2.0""", epilog="""Can output GC and GC-Skew of FASTA files in intervals for plotting""")
parser.add_argument('-f',metavar='FASTA file', dest='fasta', required=True, type=str, help='FASTA Input')
parser.add_argument('-w',metavar='window size', dest='window', required=True, type=int, help='window size used to calculate  gc and gc-skew values')
parser.add_argument('--gc', action='store_true', dest='gc', required=False, help='output gc content')
parser.add_argument('--gc-skew', action='store_true', dest='gc_skew', required=False, help='output gc-skew')
parser.add_argument('-o',metavar='Output Prefix', dest='out', required=True, help='Output prefix')
args= parser.parse_args()

#Store Date for output format:
date = datetime.datetime.now().strftime("%d%m%Y")

#For each chromosome generate data ans store in output:
for record in SeqIO.parse(args.fasta,'fasta'):
    if args.gc_skew is True:
        with open ((str(args.out) + "_" + "GC-Skew" + "_W" + str(args.window) + "_" + date + ".bed"),'a') as out_bed:
            for i in range(0, len(record.seq), (int(args.window))):
                gc_Skew = GC_skew(record.seq[i:i+(int(args.window)-1)], window=args.window)
                out_bed.write(str(record.id) + '\t' + str(i) + '\t' + str(i+(int(args.window)-1)) + '\t' +str(gc_Skew[0]) + '\n')
    if args.gc is True:
        with open ((str(args.out) + "_" + "GC-Content" + "_W" + str(args.window) + "_" + date + ".bed"),'a') as out_bed:
            for i in range(0, len(record.seq), (int(args.window))):
                gc_v = GC(record.seq[i:i+(int(args.window)-1)])
                out_bed.write(str(record.id) + '\t' + str(i) + '\t' + str(i+(int(args.window)-1)) + '\t' +str(gc_v) + '\n')

#Meriam Guellil 2021