#!/usr/bin/env python3
from collections import defaultdict
import re
import statistics
import numpy as np
import os


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        try:
            count += float(i)
        except ValueError:
            pass
    return count/len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count/len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


# *********************************************************************************************
# ************ EXTRACTING AND BINNING CONTIGS FROM THE UNICYCLER CO-ASSEMBLY ******************
# *********************************************************************************************
'''
assemblyDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
assembly = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/assembly.fasta")
assembly = fasta(assembly)
for i in assembly.keys():
    contig = (i.split(" ")[0])
    assemblyDict[contig]["header"] = i
    assemblyDict[contig]["seq"] = assembly[i]

# CONTIG AFFILIATIONS (NUMBERS CORRESPOND TO SEQ HEADERS) INFERRED FROM SPRAYANDPRAY OUTPUT BY VISUAL INSPECTION
plon1 = ["10", "5", "8", "15", "2", "13", "9", "7", "6", "3"]
plon2 = ["1", "4", "11", "12"]

out1 = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SYM.fasta", "w")
out2 = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SOD.fasta", "w")
for i in assemblyDict.keys():
    if i in plon1:
        out1.write(">" + assemblyDict[i]["header"] + "\n")
        out1.write(assemblyDict[i]["seq"] + "\n")
    elif i in plon2:
        out2.write(">" + assemblyDict[i]["header"] + "\n")
        out2.write(assemblyDict[i]["seq"] + "\n")
    else:
        pass
        #print(i)
'''


# *********************************************************************************************
# ***************************** RE-ANALYSIS OF NEW ASSEMBLY ***********************************
# COMBINING THE VARIOUS ANNOTATION PIPELINE OUTPUT FILES (PROKKA/GHOSTKOALA) WITH PSEUDOFINDER OUTPUT
# *********************************************************************************************
'''
seqs = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SYM.faa")
seqs = fasta2(seqs)

nameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
prokka = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SYM.faa")
for i in prokka:
    if re.match(r'>', i):
        locus = i.rstrip()[1:].split(" ")[0]
        annotation = (allButTheFirst(i.rstrip(), " "))
        nameDict[locus] = annotation


koDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
ko = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SYM.ghostkoala.txt")
for i in ko:
    ls = i.rstrip().split(",")
    koDict[ls[0]]["ko"] = ls[1]
    koDict[ls[0]]["cat"] = ls[2]
    koDict[ls[0]]["fam"] = ls[3]
    koDict[ls[0]]["sys"] = ls[4]
    koDict[ls[0]]["gene"] = ls[5]
    koDict[ls[0]]["taxa1"] = ls[6]
    koDict[ls[0]]["taxa2"] = ls[7]
    koDict[ls[0]]["taxa3"] = ls[8]

pseudoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
pseudos = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/pseudofinder_sym2/pseudoSYM_pseudos.gff")
for i in pseudos:
    ls = i.rstrip().split("\t")
    if not re.match(r'#', i):
        locus = (lastItem(ls[8].split(";")).split("=")[1])
        # if len(locus.split(",")) > 1:
        #     print(locus)
        # contig = ls[0].split("|")[2]
        contig = ls[0]
        start = ls[3]
        end = ls[4]
        strand = ls[6]
        pseudoDict[locus]["contig"] = contig
        pseudoDict[locus]["start"] = start
        pseudoDict[locus]["end"] = end
        pseudoDict[locus]["strand"] = strand
        pseudoDict[locus]["pg"] = "y"


pseudos = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/pseudofinder_sym2/pseudoSYM_intact.gff")
for i in pseudos:
    ls = i.rstrip().split("\t")
    if not re.match(r'#', i):
        locus = (lastItem(ls[8].split(";")).split("=")[1])
        # print(locus)
        # if len(locus.split(",")) > 1:
        #     print(locus)
        # contig = ls[0].split("|")[2]
        contig = ls[0]
        start = ls[3]
        end = ls[4]
        strand = ls[6]
        pseudoDict[locus]["contig"] = contig
        pseudoDict[locus]["start"] = start
        pseudoDict[locus]["end"] = end
        pseudoDict[locus]["strand"] = strand
        pseudoDict[locus]["pg"] = "n"


blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
blastx = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/pseudofinder_sym2/pseudoSYM_intergenic.fasta.blastX_output.tsv")
for i in blastx:
    ls = i.rstrip().split("\t")
    if ls[0] not in blastDict.keys():
        blastDict[ls[0]] = ls[12]

blastp = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/pseudofinder_sym2/pseudoSYM_proteome.faa.blastP_output.tsv")
for i in blastp:
    ls = i.rstrip().split("\t")
    if ls[0] not in blastDict.keys():
        blastDict[ls[0]] = ls[12]


dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
dnds = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/pseudofinder_sym2/pseudoSYM_dnds/dnds-summary.csv")
for i in dnds:
    ls = i.rstrip().split(",")
    dndsDict[ls[0]]["homolog"] = ls[1]
    dndsDict[ls[0]]["dn"] = ls[2]
    dndsDict[ls[0]]["ds"] = ls[3]
    dndsDict[ls[0]]["dnds"] = ls[4]


MasterDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
for i in pseudoDict.keys():
    key = (replace(i, [","], "|"))
    for j in i.split(","):

        if len(nameDict[i]) > 0:
            MasterDict[key][j]["prokka"] = nameDict[j]

        if len(blastDict[i]) > 0:
            MasterDict[key][j]["blastHit"] = blastDict[j]

        if len(koDict[j]) > 0:
            MasterDict[key][j]["ko"] = koDict[j]["ko"]
            MasterDict[key][j]["cat"] = koDict[j]["cat"]
            MasterDict[key][j]["fam"] = koDict[j]["fam"]
            MasterDict[key][j]["sys"] = koDict[j]["sys"]
            MasterDict[key][j]["gene"] = koDict[j]["gene"]
            MasterDict[key][j]["taxa1"] = koDict[j]["taxa1"]
            MasterDict[key][j]["taxa2"] = koDict[j]["taxa2"]
            MasterDict[key][j]["taxa3"] = koDict[j]["taxa3"]


        if len(dndsDict[j]) > 0:
            MasterDict[key][j]["ref"] = dndsDict[j]["homolog"]
            MasterDict[key][j]["dn"] = dndsDict[j]["dn"]
            MasterDict[key][j]["ds"] = dndsDict[j]["ds"]
            MasterDict[key][j]["dnds"] = dndsDict[j]["dnds"]

        if len(pseudoDict[i]) > 0:
            MasterDict[key][j]["contig"] = pseudoDict[i]["contig"]
            MasterDict[key][j]["start"] = pseudoDict[i]["start"]
            MasterDict[key][j]["end"] = pseudoDict[i]["end"]
            MasterDict[key][j]["strand"] = pseudoDict[i]["strand"]


MasterDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in MasterDict.keys():
    for j in MasterDict[i]:
        print(j)
        print(replace(MasterDict[i][j]["gene"], [","], ";"))
        print("")
        MasterDict2[i]["ref"] = replace(MasterDict[i][j]["ref"], [","], ";")
        MasterDict2[i]["dn"] = replace(MasterDict[i][j]["dn"], [","], ";")
        MasterDict2[i]["ds"] = replace(MasterDict[i][j]["ds"], [","], ";")
        MasterDict2[i]["dnds"] = replace(MasterDict[i][j]["dnds"], [","], ";")
        MasterDict2[i]["contig"] = replace(MasterDict[i][j]["contig"], [","], ";")
        MasterDict2[i]["start"] = replace(MasterDict[i][j]["start"], [","], ";")
        MasterDict2[i]["end"] = replace(MasterDict[i][j]["end"], [","], ";")
        MasterDict2[i]["strand"] = replace(MasterDict[i][j]["strand"], [","], ";")
        MasterDict2[i]["ko"] = replace(MasterDict[i][j]["ko"], [","], ";")
        MasterDict2[i]["cat"] = replace(MasterDict[i][j]["cat"], [","], ";")
        MasterDict2[i]["fam"] = replace(MasterDict[i][j]["fam"], [","], ";")
        MasterDict2[i]["sys"] = replace(MasterDict[i][j]["sys"], [","], ";")
        MasterDict2[i]["gene"] = replace(MasterDict[i][j]["gene"], [","], ";")
        MasterDict2[i]["taxa1"] = replace(MasterDict[i][j]["taxa1"], [","], ";")
        MasterDict2[i]["taxa2"] = replace(MasterDict[i][j]["taxa2"], [","], ";")
        MasterDict2[i]["taxa3"] = replace(MasterDict[i][j]["taxa3"], [","], ";")
        MasterDict2[i]["prokka"] = replace(MasterDict[i][j]["prokka"], [","], ";")
        MasterDict2[i]["blastHit"] = replace(MasterDict[i][j]["blastHit"], [","], ";")


out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SYM-summary.csv", "w")
out.write("locus,pseudogene,contig,start,end,strand,prokka,gene,cat,fam,sys,ko,taxa1,taxa2,taxa3,ref,ds,dn,dnds,blastHit,seq\n")
for i in MasterDict2.keys():
    seq = ''
    for j in i.split("|"):
        if len(seqs[j]) > 0:
            seq += seqs[j]

    originaLocus = replace(i, ["|"], ",")
    out.write(i + "," + pseudoDict[originaLocus]["pg"] + "," + MasterDict2[i]["contig"] + "," + MasterDict2[i]["start"] + "," +
              MasterDict2[i]["end"] + "," + MasterDict2[i]["strand"] + "," + MasterDict2[i]["prokka"] + "," +
              MasterDict2[i]["gene"] + "," + MasterDict2[i]["cat"] + "," +
              MasterDict2[i]["fam"] + "," + MasterDict2[i]["sys"] + "," +
              MasterDict2[i]["ko"] + "," + MasterDict2[i]["taxa1"] + "," + MasterDict2[i]["taxa2"] + "," +
              MasterDict2[i]["taxa3"] + "," + MasterDict2[i]["ref"] + "," + MasterDict2[i]["ds"] + "," +
              MasterDict2[i]["dn"] + "," + MasterDict2[i]["dnds"] + "," + MasterDict2[i]["blastHit"] + "," +
              str(seq) + "\n")
'''

# *********************************************************************************************
# ***************** DEFINING VARIOUS METABOLIC AND CELL STRUCTURAL PATHWAYS *******************
# *********************************************************************************************

'''
membranes = ["mscL;", "cydC;", "cydD;", "lptA;", "lptB;;", "lptC;", "lptD;", "lptE;", "lptF;", "lptG;",
       "znuA;", "znuB;", "znuC;", "pit;", "murJ;", "manX;", "manY;", "manZ;", "ptsH;", "ptsI;"]

peptidoglycan = ["murA", "murB", "murC", "murD", "murE", "murF", "murG", "murH", "murI", "murJ", "glmM", "glmS",
                  "glmU", "amiD", "dapF", "ddl", "ftsW", "glyA", "lpoB", "mepM", "metC", "mltB", "mltC", "mraY",
                  "nlpD", "pbp"]

arg = ["carA", "carB", "OAT", "argF", "argI", "argG", "argH", "argD"]
lys = ["dapA", "dapB", "dapC", "dapD", "dapE",	"dapF",	"lysA", "lysC", "asd"]
thr = ["thrA",	"thrB",	"thrC"]
phe = ["pheA", "tyrA", "AAT"]
trp = ["trpE",	"trpG",	"trpD",	"trpC",	"trpA",	"trpB"]
val = ["ilvN", "ilvB", "ilvG", "ilvE", "ilvY", "ilvY",	"ilvC",	"ilvD",	"BCA"]
his = ["hisG",	"hisI", "hisE",	"hisA",	"hisF",	"hisH",	"hisB",	"hisC", "hisN",	"hisD"]
cys = ["cysE",	"cysK", "cysM"]
met = ["CGL", "metE", "metF", "metH", "metR"]

chor = ["aroG",	"aroB",	"aroD", "aroQ", "aroH", "aroL",	"aroE",	"aroK",	"aroA",	"aroC", "tyrR"]

riboflavin = ["ribA", "ribD", "ribB", "ribE", "ribH", "ribC"]

biotin = ["bioA", "bioB", "bioC", "bioD", "bioF"]

trans = ["infA", "infB", "infC", "tufA", "fusA", "tsf", "prfA", "prfB", "rrf", "fmt", "def"]

trna = ["argS", "cysC", "glnS", "gltX", "ileS", "leuS", "valS", "alaS", "asnS", "aspS", "tyrS",
        "trpS", "thrS", "serS", "proS", "pheST", "pheS", "pheT", "metG", "lysS", "hisS", "glySQ"]

suf = ["sufA", "sufB", "sufC", "sufD", "sufS", "sufE", "iscS", "iscU", "iscA", "hscB", "hscA", "fdx"]

tRNAsuf = ["iscS", "tusA", "tusB", "tusC", "tusD", "tusE", "mnmA", "mnmG", "mnmE", "sufA", "sufB", "sufC", "sufD", "sufE", "sufS"]

sec = ["secM", "secB", "secG", "secA", "secE", "secY", "tatA", "tatB", "tatE", "tatC", "ftsY", "ffh", "lepA", "lepB"]

tat = ["tatA", "tatB", "tatC", "tatD", "torA", "torD"]

central = ["eno;", "pgk", "prsA", "pta;", "gapA", "pyk;", "zwf;", "aceE", "aceF", "nifJ", "gpmA", "glk;",
           "pgi", "pfkA", "fbaA", "ackA", "pgl;", "gnd;", "rpiA", "prsA"]

chaperones = ["dnaK", "GRPE", "djlA", "dnaJ", "ibpA", "ibpB"]
'''


# *********************************************************************************************
# ************************ CALCULATING AVERAGE DEPTH PER SYMBIONT *****************************
# *********************************************************************************************
'''
total = 0
totalLength = 0
sym = ["2", "3", "5", "6", "7", "8", "9", "10", "13"] # no plasmid
sod = ["1", "4", "12"] # no plasmid
depth = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/assembly.depth")
for i in depth:
    ls = i.rstrip().split("\t")
    print(ls)
    if ls[0] in sod:
        depth = float(ls[1]) * float(ls[2])
        total += depth
        totalLength += float(ls[1])

print(total/totalLength)
'''








