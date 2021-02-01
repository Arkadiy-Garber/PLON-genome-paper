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
    return x[0:len(x)-1]


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


# SUMMARIZING PSEUDOFINDER RESULTS INTO SUMMARY FILE FOR PLOTTING IN MULTIPLE DIMENSIONS (I.E. MAKING FIGURE 3)
# THIS SAME SCRIPT CAN BE USED FOR BOTH SYM AND SOD ENDOLONGISPINUS RESULTS, JUST CHANGE THE FILENAME DESIGNATIONS/ADDRESSES
'''
##### PROTEOME
seq = open("pseudofinder_proteome.faa")
seq = fasta2(seq)

blastDict = defaultdict(lambda: defaultdict(list))
blast = open("pseudofinder_proteome.faa.blastP_output.tsv")
for i in blast:
    ls = i.rstrip().split("\t")
    geneLength = len(seq[ls[0]])
    blastDict[ls[0]]["geneLength"].append(len(seq[ls[0]]))
    blastDict[ls[0]]["homologLength"].append(int(ls[3]))

MasterDict = defaultdict(lambda: defaultdict(lambda: "NA"))
for i in blastDict.keys():
    geneLength = blastDict[i]["geneLength"][0]
    homologLength = statistics.mean(blastDict[i]["homologLength"])
    perc = geneLength/homologLength
    MasterDict[i]["perc"] = perc

##### INTERGENIC
ig = open("pseudofinder_intergenic.fasta")
ig = fasta2(ig)


igDict = defaultdict(lambda: defaultdict(list))
blast = open("pseudofinder_intergenic.fasta.blastX_output.tsv")
for i in blast:
    ls = i.rstrip().split("\t")
    aln = int(ls[9]) - int(ls[8])
    igDict[ls[0]]["aln"].append(aln)
    igDict[ls[0]]["homologLength"].append(int(ls[3]))

for i in igDict.keys():
    if len(igDict[i]["aln"]) > 4:
        aln = (statistics.mean(igDict[i]["aln"]))
        homologLength = (statistics.mean(igDict[i]["homologLength"]))
        perc = (aln/homologLength)
        MasterDict[i]["perc"] = perc
        MasterDict[i]["ig"] = "y"


pseudos = open("pseudofinder_pseudos.gff")
for i in pseudos:
    ls = i.rstrip().split("\t")
    if re.findall(r'fragmentation', i):

        LS = (lastItem(ls[8].split(";")).split("=")[1].split(","))
        for j in LS:
            MasterDict[j]["fragment"] = "y"

dnds = open("pseudofinder_dnds/dnds-summary.csv")
for i in dnds:
    ls = i.rstrip().split(",")
    if ls[2] != "dN":
        MasterDict[ls[0]]["dn"] = float(ls[2])
        MasterDict[ls[0]]["ds"] = float(ls[3])
        MasterDict[ls[0]]["dnds"] = float(ls[4])

for i in MasterDict.keys():
    if MasterDict[i]["fragment"] != "y":
        MasterDict[i]["fragment"] = "n"


out = open("pseudofinder_summary.csv", "w")
out.write("orf,length,dn,ds,dnds,fragment,ig\n")
for i in MasterDict.keys():
    out.write(i + "," + str(MasterDict[i]["perc"]) + "," + str(MasterDict[i]["dn"]) + "," + str(MasterDict[i]["ds"]) + "," +
              str(MasterDict[i]["dnds"]) + "," + MasterDict[i]["fragment"] + "," + MasterDict[i]["ig"] + "\n")

out.close()

pseudoDict = defaultdict(list)
fragments = 0
intergenic = 0
short = 0
long = 0
dnds = 0
for i in MasterDict.keys():
    try:
        if MasterDict[i]["dnds"] > 0.3:
            dnds += 1
            pseudoDict[i].append("dnds")
    except:
        pass

    try:
        if MasterDict[i]["fragment"] == 'y':
            fragments += 1
            pseudoDict[i].append("fragments")
    except TypeError:
        pass

    try:
        if MasterDict[i]["ig"] == 'y':
            intergenic += 1
            pseudoDict[i].append("intergenic")
    except TypeError:
        pass

    try:
        if MasterDict[i]["perc"] < 0.75:
            short += 1
            pseudoDict[i].append("short")
    except TypeError:
        pass

    try:
        if MasterDict[i]["perc"] > 1.25:
            long += 1
            pseudoDict[i].append("long")
    except TypeError:
        pass

pseudoDict2 = defaultdict(list)
for i in pseudoDict.keys():
    CAT = "_".join(sorted(pseudoDict[i]))
    pseudoDict2[CAT].append(i)

total = 0
for i in sorted(pseudoDict2.keys()):
    print(i + "\t" + str(len(pseudoDict2[i])))
    total += len(pseudoDict2[i])


print(total)

print(long)
print(short)
print(intergenic)
print(fragments)
print(dnds)
'''


