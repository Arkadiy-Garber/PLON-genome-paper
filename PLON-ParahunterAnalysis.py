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


# COMBINING PSEUDOGENE, ANNOTATIONS, BLASTHITS, AND GHOSTKOALA SUMMARIES WITH THE PARAHUNTER OUTPUT
'''
sodDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
blast = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/version-3/pseudofinder_sod/pseudoSOD_proteome.faa.blastP_output.tsv")
for i in blast:
    ls = i.rstrip().split("\t")
    if ls[0] not in sodDict.keys():
        sodDict[ls[0]]["nr"] = replace(ls[12], [","], ";")

ko = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/SOD.ghostkoala.csv")
for i in ko:
    ls = i.rstrip().split(",")
    KO = (ls[5])
    sodDict[ls[0]]["ko"] = KO

prokka = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/PROKKA_10082020.faa")
for i in prokka:
    if re.match(r'>', i):
        line = (i.rstrip()[1:])
        line = allButTheFirst(line, " ")
        line = replace(line, [","], ";")
        locus = i.rstrip()[1:].split(" ")[0]
        sodDict[locus]["prokka"] = line


dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
dnds = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/version-3/pseudofinder_sod/pseudoSOD_dnds/dnds-summary.csv")
for i in dnds:
    ls = i.rstrip().split(",")
    dndsDict[ls[0]]["homolog"] = ls[1]
    dndsDict[ls[0]]["ds"] = ls[3]
    dndsDict[ls[0]]["dnds"] = ls[4]


pseudoDict = defaultdict(lambda: defaultdict(lambda: 'missing'))
sod = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/version-3/pseudofinder_sod/pseudoSOD_summary.csv")
for i in sod:
    ls = i.rstrip().split(",")
    if ls[0] != "orf":
        locus = ls[0]
        pg = "n"

        if re.findall(r'ign', locus):
            pg = "y"

        if ls[1] != "NA" and (float(ls[1]) < 0.75 or float(ls[1]) > 1.25):
            pg = "y"

        if ls[4] != "NA" and float(ls[4]) > 0.3:
            pg = "y"

        if ls[5] == "y":
            pg = "y"
        pseudoDict[locus] = pg


dupDict = defaultdict(list)
duplicates = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/parahunter/dS_summary.csv")
out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/parahunter/dS_summary-annotated.csv", "w")
for i in duplicates:
    if not re.match(r'#', i):
        ls = i.rstrip().split(",")
        print(ls)
        if ls[0] != "cluster":
            dupDict[ls[0]].append(ls)

            if float(ls[3]) > 3:
                ds = "saturated"
                dnds = "NA"
            elif float(ls[3]) < 0.001:
                ds = "low"
                dnds = "NA"
            else:
                ds = ls[3]
                dnds = str(float(ls[4]) / float(ls[3]) )

            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ds + "," + dnds + "," + str(sodDict[ls[1]]["nr"]) +
                      "," + str(sodDict[ls[1]]["ko"]) + "," + str(sodDict[ls[1]]["prokka"]) + "," +
                str(dndsDict[ls[1]]["homolog"]) + "," + str(dndsDict[ls[1]]["ds"]) + "," + str(
                    dndsDict[ls[1]]["dnds"]) + "," + str(pseudoDict[ls[1]]) + "\n")

        else:
            out.write(ls[0] + "," + ls[1] + "," + ls[2] + "," + ls[3] + ",dnds,blast,ko,prokka,HShomolog,ds,dnds,pg" + "\n")
    else:
        out.write(i.rstrip() + "\n")
out.close()


count = 0
for i in dupDict.keys():
    for j in dupDict[i]:
        print(str(j) + "\t" + str(dndsDict[j[1]]))
        count += 1
    print("")

print(count)
print(len(dupDict.keys()))
'''


# DOING THE SAME FOR GLOSSINIDIUS SOPE
'''
Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
prokka = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/SOPE/PROKKA_10122020.faa")
for i in prokka:
    if re.match(r'>', i):
        line = (i.rstrip()[1:])
        line = allButTheFirst(line, " ")
        line = replace(line, [","], ";")
        locus = i.rstrip()[1:].split(" ")[0]
        Dict[locus]["prokka"] = line


out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/SOPE/parahunter/dS_summary-annotated.csv", "w")
prokka = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/SOPE/parahunter/dS_summary.csv")
for i in prokka:
    ls = i.rstrip().split(",")
    if ls[0] != "cluster":
        if not re.match(r'#', i):
            out.write(i.rstrip() + "," + str(Dict[ls[1]]["prokka"]) + "\n")
        else:
            out.write(i.rstrip() + "\n")
    else:
        out.write(i.rstrip() + "\n")
out.close()
'''

# COMBINING PARAHUNTER OUTPUTS
'''
out = open("/Users/arkadiygarber/Desktop/manuscripts/PLON-genome-manuscript/Supplemental Files/SupplementalFile1.csv", "w")
file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/parahunter/dS_summary-annotated.csv")
for i in file:
    ls = i.rstrip().split(",")
    if ls[0] != "cluster":
        if not re.match(r'#', i):
            out.write("Sod. endolongispinus," + i.rstrip() + "\n")
        else:
            out.write("#####################################\n")
    else:
        out.write("endosymbionts" + "," + i.rstrip() + "\n")

file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sym/parahunter/dS_summary-annotated.csv")
for i in file:
    ls = i.rstrip().split(",")
    if ls[0] != "cluster":
        if not re.match(r'#', i):
            out.write("Sym. endolongispinus," + i.rstrip() + "\n")
        else:
            out.write("#####################################\n")
    else:
        pass

file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/SOPE/parahunter/dS_summary-annotated.csv")
for i in file:
    ls = i.rstrip().split(",")
    if ls[0] != "cluster":
        if not re.match(r'#', i):
            out.write("S. pierantonius str. SOPE," + i.rstrip() + "\n")
        else:
            out.write("#####################################\n")
    else:
        pass

file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/other_insects/tsetse/Sodalis_glossinidus/parahunter/dS_summary-annotated.csv")
for i in file:
    ls = i.rstrip().split(",")
    if ls[0] != "cluster":
        if not re.match(r'#', i):
            out.write("S. glossinidius," + i.rstrip() + "\n")
        else:
            out.write("#####################################\n")
    else:
        pass

out.close()
'''


# MAKING PLOTS
'''
paraDict = defaultdict(list)
parahunter = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/parahunter/dS_summary-annotated.csv")
for i in parahunter:
    ls = i.rstrip().split(",")
    if not re.match(r'#', i):
        paraDict[ls[0]].append(ls)

redundDict = defaultdict(list)
out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/prokka_sod/parahunter/dS_summary-plot.csv", "w")
out.write("locus,ds,transposase\n")
for i in paraDict.keys():
    counter = 0
    for j in paraDict[i]:
        if re.findall(r'ransposa', j[5]):
            counter += 1
    if counter > 0:
        for j in paraDict[i]:
            if j[1] not in redundDict.keys():
                redundDict[j[2]].append(j)
                if j[3] == "low":
                    ds = 0
                else:
                    ds = float(j[3])
                out.write(j[1] + "," + str(ds) + ",y" + "\n")
    else:
        for j in paraDict[i]:
            if j[1] not in redundDict.keys():
                redundDict[j[2]].append(j)
                if (j[0]) != "cluster":
                    print(j)
                    if j[3] == "low":
                        ds = 0
                    elif j[3] == "saturated":
                        ds = "NA"

                    else:
                        ds = float(j[3])
                    out.write(j[1] + "," + str(ds) + ",n" + "\n")


        print("")
'''








