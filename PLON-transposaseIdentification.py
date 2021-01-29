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
            prot.append("")
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


# MANUAL IDENTIFICATION FROM FRAGMENTED ASSEMBLY (DID NOT WORK, NOT VERY MANY STILL IDENTIFIED)
'''
assembly = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/contigs/SyPl.bgpipe.output.fasta")
assembly = fasta2(assembly)

blastDict = defaultdict(list)
blast = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/contigs/IS.SyPl.blast")
for i in blast:
    ls = i.rstrip().split("\t")
    start = sorted([int(ls[8]), int(ls[9])])[0]
    end = sorted([int(ls[8]), int(ls[9])])[1]
    blastDict[ls[1]].append([ls[0], start, end, ls[10]])

blastDict2 = defaultdict(lambda: defaultdict(list))
for i in blastDict.keys():
    # print(i)
    # print(len(assembly[i]))
    # print(blastDict[i])
    for j in blastDict[i]:
        # print(j)
        for k in range(j[1], j[2]+1):
            blastDict2[i][k].append(j[0])
    # print("")

counter = 0
for i in blastDict2.keys():
    print(i)
    print(len(assembly[i]))
    LS = []
    for j in sorted(blastDict2[i]):
        LS.append(j)

    clu = cluster(LS, 1)
    for j in clu:
        print(j)
        counter += 1


    # for j in sorted(blastDict2[i]):
    #     print(str(j) + "\t" + str(blastDict2[i][j]))
    print("")

print(counter)
'''

# REWRITING TRANSPOSASE SUMMARY FILE WITH NEW RESULTS, TAKING INTO ACCOUNT PSEUDOGENIZATION IN TRANSPOSASES
# ALL THE FILES IN THIS DIRECTORY (GFFDIR AND PSEUDODIR) ARE AVAILABLE HERE: https://doi.org/10.6084/m9.figshare.13661189
'''
isDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
gffDir = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/IS")
for i in gffDir:
    if re.findall(r'gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/IS/%s" % i)
        for j in file:
            if not re.match(r'#', j):
                ls = j.rstrip().split("\t")
                ID = (ls[8].split(";")[0].split("=")[1])
                family = (ls[8].split("product=")[1].split(";")[0].split(" ")[0])
                if family == "Putative":
                    family = lastItem(ls[8].split("product=")[1].split(";")[0].split(" "))
                isDict[ID]["family"] = family
                isDict[ID]["genome"] = i


masterDict = defaultdict(lambda: defaultdict(list))
pseudoDir = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis")
for i in pseudoDir:
    if re.findall(r'functional.gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if not re.match(r'#', j):
                ls = j.rstrip().split("\t")
                ID = (ls[8].split("gbk_locus_tags=")[1])
                if ID in isDict.keys():
                    genome = isDict[ID]["genome"]
                    family = isDict[ID]["family"]
                    masterDict[genome]["functional"].append(family)

    elif re.findall(r'intact.gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if not re.match(r'#', j):
                ls = j.rstrip().split("\t")
                ID = (ls[8].split("locus_tag=")[1])
                if ID in isDict.keys():
                    genome = isDict[ID]["genome"]
                    family = isDict[ID]["family"]
                    masterDict[genome]["functional"].append(family)

    elif re.findall(r'pseudos.gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            counter = 0
            if not re.match(r'#', j):
                ls = j.rstrip().split("\t")
                if re.findall(r'old_locus_tag', ls[8]):
                    ID = (ls[8].split("old_locus_tag=")[1])
                    for k in ID.split(","):
                        if k in isDict.keys():
                            counter += 1
                            isID = k
                    if counter > 0:
                        genome = isDict[isID]["genome"]
                        family = isDict[isID]["family"]
                        masterDict[genome]["pseudo"].append(family)

                elif re.findall(r'gbk_locus_tag', ls[8]):
                    ID = (ls[8].split("gbk_locus_tags=")[1])
                    for k in ID.split(","):
                        if k in isDict.keys():
                            counter += 1
                            isID = k

                    if counter > 0:
                        genome = isDict[isID]["genome"]
                        family = isDict[isID]["family"]
                        masterDict[genome]["pseudo"].append(family)

endos = ["IS_Pectobacterium_wasabiae.gff", "IS_Symbiopectobacteriom_endolongispinus.gff", "IS_SyHa.gff", "IS_SyDa.gff", "IS_SyCt.gff", "IS_SyCl.gff",
         "IS_Sodalis_praecaptivus.gff", "IS_Sodalis_SOPE.gff", "IS_Sodalis_glossinidius.gff", "IS_Sodalis_endolongispinus.gff",
         "IS_Sodalis_TME1.gff", "IS_Sodalis_SCIS.gff", "IS_Sodalis_PSPU.gff", "IS_Sodalis_BTRI.gff"]

cats = ["functional", "pseudo"]

collapseDict = defaultdict(list)
totalDict = defaultdict(list)
Dict2 = defaultdict(lambda: defaultdict(list))
for i in endos:
    for j in cats:
        for k in masterDict[i][j]:
            collapseDict[k].append(k)
            combo = i + "_" + j
            Dict2[combo][k].append(k)
            totalDict[combo].append(k)


summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/IS/summary.csv", "w")
out.write("genome_status,fam,count,prop,total\n")
for i in Dict2.keys():
    for j in Dict2[i]:
        total = len(totalDict[i])
        famTotal = len(Dict2[i][j])
        famProp = famTotal/total
        out.write(i + "," + j + "," + str(famTotal) + "," + str(famProp) + "," + str(total) + "\n")
        summaryDict[i] = total

summaryDict2 = defaultdict(list)
for i in summaryDict.keys():
    summaryDict2[i.split(".")[0]].append(summaryDict[i])

geneDict = {
    "IS_Pectobacterium_wasabiae": 4648,
    "IS_Symbiopectobacteriom_endolongispinus": 5878,
    "IS_SyHa": 5550,
    "IS_SyDa": 3864,
    "IS_SyCt": 1853,
    "IS_SyCl": 304,
    "IS_Sodalis_praecaptivus": 4518,
    "IS_Sodalis_SOPE": 5714,
    "IS_Sodalis_glossinidius": 6116,
    "IS_Sodalis_endolongispinus": 4412,
    "IS_Sodalis_TME1": 4019,
    "IS_Sodalis_SCIS": 3311,
    "IS_Sodalis_PSPU": 2063,
    "IS_Sodalis_BTRI": 3145,
}

# GETTING PERCENTAGES FOR TRANSPOSASES PER GENOME
for i in summaryDict2.keys():
    print(i)
    print(sum(summaryDict2[i]))
    print( (sum(summaryDict2[i]) / geneDict[i])*100 )
    print("")
'''



# CRUDE BINNING OF SOD FROM THE ILLUMINA-ONLY ASSEMBLY (TO LOOK AT HOW MANY TRANSPOSASES WILL BE PROKKA-IDENTIFIED)
'''
counter = 0
total = 0
out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/illumina_only/sod.fasta", "w")
assembly = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/PLON_gammas/illumina_only/scaffolds.fasta")
assembly = fasta(assembly)
for i in assembly.keys():
    depth = float(i.split("cov_")[1])
    length = int(i.split("length_")[1].split("_cov")[0])
    if depth < 12 and depth > 6 and length > 1000:
        out.write(">" + i + "\n")
        out.write(assembly[i])
        print(i)
        total += length
        counter += 1
    elif depth > 25 and length > 1000:
        out.write(">" + i + "\n")
        out.write(assembly[i])
        total += length
        counter += 1
print(total)
print(counter)
'''



# LOOKING INTO PSEUDOGENIZATION REASONS IN TRANPSOSASES
'''
isDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
IS = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/IS/IS_Sodalis_SOPE.gff")
for i in IS:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        ID = (ls[8].split(";")[0].split("=")[1])
        isDict[ID] = ls

pseudo = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/Sodalis_SOPE_pseudos.gff")
for i in pseudo:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        if re.findall(r'gbk_locus_tags', ls[8]):
            ID = (ls[8].split("gbk_locus_tags=")[1])
            for j in ID.split(","):
                if j in isDict.keys():
                    print(ls)
'''





















