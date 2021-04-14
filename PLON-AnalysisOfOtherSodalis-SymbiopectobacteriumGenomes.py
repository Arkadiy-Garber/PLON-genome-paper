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


# SUMMARIZING PSEUDOHFINDER AND ANNOTATION (PROKKA/GHOSTKOALA) OUTPUT INTO BIOSYNTHETIC PATHWAYS

pseudoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
ko = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/summaries/sodalis_symbionts.sod.sym.refs.othersym.ghostkoala.csv")
for i in ko:
    ls = i.rstrip().split(",")
    if ls[0] != "ORF":
        pseudoDict[ls[0]]["ko"] = ls[1]
        pseudoDict[ls[0]]["cat"] = ls[2]
        pseudoDict[ls[0]]["fam"] = ls[3]
        pseudoDict[ls[0]]["sys"] = ls[4]
        pseudoDict[ls[0]]["gene"] = ls[5]
        pseudoDict[ls[0]]["taxa1"] = ls[6]
        pseudoDict[ls[0]]["taxa2"] = ls[7]
        pseudoDict[ls[0]]["taxa3"] = ls[8]


prots = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/sodalis_prots.sod.sym.refs.othersym.faa")
prots = fasta2(prots)

nameDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
names = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/sodalis_prots.sod.sym.refs.othersym.faa")
for i in names:
    if re.match(r'>', i):
        locus = i.rstrip()[1:].split(" ")[0]
        annotation = (allButTheFirst(i.rstrip(), " "))
        pseudoDict[locus]["prokka_annotation"] = annotation


sodalisMap = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
DIR = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis")
for i in DIR:
    if re.findall(r'pseudos.gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            ls = j.rstrip().split("\t")
            if not re.match(r'#', j):
                locus = (lastItem(ls[8].split(";")).split("=")[1])
                if len(locus.split(",")) == 1:
                    try:
                        contig = ls[0].split("|")[2]
                    except IndexError:
                        contig = ls[0]
                    start = ls[3]
                    end = ls[4]
                    strand = ls[6]
                    pseudoDict[locus]["contig"] = contig
                    pseudoDict[locus]["start"] = start
                    pseudoDict[locus]["end"] = end
                    pseudoDict[locus]["strand"] = strand
                    pseudoDict[locus]["pg"] = "y"
                else:
                    for k in locus.split(","):
                        try:
                            contig = ls[0].split("|")[2]
                        except IndexError:
                            contig = ls[0]
                        start = ls[3]
                        end = ls[4]
                        strand = ls[6]
                        pseudoDict[k]["contig"] = contig
                        pseudoDict[k]["start"] = start
                        pseudoDict[k]["end"] = end
                        pseudoDict[k]["strand"] = strand
                        pseudoDict[k]["pg"] = "y"

    elif re.findall(r'functional.gff', i) or re.findall(r'_intact.gff', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            ls = j.rstrip().split("\t")
            if not re.match(r'#', j):
                locus = (lastItem(ls[8].split(";")).split("=")[1])
                try:
                    contig = ls[0].split("|")[2]
                except IndexError:
                    contig = ls[0]

                start = ls[3]
                end = ls[4]
                strand = ls[6]
                pseudoDict[locus]["contig"] = contig
                pseudoDict[locus]["start"] = start
                pseudoDict[locus]["end"] = end
                pseudoDict[locus]["strand"] = strand
                pseudoDict[locus]["pg"] = "n"

    elif re.findall(r'blast', i):
        blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        blast = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in blast:
            ls = j.rstrip().split("\t")
            if ls[0] not in blastDict.keys():
                blastDict[ls[0]] = ls[12]
                pseudoDict[ls[0]]["blast"] = ls[12]

    elif re.findall('dnds', i):
        dnds = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s/dnds-summary.csv" % i)
        for j in dnds:
            ls = j.rstrip().split(",")
            if ls[0] != "ORF":
                pseudoDict[ls[0]]["homolog"] = ls[1]
                pseudoDict[ls[0]]["dn"] = ls[2]
                pseudoDict[ls[0]]["ds"] = ls[3]
                pseudoDict[ls[0]]["dnds"] = ls[4]

    elif re.findall(r'cds.fasta', i):
        if not re.findall(r'ref', i):
            print(i)
            file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
            for j in file:
                if re.match(r'>', j):
                    if not re.findall(r'symbiont', i):
                        sodalisMap[j.split(" ")[0][1:]] = i.split("_")[0] + "_" + i.split("_")[1]
                    else:
                        sodalisMap[j.split(" ")[0][1:]] = i.split("_")[0] + "_" + i.split("_")[1] + "_" + i.split("_")[2]

    elif re.findall(r'intergenic.fasta', i):
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if re.match(r'>', j):
                if not re.findall(r'symbiont', i):
                    sodalisMap[j.split(" ")[0][1:]] = i.split("_")[0] + "_" + i.split("_")[1]
                else:
                    sodalisMap[j.split(" ")[0][1:]] = i.split("_")[0] + "_" + i.split("_")[1] + "_" + i.split("_")[2]


MasterDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'n')))
for i in pseudoDict.keys():
    for j in pseudoDict[i]:
        try:
            MasterDict[sodalisMap[i]][i][j] = pseudoDict[i][j]
        except TypeError:
            print(i)

# DEFINING PATHWAYS
met = ["metE", "metF", "metH", "metR"]
bio = ["bioA", "bioB", "bioC", "bioD", "bioF"]
arg = ["argA", "argB", "argC", "argD", "argE", "argF", "argG", "argH","carA", "carB"]
om = ["bamA", "bamB", "bamC", "bamD", "bamE", "tamA", "tamB", "ompA", "ompL", "ompF", "ompC"]
im = ["secA", "secD", "secF", "secY", "secE", "secG", "yagC", "groEL", "groES", "tatA", "tatB", "tatC", "tatE"]
chaperones = ["secB", "ftsY", "ffh", "lepA", "lepB", "hlpA", "grpE", "dnaJ", "dnaK", "ibpA", "ibpB", "surA", "ompH"]
membrane = ["norM", "pit", "murJ", "mntH", "tatE", "tatC", "tatB", "tatA", "secF", "secD", "atpI", "atpH", "atpG",
            "atpF", "atpE", "atpD", "atpC", "atpB", "atpA", "nuoN", "nuoM", "nuoL", "ndh", "nuoJ", "nuoI", "nuoH",
            "nuoG", "nuoF", "nuoE", "nuoD", "nuoC", "nuoB", "nuoA", "cydB", "cydA"]

transport = ["norM", "pit", "murJ", "mntH", "mntH", "tatE", "tatC", "tatB", "tatA", "secF", "secD", "secY", "secM",
             "secG", "secE", "secB", "secA", "lepB", "lepA", "ffh", "ftsY", "tamB", "tamA", "bamE", "bamD", "bamC",
             "bamB", "bamA", "pal", "sufC", "mlaF", "mlaE", "mlaD", "mlaC", "mlaB", "lptG", "lptF", "lptE", "lptD",
             "lptC", "lptB", "lptA", "lolD", "lolC_E", "znuC", "znuB", "znuA", "ABC.ZM.S", "ABC.ZM.P", "ABC.ZM.A",
             "msbA", "cydD", "cydC", "manZ", "manY", "manX", "ptsH", "ptsI", "mscS"]

sodalis = ["Sodalis_praecaptivus", "Sodalis_SOPE",
             "Sodalis_glossinidius", "Sodalis_endolongispinus", "Sodalis_TME1", "Sodalis_SCIS", "Sodalis_PSPU", "Sodalis_PFLU",
             "Sodalis_HHAL", "Sodalis_BTRI", "Sodalis_symbiont_GEFVIR", "Sodalis_symbiont_DEMHIR",
             "Sodalis_symbiont_HETPER", "Sodalis_symbiont_MEPCIT", "Sodalis_symbiont_MEPMAR", "Sodalis_symbiont_TRABTM"]

symbiopecto = ["Pectobacterium_wasabiae", "Symbiopectobacteriom_endolongispinus",
               "Symbiopectobacteriom_Da", "Symbiopectobacteriom_Ha",
           "Symbiopectobacteriom_Ct", "Symbiopectobacteriom_Cl",]

symbionts = ["Pectobacterium_wasabiae", "Symbiopectobacteriom_endolongispinus", "Sodalis_praecaptivus", "Sodalis_SOPE",
             "Sodalis_glossinidius", "Sodalis_endolongispinus", "Sodalis_TME1", "Sodalis_SCIS", "Sodalis_PSPU", "Sodalis_PFLU",
             "Sodalis_HHAL", "Sodalis_BTRI", "Sodalis_symbiont_GEFVIR", "Sodalis_symbiont_DEMHIR",
             "Sodalis_symbiont_HETPER", "Sodalis_symbiont_MEPCIT", "Sodalis_symbiont_MEPMAR", "Sodalis_symbiont_TRABTM"]
peptidoglycan = ["murA", "murB", "murC", "murD", "murE", "murF", "murG", "murH", "murI", "murJ", "glmM", "glmS",
                  "glmU", "amiD", "dapF", "ddl", "ftsW", "glyA", "lpoB", "mepM", "metC", "mltB", "mltC", "mraY",
                  "nlpD", "pbp"]

# WRITING SUMMARY FILE (PEPTIDOGLYCAN CAN BE SWITCHED OUT WITH OTHER PATHWAYS IN THIS SCRIPT)
out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/peptidoglycan.csv", "w")
out.write("genome,orf,gene,status\n")
for i in symbiopecto:
    print(i)
    for k in peptidoglycan:
        count = 0
        for j in MasterDict[i]:
            if re.findall(k, str(MasterDict[i][j]["gene"])):
                print(j + "\t\t" + str(MasterDict[i][j]["gene"]) + "\t\t" + str(MasterDict[i][j]["pg"]))
                out.write(i + "," + j + "," + str(MasterDict[i][j]["gene"]) + "," + str(MasterDict[i][j]["pg"]) + "\n")
                count += 1
        if count == 0:
            print("-" + "\t\t" + k + "\t\t" + "-")
            out.write(i + "," + "-" + "," + k + "," + "-" + "\n")
    print("")
    out.write("################################################################\n")

out.close()



# SUMMARIZING VARIOUS STATS ON THE ENDOSYMBIONT GENOMES (E.G. DS, DN/DS, # OF PSEUDOGENES, ETC.)
'''
pseudoDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
intactDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
DIR = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis")
for i in DIR:
    try:
        if not re.findall(r'symbiont_', i):
            key = i.split("_")[0] + "_" + i.split("_")[1]
        else:
            key = i.split("_")[0] + "_" + i.split("_")[1] + "_" + i.split("_")[2]
    except IndexError:
        pass

    if re.findall(r'functional.gff', i) or re.findall(r'_intact.gff', i):
        totalLength = 0
        intactLength = 0
        intactGenes = 0
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if not re.match(r'#', j):
                ls = j.rstrip().split("\t")
                locus = (lastItem(ls[8].split(";")).split("=")[1])
                intactDict[locus] = 'intact'

                length = int(ls[4]) - int(ls[3])
                intactLength += length
                intactGenes += 1

            else:
                if re.findall(r'region', j):
                    length = int(j.rstrip().split(" ")[3])
                    totalLength += length
        percLength = (totalLength/5159425)
        codingDensity = (intactLength/totalLength)

        pseudoDict[key]["codingDensity"] = codingDensity
        pseudoDict[key]["percLength"] = percLength
        pseudoDict[key]["intactGenes"] = intactGenes

    elif re.findall(r'cds.fasta', i):
        if not re.findall(r'ref', i):
            file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
            file = fasta(file)
            totalGenes = len(file.keys())
            pseudoDict[key]["totalGenes"] = totalGenes

    else:
        pass


for i in DIR:
    try:
        if not re.findall(r'symbiont_', i):
            key = i.split("_")[0] + "_" + i.split("_")[1]
        else:
            key = i.split("_")[0] + "_" + i.split("_")[1] + "_" + i.split("_")[2]
    except IndexError:
        pass

    if re.findall('dnds', i):
        averageDNDS = []
        intactDNDS = []
        pseudoDNDS = []
        averageDS = []
        intactDS = []
        pseudoDS = []
        dnds = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s/dnds-summary.csv" % i)
        for j in dnds:
            ls = j.rstrip().split(",")
            if ls[0] != "ORF":
                dnds = float(ls[4])
                ds = float(ls[3])

                orf = ls[0]
                averageDNDS.append(dnds)
                averageDS.append(ds)

                if orf in intactDict.keys():
                    intactDNDS.append(dnds)
                    intactDS.append(ds)
                else:
                    pseudoDNDS.append(dnds)
                    pseudoDS.append(ds)

        pseudoDict[key]["averageDNDS"] = statistics.mean(averageDNDS)
        pseudoDict[key]["intactDNDS"] = statistics.mean(intactDNDS)
        pseudoDict[key]["pseudoDNDS"] = statistics.mean(pseudoDNDS)
        pseudoDict[key]["averageDS"] = statistics.mean(averageDS)
        pseudoDict[key]["intactDS"] = statistics.mean(intactDS)
        pseudoDict[key]["pseudoDS"] = statistics.mean(pseudoDS)


for i in pseudoDict:
    percPseudo = (pseudoDict[i]["intactGenes"] / pseudoDict[i]["totalGenes"])
    percPseudo = 1 - percPseudo
    pseudoDict[i]["percPseudo"] = percPseudo

sodalis = ["Sodalis_praecaptivus", "Sodalis_SOPE",
             "Sodalis_glossinidius", "Sodalis_endolongispinus", "Sodalis_TME1", "Sodalis_SCIS", "Sodalis_PSPU", "Sodalis_PFLU",
             "Sodalis_HHAL", "Sodalis_BTRI", "Sodalis_symbiont_GEFVIR", "Sodalis_symbiont_DEMHIR",
             "Sodalis_symbiont_HETPER", "Sodalis_symbiont_MEPCIT", "Sodalis_symbiont_MEPMAR", "Sodalis_symbiont_TRABTM"]

symbiopecto = ["Pectobacterium_wasabiae", "Symbiopectobacteriom_endolongispinus",
               "Symbiopectobacteriom_Da", "Symbiopectobacteriom_Ha",
           "Symbiopectobacteriom_Ct", "Symbiopectobacteriom_Cl",]

out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/barplot-2.csv", "w")
out.write("genome,metric,value\n")
for i in symbiopecto:
    print(i)
    print(i + "\tpercLength\t" + str(pseudoDict[i]["percLength"]))
    print(i + "\tcodingDensity\t" + str(pseudoDict[i]["codingDensity"]))
    print(i + "\tpercPseudo\t" + str(pseudoDict[i]["percPseudo"]))
    print(i + "\tintactDS\t" + str(pseudoDict[i]["intactDS"]))
    print(i + "\tpseudoDS\t" + str(pseudoDict[i]["pseudoDS"]))
    print(i + "\tintactDNDS\t" + str(pseudoDict[i]["intactDNDS"]))
    print(i + "\tpseudoDNDS\t" + str(pseudoDict[i]["pseudoDNDS"]))
    print("")

    out.write(i + ",percLength," + str(pseudoDict[i]["percLength"]) + "\n")
    out.write(i + ",codingDensity," + str(pseudoDict[i]["codingDensity"]) + "\n")
    out.write(i + ",percPseudo," + str(pseudoDict[i]["percPseudo"]) + "\n")
    out.write(i + ",intactDS," + str(pseudoDict[i]["intactDS"]) + "\n")
    out.write(i + ",pseudoDS," + str(pseudoDict[i]["pseudoDS"]) + "\n")
    out.write(i + ",intactDNDS," + str(pseudoDict[i]["intactDNDS"]) + "\n")
    out.write(i + ",pseudoDNDS," + str(pseudoDict[i]["pseudoDNDS"]) + "\n")
'''



# EXPLORING DATA WITH RESPECT TO PSEUDOGENES AND TOTAL NUMBERS OF PREDICTED GENES
'''
sodalisMap = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
DIR = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis")
for i in DIR:
    if re.findall(r'functional.gff', i):
        print(i)
        total = 0
        coding = 0
        count = 0
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if re.findall(r'sequence-region', j):
                length = int(j.rstrip().split(" ")[3])
                total += length
            else:
                if not re.match(r'#', j):
                    ls = (j.rstrip().split("\t"))
                    length = int(ls[4]) - int(ls[3])
                    coding += length
                    count += 1

        print(count)
        print("")
    if re.findall(r'pseudos.gff', i):
        pseudos = 0
        file = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/OtherSodalis/%s" % i)
        for j in file:
            if not re.match(r'#', j):
                ls = (j.rstrip().split("\t"))
                pseudos += 1

        print(pseudos)
        print(count)
        print("")
'''

# CALCULTING ANI FOR SYMBIOPECTOBACTERIUM GENOMES
# THE masterDir DIRECTORY HAS MANY FOLDERS. SOME OF THESE, WITH 'prokka' IN THEIR NAMES, HAVE BLAST OUTPUT FILES.
# THESE BLAST OUTPUT FILES ARE THE RESULTS OF BLAST SEARCH OF EACH SYMBIOPECTOBACTERIUM ENDOSYMBIONT AGAINST ALL OTHER SYMBIOPECTOBACTERIUM ENDOSYMBIONT
'''
#PAIRWISE ANI
blastDict = defaultdict(lambda: defaultdict(list))
masterDir = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/contigs")
for folder in masterDir:
    if re.findall(r'prokka', folder):
        DIR = os.listdir("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/contigs/%s" % folder)
        for i in DIR:
            if re.findall(r'blast', i):
                blast = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/contigs/%s/%s" % (folder, i))
                for j in blast:
                    ls = j.rstrip().split("\t")
                    blastDict[folder][i].append(float(ls[2]))

out = open("/Users/arkadiygarber/Desktop/Ongoing_Research_Projects/Endosymbionts/otherPecto/ANI.csv", "w")
for i in blastDict.keys():
    for j in blastDict[i]:
        print(i + "\t\t" + j + "\t\t" + str(statistics.mean(blastDict[i][j])))
        out.write(i + "," + j + "," + str(statistics.mean(blastDict[i][j])) + "\n")
    print("")
out.close()
'''



















