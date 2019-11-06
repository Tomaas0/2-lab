from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils import GC

def doesQualityIsInSequence(qual, seq):
    x = True
    for c in qual:
        if (c not in seq):
            x = False
    return x

result = open("result.txt", "w")

#Nuskaitymas
file = open("reads_for_analysis.fastq")

sequences = []
qualities = []

for title, seq, qual in FastqGeneralIterator(file):
    sequences.append(seq)
    qualities.append(qual)

#Koduotės radimas
sanger = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"
solexa = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
ilumina13 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
ilumina15 = "CDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
ilumina18 = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"

isSanger = True
isSolexa = True
isIlumina13 = True
isIlumina15 = True
isIlumina18 = True

for qual in qualities:
    if (isSanger == True):
        if (not doesQualityIsInSequence(qual, sanger)):
            isSanger = False
    if (isSolexa == True):
        if (not doesQualityIsInSequence(qual, solexa)):
            isSolexa = False
    if (isIlumina13 == True):
        if (not doesQualityIsInSequence(qual, ilumina13)):
            isIlumina13 = False
    if (isIlumina15 == True):
        if (not doesQualityIsInSequence(qual, ilumina15)):
            isIlumina15 = False
    if (isIlumina18 == True):
        if (not doesQualityIsInSequence(qual, ilumina18)):
            isIlumina18 = False

result.write("Is it Sanger: ")
if (isSanger):
    result.write("Yes\n")
else:
    result.write("No\n")
    
result.write("Is it Solexa: ")
if (isSolexa):
    result.write("Yes\n")
else:
    result.write("No\n")
    
result.write("Is it Ilumina1.3: ")
if (isIlumina13):
    result.write("Yes\n")
else:
    result.write("No\n")
    
result.write("Is it Ilumina1.5: ")
if (isIlumina15):
    result.write("Yes\n")
else:
    result.write("No\n")
    
result.write("Is it Ilumina1.8: ")
if (isIlumina18):
    result.write("Yes\n")
else:
    result.write("No\n")

#Grafikas
gcsFile = open("gc.txt", "w")
GCs = []
for seq in sequences:
    gc = GC(seq)
    GCs.append(gc)
    gcsFile.write(str(gc))
    gcsFile.write("\n")

result.close()