from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def doesQualityIsInSequence(qual, seq):
    x = True
    for c in qual:
        if (c not in seq):
            x = False
    return x

result = open("result.txt", "w")

#Nuskaitymas
file = open("reads_for_analysis.fastq")

titles = []
sequences = []
qualities = []

for title, seq, qual in FastqGeneralIterator(file):
    titles.append(title)
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
GCs = []
for seq in sequences:
    gc = GC(seq)
    GCs.append(gc)
plt.hist(GCs, bins=300)
plt.show()

#Finding 5 elements from each peak
'''THIS CODE WAS USED TO FIND THE BIGGEST PEAK
maxCount = 0
maxi = 0
for i in range(1, 100, 1):
    _min = i - 0.5
    _max = i + 0.5
    count = 0
    for gc in GCs:
        if (gc > _min) and (gc < _max):
            count = count + 1
    if (count > maxCount):
        maxCount = count
        maxi = i

print(str(maxCount))
print("\n")
print(str(maxi))
'''

ids = []

#Peak1
_min = 35
_max = 36
count = 0
found = 0
i = 0
for gc in GCs:
    if (gc > _min) and (gc < _max):
        count = count + 1
        if (found < 5):
            ids.append(i)
            found = found + 1
    i = i + 1

#Peak2
_min = 53
_max = 54
count = 0
found = 0
i = 0
for gc in GCs:
    if (gc > _min) and (gc < _max):
        count = count + 1
        if (found < 5):
            ids.append(i)
            found = found + 1
    i = i + 1

#Peak3
_min = 69.5
_max = 70.5
count = 0
found = 0
i = 0
for gc in GCs:
    if (gc > _min) and (gc < _max):
        count = count + 1
        if (found < 5):
            ids.append(i)
            found = found + 1
    i = i + 1

#Blast
for i in ids:
    print("Going to Blast")
    print("Id: " + str(i))
    result_handle = NCBIWWW.qblast("blastn", "nt", sequences[i], hitlist_size=1, megablast=True, entrez_query='bacteria')
    print("Just got home")
    res = result_handle.read()
    f = open("blast" + str(i) + ".xml", "w")
    f.write(res)
    f.close()

result.close()