
"""Extract barcode, read and adpter from the source file."""

import os

SOURCEFILEPATH = "../data/ref350.fasta"
def getSeqsFromFasta(file = SOURCEFILEPATH):
    _file = open(file)
    seqList = []
    for line in _file:
        if line[0] != '>':
            seq = line.strip('\n')
            seqList.append(seq)
    _file.close()
    return seqList

def writeSeqs2Fasta(seqList, filePath):
    file = open(filePath, 'w')
    t = 0
    for seq in seqList:
        file.write('>%d\n'%t)
        file.write('%s\n'%seq)
        t += 1
    file.close()

def main():
    seqs = getSeqsFromFasta()
    reads, barcodes = [], []
    # top adater and end adapter.
    topAdapter = seqs[0][0:61]
    endAdapter = seqs[0][len(seqs[0])-61:]
    adapters = [topAdapter, endAdapter]
    # reads and barcodes
    for seq in seqs:
        read = seq[101:len(seq)-101]
        reads.append(read)
        topBarcode = seq[61:101]
        if topBarcode not in barcodes:
            barcodes.append(topBarcode)
        # endBarcode = seq[len(seq)-101:len(seq)-61]
        # if (topBarcode, endBarcode) not in barcodes:
        #     barcodes.append((topBarcode, endBarcode))

        # print((topBarcode, endBarcode))
        # break
    
    readFile = '../data/reads.fasta'
    barcodeFile = '../data/barcodes.fasta'
    adapterFile = '../data/adapter.fasta'

    writeSeqs2Fasta(reads, readFile)
    writeSeqs2Fasta(barcodes, barcodeFile)
    writeSeqs2Fasta(adapters, adapterFile)
    
def get_seq_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x != '', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    a = list(seq_list)
    return a


def get_id_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: x.split('|')[0][1:], lines)
    b = list(id_list)
    return b


if __name__ == "__main__":
    # main()

    seqs = getSeqsFromFasta()[0:10000]
    writeSeqs2Fasta(seqs, 'testGUPPY.fasta')


