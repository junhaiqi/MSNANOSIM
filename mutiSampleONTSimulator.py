"""
    The main function of this script is to simulate negative samples or normal samples
for multi-sample nanopore sequencing data. 
    The simulation scheme of negative samples is mainly: 
        1) The adapter sequence and the barcode sequence are missing. 
        2) The absence of the barcode sequence. 
        3) Incomplete barcode sequence (More than 40 percent missing).
    The normal sequence: <-----------><-------------><--------------><-------------><----------->
                          topAdapter    topBarcode         read         endBarcode     endAdapter
"""
import random
import sys
import argparse
import os
if 'script' not in sys.path:
    sys.path.append('script')
from ex_info_from_read import getSeqsFromFasta as GSFF
from ex_info_from_read import writeSeqs2Fasta as S2F
from ex_info_from_read import get_seq_list, get_id_list
if 'module/simulateNanoSigs' not in sys.path:
    sys.path.append('module/simulateNanoSigs')
if 'module/simulateNanoSigs/module/' not in sys.path:
    sys.path.append('module/simulateNanoSigs')
from generatNoiseSignal import sequence_to_true_signal
from multiprocessing import Pool

import Bio.SeqIO
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities
MISBASELINE = 0.6

# ERRORMODELIST = [1, 2, 3]  # If there is a normal sequence, change to [0, 0,...,0, 1, 2, 3], 
                       # the number of '0' controls the probability of normal sequence occurrence.
# ERRORMODELIST = [0]  # If we only want the data that was successfully barcoded, ERRORMODELIST is equal to [0].

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join([complement_dict[base] for base in reversed(dna_sequence)])
    return reverse_complement_seq

def sim_read(perfect_read, error_model = 'nanopore2023', qscore_model = 'nanopore2023', identities = '95,5,99'):
    """Simulate error into perfect read using Badread error and qscore model

        read (str): perfect reads
    """
    output=sys.stderr
    mean, sd, maxi = [float(x) for x in identities.split(',')]
    identities = Identities(mean, sd, maxi, output)

    error_model = ErrorModel(error_model, output)
    qscore_model = QScoreModel(qscore_model, output)

    seq, quals, actual_identity, identity_by_qscores = \
                    sequence_fragment(perfect_read, identities.get_identity(), 
                                        error_model, qscore_model)
    return seq

def generateRefNegseq(topAdapter = "", \
                    topBarcode = "", \
                    read = "", \
                    endAdapter = "", \
                    endBarcode = "", \
                    errorMode = 1, \
                    libaryMode = 'ONTMUTISAMPLE'):
    """
    When the libaryMode == 'ONTMUTISAMPLE':
        The correct sequence: <-----------><-------------><--------------><-------------><----------->
                                topAdapter    topBarcode       read         endBarcode     endAdapter
        The main characteristic of the failed sequence (negative sample sequence) is the absence of barcode.
    """

    if libaryMode == 'ONTMUTISAMPLE':  # No error.
        if errorMode == 0:
            refRead = topAdapter + topBarcode + read + endBarcode + endAdapter
            return refRead
        elif errorMode == 1:  # The adapter sequence and the barcode sequence are missing.
            return read
        elif errorMode == 2:  # The absence of the barcode sequence. 
            refRead = topAdapter + read + endAdapter
            return refRead
        else:  # Incomplete barcode sequence.
            topRate = random.uniform(MISBASELINE, 1)
            incompTopBarcode = topBarcode[0:int(len(topBarcode)*(1-topRate))]

            endRate = random.uniform(MISBASELINE, 1)
            incompEndBarcode = endBarcode[0:int(len(endBarcode)*(1-endRate))]

            refRead = topAdapter + incompTopBarcode + read + incompEndBarcode + endAdapter
            return refRead
    else:
        pass

def getRefNegSeqs(inReadFile, inBarFile, inAdapterFile, outFile, ERRORMODELIST = [1, 2, 3], trueBarFile = 'trueBars.csv'):
    '''generate the negtive sequences.'''
    inReads = GSFF(inReadFile)
    inBars = GSFF(inBarFile)
    inAdapters = GSFF(inAdapterFile)
    outRefReads = []
    file = open(trueBarFile, 'w')
    file.write('top barcode,end barcode\n')
    for read in inReads:
        errorMode = random.choice(ERRORMODELIST)
        barcode = random.choice(inBars)
        reverBarcode = reverse_complement(barcode)
        refRead = generateRefNegseq(topAdapter = inAdapters[0], \
                    topBarcode = barcode, \
                    read = read, \
                    endBarcode = reverBarcode, \
                    endAdapter = inAdapters[1], \
                    errorMode = errorMode)
        outRefReads.append(refRead)
        file.write('%s,%s\n'%(barcode, reverBarcode))
    S2F(outRefReads, outFile)
    file.close()

def getSigBySeq(seq, index, output_folder, sigroot):
    sequence_to_true_signal((seq, index), output_folder = output_folder, sigroot = sigroot)

def getRefNegSigs(inRefNegFile, outSigsDir, sigRootName = 'timeSeries'):
    seq_list = get_seq_list(inRefNegFile)
    id_list = get_id_list(inRefNegFile)
    # zip_id_seq = list(zip(seq_list, id_list))
    outList = [outSigsDir for i in range(len(id_list))]
    rootList = [sigRootName for i in range(len(seq_list))]
    args = list(zip(seq_list, id_list, outList, rootList))
    isExists_out = os.path.exists(outSigsDir)
    if not isExists_out:
        os.makedirs(outSigsDir)

    pool = Pool(8)
    list(pool.starmap(getSigBySeq, args))
    pool.close()
    pool.join()
    # for seq in zip_id_seq:
    #     sequence_to_true_signal(seq, output_folder = outSigsDir, sigroot = sigRootName)

def getSimNegs(inRefNegFile, outSimNegFile):
    cmd = 'bash fast_sim_reads.sh %s %s'%(inRefNegFile, outSimNegFile)
    os.system(cmd)

def mainFunction(inReadFile, inBarFile, inAdapterFile, outRefFile, outSimRefFile, outSigsDir, sigRootName = 'timeSeries', 
                ERRORMODELIST = [1, 2, 3], trueBarFile = 'trueBars.csv'):
    getRefNegSeqs(inReadFile = inReadFile, inBarFile = inBarFile, inAdapterFile = inAdapterFile, 
                    outFile = outRefFile, ERRORMODELIST = ERRORMODELIST, trueBarFile = trueBarFile)

    getSimNegs(inRefNegFile = outRefFile, outSimNegFile = outSimRefFile)

    getRefNegSigs(inRefNegFile = outRefFile, outSigsDir = outSigsDir, sigRootName = sigRootName)

def get_parameters():
    """This script generates noiseless nanopore signals based on DNA sequences"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--read-file', type=str, required=True,
                        help='It is the fasta file that contains reads.')

    parser.add_argument('--barcode-file', type=str, required=True,
                        help='It is the fasta file that contains the barcode sequences.')

    parser.add_argument('--adapter-file', type=str, required=True,
                        help='It is the fasta file that contains the adapter sequences.')

    parser.add_argument('--ref-negtive-file', type=str, required=True,
                        help='It is the fasta file containing the final reference negative sequences.')

    parser.add_argument('--sim-negtive-file', type=str, required=True,
                        help='It is the fasta file containing the final simulated negative sequences.')

    parser.add_argument('--sim-sig-floder', type=str, required=True,
                        help='It is the folder path containing the nanopore signals simulated based on the 6-mer model.')
    
    parser.add_argument('--sim-mode', type=int, required=False, default = 0,
                        help='It is the mode of simulating data. If it is set to 0, it will only \
                            simulate the negative sequence, otherwise it will simulate the sequence that is successfully inserted into the barcode.')
    
    parser.add_argument('--true-bar-file', type=str, required=False, default = "trueBars.csv",
                        help='It is a file implying the real barcode carried by the simulated sequence (top barcode, end barcode)')                  
    
    parser.add_argument('--sig-root', type=str, required=False, default = "timeSeries",
                        help='It is the prefix to the filename of the nanopore signal being simulated.')

    args = parser.parse_args()

    return args

def main():

    args = get_parameters()
    if args.sim_mode == 0:
        BAR_ADA_ALLLOSE_LIST = [1 for i in range(60)]
        BAR_LOSE_LIST = [2 for i in range(20)]
        BAR_ABSENCE_LIST = [3 for i in range(20)]
        ERRORMODELIST = BAR_ADA_ALLLOSE_LIST + BAR_LOSE_LIST + BAR_ABSENCE_LIST
        random.shuffle(ERRORMODELIST)
    else:
        ERRORMODELIST = [0]

    mainFunction(inReadFile = args.read_file, \
                inBarFile = args.barcode_file, \
                inAdapterFile = args.adapter_file, \
                outRefFile = args.ref_negtive_file, \
                outSimRefFile = args.sim_negtive_file, \
                outSigsDir = args.sim_sig_floder, \
                ERRORMODELIST = ERRORMODELIST, \
                trueBarFile = args.true_bar_file, \
                sigRootName = args.sig_root)

if __name__ == "__main__":
   main()

    
