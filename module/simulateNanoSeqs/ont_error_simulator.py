
import sys
from multiprocessing import Pool

import Bio.SeqIO
from badread.simulate import sequence_fragment, ErrorModel, QScoreModel, Identities

import argparse

def sim_reads(perfect_read_list, error_model = 'nanopore2023', qscore_model = 'nanopore2023', identities = '95,5,99', threadNum = 8):
    """Simulate error into perfect read using Badread error and qscore model

        read (str): perfect reads
    """
    output = sys.stderr
    mean, sd, maxi = [float(x) for x in identities.split(',')]
    identities = Identities(mean, sd, maxi, output)

    error_model = ErrorModel(error_model, output)
    qscore_model = QScoreModel(qscore_model, output)

    args = [(perfect_read, identities.get_identity(), error_model, qscore_model) for perfect_read in perfect_read_list]
    # seq, quals, actual_identity, identity_by_qscores = \
    #                 sequence_fragment(perfect_read, identities.get_identity(), 
    #                                     error_model, qscore_model)
    pool = Pool(threadNum)
    simList = list(pool.starmap(sequence_fragment, args))
    pool.close()
    pool.join()
    simSeqList = [item[0] for item in simList]
    return simSeqList

def smartSim(perfect_read):
    return sim_read(perfect_read = perfect_read)
    

def mutiSim(seqList = [], threadNum = 8):

    pool = Pool(threadNum)
    simList = list(pool.imap(smartSim, seqList))
    pool.close()
    pool.join()
    return simList

def getSeqsFromFasta(file = ""):
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

def mainFunction(refFile, simFile, errorModel, qsModel, identities, threadNum):
    seqList = getSeqsFromFasta(file = refFile)
    simList = sim_reads(perfect_read_list = seqList, 
                        error_model = errorModel, 
                        qscore_model = qsModel, 
                        identities = identities, 
                        threadNum = threadNum)
    writeSeqs2Fasta(simList, simFile)


def get_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read-file', type=str, required=True,
                        help='It is the fasta file that contains reads.')

    parser.add_argument('--sim-file', type=str, required=True,
                        help='It is the fasta file that contains the simulated nanopore sequences.')

    parser.add_argument('--badread-error-model', type=str, required=False, default='nanopore2023',
                        help='Error model from babread.')

    parser.add_argument('--badread-qscore-model', type=str, default='nanopore2023',
                        help='Qscore model from babread.')

    parser.add_argument('--badread-identity', type=str, default='95,5,99',
                        help='Identity/accuracy (in percentage) parameter pass to badread: format mean,st,max.')

    parser.add_argument('--thread', type=int, required=False, default=8,
                        help='Specifies the number of threads to use.')

    args = parser.parse_args()

    return args

def main():
    args = get_parameters()
    mainFunction(refFile = args.read_file, \
                simFile = args.sim_file, \
                errorModel = args.badread_error_model, \
                qsModel = args.badread_qscore_model, \
                identities = args.badread_identity, \
                threadNum = args.thread)


if __name__ == "__main__":
    main()