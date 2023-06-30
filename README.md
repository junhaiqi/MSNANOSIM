## Overview
The main functions of this project are:
* Based on the badread model (https://github.com/rrwick/Badread), simulate the nanopore sequencing data of multiple samples with barcode.
* Based on the badread model, simulate the nanopore signal of the DNA sequence and the corresponding basecalled sequence.

For any questions about this project, please contact the developer via the following email: 201911865@mail.sdu.edu.cn.
## Requirements
The main scripts runs on MacOS and Linux. It requires [Python](https://www.python.org/) 3.7 or later.

## Installation
You can install some necessary packages using pip:
```bash
pip install module/simulateNanoSeqs/Badread/
pip install h5py
pip install biopython
```

## Quick usage
You can use the following commands to quickly generate simulated multi-sample nanopore sequencing data:

```bash
python mutiSampleONTSimulator.py --read-file data/reads.fasta --barcode-file data/barcodes.fasta --adapter-file data/adapter.fasta --ref-negtive-file ref_barcoded_file.fasta --sim-negtive-file sim_barcoded_file.fasta --sim-sig-floder barcoded_sigs --sim-mode 1 --true-bar-file trueBars.csv
```
For a detailed explanation of some parameters, you can get it with the following command:

```bash
python mutiSampleONTSimulator.py --h
```

You can use the following commands to quickly generate simulated nanopore sequencing sequences:

```bash
python module/simulateNanoSeqs/ont_error_simulator.py --read-file data/reads.fasta --sim-file reads_error.fasta --badread-error-model nanopore2023 --badread-qscore-model nanopore2023 --badread-identity 95,5,99 --thread 1
```

For a detailed explanation of some parameters, you can get it with the following command:

```bash
python module/simulateNanoSeqs/ont_error_simulator.py --h
```

You can use the following commands to quickly generate simulated nanopore signals:
```bash
python module/simulateNanoSigs/generatNoiseSignal.py -i data/reads.fasta -o sigs -rootName timeSeries -seed 0
```

For a detailed explanation of some parameters, you can get it with the following command:

```bash
python module/simulateNanoSigs/generatNoiseSignal.py --h
```

We provide a shell script to quickly generate nanopore signals from DNA sequences:
```bash
bash fast_sim_6mer_sigs.sh data/reads.fasta Sigs
```
Here, "data/reads.fasta" is the file containing some reference sequences and "Sigs" is the folder path containing the corresponding simulated nanopore signals.

We also provide a shell script to quickly generate simulated nanopore sequencing sequences of reference DNA sequences:
```bash
bash fast_sim_reads.sh data/reads.fasta sim_reads.fasta
```
Here, "data/reads.fasta" is the file containing some reference sequences and "sim_reads.fasta" is the file path containing the corresponding simulated nanopore sequences.

## Detailed usage
For "mutiSampleONTSimulator.py", the specific parameters are as follows:
   | Parameters   | Description |
   |  :----:  | :----:  |
   | --h  | show this help message and exit |
   | --read-file  | It is the fasta file that contains reference reads. |
   | --barcode-file  | It is the fasta file that contains the barcode sequences. |
   | --adapter-file  | It is the fasta file that contains the adapter sequences. |
   | --ref-negtive-file  | It is the fasta file containing the final reference negative sequences. |
   | --sim-negtive-file  | It is the fasta file containing the final simulated negative sequences. |
   | --sim-sig-floder  | It is the folder path containing the nanopore signals simulated based on the 6-mer model. |
   | --sim-mode  | It is the mode of simulating data. If it is set to 0, it will only simulate the negative sequence, otherwise it will simulate the sequence that is successfully inserted into the barcode. |
   | --true-bar-file  | It is a file implying the real barcode carried by the simulated sequence (top barcode, end barcode) |
   | --sig-root   | It is the prefix to the filename of the nanopore signal being simulated. |

For "module/simulateNanoSeqs/ont_error_simulator.py", the specific parameters are as follows:
   | Parameters   | Description |
   |  :----:  | :----:  |
   | --h  | show this help message and exit |
   | --read-file  | It is the fasta file that contains reference reads. |
   | --sim-file  | It is the fasta file that contains the simulated nanopore sequences. |
   | --badread-error-model  | Error model from babread. |
   | --badread-qscore-model  | Qscore model from babread. |
   | --badread-identity  | Identity/accuracy (in percentage) parameter pass to badread: format mean,st,max. |
   | --thread  | Specifies the number of threads to use. |

For "--badread-error-model", "--badread-qscore-model", and "--badread-identity", users can refer to the parameter introduction of Badread (https://github.com/rrwick/Badread).

For "module/simulateNanoSigs/generatNoiseSignal.py", the specific parameters are as follows:
| Parameters   | Description |
   |  :----:  | :----:  |
   | --h  | show this help message and exit |
   | -i  | It is the fasta file that contains reference reads. |
   | -o  | output folder, the output simulated nanopore signals in here. |
   | -rootName  | The root name of simulated 6-mer nanpore signal files. |
   | -seed  | Random seed used to control the noise of the simulation. |

## License
No license.

   




