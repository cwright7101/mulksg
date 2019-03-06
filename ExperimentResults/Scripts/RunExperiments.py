#!/usr/bin/env python
#For keeran: 
# module load gcc-8.2.0-gcc-8.2.0-yttbr4w
# module load cmake-3.10.3-gcc-8.2.0-ccmnmdj
# module load boost-1.68.0-gcc-8.2.0-rpfwdaj
# module load openmpi-2.1.1-gcc-8.2.0-pvznauq #the openmpi-3.1.3-gcc-8.2.0-7adnbag does not work for some reason, not sure why
# pip install mpi4py
from subprocess import call
import subprocess
import gzip
import os
import sys
import tarfile
from optparse import OptionParser
import time
import datetime

def run_tests(options, args):
  print("Running tests...")
  d = options.datapath
  children = [os.path.join(options.datapath, child) for child in os.listdir(options.datapath)]
  testfolders= filter(os.path.isdir, children)

  now = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M")
  logfilename = options.outdir + "/timing_" +now + ".log"
  log = open(logfilename, "w")
  for test in testfolders:
    frag1 = test + "/raw/frag_1.fastq"
    frag2 = test + "/raw/frag_2.fastq"
    reference = test + "/reference.fa"
    if os.path.isfile(frag1) and os.path.isfile(frag2):
      print("Running tests for ", test)           

      spades_outfolder = options.outdir + "/spades_"+test.split("/")[-1]
      mulksg_outfolder = options.outdir + "/mulksg_"+test.split("/")[-1]
      #SPADES FIRST
      print("Running SPAdes-3.13.0, Putting Results in: ", spades_outfolder)
      print ("frag1: ", frag1, "\nfrag2:", frag2)
      spades_start = time.time()               
      p = subprocess.Popen(["/usr/bin/time", "-v","./SPAdes-3.13.0/spades.py","-1", frag1, "-2", frag2, "-o", 
            spades_outfolder,"-t", options.numthreads, "--careful", "-k", "21,33,45,57,69,81,93"], 
            stderr = log)
      p.wait()
      spades_end = time.time()
      log.write("\nRunning time for " + spades_outfolder+ " is: " + str(spades_end - spades_start) + " seconds")

      #MULKSG NOW
      print("Running MULKSG, Putting Results in: ", mulksg_outfolder)
      mulksg_start = time.time()    
      #,"mpirun", "-n", options.numprocs,           
      p = subprocess.Popen(["/usr/bin/time", "-v","./SPAdes-3.13.0/mulksg.py","-1", 
            frag1, "-2", frag2, "-o", mulksg_outfolder,"-t", options.numthreads, "--careful", "-k", "21,33,45,57,69,81,93"],
            stderr = log)
      p.wait()
      mulksg_end = time.time()
      log.write("\nRunning time for " + mulksg_outfolder+ " is: " + str(mulksg_end - mulksg_start) + " seconds")
      # Quast now
      if os.path.isfile(reference):
        call(["./quast-5.0.2/quast.py", "-R", reference, spades_outfolder+"/contigs.fasta", spades_outfolder+"/scaffolds.fasta",
            mulksg_outfolder+"/contigs.fasta", mulksg_outfolder+"/scaffolds.fasta",
            "-l", '"SPAdes_contigs, SPAdes_scaffolds, MULKSG_contigs, MULKSG_scaffolds"',
            "-o", options.outdir +"/quast_results_"+test.split("/")[-1]])
      else:
        call(["./quast-5.0.2/quast.py", spades_outfolder+"/contigs.fasta", spades_outfolder+"/scaffolds.fasta",
            mulksg_outfolder+"/contigs.fasta", mulksg_outfolder+"/scaffolds.fasta",
            "-l", '"SPAdes_contigs, SPAdes_scaffolds, MULKSG_contigs, MULKSG_scaffolds"',
            "-o", options.outdir +"/quast_results_"+test.split("/")[-1]])
  log.close()

def main(options, args):
  if not os.path.exists(options.datapath):
    print("Datapath of where to find the tests doesn't exist, exiting")
    sys.exit()

  if not os.path.exists(options.outdir):
    os.mkdir(options.outdir)

  run_tests(options, args)

if __name__ == '__main__':
  parser = OptionParser()
  parser.add_option("-p", "--datapath", dest="datapath", 
                    default="./data", help="path of where datasets are located, default is ./data")
  parser.add_option("-o", "--outdir", dest="outdir", 
                    default="./results", help="output directory, default is results")
  parser.add_option("-n", "--numprocs", dest="numprocs", 
                    default="1", help="number of mpi processes")
  parser.add_option("-t", "--numthreads", dest="numthreads", 
                    default="16", help="number of threads per process")

  (options, args) = parser.parse_args()

  main(options, args)
  print("")