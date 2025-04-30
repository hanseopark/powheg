#!/bin/sh
date

echo "Executing run $2 on" `hostname` "in $PWD"

eval $(/cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@POWHEG::r3964-alice2-33)

echo "print POWHEG"

#export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/afs/cern.ch/work/h/hapark/work/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/LHAPDF/lhapdf6/share/LHAPDF
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current

echo "iseed $RANDOM" >> powheg.input

pwhg_main_dijet >&pwhg.log

#pwhg_main_hvq >&pwhg.log

date

export CONFIG_SEED=$RANDOM

/cvmfs/alice.cern.ch/bin/alienv unload VO_ALICE@POWHEG::r3964-alice2-33
eval $(/cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@O2Physics::daily-20241122-0000-1)

echo "Done load POWHEG and O2Physics"

#export PYTHIA8DATA=$ALICE_ROOT/PYTHIA8/pythia8/xmldoc/

root.exe -l -b -q 'RunPythia8.C++(100000)' > ana.log 2> ana.err

date

WORKDIR="/eos/user/h/hapark/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/pp13TeVdijet_CT18nlo/$2"

mkdir -p $WORKDIR

rsync -av --exclude='_condor_stdout' --exclude='_condor_stderr' --exclude='tmp' --exclude='var' --exclude='condor_exec.exe' ./* $WORKDIR
