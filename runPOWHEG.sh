#!/bin/sh
date

echo "Executing run $2 on" `hostname` "in $PWD"

eval $(/cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@POWHEG::r3964-alice2-33)

echo "print POWHEG"

echo "remove pwgevents.lhe"
rm pwgevents.lhe

#export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/afs/cern.ch/work/h/hapark/work/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/LHAPDF/lhapdf6/share/LHAPDF
export LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH:/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current

### powheg parameter set on powheg.input
sed -i '/^iseed/d' powheg.input
echo "iseed $RANDOM" >> powheg.input

pwhg_main_dijet >&pwhg.log

#pwhg_main_hvq >&pwhg.log

date

export CONFIG_SEED=$RANDOM

/cvmfs/alice.cern.ch/bin/alienv unload VO_ALICE@POWHEG::r3964-alice2-33
eval $(/cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@O2Physics::daily-20241122-0000-1)

echo "Done load POWHEG and O2Physics"

#export PYTHIA8DATA=$ALICE_ROOT/PYTHIA8/pythia8/xmldoc/

### Set pythia
NEV=2000000000
root.exe -l -b -q "RunPythia8.C++($NEV)" > ana.log 2> ana.err

date
BASE_RESDIR="pp13TeVdijet_CT18nlo"
DATE_TAG=$(date +%y_%m_%d)
RESDIR="$BASE_RESDIR"
CHECK_PATH="/eos/user/h/hapark/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/$RESDIR" # output dir

if [ -d "$CHECK_PATH" ]; then
    i=0
    while true; do
        RESDIR="${BASE_RESDIR}_${DATE_TAG}_$i"
        CHECK_PATH="/eos/user/h/hapark/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/$RESDIR"
        if [ ! -d "$CHECK_PATH" ]; then
            break
        fi
        ((i++))
    done
fi
WORKDIR="/eos/user/h/hapark/PHD_Analysis/Run3/pp/13.6TeV/HfJets/PYTHIA/POWHEG/$RESDIR/$2"
echo "work dir: $WORKDIR"

mkdir -p $WORKDIR
rm *.so *.pcm *.d

rsync -av --exclude="*/" --exclude='_condor_stdout' --exclude='_condor_stderr' --exclude='tmp' --exclude='var' --exclude='condor_exec.exe' ./* $WORKDIR
