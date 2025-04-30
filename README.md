POWHEG.sub, this is the condor submission script, you can decide where to put the logs and how many jobs you want you want to submit.
RunPOWHEG.sh, the running scripts, you can decide where to put the result files (put them in your EOS)
powheg.input, the powheg parameters file. You can decide the number of events and PDF here.
RunPythia8.C, this is the shower MC macro. Here you can do all the physics.

condor_submit POWHEG.sub
