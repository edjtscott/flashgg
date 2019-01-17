# NB this command is specific to the configuration on lxplus and is not gaurenteed elsewhere
#outdir="/afs/cern.ch/work/s/sethzenz/ws/" # can't set absolute path on lsf because we're expecting to stage
queue="testmatch"
useAAA=0
version="reoptimisedClean"
fggRunJobs.py --load sig_jobs_2017.json -d sig_jobs_$version -x cmsRun workspaceStd.py maxEvents=-1 -n 500 -q $queue -D -P useAAA=$useAAA doHTXS=False doStage1=True doFiducial=False tthTagsOnly=False
