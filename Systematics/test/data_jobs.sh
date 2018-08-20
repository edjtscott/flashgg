queue="8nh"
LM=${CMSSW_BASE}/src/flashgg/MetaData/work/jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
version="fullNewTest"
fggRunJobs.py --load data_jobs_Legacy16_DoubleEG.json -d data_jobs_${version} -x cmsRun workspaceStd.py maxEvents=-1 -n 100 -q ${queue} -D -P useAAA=0 doStage1=True doFiducial=False tthTagsOnly=False lumiMask=${LM}
