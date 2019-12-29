export X509_USER_PROXY=~/x509up_u243093
fggRunJobs.py --load data_2016.json -d Data2016_17Dec19 --stage-to /vols/cms/es811/HggGeneral/WorkspaceTest/Pass1/2016/Data/Raw -x cmsRun workspaceStd.py maxEvents=-1 -q hepmedium.q -n 50 --no-copy-proxy dumpWorkspace=True doStageOne=True  copyInputMicroAOD=True
