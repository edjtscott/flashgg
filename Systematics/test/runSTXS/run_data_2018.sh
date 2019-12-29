export X509_USER_PROXY=~/x509up_u243093
fggRunJobs.py --load data_2018.json -d Data2018_17Dec19 --stage-to /vols/cms/es811/HggGeneral/WorkspaceTest/Pass1/2018/Data/Raw -x cmsRun workspaceStd.py maxEvents=-1 -q hepmedium.q -n 50 --no-copy-proxy dumpWorkspace=True doStageOne=True  copyInputMicroAOD=True
