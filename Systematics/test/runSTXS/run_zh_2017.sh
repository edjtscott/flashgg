export X509_USER_PROXY=~/x509up_u243093
fggRunJobs.py --load zh_2017.json -d ZH2017_17Dec19 --stage-to /vols/cms/es811/HggGeneral/WorkspaceTest/Pass1/2017/ZH/Raw -x cmsRun workspaceStd.py maxEvents=-1 -q hepmedium.q -n 50 --no-copy-proxy  dumpWorkspace=True doStageOne=True doSystematics=True useParentDataset=True copyInputMicroAOD=True
