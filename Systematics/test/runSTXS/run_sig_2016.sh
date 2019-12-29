export X509_USER_PROXY=~/x509up_u243093
fggRunJobs.py --load sig_2016.json -d Sig2016_17Dec19 --stage-to /vols/cms/es811/HggGeneral/WorkspaceTest/Pass1/2016/Sig/Raw -x cmsRun workspaceStd.py maxEvents=-1 -q hepmedium.q -n 50 --no-copy-proxy  dumpWorkspace=True doStageOne=True doSystematics=True useParentDataset=True
