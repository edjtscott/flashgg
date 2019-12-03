export X509_USER_PROXY=~/x509up_u72495
fggRunJobs.py --load sig_2017.json -d Sig2017_8Nov19 --stage-to /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass2/2017/Sig/Raw/ -x cmsRun legacy_dumper.py maxEvents=-1 runOnZee=False -q testmatch pujidWP=tight dumpJetSysTrees=False useParentDataset=True -n 200 --no-copy-proxy
