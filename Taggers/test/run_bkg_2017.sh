export X509_USER_PROXY=~/x509up_u72495
fggRunJobs.py --load bkg_2017.json -d Bkg2017_16Sep19 --stage-to /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass0/2017/Bkg/ -x cmsRun legacy_dumper.py maxEvents=-1 runOnZee=False -q testmatch pujidWP=tight dumpJetSysTrees=True -n 200 --no-copy-proxy
