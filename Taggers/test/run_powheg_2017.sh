export X509_USER_PROXY=~/x509up_u72495
fggRunJobs.py --load powheg_2017.json -d Powheg2017_19Oct19 --stage-to /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass1/2017/Powheg/Raw/ -x cmsRun legacy_dumper.py maxEvents=-1 runOnZee=False -q testmatch pujidWP=tight dumpJetSysTrees=False -n 200 --no-copy-proxy
