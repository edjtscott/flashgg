export X509_USER_PROXY=~/x509up_u72495
fggRunJobs.py --load data_2018.json -d Data2018_15Nov19 --stage-to /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass2/2018/Data/Raw/ -x cmsRun legacy_dumper.py maxEvents=-1 runOnZee=False -q testmatch pujidWP=tight dumpJetSysTrees=False -n 200 --no-copy-proxy
