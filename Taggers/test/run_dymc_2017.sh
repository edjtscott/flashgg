export X509_USER_PROXY=~/x509up_u72495
fggRunJobs.py --load dy_2017.json -d DYMC2017_3Dec19 --stage-to /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass2/2017/DYMC/Raw/ -x cmsRun legacy_dumper.py maxEvents=-1 -q testmatch pujidWP=tight dumpJetSysTrees=True -n 200 --no-copy-proxy runOnZee=True
