#!/bin/bash

export SOURCE=${CMSSW_BASE}/flashgg/Taggers/test/VBFProduction
today=`date +%F`
outdir=/vols/cms/es811/ParamTidied/testing-templates-${today}

fggRunJobs.py --load samples_param.json -d ${outdir} \
	      -x cmsRun jobs_template_maker_cfg.py QCDParam=False \
	      -q hepmedium.q --no-use-tarball useAAA=0 atIC=1 targetLumi=1.00e+3 \
	      -n 100 
