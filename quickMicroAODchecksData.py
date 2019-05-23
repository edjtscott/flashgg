#! /usr/bin/env python

import ROOT as r
import sys
from DataFormats.FWLite import Events, Handle
from numpy import sqrt
from collections import OrderedDict as od

r.gSystem.Load("libFWCoreFWLite.so")
r.gSystem.Load("libDataFormatsFWLite.so")
r.FWLiteEnabler.enable()
r.gROOT.SetBatch()

fName = 'myMicroAODOutputFile.root'
tFile = r.TFile(fName)

events = Events(fName)
jetsHandle = Handle("vector<vector<flashgg::Jet> >")
jetsLabel = ("flashggFinalJets", "", "FLASHggMicroAOD")
miniHandle = Handle("vector<pat::Jet>")
#miniLabel = ("slimmedJets", "", "RECO")
miniLabel = ("slimmedJets", "", "DQM")

for iEv,event in enumerate(events):
    event.getByLabel(jetsLabel, jetsHandle)
    jetColls = jetsHandle.product()
    event.getByLabel(miniLabel, miniHandle)
    miniJets = miniHandle.product()
    print 'Processing event %g'%iEv
    #for iColl,coll in enumerate(jetColls):
    #  for iJet,jet in enumerate(coll):
    #    if iJet==0:
    #      print 'Info on the %gth collection zeroth jet'%iColl
    #      print jet.pt()
    #      print jet.eta()
    #      print jet.phi()
    #      #print 'Leading jet pt, eta, phi is %.3f, %.3f, %.3f'%(jet.pt, jet.eta, jet.phi)
    #print 

    for iJet in range(2):
        if len(jetColls)<2: continue
        if not (len(miniJets)>=2 and len(jetColls[0])>=2 and len(jetColls[1])>=2): continue
        miniJet = miniJets[iJet]
        zerothJet = jetColls[0][iJet]
        firstJet = jetColls[1][iJet]
        print 'Comparing jet %g between mini, zeroth, and first vertex collections'%iJet
        print 'pt : %.3f, %.3f, %.3f'%(miniJet.pt(),  zerothJet.pt(),  firstJet.pt())
        print 'eta: %.3f, %.3f, %.3f'%(miniJet.eta(), zerothJet.eta(), firstJet.eta())
        print 'phi: %.3f, %.3f, %.3f'%(miniJet.phi(), zerothJet.phi(), firstJet.phi())
