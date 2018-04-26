#!/usr/bin/env python
# first attempt constucting genuine class boundaries with dipho MVA and sigmaMoM

import os
import ROOT as r
import numpy as np
from collections import OrderedDict as od
from diphoHelpers import evalSignif, evalSignifImproved, plotGraphCollection

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d', '--diphoLower', type='int', default=30, help='set lower diphoton MVA boundary for two-class optimisation')
parser.add_option('-s', '--sigmaLower', type='int', default=70, help='set lower sigmaMoM boundary for two-class optimisation')
parser.add_option('-f', '--doFancy', default=False, action='store_true', help='use the approximate mean significance metric')
(opts,args) = parser.parse_args()

r.gROOT.SetBatch(True)

name2num = {"UNKNOWN":0,
                     # Gluon fusion
                     "GG2H_FWDH" : -1, "GG2H_VBFTOPO_JET3VETO" : 1, "GG2H_VBFTOPO_JET3" : 2, "GG2H_0J" : 3,
                     "GG2H_1J_PTH_0_60" : 4, "GG2H_1J_PTH_60_120" : 5, "GG2H_1J_PTH_120_200" : 6, "GG2H_1J_PTH_GT200" : 7,
                     "GG2H_GE2J_PTH_0_60" : 8, "GG2H_GE2J_PTH_60_120" : 9, "GG2H_GE2J_PTH_120_200" : 10, "GG2H_GE2J_PTH_GT200" : 11,
                     # "VBF" / QQ2HQQ
                     "QQ2HQQ_FWDH" : -2, "QQ2HQQ_VBFTOPO_JET3VETO" : 12, "QQ2HQQ_VBFTOPO_JET3" : 13, "QQ2HQQ_VH2JET" : 14, "QQ2HQQ_REST" : 15, "QQ2HQQ_PTJET1_GT200" : 16,
                     # qq -> WH
                     "QQ2HLNU_FWDH" : -3, "QQ2HLNU_PTV_0_150" : 17, "QQ2HLNU_PTV_150_250_0J" : 18, "QQ2HLNU_PTV_150_250_GE1J" : 19, "QQ2HLNU_PTV_GT250" : 20,
                     # qq -> ZH
                     "QQ2HLL_FWDH" : -4, "QQ2HLL_PTV_0_150" : 21, "QQ2HLL_PTV_150_250_0J" : 22, "QQ2HLL_PTV_150_250_GE1J" : 23, "QQ2HLL_PTV_GT250" : 24,
                     # gg -> ZH
                     "GG2HLL_FWDH" : -5, "GG2HLL_PTV_0_150" : 25, "GG2HLL_PTV_GT150_0J" : 26, "GG2HLL_PTV_GT150_GE1J" : 27,
                     # ttH
                     "TTH_FWDH" : -6, "TTH" : 28,
                     # bbH
                     "BBH_FWDH" : -7, "BBH" : 29,
                     # tH
                     "TH_FWDH" : -8, "TH" : 30}
num2name = {v: k for k, v in name2num.iteritems()}

stage1tagMap = { 0:'LOGICERROR', 1:'RECO_0J', 2:'RECO_1J_PTH_0_60', 3:'RECO_1J_PTH_60_120', 4:'RECO_1J_PTH_120_200', 5:'RECO_1J_PTH_GT200', 
                             6:'RECO_GE2J_PTH_0_60', 7:'RECO_GE2J_PTH_60_120', 8:'RECO_GE2J_PTH_120_200', 9:'RECO_GE2J_PTH_GT200', 10:'RECO_VBFTOPO_JET3VETO', 11:'RECO_VBFTOPO_JET3', 12:'RECO_VH2JET',
                             13:'RECO_0LEP_PTV_0_150', 14:'RECO_0LEP_PTV_150_250_0J', 15:'RECO_0LEP_PTV_150_250_GE1J', 16:'RECO_0LEP_PTV_GT250', 
                             17:'RECO_1LEP_PTV_0_150', 18:'RECO_1LEP_PTV_150_250_0J', 19:'RECO_1LEP_PTV_150_250_GE1J', 20:'RECO_1LEP_PTV_GT250', 
                             21:'RECO_2LEP_PTV_0_150', 22:'RECO_2LEP_PTV_150_250_0J', 23:'RECO_2LEP_PTV_150_250_GE1J', 24:'RECO_2LEP_PTV_GT250',
                             25:'RECO_TTH_LEP', 26:'RECO_TTH_HAD' }

def main():
  print 'running with fancy metric set to %s'%str(opts.doFancy)
  sigDir = 'Systematics/test/sig_jobs_stage1_forClassification/'
  bkgDir = 'Systematics/test/data_jobs_stage1_forClassification/'

  procMap = {'GluGluHToGG':'ggh', 'VBFHToGG':'vbf','ttHJetToGG':'tth','WHToGG':'wh','ZHToGG':'zh'}
  #theProcs = ['GluGluHToGG', 'VBFHToGG','ttHJetToGG','WHToGG','ZHToGG']
  #theProcs = ['GluGluHToGG'] #ggH phase space
  theProcs = ['GluGluHToGG','VBFHToGG'] #VBF phase space

  #theTags = 'RECO_0J,RECO_1J_PTH_0_60,RECO_1J_PTH_60_120,RECO_1J_PTH_120_200,RECO_1J_PTH_GT200,RECO_GE2J_PTH_0_60,RECO_GE2J_PTH_60_120,RECO_GE2J_PTH_120_200,RECO_GE2J_PTH_GT200,RECO_VBFTOPO_JET3VETO,RECO_VBFTOPO_JET3,RECO_VH2JET,RECO_0LEP_PTV_0_150,RECO_0LEP_PTV_150_250_0J,RECO_0LEP_PTV_150_250_GE1J,RECO_0LEP_PTV_GT250,RECO_1LEP_PTV_0_150,RECO_1LEP_PTV_150_250_0J,RECO_1LEP_PTV_150_250_GE1J,RECO_1LEP_PTV_GT250,RECO_2LEP_PTV_0_150,RECO_2LEP_PTV_150_250_0J,RECO_2LEP_PTV_150_250_GE1J,RECO_2LEP_PTV_GT250,RECO_TTH_LEP,RECO_TTH_HAD'
  #theTags = 'RECO_0J,RECO_1J_PTH_0_60,RECO_1J_PTH_60_120,RECO_1J_PTH_120_200,RECO_1J_PTH_GT200,RECO_GE2J_PTH_0_60,RECO_GE2J_PTH_60_120,RECO_GE2J_PTH_120_200,RECO_GE2J_PTH_GT200' #ggH phase space
  theTags = 'RECO_VBFTOPO_JET3VETO,RECO_VBFTOPO_JET3' #VBF phase space
  theTags = theTags.split(',')

  sigTrees = {}
  bkgTrees = {}
  print 'getting trees'
  for theProc in theProcs:
    for theTag in theTags: 
      sigTrees[(theProc,theTag)] = r.TChain('tagsDumper/trees/%s_125_13TeV_%s'%(procMap[theProc],theTag))
      for root, dirs, files in os.walk(sigDir):
        for fileName in files:
          if not theProc in fileName: continue
          if not '.root' in fileName: continue
          if not 'output_' in fileName: continue
          sigTrees[(theProc,theTag)].Add(sigDir+fileName)
  for theTag in theTags: 
    bkgTrees[theTag] = r.TChain('tagsDumper/trees/Data_13TeV_%s'%(theTag))
    for root, dirs, files in os.walk(bkgDir):
      for fileName in files:
        if not '.root' in fileName: continue
        if not 'output_DoubleEG' in fileName: continue
        bkgTrees[theTag].Add(bkgDir+fileName)
  print 'got trees'

  #main code starts here
  theDiphoHists = {}
  theSigmaHists = {}
  #stage1procs = ['GG2H_0J', 'GG2H_1J_PTH_0_60', 'GG2H_1J_PTH_60_120', 'GG2H_1J_PTH_120_200', 'GG2H_1J_PTH_GT200', 'GG2H_GE2J_PTH_0_60', 'GG2H_GE2J_PTH_60_120', 'GG2H_GE2J_PTH_120_200', 'GG2H_GE2J_PTH_GT200'] #ggH phase space
  stage1procs = ['GG2H_VBFTOPO_JET3VETO', 'GG2H_VBFTOPO_JET3', 'QQ2HQQ_VBFTOPO_JET3VETO', 'QQ2HQQ_VBFTOPO_JET3'] #VBF phase space
  procNameMap = {'GG2H_0J':'0J', 'GG2H_1J_PTH_0_60':'1J low', 'GG2H_1J_PTH_60_120':'1J med', 'GG2H_1J_PTH_120_200':'1J high', 'GG2H_1J_PTH_GT200':'1J BSM', 'GG2H_GE2J_PTH_0_60':'GE2J low', 'GG2H_GE2J_PTH_60_120':'GE2J med', 'GG2H_GE2J_PTH_120_200':'GE2J high', 'GG2H_GE2J_PTH_GT200':'GE2J BSM', 'GG2H_VBFTOPO_JET3VETO':'ggH in VBF 2jet', 'GG2H_VBFTOPO_JET3':'ggH in VBF 3jet', 'QQ2HQQ_VBFTOPO_JET3VETO':'VBF in VBF 2jet', 'QQ2HQQ_VBFTOPO_JET3':'VBF in VBF 3jet'}
  #initialise hists
  nBins = 100
  for stage1proc in stage1procs:
    for theTag in theTags: 
      theDiphoHists[(stage1proc,theTag)] = r.TH2F('diphoHist_%s_%s'%(stage1proc,theTag), 'diphoHist_%s_%s'%(stage1proc,theTag), nBins, -1., 1., nBins, 100., 180.)
      theSigmaHists[(stage1proc,theTag)] = r.TH2F('sigmaHist_%s_%s'%(stage1proc,theTag), 'sigmaHist_%s_%s'%(stage1proc,theTag), nBins, 0., 0.05, nBins, 100., 180.)
  for theTag in theTags: 
    theDiphoHists[('bkg',theTag)] = r.TH2F('diphoHist_bkg_%s'%(theTag), 'diphoHist_bkg_%s'%(theTag), nBins, -1., 1., nBins, 100., 180.)
    theSigmaHists[('bkg',theTag)] = r.TH2F('sigmaHist_bkg_%s'%(theTag), 'sigmaHist_bkg_%s'%(theTag), nBins, 0., 0.05, nBins, 100., 180.)
    theDiphoHists[('all',theTag)] = r.TH2F('diphoHist_all_%s'%(theTag), 'diphoHist_all_%s'%(theTag), nBins, -1., 1., nBins, 100., 180.)
    theSigmaHists[('all',theTag)] = r.TH2F('sigmaHist_all_%s'%(theTag), 'sigmaHist_all_%s'%(theTag), nBins, 0., 0.05, nBins, 100., 180.)
  #loop over trees to fill hists
  print 'now fill hists'
  for theProc in theProcs:
    for theTag in theTags: 
      print 'processing %s, %s'%(theProc, theTag)
      for i,ev in enumerate(sigTrees[(theProc,theTag)]):
        stage1proc = getattr(ev,'stage1cat')
        stage1proc = num2name[stage1proc]
        if not stage1proc in stage1procs: continue
        weight     = getattr(ev,'weight')
        weight     = 35.9 * weight
        diphoMVA   = getattr(ev,'diphoMVA')
        sigmaMoM   = getattr(ev,'decorrSigmarv')
        mass       = getattr(ev,'CMS_hgg_mass')
        theDiphoHists[(stage1proc,theTag)].Fill( diphoMVA, mass, weight)
        theSigmaHists[(stage1proc,theTag)].Fill( sigmaMoM, mass, weight)
        theDiphoHists[('all',theTag)].Fill( diphoMVA, mass, weight)
        theSigmaHists[('all',theTag)].Fill( sigmaMoM, mass, weight)
  for theTag in theTags: 
    print 'processing data, %s'%(theTag)
    for i,ev in enumerate(bkgTrees[theTag]):
      mass = getattr(ev,'CMS_hgg_mass')
      diphoMVA   = getattr(ev,'diphoMVA')
      sigmaMoM   = getattr(ev,'decorrSigmarv')
      theDiphoHists[('bkg',theTag)].Fill( diphoMVA, mass )
      theSigmaHists[('bkg',theTag)].Fill( sigmaMoM, mass )
  print 'hists be full'

  #now calculate significance type quantities
  diphoSignifs = od()
  diphoBests   = od()
  sigmaSignifs = od()
  sigmaBests   = od()
  print 'now calculating significances'
  theTempLine = ''
  for stage1proc in stage1procs:
    print 'processing %s'%stage1proc
    diphoSignifs[stage1proc] = r.TGraph()
    sigmaSignifs[stage1proc] = r.TGraph()
    diagonalTag = 'RECO' + stage1proc.replace( stage1proc.split('_')[0], '')
    diphoMax = -1.
    diphoBest = -1.
    counter = 0
    #tempChisq = evalSignifImproved( theDiphoHists, stage1proc, diagonalTag, 50, -1 )
    #theTempLine += 'signif = %1.3f \n'%(tempChisq)
    #for iBin in range(0,1):
    for iBin in range(opts.diphoLower, nBins-3):
      signifZero = evalSignif( theDiphoHists, stage1proc, diagonalTag, opts.diphoLower, iBin, opts.doFancy )
      signifOne  = evalSignif( theDiphoHists, stage1proc, diagonalTag, iBin, -1, opts.doFancy )
      #signifZero = evalSignifImproved( theDiphoHists, stage1proc, diagonalTag, opts.diphoLower, iBin )
      #signifOne  = evalSignifImproved( theDiphoHists, stage1proc, diagonalTag, iBin, -1 )
      totSignif  = np.sqrt( signifZero*signifZero + signifOne*signifOne )
      diphoCutVal = -1. + 0.02*iBin
      diphoSignifs[stage1proc].SetPoint( counter, diphoCutVal, totSignif )
      if totSignif > diphoMax: 
        diphoMax  = totSignif
        diphoBest = diphoCutVal
      counter += 1
    diphoBests[stage1proc] = (diphoBest, diphoMax)
    sigmaMax = -1.
    sigmaBest = -1.
    counter = 0
    #for iBin in range(nBins-1,nBins):
    for iBin in range(15, opts.sigmaLower):
      signifZero = evalSignif( theSigmaHists, stage1proc, diagonalTag, iBin, opts.sigmaLower, opts.doFancy)
      signifOne  = evalSignif( theSigmaHists, stage1proc, diagonalTag, 0, iBin, opts.doFancy )
      #signifZero = evalSignifImproved( theSigmaHists, stage1proc, diagonalTag, iBin, opts.sigmaLower )
      #signifOne  = evalSignifImproved( theSigmaHists, stage1proc, diagonalTag, 0, iBin )
      totSignif  = np.sqrt( signifZero*signifZero + signifOne*signifOne )
      sigmaCutVal = 0.0005*iBin
      sigmaSignifs[stage1proc].SetPoint( counter, sigmaCutVal, totSignif )
      if totSignif > sigmaMax: 
        sigmaMax = totSignif
        sigmaBest = sigmaCutVal
      counter += 1
    sigmaBests[stage1proc] = (sigmaBest, sigmaMax)
    #for some reason this has to be done last...
    diphoSignifs[stage1proc].SetMaximum(diphoMax)
    diphoSignifs[stage1proc].SetName('diphoSignif_%s'%stage1proc)
    diphoSignifs[stage1proc].GetXaxis().SetRangeUser(-1., 1.)
    diphoSignifs[stage1proc].GetXaxis().SetTitle('Diphoton MVA cut')
    diphoSignifs[stage1proc].SetTitle('Significance for %s'%procNameMap[stage1proc])
    sigmaSignifs[stage1proc].SetMaximum(sigmaMax)
    sigmaSignifs[stage1proc].GetXaxis().SetRangeUser(0., 0.05)
    sigmaSignifs[stage1proc].SetName('sigmaSignif_%s'%stage1proc)
    sigmaSignifs[stage1proc].GetXaxis().SetTitle('#sigma_{M}/M cut')
    sigmaSignifs[stage1proc].SetTitle('Significance for %s'%procNameMap[stage1proc])
  print 'all done, write hists'

  #draw hists, send to web
  #doPlotting = False
  doPlotting = True
  if doPlotting:
    #plotGraphCollection(diphoSignifs, outdir='TwoClassPlots/DiphotonMVA/', copydir='/afs/cern.ch/user/e/escott/www/FinalFits/Stage1STXS/ClassificationStudy/TwoClass/DiphotonMVA/',
    #                    toAdd=str((-1.+0.02*opts.diphoLower)))
    #plotGraphCollection(sigmaSignifs, outdir='TwoClassPlots/SigmaMoM/', copydir='/afs/cern.ch/user/e/escott/www/FinalFits/Stage1STXS/ClassificationStudy/TwoClass/SigmaMoM/',
    #                    toAdd=str((0.0005*opts.sigmaLower)))
    plotGraphCollection(diphoSignifs, outdir='TwoClassVBF/DiphotonMVA/', copydir='/afs/cern.ch/user/e/escott/www/FinalFits/Stage1STXS/ClassificationStudy/TwoClassVBF/DiphotonMVA/',
                        toAdd=str((-1.+0.02*opts.diphoLower)))
    plotGraphCollection(sigmaSignifs, outdir='TwoClassVBF/SigmaMoM/', copydir='/afs/cern.ch/user/e/escott/www/FinalFits/Stage1STXS/ClassificationStudy/TwoClassVBF/SigmaMoM/',
                        toAdd=str((0.0005*opts.sigmaLower)))

  #save hists
  outFile = r.TFile('twoClassHists.root','RECREATE')
  for key,hist in theDiphoHists.iteritems():
    hist.Write()
  for key,hist in theSigmaHists.iteritems():
    hist.Write()
  for key,graph in diphoSignifs.iteritems():
    graph.Write()
  for key,graph in sigmaSignifs.iteritems():
    graph.Write()
  outFile.Close()

  #with open('TwoClassPlots/signifValues_%s_dipho%s_sigma%s.txt'%(str(opts.doFancy),opts.diphoLower,opts.sigmaLower),'w') as textFile:
  with open('TwoClassVBF/signifValues_%s_dipho%s_sigma%s.txt'%(str(opts.doFancy),opts.diphoLower,opts.sigmaLower),'w') as textFile:
    textFile.write('Diphoton MVA optimal cut points for lower bound of %s\n'%(-1.+ 0.02*opts.diphoLower))
    for proc,vals in diphoBests.iteritems():
      textFile.write('Proc %s: dipho cuts at %1.2f, %1.2f gives signif %1.2f\n' %(procNameMap[proc], (-1.+0.02*opts.diphoLower), vals[0], vals[1]))

    textFile.write('\nSigmaMoM optimal cut points for lower bound of %s\n'%(0.0005*opts.sigmaLower))
    for proc,vals in sigmaBests.iteritems():
      textFile.write('Proc %s: sigma cuts %1.3f, %1.3f gives signif %1.2f\n' %(procNameMap[proc], (0.0005*opts.sigmaLower), vals[0], vals[1]))
    textFile.write('\n\n')

  #print theTempLine

if __name__ == '__main__':
  main()
