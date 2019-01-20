from numpy import sqrt
import ROOT as r
r.gROOT.SetBatch(True)
#filename = 'VBF_photononly.root'
#filename = 'VBF_jetonly.root'
filename = 'VBF_photonandjet.root'
theFile = r.TFile(filename)
theDir = theFile.Get('tagsDumper/trees')
for key in theDir.GetListOfKeys():
  theTree = theDir.Get(key.GetName())
  if not ('PTH' in theTree.GetName() or 'VBF' in theTree.GetName() or 'RECO_0J' in theTree.GetName()): continue
  if theTree.GetEntries()==0:
    print 'no events for cat %s'%(theTree.GetName())
    continue
  sumN = 0
  sumW = 0.
  sumW2 = 0.
  sumP = 0.
  for i in range(theTree.GetEntries()):
    theTree.GetEntry(i)
    prob = getattr(theTree,'prefireProbability')
    w = getattr(theTree,'weight')
    sumN += 1
    sumW += w
    sumW2 += w*w
    sumP += w * prob
  ave = sumP / sumW
  neff = sumW**2 / sumW2
  frac_err = sqrt(neff) / neff
  err = frac_err * ave
  print
  #print 'average prefire probability for cat %s:  %1.3f%%  (with %.3f effective events)'%(theTree.GetName(), 100.*sumP/sumW, (sumW*sumW/sumW2))
  print '%s:  %1.3f%%  +-  %.3f%%'%(theTree.GetName(), 100.*ave, 100.*err)
