import ROOT as r
import numpy as np
import usefulStyle
from os import system


def getEffSigma(theHist, wmin=115., wmax=135., step=0.001, epsilon=0.005):
  point = wmin
  weight = 0.
  points = [] #vector<pair<double,double> > 
  thesum = theHist.Integral()
  if thesum <= 0: return 0.
  for i in range(theHist.GetNbinsX()):
    weight += theHist.GetBinContent(i)
    if weight > epsilon:
      points.append( [theHist.GetBinCenter(i),weight/thesum] )
  low = wmin
  high = wmax
  width = wmax-wmin
  for i in range(len(points)):
    for j in range(i,len(points)):
      wy = points[j][1] - points[i][1]
      if abs(wy-0.683) < epsilon:
        wx = points[j][0] - points[i][0]
        if wx < width:
          low = points[i][0]
          high = points[j][0]
          width=wx
  return 0.5*(high-low)

def getRealSigma( theHist ):
  sigma = 2.
  if theHist.GetEntries() > 0:
    theHist.Fit('gaus')
    fit = theHist.GetFunction('gaus')
    sigma = fit.GetParameter(2)
  return sigma


def computeBkg( hist, binCutLow, binCutHigh, effSigma ):
  '''Use this to go from 2D hist in score (x) and mass (y) to a background estimation via a basic exponential fit'''
  projHist = hist.ProjectionY('tempProj', binCutLow, binCutHigh)
  bkgVal = 0.
  if projHist.GetEntries() > 0:
    projHist.Fit('expo')
    fit = projHist.GetFunction('expo')
    bkgVal = fit.Integral(125. - effSigma, 125. + effSigma)
  return bkgVal

def computeSig( hist, binCutLow, binCutHigh ):
  '''Use this to go from 2D hist in score (x) and mass (y) to a signal count plus effective sigma calculation'''
  projHist = hist.ProjectionY('tempProj', binCutLow, binCutHigh)
  sigCount = projHist.Integral()
  #FIXME
  #sigmaEff = getEffSigma( projHist )
  sigmaEff = getRealSigma( projHist )
  return sigCount, sigmaEff


def evalMetric(S, B):
  val = 0.
  if S+B > 0.: val = S / np.sqrt(S+B)
  return val

def evalFancy(S, B, reg=5.):
  bkg = B + reg
  val = 0.
  if bkg > 0.:
    val = (S + bkg)*np.log(1. + (S/bkg))
    val = 2*(val - S)
    val = np.sqrt(val)
  return val

def evalSignif( hists, proc, tag, low, high, isFancy=False ):
  sig, effSigma = computeSig( hists[(proc,tag)], low, high )
  totSig, irrel = computeSig( hists[('all',tag)], low, high )
  bkg = computeBkg( hists[('bkg',tag)], low, high , effSigma )
  totBkg = totSig + bkg - sig
  metric = evalMetric(sig, totBkg)
  if isFancy: metric = evalFancy(sig, totBkg)
  return metric

def evalSignifImproved( hists, proc, tag, low, high ):
  projSigHist = hists[('all',tag)].ProjectionY('projSigHist', low, high)
  projBkgHist = hists[('bkg',tag)].ProjectionY('projBkgHist', low, high)
  if not (projSigHist.GetEntries() > 0 and projBkgHist.GetEntries() > 0): return 0.
  projSigHist.Fit('gaus')
  initSigFit  = projSigHist.GetFunction('gaus')
  projBkgHist.Fit('expo')
  initBkgFit  = projBkgHist.GetFunction('expo')
  sPlusBfunc = r.TF1('sPlusBfunc', 'gaus(0)+expo(3)', 100., 180.)
  sPlusBfunc.SetParameter(0, initSigFit.GetParameter(0))
  sPlusBfunc.SetParameter(1, initSigFit.GetParameter(1))
  sPlusBfunc.SetParameter(2, initSigFit.GetParameter(2))
  sPlusBfunc.SetParameter(3, initBkgFit.GetParameter(0))
  sPlusBfunc.SetParameter(4, initBkgFit.GetParameter(1))
  #nBins=320
  nBins=160
  sPlusBhist = r.TH1F('sPlusBhist', 'sPlusBhist', nBins, 100., 180.)
  for i in range(1,nBins):
    lowEdge  = 100.+ (i-1)*(180.-100.)/nBins
    highEdge = 100.+ i*(180.-100.)/nBins
    count = sPlusBfunc.Integral(lowEdge, highEdge)
    sPlusBhist.SetBinContent( i, count )
  sPlusBhist.Fit(sPlusBfunc)
  sPlusBfit    = sPlusBhist.GetFunction('sPlusBfunc')
  sPlusBchisq  = sPlusBfit.GetChisquare()
  signif = r.Math.gaussian_quantile_c( r.TMath.Prob(sPlusBchisq,1) / 2., 1 )
  return signif


def plotGraphCollection( coll, ytitle='S/#sqrt{S+B}', outdir='ClassificationPlots/', copydir='', toAdd='' ):
  can = usefulStyle.setCanvas()
  for key,graph in coll.iteritems():
    graph.GetYaxis().SetRangeUser(0., 1.2*graph.GetMaximum())
    graph.SetLineColor(r.kGreen+2)
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(r.kGreen+2)
    graph.Draw("AL")
    #graph.Draw("AP")
    graph.GetHistogram().SetTitle( graph.GetTitle() )
    graph.GetHistogram().GetXaxis().SetTitle( graph.GetXaxis().GetTitle() )
    graph.GetHistogram().GetYaxis().SetTitle( ytitle )
    usefulStyle.formatHisto(graph.GetHistogram())
    graph.GetHistogram().GetXaxis().SetTitleOffset(1.3)
    usefulStyle.drawCMS()
    usefulStyle.drawEnPu()
    fileName = outdir + toAdd + graph.GetName()
    can.SaveAs(fileName+'.pdf')
    can.SaveAs(fileName+'.png')
    if not copydir == '': system('cp %s* %s'%(outdir, copydir))
