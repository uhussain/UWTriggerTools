'''
Makes optionally EG efficiency, rate and resolution plots
Authors: T.M.Perry, E.K.Friis, M.Cepeda, A.G.Levine, N.Woods UW Madison
'''
from sys import argv, stdout, stderr
import ROOT

################
# Choose Plots #
################
resolutionPlots = False
efficiencyPlots = False
ratePlots = True

# which curves to draw on rate and efficiency plots
aIsoID = True    # I region + Iso<ISOTHRESHOLD
aPUa = True      # P region + ID + ISO + PU subtraction (puValA)
aPUb = True      # P region + ID + ISO + PU subtraction (puValB)
aNoIsoID = True  # R region
aLOne = True     # C current
aMC_isoID = True  # Mi MC region + ID + ISO
aMC_lOne = True   # Mc MC current

#rate plot
rateLine = False # line at recoPtVal

#resolution plot
resL1=True # cut on l1(g) pt at l1ptVal
resReco=True # cut on reco pt at recoPtVal
##################
# Set Parameters #
##################
LIso=3
LSB=50
l1ptVal=20
recoPtVal= 40
ISOTHRESHOLD=0.20
puValA = 5
puValB = 9
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0
ZEROBIAS_RATE=15000000.00
saveWhere='../plots/Tau_'

########
# File #
########
# Efficiency
eff_ntuple = '../data/LSB'+str(LSB)+'/uct_tau_efficiency.root'
eff_rlx_spot = 'rlxTauEfficiency/Ntuple'
eff_iso_spot = 'isoTauEfficiency/Ntuple'
eff_ntuple_file = ROOT.TFile(eff_ntuple)
eff_rlx_tau_ntuple = eff_ntuple_file.Get(eff_rlx_spot)
eff_iso_tau_ntuple = eff_ntuple_file.Get(eff_iso_spot)

# Rate
rate_ntuple = '../data/LSB'+str(LSB)+'/uct_rates_eic'+str(LIso)+'.root'
rate_rlx_l1_spot = 'tauL1Rate/Ntuple'
rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple'
rate_ntuple_file = ROOT.TFile(rate_ntuple)
rate_rlx_l1_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1_spot)
rate_rlx_l1g_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1g_spot)

# Rate
mc_rate_ntuple = '../data/LSB'+str(LSB)+'/uct_rates_eic'+str(LIso)+'_mc.root'
mc_rate_rlx_l1_spot = 'tauL1Rate/Ntuple'
mc_rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple'
mc_rate_ntuple_file = ROOT.TFile(mc_rate_ntuple)
mc_rate_rlx_l1_tau_ntuple = mc_rate_ntuple_file.Get(mc_rate_rlx_l1_spot)
mc_rate_rlx_l1g_tau_ntuple = mc_rate_ntuple_file.Get(mc_rate_rlx_l1g_spot)

name=''
if aIsoID: name+='I'
if aPUa or aPUb: name+='P'
if aNoIsoID: name+='N'
if aLOne: name+='C'
if aMC_isoID: name+='Mi'
if aMC_lOne: name+='Mc'
extraName=''

log = open(saveWhere+name+'_reco'+str(recoPtVal)+'_iso'+str(ISOTHRESHOLD)+extraName+'.log','w')
log.write('LIso = '+str(LIso)+'\n')
log.write('LSB = '+str(LSB)+'\n')
log.write('l1ptVal = '+str(l1ptVal)+'\n')
log.write('recoPtVal = '+str(recoPtVal)+'\n')
log.write('puValA = '+str(puValA)+'\n')
log.write('puValB = '+str(puValB)+'\n')
log.write('ISOTHRESHOLD = '+str(ISOTHRESHOLD)+'\n')
log.write('L1_CALIB_FACTOR = '+str(L1_CALIB_FACTOR)+'\n')
log.write('L1G_CALIB_FACTOR = '+str(L1G_CALIB_FACTOR)+'\n')
log.write('ZEROBIAS_RATE = '+str(ZEROBIAS_RATE)+'\n\n')

#########
# STYLE #
#########
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

tex = ROOT.TLatex()
tex.SetTextSize(0.07)
tex.SetTextAlign(11)
tex.SetNDC(True)

colorI=ROOT.EColor.kBlue
markerI=20
colorN=ROOT.EColor.kBlack
markerN=24
colorC=ROOT.EColor.kRed
markerC=21
colorPUa=ROOT.EColor.kViolet-7
markerPUa=33
colorPUb=ROOT.EColor.kYellow+3
markerPUb=23
colorMi=ROOT.EColor.kGray+3
#colorMi=ROOT.EColor.kOrange+10
markerMi=24
colorMc=ROOT.EColor.kGray+3
#colorMc=ROOT.EColor.kOrange+10
markerMc=25


canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title='',calFactor=1):
 ''' Plot a variable using draw and return the histogram '''
 draw_string = "%s * %0.2f>>htemp(%s)" % (variable,calFactor, ", ".join(str(x) for x in binning))
 print draw_string
 tree.Draw(draw_string, selection, "goff")
 output_histo = ROOT.gDirectory.Get("htemp").Clone()
 output_histo.GetXaxis().SetTitle(xaxis)
 output_histo.SetTitle(title)
 return output_histo

######################################################################
##### RESOLUTION #####################################################
######################################################################
def make_res_nrml(
ntuple,reco,l1,l1g,binning,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='',setLOG=False):
 '''variable: (reco-l1g)/reco  with cut l1g > l1pt'''
 canvas.SetLogy(setLOG)
 info = 'RES_'+str(reco)+'_'+cutPtVar+'Cut'
 
 l1gplot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+' * '+str(l1g)+'))/'+str(reco), 
  'l1gMatch&&'+cutPtVarg+'>'+str(cutPt),binning
 )
 l1gplot.SetTitle('Resolution')
 l1gplot.GetXaxis().SetTitle(reco+'-'+l1+'/'+reco)
 l1gplot.SetLineColor(ROOT.EColor.kBlue)
 l1gplot.Scale(1/l1gplot.Integral())

 l1plot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(l1)+'))/'+str(reco), 
  'l1Match&&'+cutPtVar+'>'+str(cutPt),binning
 )
 l1plot.SetLineColor(ROOT.EColor.kRed)
 l1plot.Scale(1/l1plot.Integral())

 legend = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
 legend.SetFillColor(ROOT.EColor.kWhite)
 legend.SetBorderSize(0)
 legend.AddEntry(l1gplot,'upgrade')
 legend.AddEntry(l1plot,'current')

 l1gplot.Draw()
 l1plot.Draw('sames')
 l1gplot.SetMaximum(1.1*max(l1gplot.GetMaximum(),l1plot.GetMaximum()))
 legend.Draw()
 
 canvas.SaveAs(saveWhere+info+filename+'.png')

def make_res_angl(
ntuple,reco,l1,l1g,binning,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='',setLOG=False):

 canvas.SetLogy(setLOG)
 info = 'RES_'+str(reco)+'_'+cutPtVar+'Cut'
 
 l1gplot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+' * '+str(l1g)+'))', 
  'l1gMatch && '+cutPtVarg+'>'+str(cutPt),binning
 )
 l1gplot.SetTitle('Resolution')
 l1gplot.GetXaxis().SetTitle(reco+'-'+l1)
 l1gplot.SetLineColor(ROOT.EColor.kBlue)
 l1gplot.Scale(1/l1gplot.Integral())

 l1plot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(l1)+'))', 
  'l1Match&&'+cutPtVar+'>'+str(cutPt),binning
 )
 l1plot.SetLineColor(ROOT.EColor.kRed)
 l1plot.Scale(1/l1plot.Integral())

 legend = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
 legend.SetFillColor(ROOT.EColor.kWhite)
 legend.SetBorderSize(0)
 legend.AddEntry(l1gplot,'upgrade')
 legend.AddEntry(l1plot,'current')

 l1gplot.Draw()
 l1plot.Draw('sames')
 l1gplot.SetMaximum(1.1*max(l1gplot.GetMaximum(),l1plot.GetMaximum()))
 legend.Draw()
 
 canvas.SaveAs(saveWhere+filename+info+'.png')
######################################################################
##### RESOLUTION #####################################################
######################################################################

######################################################################
##### EFFICIENCY #####################################################
######################################################################
def make_l1_efficiency(denom, num,color=ROOT.EColor.kBlue,marker=20):
 ''' Make an efficiency graph '''
 eff = ROOT.TGraphAsymmErrors(num, denom)
 eff.SetMarkerStyle(marker)
 eff.SetMarkerColor(color)
 eff.SetMarkerSize(1.5)
 eff.SetLineColor(color)
 return eff

def effi_histo(ntuple,variable,cut,binning,denom,title,leg,color,marker,logg):
 num = make_plot(ntuple,variable,cut,binning)
 efi = make_l1_efficiency(denom,num,color,marker)
 leg.AddEntry(efi,title,'pe')
 efi.Draw('pe')
 logg.write(title+'\n')
 logg.write('Cut: '+cut+'\n')
 return efi

def compare_efficiencies(
 variable,
 binning,
 ntuple,
 l1ntuple=None,
 recoPtCut='(2>1)',l1PtCut='(2>1)',l1gPtCut='(2>1)',
 isoCut='(2>1)',aPUCut='(2>1)',bPUCut='(2>1)',extraCut='(2>1)',
 isoID=False,noIsoID=False,lOne=False,PUa=False,PUb=False,
 legExtra='',
 setLOG=False
):
 ''' Returns a (L1, L1G) tuple of TGraphAsymmErrors

 If [l1ntuple] is None, use [ntuple].  If not None, separate ntuples will
 be used for L1G and L1.  The use case for this is to use Rlx UCT with Iso L1
 '''
 if l1ntuple is None:
  l1ntuple = ntuple

 cutD = recoPtCut+'&&'+extraCut
 denom = make_plot(
  ntuple,variable,
  cutD,
  binning
 )

 log.write('_____________________________\n')
 log.write('-------- Efficiency ---------\n\n')
 log.write('File: '+eff_ntuple+'\n')
 log.write('Variable: '+variable+'\n\n')
 log.write('Denominator Cut: '+cutD+'\n\n')
 
 frame = ROOT.TH1F('frame','frame',*binning)
 canvas.SetLogy(setLOG)
 frame.Draw()
 frame.SetTitle('')
 frame.GetYaxis().SetTitle('Efficiency')
 frame.SetMaximum(1.1)
 if variable is 'nPVs': frame.GetXaxis().SetTitle('Nr. Primary Vertices')
 else: frame.GetXaxis().SetTitle(variable)
 tex.DrawLatex(0.1,0.91,variable+' Efficiency')
 legend = ROOT.TLegend(0.4,0.35,0.89,0.15,'','brNDC')
 legend.SetFillColor(0)
 legend.SetBorderSize(0)
 legend.SetHeader(legExtra)
 
 info = '_EFF_'+variable+'_reco'+str(recoPtVal)+'_l1'+str(l1ptVal)
# I
#ntuple with Region and Isolation
 if isoID:
  cutI=recoPtCut+'&&'+extraCut+'&&'+isoCut +'&&'+l1gPtCut+'&&l1gMatch'
  l1g=effi_histo(ntuple,variable,cutI,binning,denom,
   'Region + Isolation < '+str(ISOTHRESHOLD),legend,
   colorI,markerI,log)
# P
#ntuple with Region and Isolation
 if PUa:
  cutPUa=recoPtCut+'&&'+extraCut+'&&'+aPUCut+'&&'+l1gPtCut+'&&l1gMatch'
  aPU=effi_histo(ntuple,variable,cutPUa,binning,denom,
   'Region + (Isolation-PU/'+str(puValA)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUa,markerPUa,log)
  info += '_PU'+str(puValA)
# P
#ntuple with Region and Isolation
 if PUb:
  cutPUb=recoPtCut+'&&'+extraCut+'&&'+bPUCut+'&&'+l1gPtCut+'&&l1gMatch'
  bPU=effi_histo(ntuple,variable,cutPUb,binning,denom,
   'Region + (Isolation-PU/'+str(puValB)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUb,markerPUb,log)
  info += '_PU'+str(puValA)
# N
#ntuple with Region and without Isolation
 if noIsoID:
  cutN=recoPtCut+'&&'+extraCut+'&&'+l1gPtCut+'&& l1gMatch'
 l1g_noiso=effi_histo(ntuple,variable,cutN,binning,denom,
  'Region',legend,
   colorN,markerN,log) 
# R
#Current trigger
 if lOne:
  cutC=recoPtCut+'&&'+extraCut+'&&'+l1PtCut+'&&l1Match'
  l1=effi_histo(l1ntuple,variable,cutC,binning,denom,
  'Current',legend,
  colorC,markerC,log)
 
 legend.Draw()
 canvas.SaveAs(saveWhere+name+info+'.png')
######################################################################
##### EFFICIENCY #####################################################
######################################################################

######################################################################
##### RATE ###########################################################
######################################################################
def make_l1_rate(pt, color=ROOT.EColor.kBlack, marker=20):
 ''' Make a rate plot out of L1Extra Pts '''
 numBins = pt.GetXaxis().GetNbins()
 rate = pt.Clone()
 for i in range(1, numBins+1):
  rate.SetBinContent(i, pt.Integral(i, numBins))
 rate.SetLineColor(color)
 rate.SetMarkerStyle(marker)
 rate.SetMarkerColor(color)
 return rate

def rate_histo(ntuple,cut,binning,calibfactor,scale,color,marker,leg,title,logg,line,ptLine,w,s):
 pt = make_plot(ntuple,'pt[0]',cut,binning,'','',calibfactor)
 rate = make_l1_rate(pt,color,marker)
 rate.Scale(scale)
 rate.Draw('phsame')
 leg.AddEntry(rate,title,'pe')
 maxx = rate.GetMaximum()
 binn = rate.GetXaxis().FindBin(ptLine)
 rateVal = rate.GetBinContent(binn)
 vert=None
 hor=None
 if line==True:
  vert=ROOT.TLine(ptLine,0,ptLine,rateVal)
  vert.SetLineWidth(w)
  vert.SetLineStyle(s)
  hor=ROOT.TLine(binning[1],rateVal,ptLine,rateVal)
  hor.SetLineWidth(w)
  hor.SetLineStyle(s)   
  vert.Draw()
  hor.Draw()
 logg.write('---------------------------------\n')
 logg.write(title+'\n')
 logg.write('---------------------------------\n')
 logg.write('Cut: '+cut+'\n\n')
 logg.write('At pT = '+str(ptLine)+', Rate = '+str(rateVal)+'\n\n')
 return rate,maxx,vert,hor

def make_rate_plot(
 l1ntuple,
 uctntuple,
 binning,
 l1mcntuple=None,
 uctmcntuple=None,
 filename='',
 setLOG=True,
 isoCut='(2>1)',puACut='(2>1)',puBCut='(2>1)',extraCut='(2>1)',
 ptLine=20,
 isoID=False,noIsoID=False,lOne=False,PUa=False,PUb=False,
 MC_isoID=False,MC_lOne=False,
 line=False
 ):

 info = '_RATE'
 scale = ZEROBIAS_RATE/uctntuple.GetEntries()
 scaleMC = ZEROBIAS_RATE/uctmcntuple.GetEntries()
 
 canvas.SetLogy(setLOG)
 frame = ROOT.TH1F('frame','frame',*binning)
 frame.Draw()
 frame.SetTitle('')
 frame.GetYaxis().SetTitle('Hz (8TeV,1E34)')
 frame.GetXaxis().SetTitle('p_{T}')
 frame.SetMaximum(10000000)
 frame.SetMinimum(100)
 tex.DrawLatex(0.1,0.91,'Tau Rate')
 legend = ROOT.TLegend(0.4,0.7,0.89,0.89,'','brNDC')
 legend.SetFillColor(0)
 legend.SetBorderSize(0)

 # line (a=width b=style)
 aI=3
 bI=3
 aN=3
 bN=3
 aR=3
 bR=3
 aC=3
 bC=3
 aPa=3
 bPa=3
 aPb=3
 bPb=3
 aM=3
 bM=3

 log.write('----------------\n')
 log.write('----- Rate -----\n\n')
 log.write('File : '+rate_ntuple+'\n')

 maxI = 1
 maxN = 1
 maxC = 1
 maxPa = 1
 maxPb = 1
 maxMi = 1
 maxMc = 1
# I
# Isolated
 if isoID:
  cutI=isoCut+'&&'+extraCut
  l1gRate,maxI,vertI,horI = rate_histo(
   uctntuple,cutI,binning,L1G_CALIB_FACTOR,
   scale,colorI,markerI,legend,
   'Region + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aI,bI)
  info+='_iso_'+str(ISOTHRESHOLD)
# P
# Isolated - PU/A
 if PUa:
  cutPUa=puACut+'&&'+extraCut
  PUaRate,maxPUa,vertPUa,horPUa = rate_histo(
   uctntuple,cutPUa,binning,L1G_CALIB_FACTOR,
   scale,colorPUa,markerPUa,legend,
   'Region + (Iso - PU/'+str(puValA)+')<'+str(ISOTHRESHOLD),
   log,line,ptLine,aPa,bPa)
  info+='_PU'+str(puValA)
# P
# Isolated - PU/A
 if PUb:
  cutPUb=puBCut+'&&'+extraCut
  PUbRate,maxPUb,vertPUb,horPUb = rate_histo(
   uctntuple,cutPUb,binning,L1G_CALIB_FACTOR,
   scale,colorPUa,markerPUa,legend,
   'Region + (Iso - PU/'+str(puValB)+')<'+str(ISOTHRESHOLD),
   log,line,ptLine,aPb,bPb)
  info+='_PU'+str(puValB)
# Non Isolated
 if noIsoID:
  cutN=extraCut
  l1gNoIsoRate,maxN,vertN,horN = rate_histo(
   uctntuple,cutN,binning,L1G_CALIB_FACTOR,
   scale,colorN,markerN,legend,
   'Region',
   log,line,ptLine,aN,bN)
# Current
 if lOne:
  l1Rate,maxC,vertC,horC = rate_histo(
   l1ntuple,'',binning,L1G_CALIB_FACTOR,
   scale,colorC,markerC,legend,
   'Current', 
   log,line,ptLine,aC,bC)
# Mi
# MC: Region + ID + ISO
 if MC_isoID:
  cutMi=isoCut+'&&'+extraCut
  l1g_mcRate,maxMi,vertMi,horMi = rate_histo(
   uctmcntuple,cutMi,binning,L1G_CALIB_FACTOR,
   scaleMC,colorMi,markerMi,legend,
   'MC: Region + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aM,bM)
# Mc
# MC: Current
 if MC_lOne:
  l1mcRate,maxMc,vertMc,horMc = rate_histo(
   l1mcntuple,'',binning,L1G_CALIB_FACTOR,
   scaleMC,colorC,markerC,legend,
   'MC: Current', 
   log,line,ptLine,aM,bM)
 print(str(maxI)+' '+str(maxN)+' '+str(maxC)+' '+str(maxMc)+' '+str(maxMi))
 frame.SetMaximum(5*max(maxI,maxN,maxC,maxMc,maxMi))
 legend.Draw()
 canvas.SaveAs(saveWhere+name+info+'.png')
######################################################################
##### RATE ###########################################################
######################################################################

######################################################################
###### DRAW PLOTS ####################################################
######################################################################
####################
# Resolution Plots #
####################
if resolutionPlots == True:
 binPt = [100,-2,2]
 binEta = [50,-0.5,0.5]

# ntuple,reco,l1,l1g,binning,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='',setLOG=False
 if resL1 == True:
  make_res_nrml(eff_rlx_tau_ntuple, 'recoPt', 'l1Pt', 'max(l1gPt,l1gRegionEt)', binPt, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoEta','l1Eta','l1gEta',binEta, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoPhi','l1Phi','l1gPhi',binEta, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')

 if resReco == True:
  make_res_nrml(eff_rlx_tau_ntuple, 'recoPt', 'l1Pt', 'max(l1gPt,l1gRegionEt)', binPt, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoEta','l1Eta','l1gEta',binEta, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoPhi','l1Phi','l1gPhi',binEta, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,filename='')
####################
# Efficiency Plots #
####################
if efficiencyPlots == True:

 binPt = [10,0,100]
 binVert=[10,0,40]
 binJetPt=[40,0,70]
 
 # variable,
 # binning,
 # ntuple,
 # l1ntuple=None,
 # recoPtCut='(2>1)',l1PtCut='(2>1)',isoCut='(2>1)',extraCut
 # isoID=False,noIsoID=False,EGT=False,lOne=False

 compare_efficiencies(
  'recoPt',
  binPt,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(max(l1gPt,l1gRegionEt) > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  isoID=aIsoID,
  noIsoID=aNoIsoID,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb
 )
 
 compare_efficiencies(
  'nPVs',
  binVert,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(max(l1gPt,l1gRegionEt) > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  isoID=aIsoID,
  noIsoID=aNoIsoID,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

##############
# Rate Plots #
##############
if ratePlots == True:
 binRate = [14,0,75]

# l1ntuple,
# uctntuple,
# binning,
# l1mcntuple=None,
# uctmcntuple=None,
# filename='',
# setLOG=True,
# isoCut='(2>1)',puACut='(2>1)',puBCut='(2>1)',extraCut='(2>1)',
# ptLine=20,
# isoID=False,noIsoID=False,lOne=False,PUa=False,PUb=False,
# MC_isoID=False,MC_lOne=False,
# line=False

 make_rate_plot(rate_rlx_l1_tau_ntuple,rate_rlx_l1g_tau_ntuple,binRate,
  l1mcntuple=mc_rate_rlx_l1_tau_ntuple,
  uctmcntuple=mc_rate_rlx_l1g_tau_ntuple,
  filename='',
  setLOG=True,
  isoCut='((jetPt[0] - max(regionPt[0], pt[0]))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  puACut='((jetPt[0] - max(regionPt[0], pt[0])-(pu[0]/'+str(puValA)+'))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  puBCut='((jetPt[0] - max(regionPt[0], pt[0])-(pu[0]/'+str(puValB)+'))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  isoID=aIsoID,
  noIsoID=aNoIsoID,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb,
  MC_isoID=aMC_isoID,
  MC_lOne=aMC_lOne,
  line = rateLine
 )
