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
efficiencyPlots = True
ratePlots = True

# which curves to draw on rate and efficiency plots
aIsoID = True      # I region + ID + ISO
aPUa = False        # P region + ID + ISO + PU subtraction (puValA)
aPUb = False        # P region + ID + ISO + PU subtraction (puValB)
aNoIsoID = True    # N region + ID 
aEGT = False        # R region (EG/Tau)
aLOne = True       # C current
aMC_isoID = False   # Mi MC region + ID + ISO
aMC_lOne = False    # Mc MC current
al1b_IsoID = False  # Bi L1b region + ID + ISO

#rate plot
rateLine = False # line at recoPtVal

#resolution plot
resL1=False # cut on l1(g) pt at l1ptVal
resReco=True # cut on reco pt at recoPtVal
##################
# Set Parameters #
##################
LIso=3
LSB=50
l1ptVal=20
recoPtVal=20
ISOTHRESHOLD=0.20
puValA = 5
puValB = 9
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0
ZEROBIAS_RATE=15000000.00
saveWhere='../plots/EG_'

########
# File #
########
#Efficiency
eff_ntuple = '../data/LSB'+str(LSB)+'/uct_eg_efficiency_eic'+str(LIso)+'.root'
eff_ntuple_file = ROOT.TFile(eff_ntuple)
# L1
eff_rlx_spot = 'rlxEGEfficiency/Ntuple'
eff_iso_spot = 'isoEGEfficiency/Ntuple'
eff_rlx_eg_ntuple = eff_ntuple_file.Get(eff_rlx_spot)
eff_iso_eg_ntuple = eff_ntuple_file.Get(eff_iso_spot)
# L1b
l1b_eff_rlx_spot = 'rlxEGEfficiencyStage1B/Ntuple'
l1b_eff_iso_spot = 'isoEGEfficiencyStage1B/Ntuple'
l1b_eff_rlx_eg_ntuple = eff_ntuple_file.Get(l1b_eff_rlx_spot)
l1b_eff_iso_eg_ntuple = eff_ntuple_file.Get(l1b_eff_iso_spot)

#Rate
rate_ntuple = '../data/LSB'+str(LSB)+'/uct_rates_eic'+str(LIso)+'.root'
rate_ntuple_file = ROOT.TFile(rate_ntuple)
# L1
rate_rlx_l1g_spot = 'rlxEGUCTRate/Ntuple'
rate_iso_l1g_spot = 'isoEGUCTRate/Ntuple'
rate_rlx_l1g_eg_ntuple = rate_ntuple_file.Get(rate_rlx_l1g_spot)
rate_iso_l1g_eg_ntuple = rate_ntuple_file.Get(rate_iso_l1g_spot)
# L1b
l1b_rate_rlx_l1g_spot = 'rlxEGUCTRateStage1B/Ntuple'
l1b_rate_iso_l1g_spot = 'isoEGUCTRateStage1B/Ntuple'
l1b_rate_rlx_l1g_eg_ntuple = rate_ntuple_file.Get(l1b_rate_rlx_l1g_spot)
l1b_rate_iso_l1g_eg_ntuple = rate_ntuple_file.Get(l1b_rate_iso_l1g_spot)

# Current
rate_rlx_l1_spot = 'rlxEGL1Rate/Ntuple'
rate_iso_l1_spot = 'isoEGL1Rate/Ntuple'
rate_rlx_l1_eg_ntuple = rate_ntuple_file.Get(rate_rlx_l1_spot)
rate_iso_l1_eg_ntuple = rate_ntuple_file.Get(rate_iso_l1_spot)

# Rate MC
mc_rate_ntuple = '../data/LSB'+str(LSB)+'/uct_rates_eic'+str(LIso)+'_mc.root'
mc_rate_rlx_l1_spot = 'rlxEGL1Rate/Ntuple'
mc_rate_rlx_l1g_spot = 'rlxEGUCTRate/Ntuple'
mc_rate_iso_l1_spot = 'isoEGL1Rate/Ntuple'
mc_rate_iso_l1g_spot = 'isoEGUCTRate/Ntuple'
mc_rate_ntuple_file = ROOT.TFile(mc_rate_ntuple)
mc_rate_rlx_l1_eg_ntuple = mc_rate_ntuple_file.Get(mc_rate_rlx_l1_spot)
mc_rate_rlx_l1g_eg_ntuple = mc_rate_ntuple_file.Get(mc_rate_rlx_l1g_spot)
mc_rate_iso_l1_eg_ntuple = mc_rate_ntuple_file.Get(mc_rate_iso_l1_spot)
mc_rate_iso_l1g_eg_ntuple = mc_rate_ntuple_file.Get(mc_rate_iso_l1g_spot)

#To Be Made
store = ROOT.TFile(saveWhere+'store.root','RECREATE')

name=''
if aIsoID: name+='I'
if aPUa or aPUb: name+='P'
if aNoIsoID: name+='N'
if aEGT: name+='R'
if aLOne: name+='C'
if aMC_isoID: name+='Mi'
if aMC_lOne: name+='Mc'
if al1b_IsoID: name+='Bi'

extraName=''
name+=extraName

log = open(saveWhere+name+'_l1pt'+str(l1ptVal)+'_reco'+str(recoPtVal)+'_iso'+str(ISOTHRESHOLD)+extraName+'.log','w')
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

colorI=ROOT.EColor.kGreen+3
markerI=20
colorN=ROOT.EColor.kBlue
markerN=22
colorR=ROOT.EColor.kBlack
markerR=24
colorC=ROOT.EColor.kViolet-5
markerC=21
colorPUa=ROOT.EColor.kRed
markerPUa=33
colorPUb=ROOT.EColor.kYellow+3
markerPUb=23
colorMi=ROOT.EColor.kGreen+3
markerMi=24
colorMc=ROOT.EColor.kViolet-7
markerMc=25
colorBi=ROOT.EColor.kBlack
markerBi=25

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

 canvas.SetLogy(setLOG)
 info = 'RES_'+str(reco)+'_'+cutPtVar+'Cut'
 
 l1gplot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+' * '+str(l1g)+'))/'+str(reco), 
  'l1gMatch&&(!l1gTauVeto&&!l1gMIP)&&'+cutPtVarg+'>'+str(cutPt),binning
 )
 l1gplot.SetTitle('Resolution')
 l1gplot.GetXaxis().SetTitle(reco+'-'+l1+'/'+reco)
 l1gplot.SetLineColor(colorN)
 l1gplot.Scale(1/l1gplot.Integral())

 l1plot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(l1)+'))/'+str(reco), 
  'l1Match&&'+cutPtVar+'>'+str(cutPt),binning
 )
 l1plot.SetLineColor(colorC)
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
  'l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& '+cutPtVarg+'>'+str(cutPt),binning
 )
 l1gplot.SetTitle('Resolution')
 l1gplot.GetXaxis().SetTitle(reco+'-'+l1)
 l1gplot.SetLineColor(colorN)
 l1gplot.Scale(1/l1gplot.Integral())

 l1plot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(l1)+'))', 
  'l1Match&&'+cutPtVar+'>'+str(cutPt),binning
 )
 l1plot.SetLineColor(colorC)
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
 logg.write('---------------------------------\n')
 logg.write(title+'\n\n')
 logg.write('Tree: '+ntuple.GetDirectory().GetName()+'\n\n')
 logg.write('Cut: '+cut+'\n\n')
 return efi

def compare_efficiencies(
 variable,
 binning,
 ntuple,
 l1ntuple=None,
 bntuple=None,
 bl1ntuple=None,
 recoPtCut='(2>1)',l1PtCut='(2>1)',l1gPtCut='(2>1)',
 isoCut='(2>1)',aPUCut='(2>1)',bPUCut='(2>1)',extraCut='(2>1)',
 isoID=False,l1b_IsoID=False,noIsoID=False,EGT=False,lOne=False,PUa=False,PUb=False,
 legExtra='',
 setLOG=False
):
 ''' Returns a (L1, L1G) tuple of TGraphAsymmErrors

 If [l1ntuple] is None, use [ntuple].  If not None, separate ntuples will
 be used for L1G and L1.  The use case for this is to use Rlx UCT with Iso L1
 '''
 if l1ntuple is None:
  l1ntuple = ntuple

 cutD = recoPtCut+'&&'+extraCut+'&&dr03CombinedEt<0.2'
 denom = make_plot(
  ntuple,variable,
  cutD,
  binning
 )
 denomB = make_plot(
  bntuple,variable,
  cutD,
  binning
 )
 
 log.write('_____________________________\n')
 log.write('-------- Efficiency ---------\n\n')
 log.write('File: '+eff_ntuple+'\n')
 log.write('Variable: '+variable+'\n\n')
 log.write('Denominator Tree: '+ntuple.GetDirectory().GetName()+'\n\n')
 log.write('Denominator Stage 1b Tree: '+bntuple.GetDirectory().GetName()+'\n\n')
 log.write('Denominator Cut: '+cutD+'\n\n')
 
 frame = ROOT.TH1F('frame','frame',*binning)
 canvas.SetLogy(setLOG)
 frame.Draw()
 frame.SetTitle('')
 frame.GetYaxis().SetTitle('Efficiency')
 frame.SetMaximum(1.1)
 if variable is 'nPVs': frame.GetXaxis().SetTitle('Nr. Primary Vertices')
 else: frame.GetXaxis().SetTitle(variable)
 tex.DrawLatex(0.1,0.91,variable+' EG Efficiency')
 tex.SetTextSize(0.03)
 tex.DrawLatex(0.1,0.87,'CMS Preliminary')
 tex.SetTextSize(0.07)
 #legend = ROOT.TLegend(0.15,0.35,0.69,0.55,'','brNDC')
 legend = ROOT.TLegend(0.35,0.35,0.89,0.55,'','brNDC')
 legend.SetFillColor(0)
 legend.SetBorderSize(0)
 legend.SetHeader(legExtra)
 
 info = '_EFF_'+variable+'_reco'+str(recoPtVal)+'_l1'+str(l1ptVal)

# I
#ntuple with our isolation and ID
 if isoID:
  cutI=recoPtCut+'&&'+isoCut+'&&'+l1gPtCut+'&&l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& %0.2f *1&&dr03CombinedEt<0.2' % (L1G_CALIB_FACTOR)
  l1g=effi_histo(ntuple,variable,cutI,binning,denom,
   'L1: Region + ID + ISO < '+str(ISOTHRESHOLD),legend,
    colorI,markerI,log)
  info += '_iso'+str(ISOTHRESHOLD)
  l1g.SetName('l1g')
  l1g.Write()
# Bi
# L1B ntuple with our isolation and ID
 if l1b_IsoID:
  cutBi=recoPtCut+'&&'+isoCut+'&&'+l1gPtCut+'&&l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& %0.2f *1&&dr03CombinedEt<0.2' % (L1G_CALIB_FACTOR)
  bl1g=effi_histo(bntuple,variable,cutBi,binning,denomB,
   'L1B: Region + ID + ISO < '+str(ISOTHRESHOLD),legend,
    colorBi,markerBi,log)
  bl1g.SetName('bl1g')
  bl1g.Write()
# P
#ntuple with ISO+ID-PU(a)
 if PUa:
  cutPUa=recoPtCut+'&&'+aPUCut+'&&'+l1gPtCut+'&&l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& %0.2f *1&&dr03CombinedEt<0.2' % (L1G_CALIB_FACTOR)
  aPU=effi_histo(ntuple,variable,cutPUa,binning,denom,
   'L1: Region + ID + (Isolation-PU/'+str(puValA)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUa,markerPUa,log)
  info += '_PU'+str(puValA)
  aPU.SetName('aPU')
  aPU.Write()
# P
#ntuple with ISO+ID-PU(a)
 if PUb:
  cutPUb=recoPtCut+'&&'+bPUCut+'&&'+l1gPtCut+'&&l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& %0.2f *1&&dr03CombinedEt<0.2' % (L1G_CALIB_FACTOR)
  bPU=effi_histo(ntuple,variable,cutPUb,binning,denom,
   'L1: Region + ID + (Isolation-PU/'+str(puValB)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUb,markerPUb,log)
  info += '_PU'+str(puValB)
  bPU.SetName('bPU')
  bPU.Write()
#N
#ntuple with ID but no ISO
 if noIsoID:
  cutN=recoPtCut+'&&'+l1gPtCut+'&&l1gMatch&&(!l1gTauVeto&&!l1gMIP)&& dr03CombinedEt < 0.2'
  l1g_noiso=effi_histo(ntuple,variable,cutN,binning,denom,
  'L1: Region + ID',legend,
   colorN,markerN,log)
  l1g_noiso.SetName('l1g_noiso')
  l1g_noiso.Write()
# R
#ntuple without tau veto
 if EGT:
  cutR=recoPtCut+'&&'+l1gPtCut+'&&l1gMatch&&dr03CombinedEt<0.2'
  l1g_egtau=effi_histo(ntuple,variable,cutR,binning,denom,
  'L1: Region (EG/Tau)',legend,
  colorR,markerR,log)
  l1g_egtau.SetName('l1g_egtau')
  l1g_egtau.Write()
# C
# Current trigger
 if lOne:
  cutC=recoPtCut+'&&'+l1PtCut+'&&l1Match&&dr03CombinedEt<0.2'
  l1=effi_histo(l1ntuple,variable,cutC,binning,denom,
  'Current',legend,
  colorC,markerC,log)
  l1.SetName('l1')
  l1.Write()

 legend.Draw()
 canvas.SaveAs(saveWhere+name+info+'.png')
######################################################################
##### EFFICIENCY #####################################################
######################################################################

######################################################################
##### RATES ##########################################################
######################################################################
def make_l1_rate(pt, color=ROOT.EColor.kBlack, marker=20):
 ''' Make a rate plot out of L1Extra Pts '''
 numBins = pt.GetXaxis().GetNbins()
 rate = pt.Clone()
 for i in range(1,numBins):
  rate.SetBinContent(i,pt.Integral(i,numBins))
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
 logg.write(title+'\n\n')
 logg.write('Tree: '+ntuple.GetDirectory().GetName()+'\n\n')
 logg.write('Cut: '+cut+'\n\n')
 logg.write('At pT = '+str(ptLine)+', Rate = '+str(rateVal)+'\n\n')
 return rate,maxx,vert,hor

def make_rate_plot(
 l1ntuple,
 uctntuple,
 binning,
 l1mcntuple=None,
 uctmcntuple=None,
 bl1ntuple=None,
 buctntuple=None,
 filename='',
 setLOG=True,
 isoCut='(2>1)',puACut='(2>1)',puBCut='(2>1)',extraCut='(2>1)',
 ptLine=20,
 isoID=False,l1b_isoID=False,
 noIsoID=False,EGT=False,lOne=False,PUa=False,PUb=False,
 MC_isoID=False,MC_lOne=False,
 line=False
 ):

 info = '_RATE_l1'+str(l1ptVal)
 scale = ZEROBIAS_RATE/uctntuple.GetEntries()
 scaleMC = ZEROBIAS_RATE/uctmcntuple.GetEntries()
 scaleB = ZEROBIAS_RATE/buctntuple.GetEntries() #same as scale = 5.844
 

 canvas.SetLogy(setLOG)
 frame = ROOT.TH1F('frame','frame',*binning)
 frame.Draw()
 frame.SetTitle('')
 frame.GetYaxis().SetTitle('Hz (8TeV,1E34)')
 frame.GetXaxis().SetTitle('p_{T}')
 tex.DrawLatex(0.1,0.91,'EG Rate')
 tex.SetTextSize(0.03)
 tex.SetTextAlign(31)
 tex.DrawLatex(0.9,0.91,'CMS Preliminary')
 tex.SetTextSize(0.07)
 tex.SetTextAlign(11)
 legend = ROOT.TLegend(0.25,0.7,0.89,0.89,'','brNDC')
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
 aBi=3
 bBi=3

 log.write('________________\n')
 log.write('----- Rate -----\n\n')
 log.write('File : '+rate_ntuple+'\n')

 maxI = 1
 maxN = 1
 maxR = 1
 maxC = 1
 maxPa = 1
 maxPb = 1
 maxMi = 1
 maxMc = 1
 maxBi = 1
# I
# Region + ID + ISO
 if isoID:
  cutI=isoCut+'&&'+'(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  l1gRate,maxI,vertI,horI = rate_histo(
   uctntuple,cutI,binning,L1G_CALIB_FACTOR,
   scale,colorI,markerI,legend,
   'Region + ID + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aI,bI)
  info+='_iso_'+str(ISOTHRESHOLD)
# Bi
# Region + ID + ISO
 if l1b_isoID:
  cutBi=isoCut+'&&'+'(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  bl1gRate,maxBi,vertBi,horBi = rate_histo(
   buctntuple,cutBi,binning,L1G_CALIB_FACTOR,
   scaleB,colorBi,markerBi,legend,
   'L1B: Region + ID + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aBi,bBi)
# Mi
# MC: Region + ID + ISO
 if MC_isoID:
  cutMi=isoCut+'&&'+'(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  l1g_mcRate,maxMi,vertMi,horMi = rate_histo(
   uctmcntuple,cutMi,binning,L1G_CALIB_FACTOR,
   scaleMC,colorMi,markerMi,legend,
   'MC: Region + ID + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aM,bM)
# P
# Region + ID + (ISO-PU/A)
 if PUa:
  cutPUa=puACut+'&&'+'(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  PUaRate,maxPUa,vertPUa,horPUa = rate_histo(
   uctntuple,cutPUa,binning,L1G_CALIB_FACTOR,
   scale,colorPUa,markerPUa,legend,
   'Region + ID + (Iso - PU/'+str(puValA)+')<'+str(ISOTHRESHOLD),
   log,line,ptLine,aPa,bPa)
  info+='_PU'+str(puValA)
# P
# Region + ID + (ISO-PU/B)
 if PUb:
  cutPUb=puBCut+'&&'+'(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  PUbRate,maxPUb,vertPUb,horPUb = rate_histo(
   uctntuple,cutPUb,binning,L1G_CALIB_FACTOR,
   scale,colorPUb,markerPUb,legend,
   'Region + ID + (Iso - PU/'+str(puValB)+')<'+str(ISOTHRESHOLD),
   log,line,ptLine,aPb,bPb)
  info+='_PU'+str(puValB)
# N
# Region + ID without ISO
 if noIsoID:
  cutN='(!tauVeto[0]&&!mipBit[0])&&'+extraCut
  l1gNoIsoRate,maxN,vertN,horN = rate_histo(
   uctntuple,cutN,binning,L1G_CALIB_FACTOR,
   scale,colorN,markerN,legend,
   'Region + ID',
   log,line,ptLine,aN,bN)
# R
# Region
 if EGT:
  l1gEGTRate,maxR,vertR,horR = rate_histo(
   uctntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scale,colorR,markerR,legend,
   'Region (EG/Tau)',
   log,line,ptLine,aR,bR)
# C
# Current
 if lOne:
  l1Rate,maxC,vertC,horC = rate_histo(
   l1ntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scale,colorC,markerC,legend,
   'Current', 
   log,line,ptLine,aC,bC)
# Mc
# MC: Current
 if MC_lOne:
  l1mcRate,maxMc,vertMc,horMc = rate_histo(
   l1mcntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scaleMC,colorMc,markerMc,legend,
   'MC: Current', 
   log,line,ptLine,aM,bM)

 frame.SetMaximum(5*max(maxI,maxPa,maxPb,maxN,maxR,maxC,maxMi,maxMc,maxBi))
 frame.SetMinimum(100)
 legend.Draw()
 canvas.SaveAs(saveWhere+name+info+'.png')
######################################################################
##### RATES ##########################################################
######################################################################

######################################################################
###### DRAW PLOTS ####################################################
######################################################################
####################
# Resolution Plots #
####################
if resolutionPlots == True:
 binPt = [100,-2,2]
 binEta = [100,-1,1]

# ntuple,reco,l1,l1g,binning,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='',setLOG=False
 if resL1 == True:
  make_res_nrml(eff_rlx_eg_ntuple, 'recoPt', 'l1Pt', 'l1gPt', binPt,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal)
  make_res_angl(eff_rlx_eg_ntuple,'recoEta','l1Eta','l1gEta',binEta,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal)
  make_res_angl(eff_rlx_eg_ntuple,'recoPhi','l1Phi','l1gPhi',binEta,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal)

 if resReco == True:
  make_res_nrml(eff_rlx_eg_ntuple, 'recoPt', 'l1Pt', 'l1gPt', binPt,cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal)
  make_res_angl(eff_rlx_eg_ntuple,'recoEta','l1Eta','l1gEta',binEta,cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal)
  make_res_angl(eff_rlx_eg_ntuple,'recoPhi','l1Phi','l1gPhi',binEta,cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal)

####################
# Efficiency Plots #
####################
if efficiencyPlots == True:
 #binPt = [10,40,80] #l120
 binPt = [35,0,140]
 binVert=[10,0,35]
 binJetPt=[40,0,70]
 
# variable,
# binning,
# ntuple,
# l1ntuple=None,
# bntuple=None,
# bl1ntuple=None,
# recoPtCut='(2>1)',l1PtCut='(2>1)',l1gPtCut='(2>1)',
# isoCut='(2>1)',aPUCut='(2>1)',bPUCut='(2>1)',extraCut='(2>1)',
# isoID=False,l1b_IsoID=False,noIsoID=False,EGT=False,lOne=False,PUa=False,PUb=False,
# legExtra='',
# setLOG=False

 compare_efficiencies(
  'recoPt',
  binPt,
  eff_rlx_eg_ntuple, l1ntuple=eff_iso_eg_ntuple,
  bntuple=l1b_eff_rlx_eg_ntuple,bl1ntuple=l1b_eff_iso_eg_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1gPt >= '+str(l1ptVal)+')',
  isoCut='(l1gPt>=63||(l1gJetPt-l1gPt)/l1gPt<'+str(ISOTHRESHOLD)+')',
  aPUCut='(l1gPt>=63||(l1gJetPt-l1gPt-(l1gPU/'+str(puValA)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  bPUCut='(l1gPt>=63||(l1gJetPt-l1gPt-(l1gPU/'+str(puValB)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  #isoCut='((l1gJetPt-l1gPt)/l1gPt<'+str(ISOTHRESHOLD)+')',
  #aPUCut='((l1gJetPt-l1gPt-(l1gPU/'+str(puValA)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  #bPUCut='((l1gJetPt-l1gPt-(l1gPU/'+str(puValB)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  isoID=aIsoID,
  l1b_IsoID=al1b_IsoID,
  noIsoID=aNoIsoID,
  EGT=aEGT,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb
 )
 
 compare_efficiencies(
  'nPVs',
  binVert,
  eff_rlx_eg_ntuple, l1ntuple=eff_iso_eg_ntuple,
  bntuple=l1b_eff_rlx_eg_ntuple,bl1ntuple=l1b_eff_iso_eg_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1gPt >= '+str(l1ptVal)+')',
  isoCut='(l1gPt>=63||(l1gJetPt-l1gPt)/l1gPt<'+str(ISOTHRESHOLD)+')',
  aPUCut='(l1gPt>=63||(l1gJetPt-l1gPt-(l1gPU/'+str(puValA)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  bPUCut='(l1gPt>=63||(l1gJetPt-l1gPt-(l1gPU/'+str(puValB)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  #isoCut='((l1gJetPt-l1gPt)/l1gPt<'+str(ISOTHRESHOLD)+')',
  #aPUCut='((l1gJetPt-l1gPt-(l1gPU/'+str(puValA)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  #bPUCut='((l1gJetPt-l1gPt-(l1gPU/'+str(puValB)+'))/l1gPt<'+str(ISOTHRESHOLD)+')',
  isoID=aIsoID,
  l1b_IsoID=al1b_IsoID,
  noIsoID=aNoIsoID,
  EGT=aEGT,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

##############
# Rate Plots #
##############
if ratePlots == True:
 binRate = [36,0,80]

# l1ntuple,
# uctntuple,
# binning,
# l1mcntuple=None,
# uctmcntuple=None,
# bl1ntuple=None,
# buctntuple=None,
# filename='',
# setLOG=True,
# isoCut='(2>1)',puACut='(2>1)',puBCut='(2>1)',extraCut='(2>1)',
# ptLine=20,
# isoID=False,l1b_isoID=False,
# noIsoID=False,EGT=False,lOne=False,PUa=False,PUb=False,
# MC_isoID=False,MC_lOne=False,
# line=False

 make_rate_plot(rate_iso_l1_eg_ntuple,rate_rlx_l1g_eg_ntuple,binRate,
  l1mcntuple=mc_rate_iso_l1_eg_ntuple,
  uctmcntuple=mc_rate_rlx_l1g_eg_ntuple,
  bl1ntuple=None,
  buctntuple=l1b_rate_rlx_l1g_eg_ntuple,
  filename='',
  setLOG=True,
  isoCut='(( pt[0]>=63 && (jetPt[0]-regionPt[0])/regionPt[0]<100)||(pt[0]<63&&(jetPt[0]-pt[0])/pt[0]<'+str(ISOTHRESHOLD)+'))',
  puACut='(( pt[0]>=63 && (jetPt[0]-regionPt[0])/regionPt[0]<100)||(pt[0]<63&&(jetPt[0]-pt[0]-(pu[0]/'+str(puValA)+'))/pt[0]<'+str(ISOTHRESHOLD)+'))',
  puBCut='(( pt[0]>=63 && (jetPt[0]-regionPt[0])/regionPt[0]<100)||(pt[0]<63&&(jetPt[0]-pt[0]-(pu[0]/'+str(puValB)+'))/pt[0]<'+str(ISOTHRESHOLD)+'))',
  #isoCut='((jetPt[0]-pt[0])/pt[0]<'+str(ISOTHRESHOLD)+')',
  #puACut='((jetPt[0]-pt[0]-(pu[0]/'+str(puValA)+'))/pt[0]<'+str(ISOTHRESHOLD)+')',
  #puBCut='((jetPt[0]-pt[0]-(pu[0]/'+str(puValB)+'))/pt[0]<'+str(ISOTHRESHOLD)+')',
  ptLine=recoPtVal,
  isoID=aIsoID,
  l1b_isoID=al1b_IsoID,
  noIsoID=aNoIsoID,
  EGT=aEGT,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb,
  MC_isoID=aMC_isoID,
  MC_lOne=aMC_lOne,
  line = rateLine
 )
