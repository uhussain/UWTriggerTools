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
aIso = True    # I region + Iso<ISOTHRESHOLD
aPUa = False      # P region + ISO + PU subtraction (puValA)
aPUb = False      # P region + ISO + PU subtraction (puValB)
aNoIso = True  # R region
aLOne = False     # C current
aLOneb = False     #L1b current
aMC_iso = False  # Mi MC region + ISO
a14MC_Noiso = False
a8MC_Noiso = False
aMC_lOne = False   # Mc MC current
al1b_Iso = True  # Bi L1b region + ISO
al1b_NoIso = True # Bn L1b region
Only_l1b = False
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
recoPtVal= 35
ISOTHRESHOLD=0.20
puValA = 5
puValB = 9
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0
ZEROBIAS_RATE=15000000.00
saveWhere='../plots/Tau_Stage1B_NewVar/'
#saveWhere=''

########
# File #
########
##The names of these files change frequently

# Efficiency
eff_ntuple = '/afs/hep.wisc.edu/cms/tperry/CMSSW_6_0_1_PostLS1v2_patch4/src/L1Trigger/UCT2015/test/data/EvanL1b/uct_tau_efficiency.root'
#L1
eff_rlx_spot = 'rlxTauEfficiency/Ntuple'
eff_iso_spot = 'isoTauEfficiency/Ntuple'
eff_ntuple_file = ROOT.TFile(eff_ntuple)
eff_rlx_tau_ntuple = eff_ntuple_file.Get(eff_rlx_spot)
eff_iso_tau_ntuple = eff_ntuple_file.Get(eff_iso_spot)
# L1b
l1b_eff_rlx_spot = 'rlxTauEfficiencyStage1B/Ntuple'
l1b_eff_iso_spot = 'isoTauEfficiencyStage1B/Ntuple'
l1b_eff_rlx_tau_ntuple = eff_ntuple_file.Get(l1b_eff_rlx_spot)
l1b_eff_iso_tau_ntuple = eff_ntuple_file.Get(l1b_eff_iso_spot)

# Rate
rate_ntupleb = '../data/uct_rates_14mc_eic3.root'

rate_ntuple = '../data/uct_rates_14mc_eic3.root'
rate_rlx_l1_spot = 'tauL1Rate/Ntuple'
rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple'
rate_ntuple_file = ROOT.TFile(rate_ntuple)
rate_rlx_l1_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1_spot)
rate_rlx_l1g_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1g_spot)
# L1b
l1b_rate_rlx_l1g_spot = 'rlxTauUCTRateStage1B/Ntuple'
l1b_rate_iso_l1g_spot = 'isoTauUCTRateStage1B/Ntuple'
l1b_rate_rlx_l1g_tau_ntuple = rate_ntuple_file.Get(l1b_rate_rlx_l1g_spot)
l1b_rate_iso_l1g_tau_ntuple = rate_ntuple_file.Get(l1b_rate_iso_l1g_spot)
#14TeV MC Rate
mc14_rate_ntuple = '../data/uct_rates_14mc_eic3.root'

mc14_rate_rlx_l1_spot = 'tauL1Rate/Ntuple'
mc14_rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple'
mc14_rate_ntuple_file = ROOT.TFile(mc14_rate_ntuple)
mc14_rate_rlx_l1_tau_ntuple = mc14_rate_ntuple_file.Get(mc14_rate_rlx_l1_spot)
mc14_rate_rlx_l1g_tau_ntuple = mc14_rate_ntuple_file.Get(mc14_rate_rlx_l1g_spot)

#14TeV MC Rate
mc8_rate_ntuple= '/afs/hep.wisc.edu/cms/tperry/CMSSW_6_0_1_PostLS1v2_patch4/src/L1Trigger/UCT2015/test/data/uct_rates_eic'+str(LIso)+'.root'
mc8_rate_rlx_l1_spot = 'tauL1Rate/Ntuple'
mc8_rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple'
mc8_rate_ntuple_file = ROOT.TFile(mc8_rate_ntuple)
mc8_rate_rlx_l1_tau_ntuple = mc8_rate_ntuple_file.Get(mc8_rate_rlx_l1_spot)
mc8_rate_rlx_l1g_tau_ntuple = mc8_rate_ntuple_file.Get(mc8_rate_rlx_l1g_spot)
store = ROOT.TFile(saveWhere+'store.root','RECREATE')

name=''
if aIso: name+='I'
if aPUa or aPUb: name+='P'
if aNoIso: name+='N'
if aLOne: name+='C'
if aLOneb: name+='Bc'
if aMC_iso: name+='Mi'
if a14MC_Noiso: name+='Mn14'
if a8MC_Noiso: name+='Mn8'
if aMC_lOne: name+='Mc'
if al1b_Iso: name+='Bi'
if al1b_NoIso: name+='Bn'
extraName=''

log = open(saveWhere+name+'_reco'+str(recoPtVal)+'_l1'+str(l1ptVal)+'_iso'+str(ISOTHRESHOLD)+extraName+'.log','w')
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
colorN=ROOT.EColor.kGreen+3
markerN=21
colorR=ROOT.EColor.kBlack
markerR=22
colorC=ROOT.EColor.kRed
markerC=23
colorPUa=ROOT.EColor.kRed
markerPUa=33
colorPUb=ROOT.EColor.kYellow+3
markerPUb=23
colorMi=ROOT.EColor.kAzure
markerMi=24
colorMc=ROOT.EColor.kViolet-7
markerMc=25
colorBi=ROOT.EColor.kBlack
markerBi=24
colorBn=ROOT.EColor.kCyan
markerBn=25
colorBc=ROOT.EColor.kMagenta+2
markerBc=23

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title='',calFactor=1):
 ''' Plot a variable using draw and return the histogram '''
 draw_string = "%s * %0.2f>>htemp(%s)" % (variable,calFactor, ", ".join(str(x) for x in binning))
 print draw_string
 print selection
 tree.Draw(draw_string, selection, "goff")
 output_histo = ROOT.gDirectory.Get("htemp").Clone()
 output_histo.GetXaxis().SetTitle(xaxis)
 output_histo.SetTitle(title)
 return output_histo

######################################################################
##### RESOLUTION #####################################################
######################################################################
def make_res_nrml(
ntuple,bntuple,reco,l1,l1g,binning,cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='',extraCut='',extraCutb='',setLOG=False):
 '''variable: (reco-l1g)/reco  with cut l1g > l1pt'''
 canvas.SetLogy(setLOG)
 info = 'RES_'+str(reco)+'_'+cutPtVar+'Cut'
 
 l1gplot = make_plot(
  ntuple, '('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+' * '+str(l1g)+'))/'+str(reco), 
  'l1gMatch&&'+cutPtVarg+'>'+str(cutPt)+extraCut,binning
 )

 l1gplot.SetTitle('Resolution')
 l1gplot.GetXaxis().SetTitle('(recoPt- l1gPt)/recoPt')
 l1gplot.SetLineColor(ROOT.EColor.kBlue)

 l1bplot = make_plot(
 bntuple, #'('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+' * '+str(l1g)+'))/'+str(reco),
  #'('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+')*(max(max(l1g2RegionEt*l1g2RegionPattern,l1g3RegionEt*l1g3RegionPattern),l1g4RegionPattern*l1g4RegionEt) + (!l1g2RegionPattern)*max(l1gRegionEt,'+str(l1g)+')))/'+str(reco), 
  '('+str(reco)+' - ('+str(L1G_CALIB_FACTOR)+')*(l1g2RegionEt*l1g2RegionPattern + (!l1g2RegionPattern)*max(l1gRegionEt,'+str(l1g)+')))/'+str(reco),
  'l1gMatch&&'+cutPtVarg+'>'+str(cutPt)+extraCutb,binning
 )
 
 l1bplot.Scale(1/l1bplot.Integral())
 l1bplot.SetLineColor(ROOT.EColor.kRed)
 legend = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
 legend.SetFillColor(ROOT.EColor.kWhite)
 legend.SetBorderSize(0)
 legend.AddEntry(l1gplot,'upgrade')
 legend.AddEntry(l1bplot,'Stage 1B upgrade')

 l1gplot.Draw()
 l1bplot.Draw('sames')
 #l1gplot.SetMaximum(1.1*max(l1gplot.GetMaximum(),l1bplot.GetMaximum()))
 l1gplot.SetMaximum(0.1)
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
 bntuple=None,
 bl1ntuple=None,
 recoPtCut='(2>1)',l1PtCut='(2>1)',l1gPtCut='(2>1)',
 isoCut='(2>1)',aPUCut='(2>1)',bPUCut='(2>1)',extraCut='(2>1)',
 Iso=False,noIso=False,lOne=False,PUa=False,PUb=False,l1b_Iso=False,l1b_noIso=False,
 legExtra='',
 setLOG=False
):
 ''' Returns a (L1, L1G) tuple of TGraphAsymmErrors

 If [l1ntuple] is None, use [ntuple].  If not None, separate ntuples will
 be used for L1G and L1.  The use case for this is to use Rlx UCT with Iso L1
 '''
 if l1ntuple is None:
  l1ntuple = ntuple
 variableb = variable
 if variable == 'l1gPt':
  variableb = '(('+str(L1G_CALIB_FACTOR)+')*(l1g2RegionEt*l1g2RegionPattern + (!l1g2RegionPattern)*max(l1gRegionEt,l1gPt)))'
 cutD = recoPtCut+'&&'+extraCut
 if not (Only_l1b):
  denom = make_plot(
   ntuple,variable,
   cutD,
   binning
  )
 denomB = make_plot(
  bntuple,variableb,
  cutD,
  binning
 )

 log.write('_____________________________\n')
 log.write('-------- Efficiency ---------\n\n')
 log.write('File: '+eff_ntuple+'\n')
 log.write('Variable: '+variable+'\n\n')
 if not (Only_l1b):
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
 tex.DrawLatex(0.1,0.91,'Tau '+variable+' Efficiency')
 legend = ROOT.TLegend(0.4,0.35,0.89,0.15,'','brNDC')
 legend.SetFillColor(0)
 legend.SetBorderSize(0)
 legend.SetHeader(legExtra)
 
 info = '_EFF_'+variable+'_reco'+str(recoPtVal)+'_l1'+str(l1ptVal)
# I
#ntuple with Region and Isolation
 if Iso:
  cutI=recoPtCut+'&&'+extraCut+'&&'+isoCut +'&&'+l1gPtCut+'&&l1gMatch'
  l1g=effi_histo(ntuple,variable,cutI,binning,denom,
   'Region + Isolation < '+str(ISOTHRESHOLD),legend,
   colorI,markerI,log)
  l1g.SetName('l1g')
  l1g.Write()
  #print 'Iso facts'
  #print cutI
  #print variable
  #print l1g.Integral()
# Bi
#L1B ntuple with Region and Isolation
 if l1b_Iso:
  cutBi=recoPtCut+'&&'+extraCut+'&&'+isoCut +'&&'+l1gPtCut+'&&l1gMatch'
  bl1g=effi_histo(bntuple,variableb,cutBi,binning,denomB,
   'L1b: Region + Isolation < '+str(ISOTHRESHOLD),legend,
   colorBi,markerBi,log)
  bl1g.SetName('bl1g')
  bl1g.Write()
  #print 'bl1g facts'
  #print cutBi
  #print bl1g.Integral()
# P
#ntuple with Region and Isolation
 if PUa:
  cutPUa=recoPtCut+'&&'+extraCut+'&&'+aPUCut+'&&'+l1gPtCut+'&&l1gMatch'
  aPU=effi_histo(ntuple,variable,cutPUa,binning,denom,
   'Region + (Isolation-PU/'+str(puValA)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUa,markerPUa,log)
  info += '_PU'+str(puValA)
  aPU.SetName('aPU')
  aPU.Write()
# P
#ntuple with Region and Isolation
 if PUb:
  cutPUb=recoPtCut+'&&'+extraCut+'&&'+bPUCut+'&&'+l1gPtCut+'&&l1gMatch'
  bPU=effi_histo(ntuple,variable,cutPUb,binning,denom,
   'Region + (Isolation-PU/'+str(puValB)+' < '+str(ISOTHRESHOLD)+')',legend,
   colorPUb,markerPUb,log)
  info += '_PU'+str(puValA)
  bPU.SetName('bPU')
  bPU.Write()
# N
#ntuple with Region and without Isolation
 if noIso:
  cutN=recoPtCut+'&&'+extraCut+'&&'+l1gPtCut+'&& l1gMatch'
  l1g_noiso=effi_histo(ntuple,variable,cutN,binning,denom,
  'Region',legend,
   colorN,markerN,log) 
  l1g_noiso.SetName('l1g_noiso')
  l1g_noiso.Write()
# Bn
# L1b: bntuple with Region and without Isolation
 if l1b_noIso:
  cutBn=recoPtCut+'&&'+extraCut+'&&'+l1gPtCut+'&& l1gMatch'
  bl1g_noiso=effi_histo(bntuple,variableb,cutBn,binning,denomB,
  'L1b: Region',legend,
   colorBn,markerBn,log) 
  bl1g_noiso.SetName('bl1g_noiso')
  bl1g_noiso.Write()
# R
#Current trigger
 if lOne:
  cutC=recoPtCut+'&&'+extraCut+'&&'+l1PtCut+'&&l1Match'
  l1=effi_histo(l1ntuple,variable,cutC,binning,denom,
  'Current',legend,
  colorC,markerC,log)
  l1.SetName('l1')
  l1.Write()
 if aLOneb:
  cutBc = recoPtCut+'&&'+extraCut+'&&'+l1PtCut+'&&l1gMatch'
  l1b = effi_histo(bl1ntuple,variableb,cutBc,binning,denomB,'L1b: Current',legend,colorBc,markerBc,log)
  l1b.SetName('l1b')
  l1b.Write() 
 
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

def rate_histo(ntuple,cut,binning,calibfactor,scale,color,marker,leg,title,logg,line,ptLine,w,s,stage1B=False):
 if not stage1B:
 	pt = make_plot(ntuple,'pt[0]',cut,binning,'','',calibfactor)
 else:
 	pt = make_plot(ntuple,'(region2Disc.totalEt[0]*region2Disc.patternPass[0] + (!region2Disc.patternPass[0])*max(regionPt[0],pt[0]))',cut,binning,'','',calibfactor)
 rate = make_l1_rate(pt,color,marker)
 rate.Scale(scale)
 rate.Draw('phsame')
 leg.AddEntry(rate,title,'pe')
 maxx = rate.GetMaximum()
 binn = rate.GetXaxis().FindBin(ptLine)
 if ptLine == 0:
	binn = 1
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
 uctmc8ntuple=None,
 uctmc14ntuple=None,
 bl1ntuple=None,
 buctntuple=None,
 filename='',
 setLOG=True,
 isoCut='(2>1)',puACut='(2>1)',puBCut='(2>1)',extraCut='(2>1)',
 ptLine=20,
 Iso=False,noIso=False,lOne=False,PUa=False,PUb=False,
 MC_Iso=False,MC_lOne=False,l1b_Iso=False,l1b_noIso=False,
 line=False
 ):

 info = '_RATE'
 scale = ZEROBIAS_RATE/uctntuple.GetEntries("lumi>40")
 scale8MC = ZEROBIAS_RATE/uctmc8ntuple.GetEntries()
 scale14MC = ZEROBIAS_RATE/uctmc14ntuple.GetEntries()
 scaleB = ZEROBIAS_RATE/buctntuple.GetEntries()
 
 canvas.SetLogy(setLOG)
 frame = ROOT.TH1F('frame','frame',*binning)
 frame.Draw()
 frame.SetTitle('')
 frame.GetYaxis().SetTitle('Hz (8TeV,1E34)')
 frame.GetXaxis().SetTitle('p_{T}')
 frame.SetMaximum(10000000)
 frame.SetMinimum(100)
 tex.DrawLatex(0.1,0.91,'Tau Rate')
 legend = ROOT.TLegend(0.5,0.7,0.89,0.89,'','brNDC')
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
 aBn=3
 bBn=3
 aBc=3
 bBc=3
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
 maxBi = 1
 maxBn = 1
 maxBc = 1
# I
# Isolated
 if Iso:
  cutI=isoCut+'&&'+extraCut
  l1gRate,maxI,vertI,horI = rate_histo(
   uctntuple,cutI,binning,L1G_CALIB_FACTOR,
   scale,colorI,markerI,legend,
   'Region + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aI,bI)
  info+='_iso_'+str(ISOTHRESHOLD)
# Bi
# L1B Region + Iso
 if l1b_Iso:
  cutBi=isoCut+'&&'+extraCut
  bl1gRate,maxBi,vertBi,horBi = rate_histo(
   buctntuple,cutBi,binning,L1G_CALIB_FACTOR,
   scaleB,colorBi,markerBi,legend,
   'L1b: Region + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aBi,bBi, stage1B=True)
  info+='_L1B_iso_'+str(ISOTHRESHOLD)
# P
# Isolated - PU/A
 #if l1b_isoID:
 # cutBi=isoCut+'&&'+extraCut
 # bl1gRate,maxBi,vertBi,horBi = rate_histo(
 #  buctntuple,cutBi,binning,L1G_CALIB_FACTOR,
 #  scaleB,colorBi,markerBi,legend,
 #  'L1B:Region + ID+ Isolation <'+str(ISOTHRESHOLD),
 #  log,line,ptLine,aBi,bBi)
 # info+='_iso_'+str(ISOTHRESHOLD)

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
 if noIso:
  cutN=extraCut
  l1gNoIsoRate,maxN,vertN,horN = rate_histo(
   uctntuple,cutN,binning,L1G_CALIB_FACTOR,
   scale,colorN,markerN,legend,
   'Region',
   log,line,ptLine,aN,bN)
# L1B: Non Isolated
 if l1b_noIso:
  cutBn=extraCut
  bl1gNoIsoRate,maxBn,vertBn,horBn = rate_histo(
   buctntuple,cutBn,binning,L1G_CALIB_FACTOR,
   scaleB,colorBn,markerBn,legend,
   'L1b: Region',
   log,line,ptLine,aBn,bBn,stage1B=True)
# Current
 if lOne:
  l1Rate,maxC,vertC,horC = rate_histo(
   l1ntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scale,colorC,markerC,legend,
   'Current', 
   log,line,ptLine,aC,bC)
 if aLOneb:
  l1bRate,maxBc,vertBc,horBc = rate_histo(
   bl1ntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scale,colorBc,markerBc,legend,
   'Stage 1B Current',
   log,line,ptLine,aBc,bBc)
# Mi
# MC: Region + ID + ISO
 if MC_Iso:
  cutMi=isoCut+'&&'+extraCut
  l1g_mcRate,maxMi,vertMi,horMi = rate_histo(
   uctmcntuple,cutMi,binning,L1G_CALIB_FACTOR,
   scaleMC,colorMi,markerMi,legend,
   'MC: Region + Isolation <'+str(ISOTHRESHOLD),
   log,line,ptLine,aM,bM)
 if a8MC_Noiso:
  cut14Mn=extraCut
  l1g_mc8Rate,maxMi,vertMi,horMi = rate_histo(
   uctmc8ntuple,cut14Mn,binning,L1G_CALIB_FACTOR,
   scale8MC,colorMi,markerMi,legend,
   '8 TeV MC: Region ',
   log,line,ptLine,aM,bM)

 if a14MC_Noiso: 
  cut8Mn = extraCut
  l1gmc14Rate,maxMi,vertMi,horMi = rate_histo(
   uctmc14ntuple,cut14Mn,binning,L1G_CALIB_FACTOR,
   scale14MC,colorMc,markerMc,legend,
   '14 TeV MC: Region ',
   log,line,ptLine,aM,bM)
 

 
   
# Mc
# MC: Current
 if MC_lOne:
  l1mcRate,maxMc,vertMc,horMc = rate_histo(
   l1mcntuple,extraCut,binning,L1G_CALIB_FACTOR,
   scaleMC,colorC,markerC,legend,
   'MC: Current', 
   log,line,ptLine,aM,bM)
 #print(str(maxI)+' '+str(maxN)+' '+str(maxC)+' '+str(maxMc)+' '+str(maxMi))
 frame.SetMaximum(5*max(maxI,maxN,maxC,maxMc,maxMi,maxBn,maxBi))
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
  make_res_nrml(eff_rlx_tau_ntuple, l1b_eff_rlx_tau_ntuple,'recoPt', 'l1Pt', 'l1gPt', binPt, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoEta','l1Eta','l1gEta',binEta, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoPhi','l1Phi','l1gPhi',binEta, cutPtVarg='l1gPt',cutPtVar='l1Pt',cutPt=l1ptVal,filename='')

 if resReco == True:
  make_res_nrml(eff_rlx_tau_ntuple, l1b_eff_rlx_tau_ntuple,'recoPt', 'l1Pt', 'l1gPt', binPt, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,extraCut = '',extraCutb = '',filename='2PatternRegionXaxis')
#extraCutb = '&&l1Pt>20&&l1g2RegionPattern==1&&l1g3RegionPattern==1&&l1g4RegionPattern==1&&l1g2RegionEt>10&&l1g3RegionEt>10&&l1g4RegionEt>10&&l1g2RegionEtEcal>10&&l1g3RegionEtEcal>10&&l1g4RegionEtEcal>10&&(l1g2RegionMips==0||l1g2RegionMips==1)&&(l1g3RegionMips==0||l1g3RegionMips==1)&&(l1g2RegionMips==0||l1g2RegionMips==1)',filename='l1Pt20l1gPatternl1gRegion10l1gRegionEcall1gMips')
  make_res_angl(eff_rlx_tau_ntuple,'recoEta','l1Eta','l1gEta',binEta, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,filename='')
  make_res_angl(eff_rlx_tau_ntuple,'recoPhi','l1Phi','l1gPhi',binEta, cutPtVarg='recoPt',cutPtVar='recoPt',cutPt=recoPtVal,filename='')
####################
# Efficiency Plots #
####################
if efficiencyPlots == True:

 binPt = [14,0,140]
 binVert=[10,0,40]
 binJetPt=[40,0,70]
 binExtra=[7,0,140]
 
# variable,
# binning,
# ntuple,
# l1ntuple=None,
# bntuple=None,
# bl1ntuple=None,
# recoPtCut='(2>1)',l1PtCut='(2>1)',l1gPtCut='(2>1)',
# isoCut='(2>1)',aPUCut='(2>1)',bPUCut='(2>1)',extraCut='(2>1)',
# Iso=False,noIso=False,lOne=False,PUa=False,PUb=False,l1b_Iso=False,l1b_noIso=False,
# legExtra='',
# setLOG=False

 compare_efficiencies(
  'l1gPt',
  binPt,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '2>1',
  l1gPtCut = '2>1',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra=''
 )
 
 compare_efficiencies(
  'nPVs',
  binVert,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1gPt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)


'''
 compare_efficiencies(
  'l1g2ndRegionEtEM',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1gEmClusterCenterEt',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1gEmClusterStripEt',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1g2ndRegionEtEM',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1gJetPtEM',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1gRegionEtEM',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)

 compare_efficiencies(
  'l1gRegionEt',
  binExtra,
  eff_rlx_tau_ntuple, eff_iso_tau_ntuple,
  bntuple=l1b_eff_rlx_tau_ntuple,bl1ntuple=l1b_eff_rlx_tau_ntuple,
  recoPtCut = '(recoPt >= '+str(recoPtVal)+')',
  l1PtCut = '(l1Pt >= '+str(l1ptVal)+')',
  l1gPtCut = '(l1Pt > '+str(l1ptVal)+')',
  isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  aPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValA)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  bPUCut='((l1gJetPt-max(l1gRegionEt,l1gPt)-(l1gPU/'+str(puValB)+'))/max(l1gRegionEt,l1gPt) <'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  PUa=aPUa,
  PUb=aPUb,
  legExtra='Reco Pt > '+str(recoPtVal)
)
'''

##############
# Rate Plots #
##############
if ratePlots == True:
 #binRate = [12,25,85]
 binRate = [14,5,140]

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
# Iso=False,noIso=False,lOne=False,PUa=False,PUb=False,
# MC_Iso=False,MC_lOne=False,l1b_Iso=False,l1b_noIso=False,
# line=False

 make_rate_plot(rate_rlx_l1_tau_ntuple,rate_rlx_l1g_tau_ntuple,binRate,
  l1mcntuple=mc14_rate_rlx_l1_tau_ntuple,
  uctmc8ntuple=mc8_rate_rlx_l1g_tau_ntuple,
  uctmc14ntuple=mc14_rate_rlx_l1g_tau_ntuple,
  bl1ntuple=l1b_rate_rlx_l1g_tau_ntuple,
  #bl1ntuple=None,
  buctntuple=l1b_rate_rlx_l1g_tau_ntuple,
  filename='',
  setLOG=True,
  isoCut='((jetPt[0] - max(regionPt[0], pt[0]))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  puACut='((jetPt[0] - max(regionPt[0], pt[0])-(pu[0]/'+str(puValA)+'))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  puBCut='((jetPt[0] - max(regionPt[0], pt[0])-(pu[0]/'+str(puValB)+'))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')',
  Iso=aIso,
  noIso=aNoIso,
  lOne=aLOne,
  PUa=aPUa,
  PUb=aPUb,
  MC_Iso=aMC_iso,
  MC_lOne=aMC_lOne,
  l1b_Iso=al1b_Iso,
  l1b_noIso=al1b_NoIso,
  line = rateLine,
  #ptLine=recoPtVal
  ptLine = 60
 )
