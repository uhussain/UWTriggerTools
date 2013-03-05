'''
Find the rate at an array of pT
Authors: T.M.Perry, M.Cepeda, E.K.Friis
'''
import ROOT
from ROOT import *
from array import array
import Tau_fitVals as fv

FixPlat=False
VarBins=True

LIso=3
LSB=50
ISOTHRESHOLD=0.20
#ZEROBIAS_RATE=15000000.00 #1e34
ZEROBIAS_RATE=30000000.00 #2e34
thresh=0.85
saveWhere='../plots/Tau%i'%(100*thresh)
if FixPlat:
 saveWhere+='_FixPlat'
if VarBins:
 saveWhere+='_BinVar'

extraName=''

cColor=ROOT.EColor.kRed
uIColor=ROOT.EColor.kBlue
uDColor=ROOT.EColor.kGreen+3
cMarker=22
uIMarker=21
uDMarker=20
cSize=1.5
uISize=1.5
uDSize=1.5

## File 
rate_file = '../data/Evan/uct_rates_eic'+str(LIso)+'.root'
rate_ntuple_file = ROOT.TFile(rate_file)
rate_rlx_l1_spot = 'tauL1Rate/Ntuple' # Current
rate_rlx_l1g_spot = 'rlxTauUCTRate/Ntuple' # Upgrade
rate_rlx_l1_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1_spot)
rate_rlx_l1g_tau_ntuple = rate_ntuple_file.Get(rate_rlx_l1g_spot)

#### WHICH NTUPLE TO USE ###
uctNtuple = rate_rlx_l1g_tau_ntuple
curNtuple = rate_rlx_l1_tau_ntuple
scale = ZEROBIAS_RATE/rate_rlx_l1g_tau_ntuple.GetEntries()

log = open(saveWhere+extraName+'.log','w')
log.write('File: %s'%rate_file)
log.write('ISOTHRESHOLD: '+str(ISOTHRESHOLD)+'\n')
log.write('ZEROBIAS_RATE: '+str(ZEROBIAS_RATE)+'\n')
log.write('Scale: '+str(scale)+'  = total entries / ZEROBIAS_RATE\n\n')

tex = ROOT.TLatex()
tex.SetNDC(True)
tex.SetTextAlign(11)
tex.SetTextSize(0.03)

# Cut Strings
uctCutStringIso='((jetPt[0] - max(regionPt[0], pt[0]))/max(regionPt[0], pt[0])<'+str(ISOTHRESHOLD)+')'
uctCutString='(2>1)'
curCutString='(2>1)'

mUIso,bUIso,mU,bU,mC,bC=fv.getVals(thresh,FixPlat,VarBins)
log.write('mUIso: %s\n'%mUIso)
log.write('bUIso: %s\n'%bUIso)
log.write('mU: %s\n'%mU)
log.write('bU: %s\n'%bU)
log.write('mC: %s\n'%mC)
log.write('bC: %s\n\n'%bC)

# y = m x + b
def extrapolate(pT95,m,b):
 L1pT = (pT95-b)/m
 return L1pT

pT95 = array('d',[30,35,40,45,50,55,60])

uctL1IsoPt = array('d',[])
uctL1Pt = array('d',[])
curL1Pt = array('d',[])
for pt in pT95:
 uctL1IsoPt.append(int(extrapolate(pt,mUIso,bUIso)))
 uctL1Pt.append(int(extrapolate(pt,mU,bU)))
 curL1Pt.append(int(extrapolate(pt,mC,bC)))

log.write('pT95 (RECO) Cut:   '+str(pT95)+'\n')
log.write('To have 95% plateau efficiency at [pT95] use L1pT cut at:\n')
log.write('UCT Iso L1pT:      '+str(uctL1IsoPt)+'\n')
log.write('UCT (no Iso) L1pT: '+str(uctL1Pt)+'\n')
log.write('Current L1pT:      '+str(curL1Pt)+'\n\n')
log.write('Exrapolated as:\n')
log.write('UCT Iso :          L1pT = (pT95 - '+str(bUIso)+')/'+str(mUIso)+'\n')
log.write('UCT (no Iso):      L1pT = (pT95 - '+str(bU)+')/'+str(mU)+'\n')
log.write('Current (no Iso):  L1pT = (pt95 - '+str(bC)+')/'+str(mC)+'\n')

uctIsoRate=array('d',[]) 
uctRate=array('d',[]) 
curRate=array('d',[])
#
log.write('-----------------------\n')
log.write('For UCT:  '+str(uctNtuple.GetDirectory().GetName())+'\n\n')
log.write('UCT CutIso: '+uctCutStringIso+'\n\n')

for pt in uctL1IsoPt:
 uIsoNr = uctNtuple.GetEntries(uctCutStringIso+'&&(pt[0]>'+str(pt)+')')
 uIsoRate = uIsoNr*scale
 uctIsoRate.append(int(uIsoRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate Iso = '+str(uIsoRate)+'\n\n')
# 
log.write('-----------------------\n')
log.write('For UCT:  '+str(uctNtuple.GetDirectory().GetName())+'\n\n')
log.write('UCT Cut (no Iso): '+uctCutString+'\n\n')

for pt in uctL1Pt:
 uNr = uctNtuple.GetEntries(uctCutString+'&&(pt[0]>'+str(pt)+')')
 uRate = uNr*scale
 uctRate.append(int(uRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate (no Iso) = '+str(uRate)+'\n\n')
# 
log.write('-----------------------\n')
log.write('For Current:  '+str(curNtuple.GetDirectory().GetName())+'\n\n')
log.write('Current Cut: '+curCutString+'\n\n')

for pt in curL1Pt:
 cNr = curNtuple.GetEntries(curCutString+'&&(pt[0]>'+str(pt)+')')
 cRate = cNr*scale
 curRate.append(int(cRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate = '+str(cRate)+'\n\n')

# Drawing the plot: pT95 vs Rate at L1 Cut
can = ROOT.TCanvas('can','can',800,800)
can.SetLogy(True)

xmin = min(pT95)-5
xmax = max(pT95)+5
ymax = max(max(uctIsoRate),max(uctRate),max(curRate))*3
nrPts = len(pT95)

frame = TH1F('frame','',1,xmin,xmax)
frame.SetMaximum(ymax)
frame.SetMinimum(100)
frame.SetStats(False)
frame.GetXaxis().SetTitle('Reco p_{T} Cut [GeV]')
#frame.GetXaxis().SetTitle('L1 threshold [GeV]')
frame.GetYaxis().SetTitleOffset(1.3)
frame.GetYaxis().SetTitle('Rate [Hz]')
frame.SetTitle('')

uctIsoGraph = ROOT.TGraph(nrPts,pT95,uctIsoRate)
uctIsoGraph.SetLineColor(ROOT.EColor.kBlack)
uctIsoGraph.SetMarkerColor(uIColor)
uctIsoGraph.SetLineWidth(2)
uctIsoGraph.SetMarkerStyle(uIMarker)
uctIsoGraph.SetMarkerSize(uISize)
uctIsoGraph.Draw('P')
#labUI0=ROOT.TText(pT95[0],uctIsoRate[0]/1.7,'%0.f' %uctL1IsoPt[0])
#labUI0.SetTextSize(0.03)
#labUI1=ROOT.TText(pT95[1],uctIsoRate[1]/1.7,'%0.f' %uctL1IsoPt[1])
#labUI1.SetTextSize(0.03)
#labUI2=ROOT.TText(pT95[2],uctIsoRate[2]/1.7,'%0.f' %uctL1IsoPt[2])
#labUI2.SetTextSize(0.03)
#labUI3=ROOT.TText(pT95[3],uctIsoRate[3]/1.7,'%0.f' %uctL1IsoPt[3])
#labUI3.SetTextSize(0.03)
#labUI4=ROOT.TText(pT95[4],uctIsoRate[4]*1.3,'%0.f' %uctL1IsoPt[4])
#labUI4.SetTextSize(0.03)
#labUI5=ROOT.TText(pT95[5],uctIsoRate[5]*1.3,'%0.f' %uctL1IsoPt[5])
#labUI5.SetTextSize(0.03)

uctGraph = ROOT.TGraph(nrPts,pT95,uctRate)
uctGraph.SetLineColor(ROOT.EColor.kBlack)
uctGraph.SetMarkerColor(uDColor)
uctGraph.SetLineWidth(2)
uctGraph.SetMarkerStyle(uDMarker)
uctGraph.SetMarkerSize(uDSize)
uctGraph.Draw('P')
#labU0=ROOT.TText(pT95[0],uctRate[0]/1.7,'%0.f'%uctL1Pt[0])
#labU0.SetTextSize(0.03)
#labU1=ROOT.TText(pT95[1],uctRate[1]/1.7,'%0.f'%uctL1Pt[1])
#labU1.SetTextSize(0.03)
#labU2=ROOT.TText(pT95[2],uctRate[2]/1.7,'%0.f'%uctL1Pt[2])
#labU2.SetTextSize(0.03)
#labU3=ROOT.TText(pT95[3],uctRate[3]/1.7,'%0.f'%uctL1Pt[3])
#labU3.SetTextSize(0.03)
#labU4=ROOT.TText(pT95[4],uctRate[4]/1.7,'%0.f'%uctL1Pt[4])
#labU4.SetTextSize(0.03)
#labU5=ROOT.TText(pT95[5],uctRate[5]/1.7,'%0.f'%uctL1Pt[5])
#labU5.SetTextSize(0.03)

curGraph = ROOT.TGraph(nrPts,pT95,curRate)
curGraph.SetLineColor(ROOT.EColor.kBlack)
curGraph.SetMarkerColor(cColor)
curGraph.SetLineWidth(2)
curGraph.SetMarkerStyle(cMarker)
curGraph.SetMarkerSize(cSize)
curGraph.Draw('P')

legend = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
legend.SetFillColor(ROOT.EColor.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(curGraph,'Current Tau','P')
legend.AddEntry(uctGraph,'Upgrade Tau','P')
legend.AddEntry(uctIsoGraph,'Upgrade Iso Tau','P')

frame.Draw()
curGraph.Draw('LP')
uctIsoGraph.Draw('LP')
uctGraph.Draw('LP')
#labUI0.Draw()
#labUI1.Draw()
#labUI2.Draw()
#labUI3.Draw()
#labUI4.Draw()
#labUI5.Draw()
#labU0.Draw()
#labU1.Draw()
#labU2.Draw()
#labU3.Draw()
#labU4.Draw()
#labU5.Draw()
legend.Draw()

## For TDR
#tex.SetTextFont(42)
#tex.SetLineWidth(2)
#tex.SetTextSize(0.04)
#tex.SetTextAlign(31) # align right
#tex.DrawLatex(0.9,0.91,'#sqrt{s} = 8 TeV')
#tex.SetTextAlign(11) # align left
#tex.DrawLatex(0.1,0.91,'CMS 2012')
#tex.DrawLatex(0.3,0.91,'L=2E34 cm^{-2}s^{-1}')

tex.SetTextAlign(11)
tex.SetTextSize(0.05)
tex.DrawLatex(0.1,0.91,'Tau Rate v Reco. p_{T} Cut')
tex.SetTextAlign(13)
tex.SetTextSize(0.03)
tex.DrawLatex(0.1,0.89,'CMS Preliminary')
#tex.DrawLatex(0.15,0.15,'Labels are L1pT cut applied')
tex.SetTextAlign(11)
tex.SetTextAlign(31) # align right
tex.DrawLatex(0.9,0.91,'L=2e34 cm^{-2}s^{-1}')
save2 = raw_input ('Press Enter to Exit (type save to save)\n')
if save2 == 'save':
 can.SaveAs(saveWhere+'_rateVreco.png')
# can.SaveAs(saveWhere+'rateVreco.pdf')

#
#### PLOT offline pT vs Rate ###
#can2 = ROOT.TCanvas('can2','can2',800,800)
#can2.SetLogy(True)
#
##xmin = min(min(uctL1IsoPt),min(uctL1Pt),min(curL1Pt))-5
##xmax = max(max(uctL1IsoPt),max(uctL1Pt),max(curL1Pt))+5
##ymax = max(max(uctIsoRate),max(uctRate),max(curRate))*3
#xmin = min(min(uctL1IsoPt),min(uctL1Pt))-5
#xmax = max(max(uctL1IsoPt),max(uctL1Pt))+5
#ymax = max(max(uctIsoRate),max(uctRate))*3
#nrPts = len(pT95)
#
#frame2 = TH1F('frame2','',1,xmin,xmax)
#frame2.SetMaximum(ymax)
#frame2.SetMinimum(100)
#frame2.SetStats(False)
#frame2.GetXaxis().SetTitle('Online p_{T} [GeV]')
#frame2.GetYaxis().SetTitle('Hz (8Tev, 1e34)')
#frame2.SetTitle('')
#
#uctIsoGraph2 = ROOT.TGraph(nrPts,uctL1IsoPt,uctIsoRate)
#uctIsoGraph2.SetLineColor(ROOT.EColor.kBlack)
#uctIsoGraph2.SetMarkerColor(uIColor)
#uctIsoGraph2.SetLineWidth(2)
#uctIsoGraph2.SetMarkerStyle(uIMarker)
#uctIsoGraph2.SetMarkerSize(uISize)
#uctIsoGraph2.Draw('P')
#
#uctGraph2 = ROOT.TGraph(nrPts,uctL1Pt,uctRate)
#uctGraph2.SetLineColor(ROOT.EColor.kBlack)
#uctGraph2.SetMarkerColor(uDColor)
#uctGraph2.SetLineWidth(2)
#uctGraph2.SetMarkerStyle(uDMarker)
#uctGraph2.SetMarkerSize(uDSize)
#uctGraph2.Draw('P')
#
##curGraph2 = ROOT.TGraph(nrPts,curL1Pt,curRate)
##curGraph2.SetLineColor(ROOT.EColor.kBlack)
##curGraph2.SetMarkerColor(cColor)
##curGraph2.SetLineWidth(2)
##curGraph2.SetMarkerStyle(cMarker)
##curGraph2.SetMarkerSize(cSize)
##curGraph2.Draw('P')
#
#legend2 = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
#legend2.SetFillColor(ROOT.EColor.kWhite)
#legend2.SetBorderSize(0)
#legend2.AddEntry(uctIsoGraph2,'UCT: Region + Iso < 0.2','P')
#legend2.AddEntry(uctGraph2,'UCT: Region (no Iso)','P')
##legend2.AddEntry(curGraph2,'Current (no Iso)','P')
#
#frame2.Draw()
##curGraph2.Draw('LP')
#uctIsoGraph2.Draw('LP')
#uctGraph2.Draw('LP')
#legend2.Draw()
#
#tex.SetTextAlign(11)
#tex.SetTextSize(0.07)
#tex.DrawLatex(0.1,0.91,'Tau Rate v p_{T}')
#tex.SetTextAlign(31)
#tex.SetTextSize(0.03)
#tex.DrawLatex(0.9,0.91,'CMS Preliminary')
#
#save = raw_input ('Press Enter for reco (pT95) v Rate (type save to save)\n')
#if save == 'save':
# can2.SaveAs(saveWhere+'_L1_Rate.png')
