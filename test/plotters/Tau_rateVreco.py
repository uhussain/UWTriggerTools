'''
Find the rate at an array of pT
Authors: T.M.Perry, M.Cepeda, E.K.Friis
'''
import ROOT
from ROOT import *
from array import array
import Tau_fitVals as fv

FixPlat=False
VarBins=False
AbsEff=True

LIso=3
LSB=50
ISOTHRESHOLD=0.20
#ZEROBIAS_RATE=15000000.00 #1e34
ZEROBIAS_RATE=30000000.00 #2e34
thresh=0.40
saveWhere='../plots/Tau%i'%(100*thresh)
if FixPlat:
 saveWhere+='_FixPlat'
if VarBins:
 saveWhere+='_BinVar'
if AbsEff:
 saveWhere+='_AbsEff'
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

mUIso,bUIso,pUIso,mU,bU,pU,mC,bC,pC=fv.getVals(thresh,FixPlat,VarBins,AbsEff)

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

#pT95 = array('d',[50,55,60,65,70,75,80])
pT95 = array('d',[30,35,40,45,50,55,60])

uctL1IsoPt = array('d',[])
uctL1Pt = array('d',[])
curL1Pt = array('d',[])
for pt in pT95:
 uctL1IsoPt.append(int(extrapolate(pt,mUIso,bUIso)))
 uctL1Pt.append(int(extrapolate(pt,mU,bU)))
 curL1Pt.append(int(extrapolate(pt,mC,bC)))

log.write('pT95 (RECO) Cut:   '+str(pT95)+'\n')
log.write('To have %i %% plateau efficiency at [pT%i] use L1pT cut at:\n'%(100*thresh,100*thresh))
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
frame.GetXaxis().SetTitle('Reco p_{T} Cut (at %i%% Plateau) [GeV]'%(100*thresh))
if AbsEff: frame.GetXaxis().SetTitle('Reco p_{T} Cut (at %i%% Efficiency) [GeV]'%(100*thresh))
frame.GetYaxis().SetTitle('Rate [Hz]')
frame.SetTitle('')

uctIsoGraph = ROOT.TGraph(nrPts,pT95,uctIsoRate)
uctIsoGraph.SetLineColor(ROOT.EColor.kBlack)
uctIsoGraph.SetMarkerColor(uIColor)
uctIsoGraph.SetLineWidth(2)
uctIsoGraph.SetMarkerStyle(uIMarker)
uctIsoGraph.SetMarkerSize(uISize)
uctIsoGraph.Draw('P')

uctGraph = ROOT.TGraph(nrPts,pT95,uctRate)
uctGraph.SetLineColor(ROOT.EColor.kBlack)
uctGraph.SetMarkerColor(uDColor)
uctGraph.SetLineWidth(2)
uctGraph.SetMarkerStyle(uDMarker)
uctGraph.SetMarkerSize(uDSize)
uctGraph.Draw('P')

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
legend.AddEntry(curGraph,'Current Tau, Plateau Efficiency: {:.2f}'.format(pC),'P')
legend.AddEntry(uctGraph,'Upgrade Tau, Plateau Efficiency: {:.2f}'.format(pU),'P')
legend.AddEntry(uctIsoGraph,'Upgrade Iso Tau, Plateau Efficiency: {:.2f}'.format(pUIso),'P')

frame.Draw()
curGraph.Draw('LP')
uctIsoGraph.Draw('LP')
uctGraph.Draw('LP')
legend.Draw()

tex.SetTextAlign(11)
tex.SetTextSize(0.05)
tex.DrawLatex(0.1,0.91,'Tau Rate v Reco. p_{T} Cut')
tex.SetTextAlign(13)
tex.SetTextSize(0.03)
tex.DrawLatex(0.1,0.89,'CMS Preliminary')
tex.SetTextAlign(31) # align right
tex.DrawLatex(0.9,0.91,'L=2e34 cm^{-2}s^{-1}')
save2 = raw_input ('Press Enter to Exit (type save to save)\n')
if save2 == 'save':
 can.SaveAs(saveWhere+'_rateVreco.png')
