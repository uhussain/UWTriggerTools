'''
Find the rate at an array of pT
Authors: T.M.Perry, M.Cepeda
'''
import ROOT
from ROOT import *
from array import array

LIso=3
LSB=50
ISOTHRESHOLD=0.20
ZEROBIAS_RATE=15000000.00
saveWhere='../plots/EG_RateThresh'
extraName=''

cColor=ROOT.EColor.kRed
uIColor=ROOT.EColor.kBlue
uDColor=ROOT.EColor.kGreen+3
cMarker=22
uIMarker=20
uDMarker=20
cSize=1
uISize=1
uDSize=1

## File 
rate_filename = '../data/uct_rates_eic3.root'
rate_file = ROOT.TFile(rate_filename)
# L1
rate_rlx_l1g_spot = 'rlxEGUCTRate/Ntuple' # USE THIS`
rate_iso_l1g_spot = 'isoEGUCTRate/Ntuple' # for testing purposes
rate_rlx_l1g_eg_ntuple = rate_file.Get(rate_rlx_l1g_spot)
rate_iso_l1g_eg_ntuple = rate_file.Get(rate_iso_l1g_spot)
# Current
rate_rlx_l1_spot = 'rlxEGL1Rate/Ntuple' # The ONE TO USE
rate_iso_l1_spot = 'isoEGL1Rate/Ntuple' # formerly used this one
rate_rlx_l1_eg_ntuple = rate_file.Get(rate_rlx_l1_spot)
rate_iso_l1_eg_ntuple = rate_file.Get(rate_iso_l1_spot)

#### WHICH NTUPLE TO USE ###
uctNtuple = rate_rlx_l1g_eg_ntuple
curNtuple = rate_rlx_l1_eg_ntuple
scale = ZEROBIAS_RATE/rate_rlx_l1g_eg_ntuple.GetEntries()

log = open(saveWhere+extraName+'.log','w')
log.write('LIso: '+str(LIso)+'\n')
log.write('LSB: '+str(LSB)+'\n')
log.write('ISOTHRESHOLD: '+str(ISOTHRESHOLD)+'\n')
log.write('ZEROBIAS_RATE: '+str(ZEROBIAS_RATE)+'\n')
log.write('Scale: '+str(scale)+'  = total entries / ZEROBIAS_RATE\n\n')

tex = ROOT.TLatex()
tex.SetNDC(True)
tex.SetTextAlign(11)
tex.SetTextSize(0.03)

#### WHICH ISOLATION + ID CUT ON UCT ####
# ID + ISO both with switch at 63GeV
#uctCutString='((pt[0]>=63&&((jetPt[0]-regionPt[0])/regionPt[0]<100))||(pt[0]<63&&((jetPt[0]-pt[0])/pt[0])<'+str(ISOTHRESHOLD)+')&&((!tauVeto[0]&&!mipBit[0])||pt[0]>=63))'

# ID w/o switch + Iso with switch at 63GeV
#uctCutString='((pt[0]>=63&&((jetPt[0]-regionPt[0])/regionPt[0]<100))||(pt[0]<63&&((jetPt[0]-pt[0])/pt[0])<'+str(ISOTHRESHOLD)+')&&((!tauVeto[0]&&!mipBit[0])))'

# Iso w/ switch at 63 + ID w/o switch
uctCutStringIsoId='(( pt[0]>=63 && (jetPt[0]-regionPt[0])/regionPt[0]<100)||(pt[0]<63&&(jetPt[0]-pt[0])/pt[0]<'+str(ISOTHRESHOLD)+'))&&(!tauVeto[0]&&!mipBit[0])'

# ID w/ switch (No Iso)
#uctCutString='(((!tauVeto[0]&&!mipBit[0])||pt[0]>=63))'

# ID w/o switch (No Iso)
uctCutStringId='(!tauVeto[0]&&!mipBit[0])'

# No Cut
#uctCutString='(2>1)'

#### CURRENT CUT STRING ###
# No Cut
curCutString='(2>1)'

# values for extrapolation
mUIsoId = 0.91
bUIsoId = 11
mUId = 1.2
bUId = 5.8
mC = 1.2
bC = 4.5

# y = m x + b
def extrapolate(l1pt,m,b):
 reco = (l1pt*m)+b
 return reco

l1CutPt = array('d',[20,25,30,35,40,45,50])

uctRecoIsoIdPt = array('d',[])
uctRecoIdPt = array('d',[])
curRecoPt = array('d',[])
for pt in l1CutPt:
 uctRecoIsoIdPt.append(int(extrapolate(pt,mUIsoId,bUIsoId)))
 uctRecoIdPt.append(int(extrapolate(pt,mUId,bUId)))
 curRecoPt.append(int(extrapolate(pt,mC,bC)))

log.write('L1 pT Cut:              '+str(l1CutPt)+'\n')
log.write('To Reach 95% plateau efficiency use reco pT cuts:\n')
log.write('UCT Iso + Id RecoPt:    '+str(uctRecoIsoIdPt)+'\n')
log.write('UCT Id (no Iso) RecoPt: '+str(uctRecoIdPt)+'\n')
log.write('Current RecoPt:         '+str(curRecoPt)+'\n\n')
log.write('Exrapolated as:\n')
log.write('UCT Iso + Id:      RecoPt = L1pT * '+str(mUIsoId)+' + '+str(bUIsoId)+'\n')
log.write('UCT Id (no Iso):   RecoPt = L1pT * '+str(mUId)+' + '+str(bUId)+'\n')
log.write('Current (no Iso):  RecoPt = L1pT * '+str(mC)+' + '+str(bC)+'\n')

uctIsoIdRate=array('d',[]) 
uctIdRate=array('d',[]) 
curRate=array('d',[])

log.write('-----------------------\n')
log.write('For UCT:  '+str(uctNtuple.GetDirectory().GetName())+'\n\n')
log.write('UCT CutIsoID: '+uctCutStringIsoId+'\n\n')

for pt in uctRecoIsoIdPt:
 uIsoIdNr = uctNtuple.GetEntries(uctCutStringIsoId+'&&(pt[0]>'+str(pt)+')')
 uIsoIdRate = uIsoIdNr*scale
 uctIsoIdRate.append(int(uIsoIdRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate Iso + Id = '+str(uIsoIdRate)+'\n\n')
 
log.write('-----------------------\n')
log.write('For UCT:  '+str(uctNtuple.GetDirectory().GetName())+'\n\n')
log.write('UCT CutID: '+uctCutStringId+'\n\n')

for pt in uctRecoIdPt:
 uIdNr = uctNtuple.GetEntries(uctCutStringId+'&&(pt[0]>'+str(pt)+')')
 uIdRate = uIdNr*scale
 uctIdRate.append(int(uIdRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate Id (no Iso) = '+str(uIdRate)+'\n\n')
 
log.write('-----------------------\n')
log.write('For Current:  '+str(curNtuple.GetDirectory().GetName())+'\n\n')
log.write('Current Cut: '+curCutString+'\n\n')

for pt in curRecoPt:
 cNr = curNtuple.GetEntries(curCutString+'&&(pt[0]>'+str(pt)+')')
 cRate = cNr*scale
 curRate.append(int(cRate))
 log.write('At pT = '+str(pt)+'\n')
 log.write('Rate = '+str(cRate)+'\n\n')

### PLOT pT vs Rate ###
can2 = ROOT.TCanvas('can2','can2',800,800)
can2.SetLogy(True)

xmin = min(min(uctRecoIsoIdPt),min(uctRecoIdPt),min(curRecoPt))-10
xmax = max(max(uctRecoIsoIdPt),max(uctRecoIdPt),max(curRecoPt))+10
ymax = max(max(uctIsoIdRate),max(uctIdRate),max(curRate))*3
nrPts = len(l1CutPt)

frame2 = TH1F('frame2','',1,xmin,xmax)
frame2.SetMaximum(ymax)
frame2.SetMinimum(100)
frame2.SetStats(False)
frame2.GetXaxis().SetTitle('p_{T}')
frame2.GetYaxis().SetTitle('Hz (8Tev, 1e34)')
frame2.SetTitle('')

uctIsoIdGraph2 = ROOT.TGraph(nrPts,uctRecoIsoIdPt,uctIsoIdRate)
uctIsoIdGraph2.SetLineColor(ROOT.EColor.kBlack)
uctIsoIdGraph2.SetMarkerColor(uIColor)
uctIsoIdGraph2.SetLineWidth(2)
uctIsoIdGraph2.SetMarkerStyle(uIMarker)
uctIsoIdGraph2.SetMarkerSize(uISize)
uctIsoIdGraph2.Draw('P')

uctIdGraph2 = ROOT.TGraph(nrPts,uctRecoIdPt,uctIdRate)
uctIdGraph2.SetLineColor(ROOT.EColor.kBlack)
uctIdGraph2.SetMarkerColor(uDColor)
uctIdGraph2.SetLineWidth(2)
uctIdGraph2.SetMarkerStyle(uDMarker)
uctIdGraph2.SetMarkerSize(uDSize)
uctIdGraph2.Draw('P')

curGraph2 = ROOT.TGraph(nrPts,curRecoPt,curRate)
curGraph2.SetLineColor(ROOT.EColor.kBlack)
curGraph2.SetMarkerColor(cColor)
curGraph2.SetLineWidth(2)
curGraph2.SetMarkerStyle(cMarker)
curGraph2.SetMarkerSize(cSize)
curGraph2.Draw('P')

legend2 = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
legend2.SetFillColor(ROOT.EColor.kWhite)
legend2.SetBorderSize(0)
legend2.AddEntry(uctIsoIdGraph2,'UCT: Region + ID + Iso < 0.2','P')
legend2.AddEntry(uctIdGraph2,'UCT: Region + ID (no Iso)','P')
legend2.AddEntry(curGraph2,'Current (no Iso)','P')

frame2.Draw()
curGraph2.Draw('LP')
uctIsoIdGraph2.Draw('LP')
uctIdGraph2.Draw('LP')
legend2.Draw()

tex.SetTextAlign(11)
tex.SetTextSize(0.07)
tex.DrawLatex(0.1,0.91,'EG Rate v p_{T}')
tex.SetTextAlign(31)
tex.SetTextSize(0.03)
tex.DrawLatex(0.9,0.91,'CMS Preliminary')

save2 = raw_input ('Press Enter for L1 pT v Rate (type save to save)\n')
if save2 == 'save':
 can2.SaveAs(saveWhere+'_L1_Rate.png')

### PLOT L1 Pt vs Rate at RecoPt corresponding to 95% plateau efficiency ###
can = ROOT.TCanvas('can','can',800,800)
can.SetLogy(True)

xmin = min(l1CutPt)-10
xmax = max(l1CutPt)+10
ymax = max(max(uctIsoIdRate),max(uctIdRate),max(curRate))*3
nrPts = len(l1CutPt)

frame = TH1F('frame','',1,xmin,xmax)
frame.SetMaximum(ymax)
frame.SetMinimum(100)
frame.SetStats(False)
frame.GetXaxis().SetTitle('p_{T}')
frame.GetYaxis().SetTitle('Hz (8Tev, 1e34)')
frame.SetTitle('')

uctIsoIdGraph = ROOT.TGraph(nrPts,l1CutPt,uctIsoIdRate)
uctIsoIdGraph.SetLineColor(ROOT.EColor.kBlack)
uctIsoIdGraph.SetMarkerColor(uIColor)
uctIsoIdGraph.SetLineWidth(2)
uctIsoIdGraph.SetMarkerStyle(uIMarker)
uctIsoIdGraph.SetMarkerSize(uISize)
uctIsoIdGraph.Draw('P')

uctIdGraph = ROOT.TGraph(nrPts,l1CutPt,uctIdRate)
uctIdGraph.SetLineColor(ROOT.EColor.kBlack)
uctIdGraph.SetMarkerColor(uDColor)
uctIdGraph.SetLineWidth(2)
uctIdGraph.SetMarkerStyle(uDMarker)
uctIdGraph.SetMarkerSize(uDSize)
uctIdGraph.Draw('P')

curGraph = ROOT.TGraph(nrPts,l1CutPt,curRate)
curGraph.SetLineColor(ROOT.EColor.kBlack)
curGraph.SetMarkerColor(cColor)
curGraph.SetLineWidth(2)
curGraph.SetMarkerStyle(cMarker)
curGraph.SetMarkerSize(cSize)
curGraph.Draw('P')

legend = ROOT.TLegend(0.6,0.7,0.89,0.89,'','brNDC')
legend.SetFillColor(ROOT.EColor.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(uctIsoIdGraph,'UCT: Region + ID + Iso < 0.2','P')
legend.AddEntry(uctIdGraph,'UCT: Region + ID (no Iso)','P')
legend.AddEntry(curGraph,'Current (no Iso)','P')

frame.Draw()
curGraph.Draw('LP')
uctIsoIdGraph.Draw('LP')
uctIdGraph.Draw('LP')
legend.Draw()

tex.SetTextAlign(11)
tex.SetTextSize(0.06)
tex.DrawLatex(0.1,0.91,'EG Rate at p_{T} for 95% Effi v L1 p_{T}')
tex.SetTextAlign(13)
tex.SetTextSize(0.03)
tex.DrawLatex(0.1,0.89,'CMS Preliminary')
tex.SetTextAlign(11)
#tex.DrawLatex(0.1,0.4,'Rate at p_{T} necessary to reach\n 95% of plateau efficiency corresponding to L1 p_{T} cut')
save = raw_input ('Press Enter to Exit (type save to save)\n')
if save == 'save':
 can.SaveAs(saveWhere+'_L1_pT95rate.png')
