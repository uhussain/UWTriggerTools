#!/usr/bin/env python
import ROOT
import random
from ROOT import *
from array import array
#eff_file = ROOT.TFile("/scratch/efriis/uct_tau_efficiency.root")
#rate_file = ROOT.TFile("/scratch/efriis/uct_rates_eic3.root")
eff_file = ROOT.TFile("../data/Evan/uct_tau_efficiency.root")

tau_eff_ntuple = eff_file.Get("rlxTauEfficiency/Ntuple")
tau_eff_ntuple_cur = eff_file.Get("isoTauEfficiency/Ntuple")
saveWhere='../plots/Tau_'

canvas = ROOT.TCanvas("asdf", "adsf", 600, 600)
rangeMax=130
def get_pt95(formula):
    for i in range(rangeMax):
        # is this actually giving pT90 instead of 95?
        if formula.Eval(i) / formula.GetParameter(0) > 0.90:
            return i
    return rangeMax

def make_turnon(ntuple,l1_cut, color, name, fix_plateau=None,save=False):
    ntuple.Draw("recoPt>>denom(20, 0, %i)"%rangeMax, "", "goff")
    ntuple.Draw("recoPt>>num(20, 0, %i)"%rangeMax, l1_cut, "goff")
    denom = ROOT.gDirectory.Get("denom")
    num = ROOT.gDirectory.Get("num")
    graph = ROOT.TGraphAsymmErrors(num, denom)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(20)
    #why does evan name the fitter with a random number? nutty
    formula = ROOT.TF1('fitter' + str(random.randint(0, 1000)),
                       "[0]*(0.5*(tanh([1]*(x-[2]))+1))", 0, rangeMax)
    formula.SetParName(0, 'plateau')
    formula.SetParName(1, 'width')
    formula.SetParName(2, 'threshold')
    formula.SetParameter(0, 0.95)
    formula.SetParLimits(0, 0.1, 1.01)
    if fix_plateau is not None:
        print 'fix plat'
        formula.FixParameter(0, fix_plateau)
    formula.SetParameter(1, 1. / 10)
    formula.SetParameter(2, 30)
    formula.SetLineColor(color)
    graph.Fit(formula)
    # refit, not using the inefficient junk at start of curve
    formula.SetRange(formula.GetParameter(2), rangeMax)
    graph.Fit(formula, "R")
    graph.Draw("ape")
    graph.GetHistogram().SetMaximum(1.1)
    graph.GetHistogram().SetMinimum(0)
    pt95 = get_pt95(formula)
    graph.GetHistogram().GetXaxis().SetTitle('pt95 = %i' % pt95)
    if save==True:
     canvas.SaveAs(saveWhere+name + ".png")
    return graph, formula, pt95

ISOTHRESHOLD=0.20

tex = ROOT.TLatex()
tex.SetNDC(True)
tex.SetTextAlign(11)
tex.SetTextSize(0.03)

isoCut='((l1gJetPt-max(l1gRegionEt,l1gPt))/max(l1gRegionEt,l1gPt) <%0.1f)'%(ISOTHRESHOLD)

l1ptVal=array('d',[20,25,30,35,40,45])
colors=[ROOT.EColor.kRed,
ROOT.EColor.kOrange,
ROOT.EColor.kYellow,
ROOT.EColor.kGreen,
ROOT.EColor.kBlue,
ROOT.EColor.kViolet]
iso95=array('d',[])
non95=array('d',[])
cur95=array('d',[])
for pt,color in zip(l1ptVal,colors):
 l1gPtCut = '(max(l1gPt,l1gRegionEt) > %i )' %pt

 cutI=isoCut +'&&'+l1gPtCut+'&&l1gMatch'
 cutN=l1gPtCut+'&&l1gMatch'
   
 resultIso = make_turnon(tau_eff_ntuple,cutI,color,'pt%i_Iso'%pt,save=True)
 resultNon = make_turnon(tau_eff_ntuple,cutN,color,'pt%i_Non'%pt,save=True)
 resultCur = make_turnon(tau_eff_ntuple_cur,cutN,color,'pt%i_Cur'%pt,save=True)
 iso95.append(resultIso[2])
 non95.append(resultNon[2])
 cur95.append(resultCur[2])
#print(l1ptVal)
#print(iso95)
#print(non95)
#print(cur95)

nrPts=len(l1ptVal)
xmin=min(l1ptVal)-5
xmax=max(l1ptVal)+10
ymin=min(min(iso95),min(non95),min(cur95))-5
ymax=max(max(iso95),max(non95),max(cur95))+5
uIColor=ROOT.EColor.kBlue
uNColor=ROOT.EColor.kGreen+3
cColor=ROOT.EColor.kRed
uIMarker=21
uNMarker=20
cMarker=22
uISize=1.5
uNSize=1.5
cSize=1.5

# but who am i to judge naming conventions
lineFit = ROOT.TF1('lineFit'+str(random.randint(0, 1000)),
   "[0]*x+[1]", 0, rangeMax)
lineFit.SetParName(0,'m')
lineFit.SetParName(1,'b')

can = ROOT.TCanvas('can','can',800,800)
can.SetLogy(False)
frame = TH1F('frame','',1,xmin,xmax)
frame.SetMinimum(ymin)
frame.SetMaximum(ymax)
frame.SetStats(False)
frame.GetXaxis().SetTitle('L1 Cut')
frame.GetYaxis().SetTitle('(RECO) pT95')
frame.SetTitle('')
frame.Draw()

uIgraph=ROOT.TGraph(nrPts,l1ptVal,iso95)
uIgraph.SetMarkerColor(uIColor)
uIgraph.SetMarkerStyle(uIMarker)
uIgraph.SetMarkerSize(uISize)
lineFit.SetLineColor(uIColor)
uIgraph.Fit(lineFit)
uIgraph.Draw('P')
tex.DrawLatex(0.5,0.3,'UCT Iso: pT95 = %0.01f L1 + %0.1f'
  %(lineFit.GetParameter(0),lineFit.GetParameter(1)))

uNgraph=ROOT.TGraph(nrPts,l1ptVal,non95)
uNgraph.SetMarkerColor(uNColor)
uNgraph.SetMarkerStyle(uNMarker)
uNgraph.SetMarkerSize(uNSize)
lineFit.SetLineColor(uNColor)
uNgraph.Fit(lineFit)
uNgraph.Draw('P')
tex.DrawLatex(0.5,0.25,'UCT NonIso: pT95 = %0.01f L1 + %0.1f'
  %(lineFit.GetParameter(0),lineFit.GetParameter(1)))
uNgraph.Draw('P')

cgraph=ROOT.TGraph(nrPts,l1ptVal,cur95)
cgraph.SetMarkerColor(cColor)
cgraph.SetMarkerStyle(cMarker)
cgraph.SetMarkerSize(cSize)
lineFit.SetLineColor(cColor)
cgraph.Fit(lineFit)
cgraph.Draw('P')
tex.DrawLatex(0.5,0.2,'Current: pT95 = %0.01f L1 + %0.1f'
  %(lineFit.GetParameter(0),lineFit.GetParameter(1)))
cgraph.Draw('P')

legend = ROOT.TLegend(0.11,0.5,0.4,0.89,'','brNDC')
legend.SetFillColor(ROOT.EColor.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(uIgraph,'UCT Iso','P')
legend.AddEntry(uNgraph,'UCT NonIso','P')
legend.AddEntry(cgraph,'Current','P')
legend.Draw()

save=raw_input('press enter to finish, type save to save\n')
if save=='save':
 can.SaveAs(saveWhere+'fits.png')
