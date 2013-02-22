'''
Usage:python jetEfficiencyPlot.py RootFile.root label[optional]

Script to make some quick efficiency plots to test ntuple generation.


Author: L. Dodd, UW Madison

'''

from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

######## File #########
if len(argv) < 2:
   print 'Usage:python jetEfficiencyPlot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

#if len(argv)>2:
   #saveWhere='~/www/'+argv[2]+'_'
#else:
   #saveWhere='~/www/'

saveWhere = './'


####### Calibration factor ####
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0


#####################################
#Get Jet Eff NTUPLE                 #
#####################################


jet_ntuple = ntuple_file.Get("jetEfficiency/Ntuple")

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_l1_efficiency(denom, num):
    ''' Make an efficiency graph with the "L1" style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerColor(ROOT.EColor.kRed)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.EColor.kBlack)
    return eff

def make_l1g_efficiency(denom, num):
    ''' Make an efficiency graph with the "UCT" style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(24)
    eff.SetMarkerColor(ROOT.EColor.kBlue)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.EColor.kBlack)
    return eff


def compare_efficiencies(ntuple, variable, l1PtCut, binning, filename,
                         title='', xaxis='', showL1=False):
    ''' Returns a (L1, L1G) tuple of TGraphAsymmErrors '''
    denom = make_plot(
        ntuple, variable,
        "", # No selection
        binning
    )
    num_l1g = make_plot(
        ntuple, variable,
        "l1gMatch && %0.2f * l1gPt > %0.2f" % (L1G_CALIB_FACTOR, l1PtCut),
        binning
    )
    num_l1 = make_plot(
        ntuple, variable,
        "l1Match && l1Pt > %0.2f" % (l1PtCut),
        binning
    )

    frame = ROOT.TH1F("frame", "frame", *binning)
    l1g = make_l1g_efficiency(denom, num_l1g)
    l1 = make_l1_efficiency(denom, num_l1)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.Draw()
    l1g.Draw('pe')
    if showL1:
        l1.Draw('pe')
    legend = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.EColor.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1g, "UCT", "pe")
    if showL1:
        legend.AddEntry(l1, "L1", "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Jet efficiency for a 20 GeV cut and 10 GeV cut on L1
################################################################################

compare_efficiencies(jet_ntuple, 'recoPt',  20, [40, 0, 200],
                     'jet_pt_trg20',
                     "Jet20 efficiency", "RECO p_{T} (GeV)")

#compare_efficiencies(jet_ntuple, 'recoPt',  10, [40, 0, 200],
                     #'jet_pt_trg10',
                     #"Jet10 efficiency", "RECO p_{T} (GeV)")

# Resolutions
jet_ntuple.Draw("(recoPt - (%0.2f * l1gPt))/recoPt>>l1gPtRes(100, -2, 2)" % L1G_CALIB_FACTOR, "l1gMatch && recoPt > 20 && l1gPt > 20 && abs(recoEta) < 3.0", "goff")
l1gPtRes = ROOT.gDirectory.Get("l1gPtRes")
l1gPtRes.SetLineColor(ROOT.EColor.kBlue)
l1gPtRes.Scale(1/l1gPtRes.Integral())
l1gPtRes.SetTitle('')
l1gPtRes.GetXaxis().SetTitle("(p_{T}^{RECO} - p_{T}^{TRG})/p_{T}^{RECO}")

#jet_ntuple.Draw("(recoPt - l1Pt)/recoPt>>l1PtRes(100, -2, 2)", "l1Match && recoPt > 20", "goff")
#l1PtRes = ROOT.gDirectory.Get("l1PtRes")
#l1PtRes.SetLineColor(ROOT.EColor.kRed)
#l1PtRes.Scale(1/l1PtRes.Integral())

#l1PtRes.Draw()
l1gPtRes.Draw()
legend1 = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
legend1.SetFillColor(ROOT.EColor.kWhite)
legend1.SetBorderSize(1)
legend1.AddEntry(l1gPtRes, "UCT")
#legend1.AddEntry(l1PtRes, "L1")
#legend1.Draw("same")
saveas=saveWhere+'jet_Pt_res.png'
canvas.SaveAs(saveas)

#>>>>>>>>Eta Res
jet_ntuple.Draw("(recoEta - l1gEta)>>l1gEtaRes(100, -2, 2)", "l1gMatch && recoPt > 20 && l1gPt > 20 && abs(recoEta) < 3.0 ", "goff")
l1gEtaRes = ROOT.gDirectory.Get("l1gEtaRes")
l1gEtaRes.SetLineColor(ROOT.EColor.kBlue)
l1gEtaRes.Scale(1/l1gEtaRes.Integral())
l1gEtaRes.SetTitle('')
l1gEtaRes.SetTitle('')
l1gEtaRes.GetXaxis().SetTitle("#eta^{RECO} - #eta^{TRG}")

#jet_ntuple.Draw("(recoEta - l1Eta)>>l1EtaRes(100, -2, 2)", "l1Match && recoPt > 20", "goff")
#l1EtaRes = ROOT.gDirectory.Get("l1EtaRes")
#l1EtaRes.SetLineColor(ROOT.EColor.kRed)
#l1EtaRes.Scale(1/l1EtaRes.Integral())

l1gEtaRes.Draw()
#l1EtaRes.Draw('same')

legend2 = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
legend2.SetFillColor(ROOT.EColor.kWhite)
legend2.SetBorderSize(1)
legend2.AddEntry(l1gEtaRes, "UCT")
#legend2.AddEntry(l1EtaRes, "L1")
#legend2.Draw("same")
saveas=saveWhere+'jet_Eta_res.png'
canvas.SaveAs(saveas)

#>>>>>>>>Phi Res
jet_ntuple.Draw("(recoPhi - l1gPhi)>>l1gPhiRes(100, -2, 2)", "l1gMatch && recoPt > 20 && l1gPt > 20 && abs(recoEta) < 3.0", "goff")
l1gPhiRes = ROOT.gDirectory.Get("l1gPhiRes")
l1gPhiRes.SetLineColor(ROOT.EColor.kBlue)
l1gPhiRes.Scale(1/l1gPhiRes.Integral())
l1gPhiRes.SetTitle('')
l1gPhiRes.GetXaxis().SetTitle("#phi^{RECO} - #phi^{TRG}")

#jet_ntuple.Draw("(recoPhi - l1Phi)>>l1PhiRes(100, -2, 2)", "l1Match && recoPt > 20", "goff")
#l1PhiRes = ROOT.gDirectory.Get("l1PhiRes")
#l1PhiRes.SetLineColor(ROOT.EColor.kRed)
#l1PhiRes.Scale(1/l1PhiRes.Integral())

l1gPhiRes.Draw()
#l1PhiRes.Draw('same')

legend3 = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
legend3.SetFillColor(ROOT.EColor.kWhite)
legend3.SetBorderSize(1)
legend3.AddEntry(l1gPhiRes, "UCT")
#legend3.AddEntry(l1PhiRes, "L1")
#legend3.Draw("same")
saveas = saveWhere+'jet_Phi_res.png'
canvas.SaveAs(saveas)




