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
if len(argv) < 3:
   print 'Usage:python jetEfficiencyPlot.py RootFile.root RootFile_new.root label[optional]'
   exit()

infile = argv[1]
infile_new = argv[2]
ntuple_file = ROOT.TFile(infile)
ntuple_file_new = ROOT.TFile(infile_new)

######## LABEL & SAVE WHERE #########

if len(argv)>3:
   saveWhere='~/www/'+argv[3]+'_'
else:
   saveWhere='~/www/Research/CorrvsNEw/'



####### Calibration factor ####
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0


#####################################
#Get Jet Eff NTUPLE                 #
#####################################


#jet_ntuple = ntuple_file.Get("corrjetEfficiency/Ntuple")
jet_ntuple_new = ntuple_file_new.Get("jetEfficiency/Ntuple")
jet_ntuple = ntuple_file.Get("corrjetEfficiency/Ntuple")

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


def compare_efficiencies(ntuple, ntuple_new, variable, l1PtCut, binning, filename,
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
#    num_l1 = make_plot(
#        ntuple, variable,
#        "l1Match && l1Pt > %0.2f" % (l1PtCut),
#        binning
#    )
    num_l1g_new = make_plot(
        ntuple_new, variable,
        "l1gMatch && l1gPt > %0.2f" % (l1PtCut),
        binning
    )
    frame = ROOT.TH1F("frame", "frame", *binning)
    l1g = make_l1g_efficiency(denom, num_l1g)
    #    l1 = make_l1_efficiency(denom, num_l1)
    l1gnew = make_l1_efficiency(denom, num_l1g_new)
    frame.SetMaximum(1.2)
    frame.SetTitle(title)
    frame.GetXaxis().SetTitle(xaxis)
    frame.Draw()
    l1g.Draw('pe')
    l1gnew.Draw('pe')
    #if showL1:
    #    l1.Draw('pe')
    legend = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
    legend.SetFillColor(ROOT.EColor.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1g, "UCT", "pe")
    legend.AddEntry(l1gnew, "UCT pile up subtracted", "pe")
    #if showL1:
    #    legend.AddEntry(l1, "Current", "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

################################################################################
# Jet efficiency for a 20 GeV cut and 10 GeV cut on L1
################################################################################

compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  100, [40, 0, 300],
                    'jet_pt_trg100',
                     "Jet 100 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  150, [40, 0, 300],
                     'jet_pt_trg150',
                     "Jet 150 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  200, [40, 0, 300],
                     'jet_pt_trg200',
                     "Jet 200 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  30, [40, 0, 200],
                     'jet_pt_trg30',
                     "Jet 30 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  40, [40, 0, 200],
                     'jet_pt_trg40',
                     "Jet 40 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  50, [40, 0, 200],
                     'jet_pt_trg50',
                     "Jet 50 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple, jet_ntuple_new,'recoPt',  60, [40, 0, 200],
                     'jet_pt_trg60',
                     "Jet 60 efficiency", "RECO p_{T} (GeV)",True)
compare_efficiencies(jet_ntuple,jet_ntuple_new, 'recoPt',  70, [40, 0, 200],
                    'jet_pt_trg70',
                   "Jet 70 efficiency", "RECO p_{T} (GeV)",True)
# Resolutions
jet_ntuple.Draw("(recoPt - (%0.2f * l1gPt))/recoPt>>l1gPtRes(150, -2, 1)" % L1G_CALIB_FACTOR, "l1gMatch && recoPt > 30", "goff")
l1gPtRes = ROOT.gDirectory.Get("l1gPtRes")
l1gPtRes.SetLineColor(ROOT.EColor.kBlue)
l1gPtRes.Scale(1/l1gPtRes.Integral())
l1gPtRes.SetTitle('')
l1gPtRes.GetXaxis().SetTitle("(p_{T}^{RECO} - p_{T}^{TRG})/p_{T}^{RECO}")

jet_ntuple_new.Draw("(recoPt - l1gPt)/recoPt>>l1gPtResNew(150, -2, 1)", "l1Match && recoPt > 30", "goff")
l1PtRes = ROOT.gDirectory.Get("l1gPtResNew")
l1PtRes.SetLineColor(ROOT.EColor.kRed)
l1PtRes.Scale(1/l1PtRes.Integral())
l1gPtRes.Draw()
l1PtRes.Draw("same")
legend1 = ROOT.TLegend(0.2, 0.2, 0.39, 0.4, "", "brNDC")
legend1.SetFillColor(ROOT.EColor.kWhite)
legend1.SetBorderSize(1)
legend1.AddEntry(l1gPtRes, "UCT")
legend1.AddEntry(l1PtRes, "New UCT Calibration")
legend1.Draw("same")
saveas=saveWhere+'jet_Pt_res_corr_cut_30_new.png'
canvas.SaveAs(saveas)

#
# Resolutions
jet_ntuple.Draw("(recoPt - (%0.2f * l1gPt))/recoPt>>l1gPtRes(150, -2, 1)" % L1G_CALIB_FACTOR, "l1gMatch && recoPt > 30 && abs(recoEta) < 3.0", "goff")
l1gPtRes = ROOT.gDirectory.Get("l1gPtRes")
l1gPtRes.SetLineColor(ROOT.EColor.kBlue)
l1gPtRes.Scale(1/l1gPtRes.Integral())
l1gPtRes.SetTitle('')
l1gPtRes.GetXaxis().SetTitle("(p_{T}^{RECO} - p_{T}^{TRG})/p_{T}^{RECO}")

jet_ntuple.Draw("(recoPt - l1Pt)/recoPt>>l1PtRes(150, -2, 1)", "l1Match && recoPt > 30 && abs(recoEta)<3.0", "goff")
l1PtRes = ROOT.gDirectory.Get("l1PtRes")
l1PtRes.SetLineColor(ROOT.EColor.kRed)
l1PtRes.Scale(1/l1PtRes.Integral())
l1gPtRes.Draw()
l1PtRes.Draw("same")
legend1 = ROOT.TLegend(0.2, 0.2, 0.39, 0.4, "", "brNDC")
legend1.SetFillColor(ROOT.EColor.kWhite)
legend1.SetBorderSize(1)
legend1.AddEntry(l1gPtRes, "UCT")
legend1.AddEntry(l1PtRes, "Current")
legend1.Draw("same")
saveas=saveWhere+'jet_Pt_res_corr_cut_30.png'
canvas.SaveAs(saveas)

# Pt difference vs pilup
########################################
# 30<Reco PT<40
########################################


jet_ntuple.Draw("recoPt-l1gPt:nPVs>>l1g","recoPt>30 &&recoPt<40 && recoEta<3.0 && l1gMatch","BOX")
l1g = ROOT.gDirectory.Get("l1g")
l1g.SetLineColor(ROOT.EColor.kBlue)
l1g.SetTitle('Upgrade Corr Trigger Pt vs PU')
#l1g.GetXaxis().SetTitle("nPVs")
profilel1g = l1g.ProfileX("_profilel1g")
profilel1g.SetMarkerStyle(23)
l1g.GetXaxis().SetTitle("nPVs")
l1g.GetYaxis().SetTitle("Reco Pt-Trigger Pt")
profilel1g.SetLineWidth(2)
profilel1g.Draw("same")

saveas=saveWhere+'jet_PU_comp_uct.png'
canvas.SaveAs(saveas)

jet_ntuple.Draw("recoPt-l1Pt:nPVs>>l1","recoPt>30 &&recoPt<40 && recoEta<3.0 && l1Match","BOX")
l1 = ROOT.gDirectory.Get("l1")
l1.SetLineColor(ROOT.EColor.kRed)
l1.SetTitle('Current Trigger Pt vs PU')
#l1g.GetXaxis().SetTitle("nPVs")
profilel1 = l1.ProfileX("_profilel1")
profilel1.SetLineWidth(2)
profilel1.SetMarkerStyle(23)
l1.GetXaxis().SetTitle("nPVs")
l1.GetYaxis().SetTitle("Reco Pt-Trigger Pt")
profilel1.Draw("same")
saveas=saveWhere+'jet_PU_comp_l1.png'
canvas.SaveAs(saveas)

jet_ntuple_uncorr.Draw("recoPt-l1gPt:nPVs>>l1gu","recoPt>30 &&recoPt<40 && recoEta<3.0 && l1gMatch","BOX")
l1gu = ROOT.gDirectory.Get("l1gu")
l1gu.SetLineColor(ROOT.EColor.kMagenta)
l1gu.SetTitle('Uncorr Upgrade Trigger Pt vs PU')
#l1g.GetXaxis().SetTitle("nPVs")
profilel1gu = l1gu.ProfileX("_profilel1gu")
profilel1gu.SetLineWidth(2)
profilel1gu.SetMarkerStyle(23)
l1gu.GetXaxis().SetTitle("nPVs")
l1gu.GetYaxis().SetTitle("Reco Pt-Trigger Pt")
profilel1gu.Draw("same")
saveas=saveWhere+'jet_PU_comp_uct_uncor.png'
canvas.SaveAs(saveas)

profilel1gu.Draw()
profilel1gu.SetTitle('Average (Reco Pt)-(Trigger Pt) Difference vs PU (30<recoPT<40)')
profilel1gu.GetYaxis().SetRangeUser(-10.0,25.0)
profilel1gu.GetXaxis().SetTitle('nPVs')
profilel1gu.GetYaxis().SetTitle('Average RecoPt-TriggerPt')
profilel1.Draw("same")
profilel1g.Draw("same")
legend1 = ROOT.TLegend(0.7, 0.82, 0.99, 0.99, "", "brNDC")
legend1.SetFillColor(ROOT.EColor.kWhite)
legend1.SetBorderSize(1)
legend1.AddEntry(profilel1, "Current Jet")
legend1.AddEntry(profilel1g, "Corrected UCT Jet")
legend1.AddEntry(profilel1gu, "Unorrected UCT Jet")
legend1.Draw("same")
saveas=saveWhere+'profilePU.png'
canvas.SaveAs(saveas)


jet_ntuple.Draw("(recoPt - l1Pt)/recoPt>>l1PtRes(150, -2, 1)", "l1Match && recoPt > 30 && abs(recoEta)<3.0", "goff")
l1PtRes = ROOT.gDirectory.Get("l1PtRes")
l1PtRes.SetLineColor(ROOT.EColor.kRed)
l1PtRes.Scale(1/l1PtRes.Integral())
l1gPtRes.Draw()
l1PtRes.Draw("same")
legend1 = ROOT.TLegend(0.2, 0.2, 0.39, 0.4, "", "brNDC")
legend1.SetFillColor(ROOT.EColor.kWhite)
legend1.SetBorderSize(1)
legend1.AddEntry(l1gPtRes, "UCT")
legend1.AddEntry(l1PtRes, "Current")
legend1.Draw("same")
saveas=saveWhere+'jet_Pt_res_corr_cut_30.png'
canvas.SaveAs(saveas)



#>>>>>>>>Eta Res
jet_ntuple.Draw("(recoEta - l1gEta)>>l1gEtaRes(100, -2, 2)", "l1gMatch && recoPt > 20 && l1gPt > 20 && abs(recoEta) < 3.0 ", "goff")
l1gEtaRes = ROOT.gDirectory.Get("l1gEtaRes")
l1gEtaRes.SetLineColor(ROOT.EColor.kBlue)
l1gEtaRes.Scale(1/l1gEtaRes.Integral())
l1gEtaRes.SetTitle('')
l1gEtaRes.SetTitle('')
l1gEtaRes.GetXaxis().SetTitle("#eta^{RECO} - #eta^{TRG}")
jet_ntuple.Draw("(recoEta - l1Eta)>>l1EtaRes(100, -2, 2)", "l1Match && recoPt > 20", "goff")
l1EtaRes = ROOT.gDirectory.Get("l1EtaRes")
l1EtaRes.SetLineColor(ROOT.EColor.kRed)
l1EtaRes.Scale(1/l1EtaRes.Integral())
l1gEtaRes.Draw()
l1EtaRes.Draw('same')
legend2 = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
legend2.SetFillColor(ROOT.EColor.kWhite)
legend2.SetBorderSize(1)
legend2.AddEntry(l1gEtaRes, "UCT")
legend2.AddEntry(l1EtaRes, "L1")
legend2.Draw("same")
saveas=saveWhere+'jet_Eta_res.png'
canvas.SaveAs(saveas)

#>>>>>>>>Phi Res
jet_ntuple.Draw("(recoPhi - l1gPhi)>>l1gPhiRes(100, -2, 2)", "l1gMatch && recoPt > 20 && l1gPt > 20 && abs(recoEta) < 3.0", "goff")
l1gPhiRes = ROOT.gDirectory.Get("l1gPhiRes")
l1gPhiRes.SetLineColor(ROOT.EColor.kBlue)
l1gPhiRes.Scale(1/l1gPhiRes.Integral())
l1gPhiRes.SetTitle('')
l1gPhiRes.GetXaxis().SetTitle("#phi^{RECO} - #phi^{TRG}")

jet_ntuple.Draw("(recoPhi - l1Phi)>>l1PhiRes(100, -2, 2)", "l1Match && recoPt > 20", "goff")
l1PhiRes = ROOT.gDirectory.Get("l1PhiRes")
l1PhiRes.SetLineColor(ROOT.EColor.kRed)
l1PhiRes.Scale(1/l1PhiRes.Integral())

l1gPhiRes.Draw()
l1PhiRes.Draw('same')

legend3 = ROOT.TLegend(0.7, 0.2, 0.89, 0.4, "", "brNDC")
legend3.SetFillColor(ROOT.EColor.kWhite)
legend3.SetBorderSize(1)
legend3.AddEntry(l1gPhiRes, "UCT")
legend3.AddEntry(l1PhiRes, "L1")
legend3.Draw("same")
saveas = saveWhere+'jet_Phi_res.png'
canvas.SaveAs(saveas)




