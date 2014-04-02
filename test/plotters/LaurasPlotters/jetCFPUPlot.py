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
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit()

######## File #########
if len(argv) < 2:
   print 'Usage:python jetEfficiencyPlot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)

######## LABEL & SAVE WHERE #########

if len(argv)>2:
   saveWhere='~/www/'+argv[2]+'_'
else:
   saveWhere='~/www/'

#saveWhere = './'


####### Calibration factor ####
L1_CALIB_FACTOR = 1.0
L1G_CALIB_FACTOR = 1.0


#####################################
#Get Jet Eff NTUPLE                 #
#####################################


jet_ntuple = ntuple_file.Get("corrjetEfficiency/Ntuple")
jet_ntuple_uncorr = ntuple_file.Get("jetEfficiency/Ntuple")

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
        legend.AddEntry(l1, "Current", "pe")
    legend.Draw()
    saveas = saveWhere+filename+'.png'
    print saveas
    canvas.SaveAs(saveas)

# Pt difference vs pilup
canvas
#########################################
# 30<Reco PT
########################################


jet_ntuple.Draw("recoPt-l1gPt:nPVs>>l1g","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
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

saveas=saveWhere+'jet_PU_comp_uct30all.png'
canvas.SaveAs(saveas)

jet_ntuple.Draw("recoPt-l1Pt:nPVs>>l1","recoPt>30 && recoEta<3.0 && l1Match","BOX")
l1 = ROOT.gDirectory.Get("l1")
l1.SetLineColor(ROOT.EColor.kRed)
l1.SetTitle('Current Trigger Pt vs PU')
#l1g.GetXaxis().SetTitle("nPVs")
profilel1 = l1.ProfileX("_profilel1")
#printf('Current Trigger Pt vs PU:30,40')
#profilel1.Fit("pol1")
profilel1.SetLineWidth(2)
profilel1.SetMarkerStyle(23)
l1.GetXaxis().SetTitle("nPVs")
l1.GetYaxis().SetTitle("Reco Pt-Trigger Pt")
profilel1.Draw("same")
saveas=saveWhere+'jet_PU_comp_l130all.png'
canvas.SaveAs(saveas)

jet_ntuple_uncorr.Draw("recoPt-l1gPt:nPVs>>l1gu","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
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
saveas=saveWhere+'jet_PU_comp_uct_uncor30all.png'
canvas.SaveAs(saveas)

print 'Average (Reco Pt)-(Trigger Pt) Difference vs PU (30<recoPT<all)'
print 'Current'
profilel1.Fit("pol1")
print 'Uncorrected'
profilel1gu.Fit("pol1")
print 'Corrected'
profilel1g.Fit("pol1")

profilel1gu.Draw()
profilel1gu.SetTitle('Average (Reco Pt)-(Trigger Pt) Difference vs PU (30<recoPT)')
profilel1gu.GetYaxis().SetRangeUser(-10.0,35.0)
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
saveas=saveWhere+'profilePU30al.png'
canvas.SaveAs(saveas)

canvas
########################
# CF Correction FActor #
########################

jet_ntuple_uncorr.Draw("recoPt:nPVs>>CFr","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
jet_ntuple_uncorr.Draw("l1gPt:nPVs>>CFu","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
CFr = ROOT.gDirectory.Get("CFr")
CFu = ROOT.gDirectory.Get("CFu")
CFr.SetLineColor(ROOT.EColor.kRed)
CFu.SetLineColor(ROOT.EColor.kMagenta)
profileCFr = CFr.ProfileX("_profileCFr")
profileCFu = CFu.ProfileX("_profileCFu")
profileCFr.SetMarkerStyle(23)
profileCFu.SetMarkerStyle(23)
profileCFu.GetXaxis().SetTitle("nPVs")
profileCFu.GetYaxis().SetTitle("Uncorrected PT^{Trigger}- PT^{Reco}")
profileCFu.SetLineWidth(2)
profileCFu.Add(profileCFr,-1)
profileCFu.Draw()
profileCFu.Fit("pol1")
saveas=saveWhere+'profilePUUR.png'
canvas.SaveAs(saveas)

########################
# CF FActor #
########################
canvas


jet_ntuple.Draw("l1gPt:nPVs>>CF","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
jet_ntuple_uncorr.Draw("l1gPt:nPVs>>CFu","recoPt>30 && recoEta<3.0 && l1gMatch","BOX")
CF = ROOT.gDirectory.Get("CF")
CFu = ROOT.gDirectory.Get("CFu")
CF.SetLineColor(ROOT.EColor.kBlue)
CFu.SetLineColor(ROOT.EColor.kMagenta)
profileCF = CF.ProfileX("_profileCF")
profileCFu = CFu.ProfileX("_profileCFu")
profileCF.SetMarkerStyle(23)
profileCFu.SetMarkerStyle(23)
profileCFu.GetXaxis().SetTitle("nPVs")
profileCFu.GetYaxis().SetTitle("Uncorrected PT^{Trigger}- Corrected PT^{Trigger}")
profileCFu.SetLineWidth(2)
profileCFu.Add(profileCF,-1)
#profileCF.SetOptStats(111111)
profileCFu.Draw()
profileCFu.Fit("pol1")
saveas=saveWhere+'profilePUCF.png'
canvas.SaveAs(saveas)




