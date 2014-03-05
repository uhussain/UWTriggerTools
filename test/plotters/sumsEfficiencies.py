'''

Script to make plots of MET/MHT efficiencies

'''
from subprocess import Popen
from sys import argv, exit, stdout, stderr

import ROOT

# So things don't look like crap.
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)

def make_plot(tree, variable, selection, binning, xaxis='', title=''):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s>>htemp(%s)" % (variable, ", ".join(str(x) for x in binning))
    print draw_string
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_l1_efficiency(denom, num):
    ''' Make an efficiency graph with the "L1" style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(25)
    eff.SetMarkerColor(ROOT.EColor.kRed)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.EColor.kBlack)
    return eff

def make_l1g_efficiency(denom, num):
    ''' Make an efficiency graph with the "UCT" style '''
    eff = ROOT.TGraphAsymmErrors(num, denom)
    eff.SetMarkerStyle(20)
    eff.SetMarkerColor(ROOT.EColor.kBlue)
    eff.SetMarkerSize(1.5)
    eff.SetLineColor(ROOT.EColor.kBlack)
    return eff

def compare_efficiencies(oct_ntuple, uct_ntuple, variable='recoMET',
                         title='', xaxis=' (GeV)',
                         uct_cut=None, oct_cut=None,
                         legend_label='',
                         uct_legend = 'UCT',
                         l1_legend = 'Current',
                         filename=None,
                         binning = (40, 0, 400)
                        ):
    ''' Returns a (L1, L1G) tuple of TGraphAsymmErrors '''
    denom = make_plot(
        oct_ntuple, variable,
        "", # No selection
        binning
    )
    num_l1g = make_plot(
        uct_ntuple, variable,
        uct_cut,
        binning
    )
    num_l1 = make_plot(
        oct_ntuple, variable,
        oct_cut,
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
    l1g.SetMaximum(1.1)
    useL1 = True
    if useL1:
        l1.Draw('pe')
    legend = ROOT.TLegend(0.65, 0.2, 0.89, 0.35, "", "brNDC")
    if legend_label:
        legend.AddEntry("NULL", legend_label,"")
    legend.SetFillColor(ROOT.EColor.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1g, uct_legend, "pe")
    if useL1:
        legend.AddEntry(l1, l1_legend, "pe")
    legend.Draw()
    canvas.SetGridy()
    canvas.SetGridx()
    canvas.SaveAs(filename+".png")

######## File #########
if len(argv) < 2:
   print 'Usage:python jetEfficiencyPlot.py RootFile.root label[optional]'
   exit()

infile = argv[1]
ntuple_file = ROOT.TFile(infile)


#ntuple_file = ROOT.TFile("/scratch/efriis/uct_jet_efficiency.root")
uct_sums_ntuple = ntuple_file.Get("uctSumsEfficiency/Ntuple")
l1_sums_ntuple = ntuple_file.Get("l1SumsEfficiency/Ntuple")
#################
compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMET',
                     uct_cut = 'l1MET > 20',
                     oct_cut = 'l1MET > 20',
                     filename = 'l1_met_20'
                    )
compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMET',
                     uct_cut = 'l1MET > 30',
                     oct_cut = 'l1MET > 30',
                     filename = 'l1_met_30'
                    )


compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMET',
                     uct_cut = 'l1MET > 50',
                     oct_cut = 'l1MET > 50',
                     filename = 'l1_met_uct_50_l1_50'
                    )

compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMHT',
                     uct_cut = 'l1MHT > 20',
                     oct_cut = 'l1MHT > 20',
                     filename = 'l1_mht_20'
                    )

compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMHT',
                     uct_cut = 'l1SHT > 30',
                     oct_cut = 'l1SHT > 30',
                     filename = 'l1_sht_30'
                    )


compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMHT',
                     uct_cut = 'l1SHT > 175',
                     oct_cut = 'l1SHT > 175',
                     filename = 'l1_sht_175'
                    )


compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMHT',
                     uct_cut = 'l1SHT > 50',
                     oct_cut = 'l1SHT > 50',
                     filename = 'l1_sht_50'
                    )


compare_efficiencies(l1_sums_ntuple, uct_sums_ntuple,
                     variable='recoMHT',
                     uct_cut = 'l1MHT > 40',
                     oct_cut = 'l1MHT > 40',
                     filename = 'l1_mht_40'
                    )

# Resolutions
uct_sums_ntuple.Draw("(recoMET - l1MET)/recoMET>>l1gRes(100, -2, 2)", "recoMET > 30", "goff")
l1gRes = ROOT.gDirectory.Get("l1gRes")
l1gRes.SetLineColor(ROOT.EColor.kBlue)
l1gRes.Scale(1/l1gRes.Integral())

l1gRes.SetTitle('')
l1gRes.GetXaxis().SetTitle("(ME_{T}^{RECO} - ME_{T}^{TRG})/ME_{T}^{RECO}")
l1gRes.Draw()
canvas.SaveAs("reco_met30_uct_res.png")

l1_sums_ntuple.Draw("(recoMET - l1MET)/recoMET>>l1Res(100, -2, 2)", "recoMET > 30", "goff")
l1Res = ROOT.gDirectory.Get("l1Res")
l1Res.SetLineColor(ROOT.EColor.kRed)
l1Res.Scale(1/l1Res.Integral())

l1Res.Draw('same')
canvas.SaveAs("reco_met30_res_cmp.png")
