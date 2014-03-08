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

####
if len(argv) < 3:
   print 'Usage:python jetEfficiencyPlot.py EfficiencyRootFile.root EfficiencyRootFile.root label[optional]'
   exit()
#####
eff_infile = argv[1]
rate_infile = argv[2]
ntuple_file = ROOT.TFile(eff_infile)
rate_ntuple_file = ROOT.TFile(rate_infile)


#ntuple_file = ROOT.TFile("/scratch/efriis/uct_jet_efficiency.root")
uct_sums_ntuple = ntuple_file.Get("uctSumsEfficiency/Ntuple")
l1_sums_ntuple = ntuple_file.Get("l1SumsEfficiency/Ntuple")
#################

sums_l1_ntuple = rate_ntuple_file.Get("sumsL1Rates/Ntuple")
sums_uct_ntuple = rate_ntuple_file.Get("sumsUCTRates/Ntuple")
#################


canvas = ROOT.TCanvas("asdf", "adsf", 800, 800)
ZEROBIAS_RATE=15000000.00
#ZEROBIAS_RATE=1.00


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
#End EFF


######## File #########
def make_plot(tree, variable, selection, binning, xaxis='', title='', calFactor=1):
    ''' Plot a variable using draw and return the histogram '''
    #draw_string = "MaxIf$(%s * %0.2f, %s)>>htemp(%s)" % (variable, calFactor, selection, ", ".join(str(x) for x in binning))
    draw_string = "%s * %0.2f>>htemp(%s)" % (variable, calFactor, ", ".join(str(x) for x in binning))
    print draw_string
    print tree
    #tree.Draw(draw_string, "", "goff")
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

def make_l1_rate(pt):
    ''' Make a rate plot out of L1Extra Pts '''
    numBins = pt.GetNbinsX()
    rate = pt.Clone()

    for i in range(numBins):
        rate.SetBinContent(i+1, pt.Integral(i+1, numBins))

    rate.SetLineColor(ROOT.EColor.kRed)
    rate.SetMarkerStyle(25)
    rate.SetMarkerColor(ROOT.EColor.kRed)

    return rate

def make_l1g_rate(pt):
    ''' Return a rate plot out of L1Extra Pts '''
    numBins = pt.GetNbinsX()
    rate = pt.Clone()

    for i in range(1,numBins+1):
        rate.SetBinContent(i, pt.Integral(i, numBins))

    rate.SetLineColor(ROOT.EColor.kBlue)
    rate.SetMarkerStyle(22)
    rate.SetMarkerColor(ROOT.EColor.kBlue)

    return rate

def plotRates(l1Ntuple, l1gNtuple, binning,
              l1Variable='', l1gVariable='',
              filename='', title='', xaxis='',
              l1Cut = '', l1gCut = '',
              showL1=True, showMC=False, l1gLabel='UCT', l1gColor=ROOT.EColor.kBlack,
              l1gStyle = 22
             ):
    ''' Save a rate Plot '''
    scale = ZEROBIAS_RATE/l1Ntuple.GetEntries()
    scaleg = ZEROBIAS_RATE/l1gNtuple.GetEntries()
	
    l1_pt = make_plot(
        l1Ntuple, l1Variable,
        l1Cut,
        binning,
        '','', # Don't bother with titles
        )
    l1g_pt = make_plot(
        l1gNtuple, l1gVariable,
        l1gCut,
        binning,
        '','', # Don't bother with titles
        )
    l1Rate = make_l1_rate(l1_pt)
    l1gRate = make_l1g_rate(l1g_pt)
    l1gRate.Sumw2()
    l1Rate.Sumw2()
#    l1Rate.Scale(scale)
#    l1gRate.Scale(scaleg)
    l1gRate.SetLineColor(l1gColor)
    l1gRate.SetMarkerStyle(l1gStyle)
    l1gRate.SetMarkerColor(l1gColor)

    l1gRate.SetMaximum(l1Rate.GetMaximum()*20)
    l1gRate.SetMinimum(
#            l1Rate.GetBinContent(l1Rate.FindBin(200))
	10
    )

    canvas.SetLogy()
    l1gRate.SetTitle(title)
    l1gRate.GetXaxis().SetTitle(xaxis)
    l1gRate.GetYaxis().SetTitle("")
    l1gRate.GetYaxis().SetTitleOffset(1.29)
    l1gRate.Draw('ph')
    if showL1:
        l1Rate.Draw('phsame')
    if showMC:
  	l1gMCRate.SetLineColor(ROOT.EColor.kRed)
	l1gMCRate.Draw('phsame')
    legend = ROOT.TLegend(0.3, 0.89, 0.89, 0.7, "", "brNDC")
    legend.SetFillColor(ROOT.EColor.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1gRate, l1gLabel, "p")
    if showL1:
        legend.AddEntry(l1Rate, "Current", "p")
    legend.Draw()
    print 'rate at 72:'
    binn = l1gRate.GetXaxis().FindBin(72)
    rateVal = l1gRate.GetBinContent(binn)
    print rateVal
    print 'rate at 84:'
    binn = l1gRate.GetXaxis().FindBin(84)
    rateVal = l1gRate.GetBinContent(binn)
    print rateVal
    canvas.SaveAs(filename+".png")

###############################################################################
# Draw Rates                                                                  #
###############################################################################

plotRates(sums_l1_ntuple, sums_uct_ntuple,
          [39, 5, 200],
          'l1MET',
          'l1MET',
          'met_rate_cmp_nocut',
          "MET Rate NoCut", " ME_{T} Threshold (GeV)",
          l1gLabel="",
          l1gColor = ROOT.EColor.kBlue,
          l1gStyle = 20,
          l1gCut = '',
          l1Cut = '',
          showL1 = True,
         )



################

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



# Resolutions
uct_sums_ntuple.Draw("(recoMET - l1MET)/recoMET>>l1gRes(100, -2, 2)", "recoMET > 30", "goff")
l1gRes = ROOT.gDirectory.Get("l1gRes")
l1gRes.SetLineColor(ROOT.EColor.kBlue)
#l1gRes.Scale(1/l1gRes.Integral())

l1gRes.SetTitle('')
l1gRes.GetXaxis().SetTitle("(ME_{T}^{RECO} - ME_{T}^{TRG})/ME_{T}^{RECO}")
l1gRes.Draw()
canvas.SaveAs("reco_met30_uct_res.png")

l1_sums_ntuple.Draw("(recoMET - l1MET)/recoMET>>l1Res(100, -2, 2)", "recoMET > 30", "goff")
l1Res = ROOT.gDirectory.Get("l1Res")
l1Res.SetLineColor(ROOT.EColor.kRed)
#l1Res.Scale(1/l1Res.Integral())

l1Res.Draw('same')
canvas.SaveAs("reco_met30_res_cmp.png")
