'''

Script to make some rate plots

Author: Cepeda, Dasu, Dodd, Friis, Woods, UW Madison

# Please provide an input file (argv[1]) with the rate root file
# names. 
# To Run: python plotRates.py rates.root L1_CAL_Factor UCT_CAL_Factor 
if no calibration is specified default cal factor is 1

'''

########
#get everything
########
from sys import argv, stdout, stderr

import ROOT


######## STYLE ########
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


######## File #########
infile = argv[1]
ntuple_file = ROOT.TFile(infile)

####### Calibration factor ####
if(len(argv) < 3):
  L1_CALIB_FACTOR = 1.0
else:
  L1_CALIB_FACTOR = float(argv[2])

if(len(argv) < 4):
  L1G_CALIB_FACTOR = 1.0
else:
  L1G_CALIB_FACTOR = float(argv[3])


print 'L1_CALIB_FACTOR: %s' %L1_CALIB_FACTOR 
####### Get NTUPLES #########
#jet_l1_ntuple = ntuple_file.Get("jetL1Rate/Ntuple")
jet_l1_ntuple = ntuple_file.Get("jetL1Rate/Ntuple")
jet_uct_ntuple = ntuple_file.Get("jetUCTRate/Ntuple")
#jet_uct_ntuple_corr = ntuple_file.Get("corrjetUCTRate/Ntuple")

#rlx_tau_L1_ntuple = ntuple_file.Get("tauL1Rate/Ntuple")
#rlx_tau_UCT_ntuple = ntuple_file.Get("rlxTauUCTRate/Ntuple")
#iso_tau_L1_ntuple = ntuple_file.Get("isoTauEfficiency/Ntuple")
#iso_tau_UCT_ntuple = ntuple_file.Get("isoTauEfficiency/Ntuple")
#iso_eg_L1_ntuple = ntuple_file.Get("isoTauEfficiency/Ntuple")
#iso_eg_UCT_ntuple = ntuple_file.Get("isoTauEfficiency/Ntuple")
#rlx_eg_L1_ntuple = ntuple_file.Get("rlxEGEfficiency/Ntuple")
#rlx_eg_UCT_ntuple = ntuple_file.Get("rlxEGEfficiency/Ntuple")


###### PLOTTING #######
canvas = ROOT.TCanvas("asdf", "adsf", 800, 600)


###### Make a plot by drawing #########

def make_plot(tree, variable, selection, binning, xaxis='', title='', calFactor=1):
    ''' Plot a variable using draw and return the histogram '''
    draw_string = "%s * %0.2f>>htemp(%s)" % (variable, calFactor, ", ".join(str(x) for x in binning))
    print draw_string
    tree.Draw(draw_string, selection, "goff")
    output_histo = ROOT.gDirectory.Get("htemp").Clone()
    output_histo.GetXaxis().SetTitle(xaxis)
    output_histo.SetTitle(title)
    return output_histo

####### Get L1 Rate
def make_l1_rate(pt):
    ''' Make a rate plot out of L1Extra Pts '''
    numBins = pt.GetXaxis().GetNbins()
    rate = pt.Clone()

    for i in range(1, numBins+1):
        rate.SetBinContent(i, pt.Integral(i, numBins))

    rate.SetLineColor(ROOT.EColor.kRed)
    rate.SetMarkerStyle(20)
    rate.SetMarkerColor(ROOT.EColor.kRed)

    return rate


def make_l1g_rate(pt):
    ''' Return a rate plot out of UCT Pts '''
    #numBins = pt.GetNBinsX()
    numBins = pt.GetXaxis().GetNbins()
    rate = pt.Clone()

    for i in range(1, numBins+1):
        rate.SetBinContent(i, pt.Integral(i, numBins))

    rate.SetLineColor(ROOT.EColor.kBlue)
    rate.SetMarkerStyle(22)
    rate.SetMarkerColor(ROOT.EColor.kBlue)

    return rate


#### PLOT the RAtes
print 'here'  

def plotRates(l1ntuple, uctntuple, binning, filename, title='', xaxis='', isIso = False):
    ''' Save a rate Plot '''
    
    l1_pt = make_plot(
        l1ntuple, 'pt[0]',
        "pt[0]>30 && abs(eta[0])<3", # No selection
       #"pt[0]>30", # No selection
        #"abs(eta[0])<3", # No selection
        binning,
        '','', # Don't bother with titles
        L1_CALIB_FACTOR
        )
    l1g_pt = make_plot(
        uctntuple, 'pt[0]',
        "pt[0]>30 && abs(eta[0])<3", # No selection
       # "pt[0]>30", # No selection
       # "abs(eta[0])<3", # No selection
        binning,
        '','', # Don't bother with titles
        L1G_CALIB_FACTOR
        )

    l1Rate = make_l1_rate(l1_pt)
    l1gRate = make_l1g_rate(l1g_pt)
        


    canvas.SetLogy()
    l1gRate.SetTitle(title)
    l1gRate.GetXaxis().SetTitle(xaxis)
    l1gRate.Draw('ph')
    l1Rate.Draw('phsame')
    legend = ROOT.TLegend(0.7, 0.5, 0.89, 0.7, "", "brNDC")
    legend.SetFillColor(ROOT.EColor.kWhite)
    legend.SetBorderSize(1)
    legend.AddEntry(l1gRate, "UCT Corrected", "p")
    legend.AddEntry(l1Rate, "Current", "p")
    legend.Draw("same")
    canvas.SaveAs(filename)
##################################
#MakePLOTS

#################################
#plotRates(rlx_tau_ntuple, 20, [40, 0, 200],
#          '~/www/UCT2015_F/rlx_tau_rate.png',
#          "Tau Rate", "L1 p_{T} (GeV)")
#plotRates(iso_tau_ntuple, 20, [40, 0, 200],
#          '~/www/UCT2015_F/iso_tau_rate.png',
#          "IsoTau Rate", "L1 p_{T} (GeV)")
#plotRates(rlx_ntuple, 20, [40, 0, 200],
#          '~/www/UCT2015_F/rlx_tau_rate.png',
#          "Rlx EG Rate", "L1 p_{T} (GeV)")
plotRates(jet_l1_ntuple,jet_uct_ntuple, [40, 0, 200],
          '~/www/Research/PUMF/jet_rate_corr_etacut.png',
          "Jet Rate Pt>30 & abs(eta)<3", "P_{T} (GeV)")





