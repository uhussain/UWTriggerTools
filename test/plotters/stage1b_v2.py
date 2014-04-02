import os
import hashlib
import glob
import math

import ROOT
#ROOT.gSystem.Load("libL1TriggerUCT2015")

colors = ROOT.EColor

basedir = "/scratch/efriis/data/uct/"
standard_rate_files = glob.glob(
    basedir + "2013-06-04-Stage1B-Norm-makeRateTrees_cfg/*root")

hcal_rate_files = glob.glob(
    basedir + "2013-06-04-Stage1B-HCAL-makeRateTrees_cfg/*root")

standard_eff_files = glob.glob(
    basedir + "2013-06-04-Stage1B-PU50-makeEfficiencyTree_cfg/*root")

hcal_eff_files = glob.glob(
    basedir + "2013-06-04-Stage1B-PU50-HCALONLY-makeEfficiencyTree_cfg/*root")

colors = ROOT.EColor
basedir = "/scratch/efriis/data/uct/"

standard_rate_files = glob.glob(
    basedir + "2013-06-04-Stage1B-Norm-makeRateTrees_cfg/*root")

hcal_rate_files = glob.glob(
    basedir + "2013-06-04-Stage1B-HCAL-makeRateTrees_cfg/*root")

standard_eff_files = glob.glob(
    basedir + "2013-06-04-Stage1B-PU50-makeEfficiencyTree_cfg/*root")

hcal_eff_files = glob.glob(
    basedir + "2013-06-04-Stage1B-PU50-HCALONLY-makeEfficiencyTree_cfg/*root")


def chain(name, files):
    chain = ROOT.TChain(name)
    for file in files:
        chain.Add(file)
    return chain

trees = {
    'tau': {
        'eff': {
            "stage1": chain("rlxTauEfficiency/Ntuple", standard_eff_files),
            "stage1b": chain("rlxTauEfficiencyStage1B/Ntuple",
                             hcal_eff_files),
        },
        'rate': {
            "current": chain("tauL1Rate/Ntuple", standard_rate_files),
            "stage1": chain("rlxTauUCTRate/Ntuple", standard_rate_files),
            "stage1b": chain("rlxTauUCTRateStage1B/Ntuple", hcal_rate_files),
        }
    },
    'eg': {
        'eff': {
            "stage1": chain("rlxEGEfficiency/Ntuple", standard_eff_files),
            "stage1b": chain("rlxEGEfficiencyStage1B/Ntuple", hcal_eff_files),
        },
        'rate': {
            "current": chain("rlxEGL1Rate/Ntuple", standard_rate_files),
            "stage1": chain("rlxEGUCTRate/Ntuple", standard_rate_files),
            "stage1b": chain("rlxEGUCTRateStage1B/Ntuple", hcal_rate_files),
        }
    }
}

for obj in ['eg', 'tau']:
    trees[obj]['rate']['stage1b'].SetAlias(
        #"ecalPt", "emClusterEt[0]")
        "ecalPt", "emClusterEt")

    trees[obj]['rate']['stage1b'].SetAlias(
        "ecalIso", "max(l1gJetPtEM[0] - pt[0] - puUICEM[0]*0.13*513627, 0)")

    trees[obj]['rate']['stage1b'].SetAlias(
        "hcalIso",
        "max(l1gJetPt[0] - (7./9)*puUICEM[0]*0.13*953311 "
        "- region2Disc[0].totalEt, 0)")

    trees[obj]['rate']['stage1b'].SetAlias(
        "hcalCore",
        "max(-(2./9)*puUICEM[0]*0.13*953311 + region2Disc[0].totalEt, 0)")

    trees[obj]['eff']['stage1b'].SetAlias(
        "ecalPt", "l1gEmClusterEt")

    trees[obj]['eff']['stage1b'].SetAlias(
        "ecalIso",
        "max(l1gJetPtEM - l1gEmClusterEt - (8./9)*l1gPUUICEM*0.13*513627, 0)")

    trees[obj]['eff']['stage1b'].SetAlias(
        "hcalIso",
        "max(l1gJetPt - (7./9)*l1gPUUICEM*0.13*953311 - l1g2RegionEt, 0)")

    trees[obj]['eff']['stage1b'].SetAlias(
        "hcalCore",
        "max(-(2./9)*l1gPUUICEM*0.13*953311 + l1g2RegionEt, 0)")

fout = None
signal = None
bkg = None
if True or not os.path.exists("test.root"):
    fout = ROOT.TFile("test.root", "RECREATE")
    signal = ROOT.TNtuple(
        "signal", "signal",
        ":".join(["ecalPt", "ecalIso", "hcalIso"])
    )
    bkg = ROOT.TNtuple(
        "bkg", "bkg",
        ":".join(["ecalPt", "ecalIso", "hcalIso"])
    )

    for row in trees['eg']['eff']['stage1b']:
        if row.l1gPt < 25:
            continue
        signal.Fill(
            row.l1gPt,
            max(row.l1gJetPtEM - row.l1gPt
                - row.l1gPUUICEM * 0.13 * 513627, 0),
            max(row.l1gJetPtEM - row.l1gPt
                - row.l1gPUUICEM * 0.13 * 513627, 0),

        )
    for row in trees['eg']['rate']['stage1b']:
        if row.pt[0] < 25:
            continue
        bkg.Fill(
            row.pt[0],
            max(row.jetPtEM[0] - row.pt[0]
                - row.puUICEM[0] * 0.13 * 513627, 0)
        )


# Do BDT
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory(
    "TMVAClassification", fout,
    ":".join([
        "!V",
        "!Silent",
        "Color",
        "DrawProgressBar",
        "AnalysisType=Classification"]
    ))

factory.AddVariable("ecalPt", "F")
#factory.AddVariable("ecalPt", "F")
factory.AddVariable("ecalIso", "F")
#factory.AddVariable("hcalIso", "F")
#factory.AddVariable("hcalCore", "F")

factory.AddSignalTree(signal)
factory.AddBackgroundTree(bkg)

sigCut = ROOT.TCut("1")
bgCut = ROOT.TCut("1")

factory.PrepareTrainingAndTestTree(
    sigCut,   # signal events
    bgCut,    # background events
    ":".join([
        "nTrain_Signal=0",
        "nTrain_Background=0",
        "SplitMode=Random",
        "NormMode=NumEvents",
        "!V"
    ]))

method = factory.BookMethod(
    ROOT.TMVA.Types.kBDT, "BDT",
    ":".join([
        "!H",
        "!V",
        "NTrees=5",
        "nEventsMin=150",
        "MaxDepth=3",
        "BoostType=AdaBoost",
        "AdaBoostBeta=0.5",
        "SeparationType=GiniIndex",
        "nCuts=20",
        "PruneMethod=NoPruning",
    ]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
