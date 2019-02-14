import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.rerecoGen_cff import *
from HeavyIonsAnalysis.JetAnalysis.rerecoRho_cff import *
from HeavyIonsAnalysis.JetAnalysis.rerecoJets_cff import *
from HeavyIonsAnalysis.JetAnalysis.rerecoTracks_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.akPu3CaloJetSequence_pponPbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_pponPbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs3PFJetSequence_pponPbPb_mc_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_pponPbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_pponPbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_pponPbPb_mc_cff import *

genSignalSequence = cms.Sequence(
    genParticlesForJets +

    hiSignalGenParticles +
    genParticlesForJetsSignal +

    #ak3HiGenJets +
    ak4HiGenJets +

    signalPartons +

    #ak3HiSignalGenJets +
    ak4HiSignalGenJets +

    #ak3HiGenNjettiness +
    ak4HiGenNjettiness
)

genCleanedSequence = cms.Sequence(
    genParticlesForJets +

    #ak3HiGenJets +
    ak4HiGenJets +

    myPartons +
    cleanedPartons +

    #ak3HiCleanedGenJets +
    ak4HiCleanedGenJets
)

## This jet sequence has been horribly hacked for the onia-jet analysis

# cluster AK including J/psi, but with no UE subtraction
ak4PFJetsWithJPsi = ak4PFJets.clone(
    src = 'pfCandComposites',
    jetCollInstanceName = cms.string('pfParticlesWithJPsi'),
    writeJetsWithConst = cms.bool(True)
)

#filter out the jet containing the jpsi
ak4PFXJetsWithJPsi = cms.EDFilter("PFJetXSelector",
                            src = cms.InputTag("ak4PFJetsWithJPsi"),
                            cut = cms.string("pt > 0.0 && abs(rapidity()) < 3.0")
)

# Now remove the J/Psi
csCandsNoJPsi = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
                             src = cms.InputTag("ak4PFXJetsWithJPsi","constituents"),
                             pdgId = cms.vint32(211,11,13,-11,-13,22,130)
                         )

# run the Cs subtraction on the J/Psi-jet, ignoring the Jpsi 
akCs4PFJets.src = cms.InputTag("csCandsNoJPsi")
akCs4PFJets.jetCollInstanceName = cms.string('pfParticlesCsNoJPsi')

### put the J/Psi back with the candidates

#here's the jpsi
pfCandJPsi = cms.EDFilter("PdgIdPFCandidateSelector",
                          src = cms.InputTag("pfCandComposites"),
                          pdgId = cms.vint32(1)
)


#pfCandsCs = cms.EDProducer("PFCandidateFromFwdPtrProducer",
#                          src = cms.InputTag("akCs4PFJets","pfParticlesCs")
#)

pfCandsCsPlusJPsi = cms.EDProducer(
    "PFCandidateListMerger",
    src = cms.VInputTag(cms.InputTag("akCs4PFJets","pfParticlesCsNoJPsi"),
                        cms.InputTag("pfCandJPsi"))    
)



akCs4PFJetsWithJPsi = ak4PFJets.clone(
    src = cms.InputTag("pfCandsCsPlusJPsi"),
    jetCollInstanceName = cms.string('pfParticlesCsWithJPsi'),
    writeJetsWithConst = cms.bool(True)
)

#akCs4PFJetAnalyzer.jetTag = "akCs4PFJetsWithJPsi"
akCs4PFmatch.src = "akCs4PFJetsWithJPsi"
akCs4PFparton.src = "akCs4PFJetsWithJPsi"
akCs4PFcorr.src = "akCs4PFJetsWithJPsi"
akCs4PFNjettiness.src = "akCs4PFJetsWithJPsi"

akCs4PFpatJetsWithBtagging.jetSource = "akCs4PFJetsWithJPsi"
akCs4PFpatJets.jetSource = "akCs4PFJetsWithJPsi"
akCs4PFPatJetPartonAssociationLegacy.jets = "akCs4PFJetsWithJPsi"

akCs4PFJetAnalyzer.pfCandidateLabel = cms.untracked.InputTag('pfCandsCsPlusJPsi')

jetSequence = cms.Sequence(
    rhoSequence +

    highPurityTracks +

    #akPu3CaloJets +
    #akPu3PFJets +
    #akCs3PFJets +

    #akPu4CaloJets +
    #akPu4PFJets +
    ak4PFJetsWithJPsi +
    ak4PFXJetsWithJPsi +
    csCandsNoJPsi +
    akCs4PFJets +
    pfCandJPsi +
    pfCandsCsPlusJPsi +
    akCs4PFJetsWithJPsi +
    akCs4PFJetSequence
)
