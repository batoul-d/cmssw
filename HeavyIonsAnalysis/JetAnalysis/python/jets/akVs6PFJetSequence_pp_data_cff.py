

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akVs6PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akVs6PFJets"),
    matched = cms.InputTag("ak6GenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.6
    )

akVs6PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("ak6GenJets"),
    matched = cms.InputTag("ak6GenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.6
    )

akVs6PFparton = patJetPartonMatch.clone(src = cms.InputTag("akVs6PFJets")
                                                        )

akVs6PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akVs6PFJets"),
    payload = "AK6PF_offline"
    )

akVs6PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akVs6CaloJets'))

#akVs6PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak6GenJets'))

akVs6PFbTagger = bTaggers("akVs6PF",0.6,True,False)

#create objects locally since they dont load properly otherwise
#akVs6PFmatch = akVs6PFbTagger.match
akVs6PFparton = patJetPartonMatch.clone(src = cms.InputTag("akVs6PFJets"), matched = cms.InputTag("genParticles"))
akVs6PFPatJetFlavourAssociationLegacy = akVs6PFbTagger.PatJetFlavourAssociationLegacy
akVs6PFPatJetPartons = akVs6PFbTagger.PatJetPartons
akVs6PFJetTracksAssociatorAtVertex = akVs6PFbTagger.JetTracksAssociatorAtVertex
akVs6PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akVs6PFSimpleSecondaryVertexHighEffBJetTags = akVs6PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akVs6PFSimpleSecondaryVertexHighPurBJetTags = akVs6PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akVs6PFCombinedSecondaryVertexBJetTags = akVs6PFbTagger.CombinedSecondaryVertexBJetTags
akVs6PFCombinedSecondaryVertexV2BJetTags = akVs6PFbTagger.CombinedSecondaryVertexV2BJetTags
akVs6PFJetBProbabilityBJetTags = akVs6PFbTagger.JetBProbabilityBJetTags
akVs6PFSoftPFMuonByPtBJetTags = akVs6PFbTagger.SoftPFMuonByPtBJetTags
akVs6PFSoftPFMuonByIP3dBJetTags = akVs6PFbTagger.SoftPFMuonByIP3dBJetTags
akVs6PFTrackCountingHighEffBJetTags = akVs6PFbTagger.TrackCountingHighEffBJetTags
akVs6PFTrackCountingHighPurBJetTags = akVs6PFbTagger.TrackCountingHighPurBJetTags
akVs6PFPatJetPartonAssociationLegacy = akVs6PFbTagger.PatJetPartonAssociationLegacy

akVs6PFImpactParameterTagInfos = akVs6PFbTagger.ImpactParameterTagInfos
akVs6PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akVs6PFJetProbabilityBJetTags = akVs6PFbTagger.JetProbabilityBJetTags

akVs6PFSecondaryVertexTagInfos = akVs6PFbTagger.SecondaryVertexTagInfos
akVs6PFSimpleSecondaryVertexHighEffBJetTags = akVs6PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akVs6PFSimpleSecondaryVertexHighPurBJetTags = akVs6PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akVs6PFCombinedSecondaryVertexBJetTags = akVs6PFbTagger.CombinedSecondaryVertexBJetTags
akVs6PFCombinedSecondaryVertexV2BJetTags = akVs6PFbTagger.CombinedSecondaryVertexV2BJetTags

akVs6PFSecondaryVertexNegativeTagInfos = akVs6PFbTagger.SecondaryVertexNegativeTagInfos
akVs6PFNegativeSimpleSecondaryVertexHighEffBJetTags = akVs6PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akVs6PFNegativeSimpleSecondaryVertexHighPurBJetTags = akVs6PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akVs6PFNegativeCombinedSecondaryVertexBJetTags = akVs6PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akVs6PFPositiveCombinedSecondaryVertexBJetTags = akVs6PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akVs6PFNegativeCombinedSecondaryVertexV2BJetTags = akVs6PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akVs6PFPositiveCombinedSecondaryVertexV2BJetTags = akVs6PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akVs6PFSoftPFMuonsTagInfos = akVs6PFbTagger.SoftPFMuonsTagInfos
akVs6PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akVs6PFSoftPFMuonBJetTags = akVs6PFbTagger.SoftPFMuonBJetTags
akVs6PFSoftPFMuonByIP3dBJetTags = akVs6PFbTagger.SoftPFMuonByIP3dBJetTags
akVs6PFSoftPFMuonByPtBJetTags = akVs6PFbTagger.SoftPFMuonByPtBJetTags
akVs6PFNegativeSoftPFMuonByPtBJetTags = akVs6PFbTagger.NegativeSoftPFMuonByPtBJetTags
akVs6PFPositiveSoftPFMuonByPtBJetTags = akVs6PFbTagger.PositiveSoftPFMuonByPtBJetTags
akVs6PFPatJetFlavourIdLegacy = cms.Sequence(akVs6PFPatJetPartonAssociationLegacy*akVs6PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub
akVs6PFPatJetFlavourAssociation = akVs6PFbTagger.PatJetFlavourAssociation
akVs6PFPatJetFlavourId = cms.Sequence(akVs6PFPatJetPartons*akVs6PFPatJetFlavourAssociation)

#adding the subjet taggers
#SUBJETDUMMY_akVs6PFSubjetImpactParameterTagInfos = akVs6PFbTagger.SubjetImpactParameterTagInfos
#SUBJETDUMMY_akVs6PFSubjetJetProbabilityBJetTags = akVs6PFbTagger.SubjetJetProbabilityBJetTags
#SUBJETDUMMY_akVs6PFSubjetSecondaryVertexTagInfos = akVs6PFbTagger.SubjetSecondaryVertexTagInfos
#SUBJETDUMMY_akVs6PFSubjetSecondaryVertexNegativeTagInfos = akVs6PFbTagger.SubjetSecondaryVertexNegativeTagInfos
#SUBJETDUMMY_akVs6PFSubjetJetTracksAssociatorAtVertex = akVs6PFbTagger.SubjetJetTracksAssociatorAtVertex
#SUBJETDUMMY_akVs6PFCombinedSubjetSecondaryVertexBJetTags = akVs6PFbTagger.CombinedSubjetSecondaryVertexBJetTags
#SUBJETDUMMY_akVs6PFCombinedSubjetSecondaryVertexV2BJetTags = akVs6PFbTagger.CombinedSubjetSecondaryVertexV2BJetTags
#SUBJETDUMMY_akVs6PFCombinedSubjetNegativeSecondaryVertexV2BJetTags = akVs6PFbTagger.CombinedSubjetNegativeSecondaryVertexV2BJetTags

akVs6PFJetBtaggingIP       = cms.Sequence(akVs6PFImpactParameterTagInfos *
            (akVs6PFTrackCountingHighEffBJetTags +
             akVs6PFTrackCountingHighPurBJetTags +
             akVs6PFJetProbabilityBJetTags +
             akVs6PFJetBProbabilityBJetTags 
            )
            )

akVs6PFJetBtaggingSV = cms.Sequence(akVs6PFImpactParameterTagInfos
            *
            akVs6PFSecondaryVertexTagInfos
            * (akVs6PFSimpleSecondaryVertexHighEffBJetTags+
                akVs6PFSimpleSecondaryVertexHighPurBJetTags+
                akVs6PFCombinedSecondaryVertexBJetTags+
                akVs6PFCombinedSecondaryVertexV2BJetTags
              )
            )

akVs6PFJetBtaggingNegSV = cms.Sequence(akVs6PFImpactParameterTagInfos
            *
            akVs6PFSecondaryVertexNegativeTagInfos
            * (akVs6PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akVs6PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akVs6PFNegativeCombinedSecondaryVertexBJetTags+
                akVs6PFPositiveCombinedSecondaryVertexBJetTags+
                akVs6PFNegativeCombinedSecondaryVertexV2BJetTags+
                akVs6PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akVs6PFJetBtaggingMu = cms.Sequence(akVs6PFSoftPFMuonsTagInfos * (akVs6PFSoftPFMuonBJetTags
                +
                akVs6PFSoftPFMuonByIP3dBJetTags
                +
                akVs6PFSoftPFMuonByPtBJetTags
                +
                akVs6PFNegativeSoftPFMuonByPtBJetTags
                +
                akVs6PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akVs6PFJetBtagging = cms.Sequence(akVs6PFJetBtaggingIP
            *akVs6PFJetBtaggingSV
            *akVs6PFJetBtaggingNegSV
#            *akVs6PFJetBtaggingMu
            )

akVs6PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akVs6PFJets"),
        genJetMatch          = cms.InputTag("akVs6PFmatch"),
        genPartonMatch       = cms.InputTag("akVs6PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akVs6PFcorr")),
        #JetPartonMapSource   = cms.InputTag("akVs6PFPatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akVs6PFPatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akVs6PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akVs6PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akVs6PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akVs6PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akVs6PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akVs6PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akVs6PFJetBProbabilityBJetTags"),
            cms.InputTag("akVs6PFJetProbabilityBJetTags"),
            #cms.InputTag("akVs6PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akVs6PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akVs6PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akVs6PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akVs6PFJetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = False,
        addGenPartonMatch = False,
        addGenJetMatch = False,
        embedGenJetMatch = False,
        embedGenPartonMatch = False,
        # embedCaloTowers = False,
        # embedPFCandidates = True
        )

akVs6PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akVs6PFJets"),
           	    R0  = cms.double( 0.6)
)
#ppDataDummy_akVs6PFpatJetsWithBtagging.userData.userFloats.src += ['akVs6PFNjettiness:tau1','akVs6PFNjettiness:tau2','akVs6PFNjettiness:tau3']

akVs6PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akVs6PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak6GenJets',
                                                             rParam = 0.6,
                                                             matchJets = cms.untracked.bool(False),
                                                             matchTag = 'patJetsWithBtagging',
                                                             pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
                                                             trackTag = cms.InputTag("generalTracks"),
                                                             fillGenJets = False,
                                                             isMC = False,
							     doSubEvent = False,
                                                             useHepMC = cms.untracked.bool(False),
							     genParticles = cms.untracked.InputTag("genParticles"),
							     eventInfoTag = cms.InputTag("generator"),
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akVs6PF"),
                                                             jetName = cms.untracked.string("akVs6PF"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(False),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("ak6GenJets"),
							     doExtendedFlavorTagging = cms.untracked.bool(False),
							     jetFlavourInfos = cms.InputTag("akVs6PFPatJetFlavourAssociation"),
							     subjetFlavourInfos = cms.InputTag("akVs6PFPatJetFlavourAssociation","SubJets"),
							     groomedJets = cms.InputTag("akVs6PFJets"),
							     isPythia6 = cms.untracked.bool(False),
                                                             doGenTaus = False
                                                            )

akVs6PFJetSequence_mc = cms.Sequence(
                                                  #akVs6PFclean
                                                  #*
                                                  akVs6PFmatch
                                                  #*
                                                  #akVs6PFmatchGroomed
                                                  *
                                                  akVs6PFparton
                                                  *
                                                  akVs6PFcorr
                                                  *
                                                  #akVs6PFJetID
                                                  #*
                                                  #akVs6PFPatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akVs6PFPatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akVs6PFJetTracksAssociatorAtVertex
                                                  *
                                                  akVs6PFJetBtagging
                                                  *
                                                  #ppDataDummy_akVs6PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  #ppDataDummy_*
                                                  akVs6PFpatJetsWithBtagging
                                                  *
                                                  akVs6PFJetAnalyzer
                                                  )

akVs6PFJetSequence_data = cms.Sequence(akVs6PFcorr
                                                    *
                                                    #akVs6PFJetID
                                                    #*
                                                    akVs6PFJetTracksAssociatorAtVertex
                                                    *
                                                    akVs6PFJetBtagging
                                                    *
                                                    #ppDataDummy_akVs6PFNjettiness 
                                                    #ppDataDummy_*
                                                    akVs6PFpatJetsWithBtagging
                                                    *
                                                    akVs6PFJetAnalyzer
                                                    )

akVs6PFJetSequence_jec = cms.Sequence(akVs6PFJetSequence_mc)
akVs6PFJetSequence_mb = cms.Sequence(akVs6PFJetSequence_mc)

akVs6PFJetSequence = cms.Sequence(akVs6PFJetSequence_data)
