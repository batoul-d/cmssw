

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akCs5PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akCs5PFJets"),
    matched = cms.InputTag("ak5GenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.5
    )

akCs5PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("ak5GenJets"),
    matched = cms.InputTag("ak5GenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.5
    )

akCs5PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCs5PFJets")
                                                        )

akCs5PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akCs5PFJets"),
    payload = "AK5PF_offline"
    )

akCs5PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akCs5CaloJets'))

#akCs5PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak5GenJets'))

akCs5PFbTagger = bTaggers("akCs5PF",0.5,True,False)

#create objects locally since they dont load properly otherwise
#akCs5PFmatch = akCs5PFbTagger.match
akCs5PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCs5PFJets"), matched = cms.InputTag("genParticles"))
akCs5PFPatJetFlavourAssociationLegacy = akCs5PFbTagger.PatJetFlavourAssociationLegacy
akCs5PFPatJetPartons = akCs5PFbTagger.PatJetPartons
akCs5PFJetTracksAssociatorAtVertex = akCs5PFbTagger.JetTracksAssociatorAtVertex
akCs5PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akCs5PFSimpleSecondaryVertexHighEffBJetTags = akCs5PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCs5PFSimpleSecondaryVertexHighPurBJetTags = akCs5PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCs5PFCombinedSecondaryVertexBJetTags = akCs5PFbTagger.CombinedSecondaryVertexBJetTags
akCs5PFCombinedSecondaryVertexV2BJetTags = akCs5PFbTagger.CombinedSecondaryVertexV2BJetTags
akCs5PFJetBProbabilityBJetTags = akCs5PFbTagger.JetBProbabilityBJetTags
akCs5PFSoftPFMuonByPtBJetTags = akCs5PFbTagger.SoftPFMuonByPtBJetTags
akCs5PFSoftPFMuonByIP3dBJetTags = akCs5PFbTagger.SoftPFMuonByIP3dBJetTags
akCs5PFTrackCountingHighEffBJetTags = akCs5PFbTagger.TrackCountingHighEffBJetTags
akCs5PFTrackCountingHighPurBJetTags = akCs5PFbTagger.TrackCountingHighPurBJetTags
akCs5PFPatJetPartonAssociationLegacy = akCs5PFbTagger.PatJetPartonAssociationLegacy

akCs5PFImpactParameterTagInfos = akCs5PFbTagger.ImpactParameterTagInfos
akCs5PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCs5PFJetProbabilityBJetTags = akCs5PFbTagger.JetProbabilityBJetTags

akCs5PFSecondaryVertexTagInfos = akCs5PFbTagger.SecondaryVertexTagInfos
akCs5PFSimpleSecondaryVertexHighEffBJetTags = akCs5PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCs5PFSimpleSecondaryVertexHighPurBJetTags = akCs5PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCs5PFCombinedSecondaryVertexBJetTags = akCs5PFbTagger.CombinedSecondaryVertexBJetTags
akCs5PFCombinedSecondaryVertexV2BJetTags = akCs5PFbTagger.CombinedSecondaryVertexV2BJetTags

akCs5PFSecondaryVertexNegativeTagInfos = akCs5PFbTagger.SecondaryVertexNegativeTagInfos
akCs5PFNegativeSimpleSecondaryVertexHighEffBJetTags = akCs5PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akCs5PFNegativeSimpleSecondaryVertexHighPurBJetTags = akCs5PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akCs5PFNegativeCombinedSecondaryVertexBJetTags = akCs5PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akCs5PFPositiveCombinedSecondaryVertexBJetTags = akCs5PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akCs5PFNegativeCombinedSecondaryVertexV2BJetTags = akCs5PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akCs5PFPositiveCombinedSecondaryVertexV2BJetTags = akCs5PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akCs5PFSoftPFMuonsTagInfos = akCs5PFbTagger.SoftPFMuonsTagInfos
akCs5PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCs5PFSoftPFMuonBJetTags = akCs5PFbTagger.SoftPFMuonBJetTags
akCs5PFSoftPFMuonByIP3dBJetTags = akCs5PFbTagger.SoftPFMuonByIP3dBJetTags
akCs5PFSoftPFMuonByPtBJetTags = akCs5PFbTagger.SoftPFMuonByPtBJetTags
akCs5PFNegativeSoftPFMuonByPtBJetTags = akCs5PFbTagger.NegativeSoftPFMuonByPtBJetTags
akCs5PFPositiveSoftPFMuonByPtBJetTags = akCs5PFbTagger.PositiveSoftPFMuonByPtBJetTags
akCs5PFPatJetFlavourIdLegacy = cms.Sequence(akCs5PFPatJetPartonAssociationLegacy*akCs5PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub
akCs5PFPatJetFlavourAssociation = akCs5PFbTagger.PatJetFlavourAssociation
akCs5PFPatJetFlavourId = cms.Sequence(akCs5PFPatJetPartons*akCs5PFPatJetFlavourAssociation)

#adding the subjet taggers
#SUBJETDUMMY_akCs5PFSubjetImpactParameterTagInfos = akCs5PFbTagger.SubjetImpactParameterTagInfos
#SUBJETDUMMY_akCs5PFSubjetJetProbabilityBJetTags = akCs5PFbTagger.SubjetJetProbabilityBJetTags
#SUBJETDUMMY_akCs5PFSubjetSecondaryVertexTagInfos = akCs5PFbTagger.SubjetSecondaryVertexTagInfos
#SUBJETDUMMY_akCs5PFSubjetSecondaryVertexNegativeTagInfos = akCs5PFbTagger.SubjetSecondaryVertexNegativeTagInfos
#SUBJETDUMMY_akCs5PFSubjetJetTracksAssociatorAtVertex = akCs5PFbTagger.SubjetJetTracksAssociatorAtVertex
#SUBJETDUMMY_akCs5PFCombinedSubjetSecondaryVertexBJetTags = akCs5PFbTagger.CombinedSubjetSecondaryVertexBJetTags
#SUBJETDUMMY_akCs5PFCombinedSubjetSecondaryVertexV2BJetTags = akCs5PFbTagger.CombinedSubjetSecondaryVertexV2BJetTags
#SUBJETDUMMY_akCs5PFCombinedSubjetNegativeSecondaryVertexV2BJetTags = akCs5PFbTagger.CombinedSubjetNegativeSecondaryVertexV2BJetTags

akCs5PFJetBtaggingIP       = cms.Sequence(akCs5PFImpactParameterTagInfos *
            (akCs5PFTrackCountingHighEffBJetTags +
             akCs5PFTrackCountingHighPurBJetTags +
             akCs5PFJetProbabilityBJetTags +
             akCs5PFJetBProbabilityBJetTags 
            )
            )

akCs5PFJetBtaggingSV = cms.Sequence(akCs5PFImpactParameterTagInfos
            *
            akCs5PFSecondaryVertexTagInfos
            * (akCs5PFSimpleSecondaryVertexHighEffBJetTags+
                akCs5PFSimpleSecondaryVertexHighPurBJetTags+
                akCs5PFCombinedSecondaryVertexBJetTags+
                akCs5PFCombinedSecondaryVertexV2BJetTags
              )
            )

akCs5PFJetBtaggingNegSV = cms.Sequence(akCs5PFImpactParameterTagInfos
            *
            akCs5PFSecondaryVertexNegativeTagInfos
            * (akCs5PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCs5PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCs5PFNegativeCombinedSecondaryVertexBJetTags+
                akCs5PFPositiveCombinedSecondaryVertexBJetTags+
                akCs5PFNegativeCombinedSecondaryVertexV2BJetTags+
                akCs5PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akCs5PFJetBtaggingMu = cms.Sequence(akCs5PFSoftPFMuonsTagInfos * (akCs5PFSoftPFMuonBJetTags
                +
                akCs5PFSoftPFMuonByIP3dBJetTags
                +
                akCs5PFSoftPFMuonByPtBJetTags
                +
                akCs5PFNegativeSoftPFMuonByPtBJetTags
                +
                akCs5PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akCs5PFJetBtagging = cms.Sequence(akCs5PFJetBtaggingIP
            *akCs5PFJetBtaggingSV
            *akCs5PFJetBtaggingNegSV
#            *akCs5PFJetBtaggingMu
            )

akCs5PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCs5PFJets"),
        genJetMatch          = cms.InputTag("akCs5PFmatch"),
        genPartonMatch       = cms.InputTag("akCs5PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akCs5PFcorr")),
        #JetPartonMapSource   = cms.InputTag("akCs5PFPatJetFlavourAssociationLegacy"),
        JetPartonMapSource   = cms.InputTag("akCs5PFPatJetFlavourAssociation"),
	JetFlavourInfoSource   = cms.InputTag("akCs5PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCs5PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = False,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCs5PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCs5PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCs5PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCs5PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCs5PFJetBProbabilityBJetTags"),
            cms.InputTag("akCs5PFJetProbabilityBJetTags"),
            #cms.InputTag("akCs5PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akCs5PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akCs5PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akCs5PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCs5PFJetID"),
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

akCs5PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akCs5PFJets"),
           	    R0  = cms.double( 0.5)
)
#ppDataDummy_akCs5PFpatJetsWithBtagging.userData.userFloats.src += ['akCs5PFNjettiness:tau1','akCs5PFNjettiness:tau2','akCs5PFNjettiness:tau3']

akCs5PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akCs5PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak5GenJets',
                                                             rParam = 0.5,
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
                                                             bTagJetName = cms.untracked.string("akCs5PF"),
                                                             jetName = cms.untracked.string("akCs5PF"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(False),
							     doSubJets = cms.untracked.bool(False),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("ak5GenJets"),
							     doExtendedFlavorTagging = cms.untracked.bool(False),
							     jetFlavourInfos = cms.InputTag("akCs5PFPatJetFlavourAssociation"),
							     subjetFlavourInfos = cms.InputTag("akCs5PFPatJetFlavourAssociation","SubJets"),
							     groomedJets = cms.InputTag("akCs5PFJets"),
							     isPythia6 = cms.untracked.bool(False),
                                                             doGenTaus = False
                                                            )

akCs5PFJetSequence_mc = cms.Sequence(
                                                  #akCs5PFclean
                                                  #*
                                                  akCs5PFmatch
                                                  #*
                                                  #akCs5PFmatchGroomed
                                                  *
                                                  akCs5PFparton
                                                  *
                                                  akCs5PFcorr
                                                  *
                                                  #akCs5PFJetID
                                                  #*
                                                  #akCs5PFPatJetFlavourIdLegacy  # works for PbPb
                                                  #*
			                          akCs5PFPatJetFlavourId  # doesn't work for PbPb yet
                                                  *
                                                  akCs5PFJetTracksAssociatorAtVertex
                                                  *
                                                  akCs5PFJetBtagging
                                                  *
                                                  #ppDataDummy_akCs5PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  #ppDataDummy_*
                                                  akCs5PFpatJetsWithBtagging
                                                  *
                                                  akCs5PFJetAnalyzer
                                                  )

akCs5PFJetSequence_data = cms.Sequence(akCs5PFcorr
                                                    *
                                                    #akCs5PFJetID
                                                    #*
                                                    akCs5PFJetTracksAssociatorAtVertex
                                                    *
                                                    akCs5PFJetBtagging
                                                    *
                                                    #ppDataDummy_akCs5PFNjettiness 
                                                    #ppDataDummy_*
                                                    akCs5PFpatJetsWithBtagging
                                                    *
                                                    akCs5PFJetAnalyzer
                                                    )

akCs5PFJetSequence_jec = cms.Sequence(akCs5PFJetSequence_mc)
akCs5PFJetSequence_mb = cms.Sequence(akCs5PFJetSequence_mc)

akCs5PFJetSequence = cms.Sequence(akCs5PFJetSequence_data)
