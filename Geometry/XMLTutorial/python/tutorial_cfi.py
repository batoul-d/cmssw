import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/materials.xml',
                               'Geometry/CMSCommonData/data/rotations.xml',
                               'Geometry/CMSCommonData/data/extend/cmsextent.xml',
                               'Geometry/CMSCommonData/data/cms.xml',
                               'Geometry/CMSCommonData/data/cmsMother.xml',
                               'Geometry/XMLTutorial/data/main.xml'
                               ),
    rootNodeName = cms.string('cms:OCMS')
)
