import re

import FWCore.ParameterSet.Config as cms

def getMicroAODHLTFilter(datasetName, options):
    """
    Return a TriggerResultsFilter customized accordingly to the data-taking period.
    The returned module is inteded to be used as first step in the MicroAOD production sequence
    """

    for tag, trgs in options["TriggerPaths"].items():
        triggerConditions = cms.vstring()
        print 'ED DEBUG a'
        print tag
        print datasetName
        print re.search(tag, datasetName)
        for trg in trgs:
            triggerConditions.append(str(trg))
        print 'ED DEBUG b'
        #if re.search(tag, datasetName):
        if tag.count(datasetName):
            print 'Only events from these path will be processed: '
            print triggerConditions
            triggerSelectionModule = cms.EDFilter("TriggerResultsFilter",
                                                  triggerConditions = triggerConditions,
                                                  hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                                  l1tResults = cms.InputTag( "" ),
                                                  throw = cms.bool(False)
                                              )
            print 'about to return'
            return triggerSelectionModule
        print 'ED DEBUG c'

    return None
