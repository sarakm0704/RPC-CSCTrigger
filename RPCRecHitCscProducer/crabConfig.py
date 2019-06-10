from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.requestName = 'CSCTriggerPrimitives_withRPCRecHitGEN'

config.section_('JobType')
config.JobType.psetName = 'RPCRecHitCscProducer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['RPCoutput_withRecHitLCT_GEN.root']

config.section_('Data')
config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
config.Data.publication = False
config.Data.unitsPerJob = 180

config.section_('Site')
#config.Site.blacklist = ['T_IT_Legnaro']
#config.Site.whitelist = ['T2_IT_Bari']
config.Site.storageSite = 'T2_KR_KNU'

