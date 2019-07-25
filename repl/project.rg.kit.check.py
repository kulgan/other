from collections import Counter as table



rgs = g.nodes(ReadGroup).props(library_strategy='WXS', sequencing_center='BI').all()


kits = []

for rg in rgs:
  kits.append((rg.sequencing_center, rg.project_id, rg.target_capture_kit))



{(u'BI', u'CPTAC-3', None),
 (u'BI', u'CPTAC-3', u'Nextera Rapid Capture Exome v1.2'),
 (u'BI', u'HCMI-CMDC', u'TruSeq Exome Enrichment - 62 Mb'),
 (u'BI',
  u'TARGET-NBL',
  u'Custom SureSelect Human All Exon v1.1 Plus 3 Boosters')}
