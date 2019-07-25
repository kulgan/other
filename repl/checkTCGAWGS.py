import pandas as pd
tcga = [
 u'TCGA-ACC',
 u'TCGA-BLCA',
 u'TCGA-BRCA',
 u'TCGA-CESC',
 u'TCGA-CHOL',
 u'TCGA-COAD',
 u'TCGA-DLBC',
 u'TCGA-ESCA',
 u'TCGA-GBM',
 u'TCGA-HNSC',
 u'TCGA-KICH',
 u'TCGA-KIRC',
 u'TCGA-KIRP',
 u'TCGA-LAML',
 u'TCGA-LGG',
 u'TCGA-LIHC',
 u'TCGA-LUAD',
 u'TCGA-LUSC',
 u'TCGA-MESO',
 u'TCGA-OV',
 u'TCGA-PAAD',
 u'TCGA-PCPG',
 u'TCGA-PRAD',
 u'TCGA-READ',
 u'TCGA-SARC',
 u'TCGA-SKCM',
 u'TCGA-STAD',
 u'TCGA-TGCT',
 u'TCGA-THCA',
 u'TCGA-THYM',
 u'TCGA-UCEC',
 u'TCGA-UCS',
 u'TCGA-UVM']

# Get TCGA ATAC-Seq aliquot
surs = g.nodes(SubmittedUnalignedReads).props(experimental_strategy = 'ATAC-Seq').all()
a = set()
for sur in surs:
  a.add(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst._ReadGroupDerivedFromAliquot_out[0].dst.node_id)
with open('atac.aliquot.txt', 'w') as f:
  for item in a:
  f.write("%s\n" % item)
    
a = set()
for sur in surs:
  a.add(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst._ReadGroupDerivedFromAliquot_out[0].dst)

d = set()
for b in a:
  for rg in b._ReadGroupDerivedFromAliquot_in:
    if rg.src.library_strategy == 'WGS': 
      d = d.add(b)

atacsamples = set()
for sur in g.nodes(SubmittedUnalignedReads).props(experimental_strategy = 'ATAC-Seq').all():
  atacsamples.add(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst._ReadGroupDerivedFromAliquot_out[0].dst.edges_out[0].dst.edges_out[0].dst.edges_out[0].dst)

wgssamples = set()
for sur in g.nodes(SubmittedUnalignedReads).props(experimental_strategy = 'WGS').all():
  wgssamples.add(sur._SubmittedUnalignedReadsDataFromReadGroup_out[0].dst._ReadGroupDerivedFromAliquot_out[0].dst.edges_out[0].dst.edges_out[0].dst.edges_out[0].dst)

s = wgssamples.intersection(atacsamples)


# Get TCGA Read Groups
rgData = pd.DataFrame(columns = ['project', \
                                  'rg', \
                                  'rg.experiment_name', \
                                  'rg.sequencing_center', \
                                  'rg_name', \
                                  'rg.is_paired_end', \
                                  'aliquot', \
                                  'aliquot.submitter_id', \
                                  'sr_type', \
                                  'sr'])
i = 0
for project in tcga:
    rgs = g.nodes(ReadGroup).props(library_strategy = 'WGS', project_id = project).all()
    for rg in rgs:
        aliquot = rg._ReadGroupDerivedFromAliquot_out[0].dst
        sars = rg._SubmittedAlignedReadsDataFromReadGroup_in
        surs = rg._SubmittedUnalignedReadsDataFromReadGroup_in
        for sar in sars:
            sar = sar.src
            rgData.loc[i] = [project, rg.node_id, rg.experiment_name, rg.sequencing_center, rg.read_group_name, rg.is_paired_end, aliquot.node_id, aliquot.submitter_id, 'SAR', sar.node_id]
            i = i + 1
        for sur in surs:
            sur = sur.src
            rgData.loc[i] = [project, rg.node_id, rg.experiment_name, rg.sequencing_center, rg.read_group_name, rg.is_paired_end, aliquot.node_id, aliquot.submitter_id, 'SUR', sur.node_id]
            i = i + 1
rgData.to_csv('tcga.wgs.rg.tsv', sep='\t', encoding='utf-8')




# Get TARGET SARs and SURs
srData = pd.DataFrame(columns =  ['sr_type', \
                                  'sr', \
                                  'sr.batch_id', \
                                  'sr.file_name', \
                                  'sr.data_format', \
                                  'sr.file_state', \
                                  'sr.state', \
                                  'sr.submitter_id', \
                                  'sr.md5sum', \
                                  'sr.file_size', \
                                  'cghub_id', \
                                  'align_wf'])
i = 0
for project in tcga:
    sars = g.nodes(SubmittedAlignedReads).props(experimental_strategy='WGS', project_id = project).all()
    for sar in sars:
        file = g.nodes(File).get(sar.submitter_id)
        cghub_id = None if file is None else file.submitter_id
        wfs = sar._AlignmentWorkflowPerformedOnSubmittedAlignedReads_in
        if len(wfs) == 0: wfs = [None]
        for wf in wfs:
            wf = None if wf is None else wf.src
            srData.loc[i] = ['SAR', sar.node_id, sar.batch_id, sar.file_name, sar.data_format, sar.file_state, sar.state, sar.submitter_id, sar.md5sum, sar.file_size, cghub_id, None if wf is None else wf.node_id]
            i = i + 1
    surs = g.nodes(SubmittedUnalignedReads).props(experimental_strategy='WGS', project_id = project).all()    
    for sur in surs:
        file = g.nodes(File).get(sur.submitter_id)
        cghub_id = None if file is None else file.submitter_id
        wfs = sur._AlignmentWorkflowPerformedOnSubmittedUnalignedReads_in
        if len(wfs) == 0: wfs = [None]
        for wf in wfs:
            wf = None if wf is None else wf.src
            srData.loc[i] = ['SUR', sur.node_id, sur.batch_id, sur.file_name, sur.data_format, sur.file_state, sur.state, sur.submitter_id, sur.md5sum, sur.file_size, cghub_id, None if wf is None else wf.node_id]
            i = i + 1
srData.to_csv('tcga.wgs.sr.tsv', sep='\t', encoding='utf-8')


# Get TARGET ARs
arData = pd.DataFrame(columns =  ['ar', \
                                  'ar.batch_id', \
                                  'ar.file_name', \
                                  'ar.file_state', \
                                  'ar.state', \
                                  'ar.submitter_id', \
                                  'ar.md5sum', \
                                  'ar.file_size', \
                                  'ar.num.align_wf', \
                                  'align_wf'])
i = 0
for project in tcga:
    ars = g.nodes(AlignedReads).props(experimental_strategy='WGS', project_id = project).all()
    for ar in ars:
        wfs = ar._AlignedReadsDataFromAlignmentWorkflow_out
        wf = wfs[0].dst
        arData.loc[i] = [ ar.node_id, ar.batch_id, ar.file_name, sur.file_state, sur.state, ar.submitter_id, ar.md5sum, ar.file_size, len(wfs), wf.node_id]
        i = i + 1
arData.to_csv('tcga.wgs.ar.tsv', sep='\t', encoding='utf-8')





