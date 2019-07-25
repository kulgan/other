import pandas as pd

target = ['TARGET-AML', 'TARGET-WT', 'TARGET-RT', 'TARGET-ALL-P1', 'TARGET-ALL-P2', 'TARGET-ALL-P3', 'TARGET-CCSK', 'TARGET-NBL']

# from collections import Counter as table

# Get TARGET Read Groups
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
for project in target:
    rgs = g.nodes(ReadGroup).props(library_strategy = 'RNA-Seq', project_id = project).all()
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
rgData.to_csv('target.rg.tsv', sep='\t', encoding='utf-8')




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
for project in target:
    sars = g.nodes(SubmittedAlignedReads).props(experimental_strategy='RNA-Seq', project_id = project).all()
    for sar in sars:
        file = g.nodes(File).get(sar.submitter_id)
        cghub_id = None if file is None else file.submitter_id
        wfs = sar._AlignmentWorkflowPerformedOnSubmittedAlignedReads_in
        if len(wfs) == 0: wfs = [None]
        for wf in wfs:
            wf = None if wf is None else wf.src
            srData.loc[i] = ['SAR', sar.node_id, sar.batch_id, sar.file_name, sar.data_format, sar.file_state, sar.state, sar.submitter_id, sar.md5sum, sar.file_size, cghub_id, None if wf is None else wf.node_id]
            i = i + 1
    surs = g.nodes(SubmittedUnalignedReads).props(experimental_strategy='RNA-Seq', project_id = project).all()    
    for sur in surs:
        file = g.nodes(File).get(sur.submitter_id)
        cghub_id = None if file is None else file.submitter_id
        wfs = sur._AlignmentWorkflowPerformedOnSubmittedUnalignedReads_in
        if len(wfs) == 0: wfs = [None]
        for wf in wfs:
            wf = None if wf is None else wf.src
            srData.loc[i] = ['SUR', sur.node_id, sur.batch_id, sur.file_name, sur.data_format, sur.file_state, sur.state, sur.submitter_id, sur.md5sum, sur.file_size, cghub_id, None if wf is None else wf.node_id]
            i = i + 1
srData.to_csv('target.sr.tsv', sep='\t', encoding='utf-8')


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
for project in target:
    ars = g.nodes(AlignedReads).props(experimental_strategy='RNA-Seq', project_id = project).all()
    for ar in ars:
        wfs = ar._AlignedReadsDataFromAlignmentWorkflow_out
        wf = wfs[0].dst
        arData.loc[i] = [ ar.node_id, ar.batch_id, ar.file_name, sur.file_state, sur.state, ar.submitter_id, ar.md5sum, ar.file_size, len(wfs), wf.node_id]
        i = i + 1
arData.to_csv('target.ar.tsv', sep='\t', encoding='utf-8')





