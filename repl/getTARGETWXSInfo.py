import pandas as pd

projects = ['TARGET-NBL','TARGET-WT','TARGET-AML']

# from collections import Counter as table


# Get SARs
sarData = pd.DataFrame(columns = ['project', \
                                  'sar', \
                                  'sar.batch_id', \
                                  'sar.file_state', \
                                  'sar.state', \
                                  'sar.submitter_id', \
                                  'sar.md5sum', \
                                  'sar.cghub_id', \
                                  'aliquot', \
                                  'aliquot.submitter_id', \
                                  'aliquot.selected_normal_wxs', \
                                  'aliquot.no_matched_normal_wxs', \
                                  'sample', \
                                  'sample.is_ffpe', \
                                  'sample.sample_type', \
                                  'sample.sample_type_id', \
                                  'sample.tissue_type', \
                                  'case', \
                                  'case.submitter_id'])
i = 0
for project in projects:
    sars = g.nodes(SubmittedAlignedReads).props(experimental_strategy='WXS', project_id = project).all()
    for sar in sars:
        rg0 = sar._SubmittedAlignedReadsDataFromReadGroup_out[0].dst
        aliquot = rg0._ReadGroupDerivedFromAliquot_out[0].dst
        sample = aliquot._AliquotDerivedFromSample_out[0].dst
        case = sample._SampleDerivedFromCase_out[0].dst
        file = g.nodes(File).get(sar.submitter_id)
        cghub_id = None if file is None else file.submitter_id
        sarData.loc[i] = [project, \
                         sar.node_id, \
                         sar.batch_id, \
                         sar.file_state, \
                         sar.state, \
                         sar.submitter_id, \
                         sar.md5sum, \
                         cghub_id, \
                         aliquot.node_id, \
                         aliquot.submitter_id, \
                         aliquot.selected_normal_wxs, \
                         aliquot.no_matched_normal_wxs, \
                         sample.node_id, \
                         sample.is_ffpe, \
                         sample.sample_type, \
                         sample.sample_type_id, \
                         sample.tissue_type, \
                         case.node_id, \
                         case.submitter_id]
        i = i + 1
sarData.to_csv('target.wxs.sar.tsv', sep='\t', encoding='utf-8')


# Get TCGA ARs
arData = pd.DataFrame(columns = ['ar', \
                                  'ar.batch_id', \
                                  'ar.file_state', \
                                  'ar.state', \
                                  'ar.submitter_id', \
                                  'ar.md5sum', \
                                  'ar.index', \
                                  'sar', \
                                  'coclean_wf'])
i = 0
for project in projects:
    ars = g.nodes(AlignedReads).props(experimental_strategy='WXS', project_id = project).all()
    for ar in ars:
        arIndex = ar._AlignedReadsIndexDerivedFromAlignedReads_in[0].src
        sar = ar._AlignedReadsMatchedToSubmittedAlignedReads_out[0].dst.node_id if len(ar._AlignedReadsMatchedToSubmittedAlignedReads_out) > 0 else None
        cocleanWf = ar._AlignedReadsDataFromAlignmentCocleaningWorkflow_out[0].dst
        arData.loc[i] = [ar.node_id, \
                         ar.batch_id, \
                         ar.file_state, \
                         ar.state, \
                         ar.submitter_id, \
                         ar.md5sum, \
                         arIndex.node_id, \
                         sar, \
                         cocleanWf.node_id]
        i = i + 1
arData.to_csv('target.wxs.ar.tsv', sep='\t', encoding='utf-8')


# Get TCGA Raw VCFs
ssmData = pd.DataFrame(columns = ['vcf_wf', \
                                  'vcf_wf.type', \
                                  'vcf_wf.version', \
                                  'vcf', \
                                  'vcf.batch_id', \
                                  'vcf.file_state', \
                                  'vcf.state', \
                                  'vcf.submitter_id', \
                                  'vcf.md5sum', \
                                  'vcf_index', \
                                  'ar1', \
                                  'ar2'])
i = 0
for project in projects:
    ssms = g.nodes(SimpleSomaticMutation).props(experimental_strategy='WXS', project_id = project).all()
    for ssm in ssms:
        ssmIndex = ssm._SomaticMutationIndexDerivedFromSimpleSomaticMutation_in[0].src.node_id if len(ssm._SomaticMutationIndexDerivedFromSimpleSomaticMutation_in) > 0 else None
        ssmWf = ssm._SimpleSomaticMutationDataFromSomaticMutationCallingWorkflow_out[0].dst
        ar1 = ssmWf._SomaticMutationCallingWorkflowPerformedOnAlignedReads_out[0].dst
        ar2 = ssmWf._SomaticMutationCallingWorkflowPerformedOnAlignedReads_out[1].dst
        ssmData.loc[i] = [ssmWf.node_id, \
                          ssmWf.workflow_type, \
                          ssmWf.workflow_version, \
                          ssm.node_id, \
                          ssm.batch_id, \
                          ssm.file_state, \
                          ssm.state, \
                          ssm.submitter_id, \
                          ssm.md5sum, \
                          ssmIndex, \
                          ar1.node_id, \
                          ar2.node_id]
        i = i + 1
ssmData.to_csv('target.wxs.vcf.tsv', sep='\t', encoding='utf-8')


