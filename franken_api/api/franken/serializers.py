from flask_restx import fields
from franken_api.api.restplus import api

status_result = api.model('server status', {'server_status': fields.Boolean(default=False, required=True)})

list = fields.String()
dropdownlist = api.model('Dropdown list for sample ids ',{'sidis': fields.List(list), 'status':fields.Boolean(required=True)})
dropdownlist_capture = api.model('Dropdown list for capture ids ',{'sample_capture': fields.List(list), 'status':fields.Boolean(required=True)})

ploturl_list = api.model('List the rul for static franken plots', {'image_url': fields.List(list), 'status':fields.Boolean(required=True)})

probio_ref_data = api.model('ProbioBloodReferral', {
    'crid': fields.String(description='referral id for reach record'),
    'pnr': fields.String(description='personal number of patient'),
    'rid': fields.String(description='referral id'),
    'datum': fields.String(description='date and of entry'),
    'tid': fields.String(description='tid of sample'),
    'sign': fields.String(description='binary stautus 1 or 0 '),
    'countyletter': fields.String(description='hospital code'),
    'new': fields.String(description='new'),
    'progression': fields.String(description='state of progression'),
    'follow_up': fields.String(description='status to follow_up'),
    'cf_dna1': fields.String(description='proof read 1'),
    'cf_dna2': fields.String(description='proof read 2'),
    'cf_dna3': fields.String(description='proof read 3'),
    'blod': fields.String(description='proof read 3'),
    'kommentar': fields.String(description='comments of each sample'),
    'filnamn': fields.String(description='path to sample report in pdf '),
    'remisstyp': fields.String(description='Remise type'),
    'studieid': fields.String(description='Study ID'),
    'hormonkänslig': fields.String(description='hormonkänslig'),
    'kastrationsresistent': fields.String(description='kastrationsresistent'),
    'cdk': fields.String(description='cdk'),
    'endofstudy': fields.String(description='end of study')
})

psff_ref_data = api.model('PsffBloodReferral', {
    'crid': fields.String(description='referral id for reach record'),
    'rid': fields.String(description='referral id'),
    'datum': fields.String(description='date and of entry'),
    'tid': fields.String(description='tid of sample'),
    'sign': fields.String(description='binary stautus 1 or 0 '),
    'blood1': fields.String(description='proof read 1'),
    'blood2': fields.String(description='proof read 2'),
    'blood3': fields.String(description='proof read 3'),
    'blood4': fields.String(description='proof read 4'),
    'comment': fields.String(description='comments on case'),
    'filnamn': fields.String(description='path to sample report in pdf '),
    'cdk': fields.String(description='cdk id')
})

header = api.model('Probio header', {'key':fields.String(description='Column name') , 'title': fields.String(description='Column display name')})
#probio_ref_data_list = api.model('Referral DB Data', { 'status':fields.Boolean(required=True), 'data': fields.List(fields.Nested(probio_ref_data)), 'header': fields.List(fields.Nested(header)), 'error':fields.String() }
probio_ref_data_list = api.model('Referral DB Data', { 'status':fields.Boolean(required=True), 'data': fields.List(fields.Nested(probio_ref_data)), 'header': fields.List(fields.Nested(header)), 'error': fields.String()})


header_psff = api.model('Psff header', {'key':fields.String(description='Column name') , 'title': fields.String(description='Column display name')})
psff_ref_data_list = api.model('Referral DB Data', { 'status':fields.Boolean(required=True), 'data': fields.List(fields.Nested(psff_ref_data)), 'header': fields.List(fields.Nested(header_psff)), 'error': fields.String()})


referral_db_out =  api.model('Referral DB update', { 'status':fields.Boolean(required=True), 'error': fields.String()})

curation_germline = api.model('IGV Germline', {
        'id': fields.String(description='Project ID'),
        'PROJECT_ID' : fields.String(description='Project ID'),
        'SDID' : fields.String(description='Sample ID'),
        'CAPTURE_ID' : fields.String(description='Capture ID'),
        'CHROM': fields.String(description='Chromosome '),
        'START':  fields.String(description='Start position'),
        'END': fields.String(description='End'),
        'REF': fields.String(description=''),
        'ALT': fields.String(description=''),
        'CALL':fields.String(description=''),
        'TAG' : fields.String(description=''),
        'ASSESSMENT': fields.String(description=''),
        'NOTES': fields.String(description=''),
        'GENE': fields.String(description=''),
        'IMPACT': fields.String(description=''),
        'CONSEQUENCE' : fields.String(description=''),
        'HGVSp': fields.String(description=''),
        'N_DP' : fields.String(description=''),
        'N_ALT' : fields.String(description=''),
        'N_VAF' : fields.String(description=''),
        'CLIN_SIG' : fields.String(description=''),
        'gnomAD': fields.String(description=''),
        'BRCAEx' : fields.String(description=''),
        'OncoKB' : fields.String(description=''),
        'purecn_probability' : fields.String(description=''),
        'purecn_status' : fields.String(description=''),
        'purecn_tot_copies' : fields.String(description=''),
        'include_variant_report_pdf' : fields.String(description=''),
        'user_name' : fields.String(description='')
})
curation_somatic = api.model('IGV Somatic', {
        'id': fields.String(description='Project ID'),
        'PROJECT_ID' : fields.String(description='Project ID'),
        'SDID' : fields.String(description='Sample ID'),
        'CAPTURE_ID' : fields.String(description='Capture ID'),
        'CHROM': fields.String(description='Chromosome '),
        'START':  fields.String(description='Start position'),
        'END': fields.String(description='End'),
        'REF': fields.String(description=''),
        'ALT': fields.String(description=''),
        'CALL': fields.String(description=''),
        'TAG': fields.String(description=''),
        'ASSESSMENT': fields.String(description=''),
        'CLONALITY': fields.String(description=''),
        'NOTES': fields.String(description=''),
        'GENE': fields.String(description=''),
        'IMPACT': fields.String(description=''),
        'CONSEQUENCE' : fields.String(description=''),
        'HGVSp': fields.String(description=''),
        'T_DP': fields.String(description=''),
        'T_ALT': fields.String(description=''),
        'T_VAF' : fields.String(description=''),
        'N_DP' : fields.String(description=''),
        'N_ALT' : fields.String(description=''),
        'N_VAF' : fields.String(description=''),
        'CLIN_SIG' : fields.String(description=''),
        'gnomAD': fields.String(description=''),
        'BRCAEx' : fields.String(description=''),
        'OncoKB' : fields.String(description=''),
        'purecn_probability' : fields.String(description=''),
        'purecn_status' : fields.String(description=''),
        'purecn_tot_copies' : fields.String(description=''),
        'include_variant_report_pdf' : fields.String(description=''),
        'user_name' : fields.String(description='')
})

curation_svs = api.model('SVS', {
        'id': fields.String(description=' ID'),
        'PROJECT_ID' : fields.String(description='Project ID'),
        'SDID' : fields.String(description='Sample ID'),
        'CAPTURE_ID' : fields.String(description='Capture ID'),
        'CHROM_A': fields.String(description='Chromosome '),
        'START_A':  fields.String(description='Start position'),
        'END_A': fields.String(description='End'),
        'CHROM_B': fields.String(description=''),
        'START_B': fields.String(description=''),
        'END_B':fields.String(description=''),
        'SVTYPE' : fields.String(description=''),
        'SV_LENGTH': fields.String(description=''),
        'SUPPORT_READS': fields.String(description=''),
        'TOOL': fields.String(description=''),
        'SAMPLE' : fields.String(description=''),
        'GENE_A': fields.String(description=''),
        'IN_DESIGN_A': fields.String(description=''),
        'GENE_B': fields.String(description=''),
        'IN_DESIGN_B' : fields.String(description=''),
        'GENE_A-GENE_B-sorted' : fields.String(description=''),
        'CALL' : fields.String(description=''),
        'TYPE' : fields.String(description=''),
        'SECONDHIT' : fields.String(description=''),
        'COMMENT': fields.String(description=''),
        'ASSESSMENT': fields.String(description=''),
        'CLONALITY': fields.String(description=''),
        'CONSEQUENCE': fields.String(description=''),
        'FUNCTIONAL_TYPE': fields.String(description=''),
        'VARIANT_STRING': fields.String(description=''),
        'include_variant_report_pdf' : fields.String(description=''),
        'user_name' : fields.String(description='')
})

curation_hotspot = api.model('IGV Hotspot', {
        'id': fields.String(description=''),
        'gene': fields.String(description=''),
        'aapos': fields.String(description=''),
        'nmut': fields.String(description=''),
        'protmut': fields.String(description=''),
        'prot2mut': fields.String(description=''),
        'dnamut': fields.String(description=''),
        'canmut': fields.String(description=''),
        'conseqmut': fields.String(description=''),
        'transcript': fields.String(description=''),
        'dn_ds': fields.String(description=''), 
        'community_notes': fields.String(description='')
})

curation_cancer_hotspot = api.model('IGV Cancer Hotspot', {
        'h_id': fields.String(description=''),
        'gene': fields.String(description=''),
        'hgvsp': fields.String(description=''),
        'amino_acid_position': fields.String(description=''),
        'start_aa': fields.String(description=''),
        'end_aa': fields.String(description='')
})

# curation_cancer_hotspot = api.model('IGV Cancer Hotspot', {
#         'id': fields.String(description=''),
#         'gene': fields.String(description=''),
#         'residue': fields.String(description=''),
#         'res_type': fields.String(description=''),        
#         'variants': fields.String(description=''),
#         'variant_arr': fields.String(description='')
# })


curation_warmspot = api.model('IGV Warmspot', {
        'id': fields.String(description=''),
        'gene': fields.String(description=''),
        'aapos': fields.String(description=''),
        'nmut': fields.String(description=''),
        'protmut': fields.String(description=''),
        'prot2mut': fields.String(description=''),
        'dnamut': fields.String(description=''),
        'canmut': fields.String(description=''),
        'conseqmut': fields.String(description=''),
        'transcript': fields.String(description=''),
        'dn_ds': fields.String(description=''), 
        'community_notes': fields.String(description='')
})

curation_psff_profile = api.model('PSFF Profile', {
        'id': fields.String(description=''),
        'sample_id': fields.String(description=''),
        'capture_id': fields.String(description=''),
        'basic_qc': fields.String(description=''),
        'franken_plot': fields.String(description=''),
        'ploidy': fields.String(description=''),
        'ctdna_fraction': fields.String(description=''),
        'ctdna_param': fields.String(description=''),
        'ctdna_method': fields.String(description=''),
        'study_code': fields.String(description=''),
        'study_site': fields.String(description=''),
        'somatic_mutations': fields.String(description=''), 
        'germline_alterations': fields.String(description=''),
        'structural_variants': fields.String(description=''),
        'cnvs': fields.String(description=''),
        'comment_info': fields.String(description=''),
})

curation_probio_profile = api.model('PROBIO Profile', {
        'id': fields.String(description=''),
        'sample_id': fields.String(description=''),
        'capture_id': fields.String(description=''),
        'basic_qc': fields.String(description=''),
        'franken_plot': fields.String(description=''),
        'ploidy': fields.String(description=''),
        'ctdna_fraction': fields.String(description=''),
        'ctdna_param': fields.String(description=''),
        'ctdna_method': fields.String(description=''),
        'ctdna_category': fields.String(description=''),
        'study_code': fields.String(description=''),
        'study_site': fields.String(description=''),
        'somatic_mutations': fields.String(description=''), 
        'germline_alterations': fields.String(description=''),
        'structural_variants': fields.String(description=''),
        'cnvs': fields.String(description=''),
        'comment_info': fields.String(description=''),
})

curation_genomic_profile = api.model('Genomic Profile Summary', {
        'id': fields.String(description=''),
        'project_name': fields.String(description=''),
        'sample_id': fields.String(description=''),
        'capture_id': fields.String(description=''),
        'study_code': fields.String(description=''),
        'study_site': fields.String(description=''),
        'dob': fields.String(description=''),
        'disease': fields.String(description=''),
        'specimen_assay': fields.String(description=''),
        'ctdna_param': fields.String(description=''),
        'ctdna_method': fields.String(description=''),
        'genome_wide': fields.String(description=''),
        'somatic_mutations': fields.String(description=''), 
        'germline_alterations': fields.String(description=''),
        'structural_variants': fields.String(description=''),
        'cnvs': fields.String(description=''),
        'summary_txt': fields.String(description='')
})


header_curation = api.model('Curation Header', {'key':fields.String(description='Column name') ,
                                        'title': fields.String(description='Column display name')})

germline_data_list = api.model('Curation Germline DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_germline)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

somatic_data_list = api.model('Curation Somatic DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_somatic)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

svs_data_list = api.model('Curation SVS DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_svs)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

hotspot_data_list = api.model('Curation Hotspot DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_hotspot)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

cancer_hotspot_data_list = api.model('Curation Cancer Hotspot DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_cancer_hotspot)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

warmspot_data_list = api.model('Curation Warmspot DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_warmspot)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

psff_profile_data_list = api.model('Curation PSFF Profile DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_psff_profile)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})                                                     

probio_profile_data_list = api.model('Curation PROBIO Profile DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_probio_profile)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})

genomic_profile_data_list = api.model('Curation Genomic Profile DB Data', { 'status':fields.Boolean(required=True),
                                                     'data': fields.List(fields.Nested(curation_genomic_profile)),
                                                     'header': fields.List(fields.Nested(header_curation)),
                                                     'error': fields.String()})                                                     
