from flask_restplus import reqparse

project_arguments = reqparse.RequestParser()
project_arguments.add_argument('project_name', choices=('PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', 'msHSPC', 'RESBIO'), required=True,
                              help="valid Project names: 'PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', msHSPC', 'RESBIO'")

capture_arguments = reqparse.RequestParser()
capture_arguments.add_argument('project_name', choices=('PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', 'msHSPC', 'RESBIO'), required=True,
                              help="valid Project names: 'PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', 'msHSPC', 'RESBIO'")
capture_arguments.add_argument('sdid', type=str, required=True,  help='sdid example : P-00360714')


common_arguments = reqparse.RequestParser()
common_arguments.add_argument('project_name', choices=('PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', 'msHSPC', 'RESBIO'), required=True,
                              help="valid Project names: 'PROBIO', 'PSFF', 'LARS_DUCTAL', 'AZ_RINGTRIAL', 'msHSPC', 'RESBIO'")
common_arguments.add_argument('sdid', type=str, required=True,  help='sdid example : P-00360714')
common_arguments.add_argument('capture_id', type=str, required=True,  help='capture id')


search_arguments = common_arguments.copy()
search_arguments.add_argument('pname', choices=('normal', 'variant_allelic', 'plot_AR', 'snpratio', 'scatter_plot', 'snpratio_density', 'scatter_density_plot' ), required=True,
                              help="valid plot names: 'normal', 'variant_allelic', 'plot_AR', 'snpratio', 'scatter_plot','scatter_density_plot', 'snpratio_density.json'")


ploturls_arguments = common_arguments.copy()

purecn_arguments = common_arguments.copy()

purecn_max_val_arguments = common_arguments.copy()

purecn_max_val_arguments.add_argument('variant_type', choices=('somatic', 'germline'), required=True,
                              help="valid Variant Type: 'somatic', 'germline'")


staticplot_arguments = common_arguments.copy()
staticplot_arguments.add_argument('imagename', type=str, required=True,  help='static franken plot name')


igv_arguments =  common_arguments.copy()
igv_arguments.add_argument('filename', type=str, required=True,  help='file required for igv track creation')


table_svs_arguments =  common_arguments.copy()
table_svs_arguments.add_argument('header', choices=('true', 'false'), default='true', required=True,  help='Boolen True or False : True returns only column header')


table_igvnav_arguments = table_svs_arguments.copy()

table_qc_arguments = table_svs_arguments.copy()

table_cnv_arguments = table_svs_arguments.copy()

igv_save_file_arguments = reqparse.RequestParser()
igv_save_file_arguments.add_argument('file_name',  required=True, help="igvnav-input.txt as File Name  ")
igv_save_file_arguments.add_argument('data',  type=str, required=True, help="Jsondata")


pdf_arguments = common_arguments.copy()

referral_update_arguments = reqparse.RequestParser()
referral_update_arguments.add_argument('db_name', choices=('probio', 'psff'), required=True,
                              help="provide db name to updae the referral db for probio / psff")

curation_germline_arguments = reqparse.RequestParser()
curation_germline_arguments.add_argument('PROJECT_ID', type=str, required=True, help="")
curation_germline_arguments.add_argument('SDID' , type=str, required=True, help="")
curation_germline_arguments.add_argument('CAPTURE_ID' , type=str, required=True, help="")
curation_germline_arguments.add_argument('CHROM', type=str, required=True, help="")
curation_germline_arguments.add_argument('START', type=str, required=True, help="")
curation_germline_arguments.add_argument('END', type=str, required=True, help="")
curation_germline_arguments.add_argument('REF', type=str, required=True, help="")
curation_germline_arguments.add_argument('ALT', type=str, required=True, help="")
curation_germline_arguments.add_argument('CALL', type=str, required=True, help="")
curation_germline_arguments.add_argument('TAG' , type=str,  help="")
curation_germline_arguments.add_argument('NOTES', type=str,  help="")
curation_germline_arguments.add_argument('GENE', type=str,  help="")
curation_germline_arguments.add_argument('IMPACT', type=str,  help="")
curation_germline_arguments.add_argument('CONSEQUENCE', type=str,  help="")
curation_germline_arguments.add_argument('HGVSp', type=str,  help="")
curation_germline_arguments.add_argument('N_DP', type=str,  help="")
curation_germline_arguments.add_argument('N_ALT', type=str,  help="")
curation_germline_arguments.add_argument('N_VAF', type=str,  help="")
curation_germline_arguments.add_argument('CLIN_SIG', type=str,  help="")
curation_germline_arguments.add_argument('gnomAD', type=str,  help="")
curation_germline_arguments.add_argument('BRCAEx' , type=str,  help="")
curation_germline_arguments.add_argument('OncoKB' , type=str,  help="")
curation_germline_arguments.add_argument('purecn_probability' , type=str,  help="")
curation_germline_arguments.add_argument('purecn_status' , type=str,  help="")
curation_germline_arguments.add_argument('purecn_tot_copies' , type=str,  help="")


curation_somatic_arguments = reqparse.RequestParser()
curation_somatic_arguments.add_argument('PROJECT_ID', type=str, required=True, help="")
curation_somatic_arguments.add_argument('SDID' , type=str, required=True, help="")
curation_somatic_arguments.add_argument('CAPTURE_ID' , type=str, required=True, help="")
curation_somatic_arguments.add_argument('CHROM', type=str, required=True, help="")
curation_somatic_arguments.add_argument('START', type=str, required=True, help="")
curation_somatic_arguments.add_argument('END', type=str, required=True, help="")
curation_somatic_arguments.add_argument('REF', type=str, required=True, help="")
curation_somatic_arguments.add_argument('ALT', type=str, required=True, help="")
curation_somatic_arguments.add_argument('CALL', type=str,  help="")
curation_somatic_arguments.add_argument('TAG' , type=str,  help="")
curation_somatic_arguments.add_argument('NOTES', type=str,  help="")
curation_somatic_arguments.add_argument('GENE', type=str,  help="")
curation_somatic_arguments.add_argument('IMPACT', type=str,  help="")
curation_somatic_arguments.add_argument('CONSEQUENCE', type=str,  help="")
curation_somatic_arguments.add_argument('HGVSp', type=str,  help="")
curation_somatic_arguments.add_argument('T_DP', type=str,  help="")
curation_somatic_arguments.add_argument('T_ALT', type=str,  help="")
curation_somatic_arguments.add_argument('T_VAF', type=str,  help="")
curation_somatic_arguments.add_argument('N_DP', type=str,  help="")
curation_somatic_arguments.add_argument('N_ALT', type=str,  help="")
curation_somatic_arguments.add_argument('N_VAF', type=str,  help="")
curation_somatic_arguments.add_argument('CLIN_SIG', type=str,  help="")
curation_somatic_arguments.add_argument('gnomAD', type=str,  help="")
curation_somatic_arguments.add_argument('BRCAEx' , type=str,  help="")
curation_somatic_arguments.add_argument('OncoKB' , type=str,  help="")
curation_somatic_arguments.add_argument('ASSESSMENT' , type=str,  help="")
curation_somatic_arguments.add_argument('CLONALITY' , type=str,  help="")
curation_somatic_arguments.add_argument('purecn_probability' , type=str,  help="")
curation_somatic_arguments.add_argument('purecn_status' , type=str,  help="")
curation_somatic_arguments.add_argument('purecn_tot_copies' , type=str,  help="")

curation_svs_arguments = reqparse.RequestParser()
curation_svs_arguments.add_argument('PROJECT_ID', type=str, required=True, help="")
curation_svs_arguments.add_argument('SDID' , type=str, required=True, help="")
curation_svs_arguments.add_argument('CAPTURE_ID' , type=str, required=True, help="")
curation_svs_arguments.add_argument('CHROM_A', type=str, required=True, help="")
curation_svs_arguments.add_argument('START_A', type=str, required=True, help="")
curation_svs_arguments.add_argument('END_A', type=str, required=True, help="")
curation_svs_arguments.add_argument('CHROM_B', type=str, required=True, help="")
curation_svs_arguments.add_argument('START_B', type=str, required=True, help="")
curation_svs_arguments.add_argument('END_B', type=str, required=True, help="")
curation_svs_arguments.add_argument('SVTYPE' , type=str,  help="")
curation_svs_arguments.add_argument('SV_LENGTH', type=str,  help="")
curation_svs_arguments.add_argument('SUPPORT_READS', type=str,  help="")
curation_svs_arguments.add_argument('TOOL', type=str,  help="")
curation_svs_arguments.add_argument('SAMPLE', type=str,  help="")
curation_svs_arguments.add_argument('GENE_A', type=str,  help="")
curation_svs_arguments.add_argument('IN_DESIGN_A', type=str,  help="")
curation_svs_arguments.add_argument('GENE_B', type=str,  help="")
curation_svs_arguments.add_argument('IN_DESIGN_B', type=str,  help="")
curation_svs_arguments.add_argument('GENE_A-GENE_B-sorted', type=str,  help="")
curation_svs_arguments.add_argument('CALL', type=str,  help="")
curation_svs_arguments.add_argument('TYPE' , type=str,  help="")
curation_svs_arguments.add_argument('SECONDHIT' , type=str,  help="")
curation_svs_arguments.add_argument('COMMENT' , type=str,  help="")


curation_psff_profile_arguments = reqparse.RequestParser()
curation_psff_profile_arguments.add_argument('sample_id' , type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('capture_id' , type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('basic_qc', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('franken_plot', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('ploidy', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('ctdna_fraction', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('ctdna_param', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('ctdna_method', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('study_code', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('study_site' , type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('somatic_mutations' , type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('germline_alterations', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('structural_variants', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('cnvs', type=str, required=True, help="")
curation_psff_profile_arguments.add_argument('comment_info', type=str, required=True, help="")

curation_psff_profile_id_arguments = reqparse.RequestParser()
curation_psff_profile_id_arguments.add_argument('sample_id' , type=str, required=True, help="")


curation_probio_profile_arguments = reqparse.RequestParser()
curation_probio_profile_arguments.add_argument('sample_id' , type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('capture_id' , type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('basic_qc', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('franken_plot', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('ploidy', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('ctdna_fraction', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('ctdna_param', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('ctdna_method', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('ctdna_category', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('study_code', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('study_site' , type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('somatic_mutations' , type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('germline_alterations', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('structural_variants', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('cnvs', type=str, required=True, help="")
curation_probio_profile_arguments.add_argument('comment_info', type=str, required=True, help="")

curation_probio_profile_id_arguments = reqparse.RequestParser()
curation_probio_profile_id_arguments.add_argument('sample_id' , type=str, required=True, help="")