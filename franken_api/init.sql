CREATE TYPE cur_user_status AS ENUM ('0','1', '2', '-1');
CREATE TYPE cur_user_status_desc AS ENUM ('0 -New user','1 - Active User', '2 - Suspended User', '-1 Delete User');

CREATE TYPE cur_role_status AS ENUM ('0', '1', '-1', '2');
CREATE TYPE cur_role_status_desc AS ENUM('0 - Normal User', '1 - Admin', '2 Course User', '-1 Super Admin');

CREATE TABLE cur_users_t (
    u_id SERIAL PRIMARY KEY, 
    first_name VARCHAR(150), 
    last_name VARCHAR(150), 
    email_id VARCHAR(150), 
    pwd Text, 
    user_status cur_user_status  default '0',
    role_id cur_role_status  default '0', 
    project_access text,
    created_on timestamptz
);

CREATE TYPE proj_status AS ENUM ('0', '1');
CREATE TYPE proj_status_desc AS ENUM('0 - Disable', '1 - Enable');

CREATE TABLE cur_projects_t (
    p_id SERIAL PRIMARY KEY, 
    project_name VARCHAR(150), 
    prefix_name VARCHAR(150), 
    nfs_path Text, 
    proj_status proj_status  default '0',
    mtbp_json proj_status NOT NULL DEFAULT '0',
    mtbp_report proj_status NOT NULL DEFAULT '0',
    pdf_report proj_status NOT NULL DEFAULT '0',
    sort_order integer NOT NULL DEFAULT 0,
    created_on timestamptz
);

/* INSERT INTO public.cur_projects_t( p_id, project_name, prefix_name, nfs_path, proj_status, created_on) 	VALUES (DEFAULT, 'IPCM', 'iPCM', '/sdata/IPCM/autoseq-output', '1', NOW()); */


==================================

CREATE SEQUENCE table_svs_id_seq  AS integer  START WITH 1  INCREMENT BY 1 NO MINVALUE  NO MAXVALUE  CACHE 1;

CREATE TABLE table_svs (
    id integer NOT NULL DEFAULT nextval('table_svs_id_seq'::regclass),
    "PROJECT_ID" character varying NOT NULL,
    "SDID" character varying NOT NULL,
    "CAPTURE_ID" character varying NOT NULL,
    "CHROM_A" character varying,
    "START_A" character varying,
    "END_A" character varying,
    "CHROM_B" character varying,
    "START_B" character varying,
    "END_B" character varying,
    "SVTYPE" character varying,
    "SV_LENGTH" character varying,
    "SUPPORT_READS" character varying,
    "TOOL" character varying,
    "SAMPLE" character varying,
    "GENE_A" character varying,
    "IN_DESIGN_A" character varying,
    "GENE_B" character varying,
    "IN_DESIGN_B" character varying,
    "GENE_A-GENE_B-sorted" character varying,
    "CALL" character varying,
    "TYPE" character varying,
    "SECONDHIT" character varying,
    "COMMENT" character varying,
    "ASSESSMENT" character varying,
    "CLONALITY" character varying,
    "CONSEQUENCE" character varying,
    "FUNCTIONAL_TYPE" character varying,
    "VARIANT_STRING" text,
    include_variant_report_pdf character varying(150)
);

ALTER SEQUENCE table_svs_id_seq OWNED BY table_svs.id;

==============================================================


CREATE SEQUENCE table_igv_somatic_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


CREATE TABLE public.table_igv_somatic (
    id integer NOT NULL DEFAULT nextval('table_igv_somatic_id_seq'::regclass),
    "PROJECT_ID" character varying NOT NULL,
    "SDID" character varying NOT NULL,
    "CAPTURE_ID" character varying NOT NULL,
    "CHROM" character varying,
    "START" character varying,
    "END" character varying,
    "REF" character varying,
    "ALT" character varying,
    "CALL" character varying,
    "TAG" character varying,
    "NOTES" character varying,
    "GENE" character varying,
    "IMPACT" character varying,
    "CONSEQUENCE" character varying,
    "HGVSp" character varying,
    "T_DP" character varying,
    "T_ALT" character varying,
    "T_VAF" character varying,
    "N_DP" character varying,
    "N_ALT" character varying,
    "N_VAF" character varying,
    "CLIN_SIG" character varying,
    "gnomAD" character varying,
    "BRCAEx" character varying,
    "OncoKB" character varying,
    "ASSESSMENT" character varying,
    "CLONALITY" character varying,
    "SECONDHIT" character varying,
    purecn_probability character varying,
    purecn_status character varying,
    purecn_tot_copies character varying,
    include_variant_report_pdf character varying(150)
);


ALTER SEQUENCE table_igv_somatic_id_seq OWNED BY table_igv_somatic.id;

==========================================================

CREATE SEQUENCE table_igv_germline_id_seq AS integer START WITH 1 INCREMENT BY 1 NO MINVALUE NO MAXVALUE CACHE 1;


CREATE TABLE table_igv_germline (
    id integer NOT NULL DEFAULT nextval('table_igv_germline_id_seq'::regclass),
    "PROJECT_ID" character varying NOT NULL,
    "SDID" character varying NOT NULL,
    "CAPTURE_ID" character varying NOT NULL,
    "CHROM" character varying,
    "START" character varying,
    "END" character varying,
    "REF" character varying,
    "ALT" character varying,
    "CALL" character varying,
    "TAG" character varying,
    "NOTES" character varying,
    "GENE" character varying,
    "IMPACT" character varying,
    "CONSEQUENCE" character varying,
    "HGVSp" character varying,
    "N_DP" character varying,
    "N_ALT" character varying,
    "N_VAF" character varying,
    "CLIN_SIG" character varying,
    "gnomAD" character varying,
    "BRCAEx" character varying,
    "OncoKB" character varying,
    purecn_probability character varying,
    purecn_status character varying,
    purecn_tot_copies character varying,
    include_variant_report_pdf character varying(150)
);

ALTER SEQUENCE table_igv_germline_id_seq OWNED BY table_igv_germline.id;

===========================================

CREATE SEQUENCE genomic_profile_summary_id_seq AS integer START WITH 1 INCREMENT BY 1 NO MINVALUE NO MAXVALUE CACHE 1;

CREATE TABLE genomic_profile_summary (
    id integer NOT NULL DEFAULT nextval('genomic_profile_summary_id_seq'::regclass),
    project_name text,
    sample_id text,
    capture_id text,
    study_code character varying(100),
    study_site character varying(250),
    dob character varying(100),
    disease character varying(250),
    specimen_assay json DEFAULT '{}'::json NOT NULL,
    ctdna_param text,
    ctdna_method character varying(250),
    genome_wide json DEFAULT '{}'::json NOT NULL,
    structural_variants json DEFAULT '{}'::json NOT NULL,
    somatic_mutations json DEFAULT '{}'::json NOT NULL,
    germline_alterations json DEFAULT '{}'::json NOT NULL,
    cnvs json DEFAULT '{}'::json NOT NULL,
    summary_txt text
);

ALTER SEQUENCE genomic_profile_summary_id_seq OWNED BY genomic_profile_summary.id;
