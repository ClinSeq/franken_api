from franken_api.database import db as pssql
from sqlalchemy.dialects.postgresql import JSON


class ProbioBloodReferral(pssql.Model):
    __tablename__ = "probio_bloodreferrals"
    crid = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    pnr  = pssql.Column(pssql.String)
    rid  = pssql.Column(pssql.String, nullable=False)
    datum = pssql.Column(pssql.Date, nullable=False)
    tid  = pssql.Column(pssql.String, nullable=False)
    sign = pssql.Column(pssql.Integer)
    countyletter  = pssql.Column(pssql.String)
    new = pssql.Column(pssql.String)
    progression = pssql.Column(pssql.String)
    follow_up    = pssql.Column(pssql.String)
    cf_dna1 = pssql.Column(pssql.String, nullable=False)
    cf_dna2 = pssql.Column(pssql.String, nullable=False)
    cf_dna3 = pssql.Column(pssql.String)
    blod = pssql.Column(pssql.String)
    kommentar = pssql.Column(pssql.String, nullable=False)
    filnamn = pssql.Column(pssql.String, nullable=False)
    remisstyp = pssql.Column(pssql.String)
    studieid = pssql.Column(pssql.String)
    hormonkänslig = pssql.Column(pssql.String)
    kastrationsresistent = pssql.Column(pssql.String)
    cdk = pssql.Column(pssql.String)
    endofstudy = pssql.Column(pssql.String)

    def __repr__(self): 
        return "<ProbioReferral (crid='%s', pnr='%s', rid='%s', date='%s', filename='%s')>" % (self.crid,self.pnr,self.rid,self.datum,self.filnamn)


class PSFFBloodReferral(pssql.Model):
    __tablename__ = "psff_bloodreferrals"
    crid = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    rid  = pssql.Column(pssql.String, nullable=False)
    datum = pssql.Column(pssql.Date, nullable=False)
    tid  = pssql.Column(pssql.String, nullable=False)
    sign = pssql.Column(pssql.Integer)
    blood1 = pssql.Column(pssql.String, nullable=False)
    blood2 = pssql.Column(pssql.String, nullable=False)
    blood3 = pssql.Column(pssql.String, nullable=False)
    blood4 = pssql.Column(pssql.String, nullable=False)
    comment = pssql.Column(pssql.String, nullable=False)
    filnamn = pssql.Column(pssql.String, nullable=False)
    cdk = pssql.Column(pssql.String, nullable=False)
    def __repr__(self):
        return "<PSFFReferral (cdk='%s', rid='%s', date='%s', filename='%s')>" % (self.cdk,self.rid,self.datum,self.filnamn)

class TableIgvGermline(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "table_igv_germline"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    PROJECT_ID = pssql.Column(pssql.String, nullable=False)
    SDID = pssql.Column(pssql.String, nullable=False)
    CAPTURE_ID = pssql.Column(pssql.String, nullable=False)
    CHROM = pssql.Column(pssql.String)
    START = pssql.Column(pssql.String)
    END = pssql.Column(pssql.String)
    REF = pssql.Column(pssql.String)
    ALT = pssql.Column(pssql.String)
    CALL = pssql.Column(pssql.String)
    TAG = pssql.Column(pssql.String)
    ASSESSMENT = pssql.Column(pssql.String)
    NOTES = pssql.Column(pssql.String)
    GENE  = pssql.Column(pssql.String)
    IMPACT = pssql.Column(pssql.String)
    CONSEQUENCE = pssql.Column(pssql.String)
    HGVSp = pssql.Column(pssql.String)
    N_DP = pssql.Column(pssql.String)
    N_ALT = pssql.Column(pssql.String)
    N_VAF = pssql.Column(pssql.String)
    CLIN_SIG = pssql.Column(pssql.String)
    gnomAD = pssql.Column(pssql.String)
    BRCAEx = pssql.Column(pssql.String)
    OncoKB = pssql.Column(pssql.String)
    purecn_probability = pssql.Column(pssql.String)
    purecn_status = pssql.Column(pssql.String)
    purecn_tot_copies = pssql.Column(pssql.String)
    include_variant_report_pdf = pssql.Column(pssql.String)
    user_name = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.PROJECT_ID = row_dict.get('PROJECT_ID', None)
        self.SDID = row_dict.get('SDID', None)
        self.CAPTURE_ID = row_dict.get('CAPTURE_ID', None)
        self.CHROM = row_dict.get('CHROM', None)
        self.START = row_dict.get('START', None)
        self.END = row_dict.get('END', None)
        self.REF = row_dict.get('REF', None)
        self.ALT = row_dict.get('ALT', None)
        self.CALL = row_dict.get('CALL', None)
        self.TAG = row_dict.get('TAG', None)
        self.ASSESSMENT = row_dict.get('ASSESSMENT', None)
        self.NOTES = row_dict.get('NOTES', None)
        self.GENE = row_dict.get('GENE', None)
        self.IMPACT = row_dict.get('IMPACT', None)
        self.CONSEQUENCE = row_dict.get('CONSEQUENCE', None)
        self.HGVSp = row_dict.get('HGVSp', None)
        self.N_DP = row_dict.get('N_DP', None)
        self.N_ALT = row_dict.get('N_ALT', None)
        self.N_VAF = row_dict.get('N_VAF', None)
        self.CLIN_SIG = row_dict.get('CLIN_SIG', None)
        self.gnomAD = row_dict.get('gnomAD', None)
        self.BRCAEx = row_dict.get('BRCAEx', None)
        self.OncoKB = row_dict.get('OncoKB', None)
        self.purecn_probability = row_dict.get('purecn_probability', None)
        self.purecn_status = row_dict.get('purecn_status', None)
        self.purecn_tot_copies = row_dict.get('purecn_tot_copies', None)
        self.include_variant_report_pdf = row_dict.get('include_variant_report_pdf', None)
        self.user_name = row_dict.get('user_name', None)

    def __repr__(self):
        return "<TableIgvGermline (id='%s', projectid='%s', sdid='%s', captureid='%s')>" % (self.id,self.PROJECT_ID,
                                                                                            self.SDID,self.CAPTURE_ID)

class TableIgvSomatic(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "table_igv_somatic"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    PROJECT_ID = pssql.Column(pssql.String, nullable=False)
    SDID = pssql.Column(pssql.String, nullable=False)
    CAPTURE_ID = pssql.Column(pssql.String, nullable=False)
    CHROM = pssql.Column(pssql.String)
    START = pssql.Column(pssql.String)
    END = pssql.Column(pssql.String)
    REF = pssql.Column(pssql.String)
    ALT = pssql.Column(pssql.String)
    CALL = pssql.Column(pssql.String)
    TAG = pssql.Column(pssql.String)
    ASSESSMENT = pssql.Column(pssql.String)
    NOTES = pssql.Column(pssql.String)
    GENE = pssql.Column(pssql.String)
    IMPACT = pssql.Column(pssql.String)
    CONSEQUENCE = pssql.Column(pssql.String)
    HGVSp = pssql.Column(pssql.String)
    T_DP = pssql.Column(pssql.String)
    T_ALT = pssql.Column(pssql.String)
    T_VAF = pssql.Column(pssql.String)
    N_DP = pssql.Column(pssql.String)
    N_ALT = pssql.Column(pssql.String)
    N_VAF = pssql.Column(pssql.String)
    CLIN_SIG = pssql.Column(pssql.String)
    gnomAD = pssql.Column(pssql.String)
    BRCAEx = pssql.Column(pssql.String)
    OncoKB = pssql.Column(pssql.String)
    purecn_probability = pssql.Column(pssql.String)
    purecn_status = pssql.Column(pssql.String)
    purecn_tot_copies = pssql.Column(pssql.String)
    include_variant_report_pdf = pssql.Column(pssql.String)
    user_name = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.PROJECT_ID = row_dict.get('PROJECT_ID', None)
        self.SDID = row_dict.get('SDID', None)
        self.CAPTURE_ID = row_dict.get('CAPTURE_ID', None)
        self.CHROM = row_dict.get('CHROM', None)
        self.START = row_dict.get('START', None)
        self.END = row_dict.get('END', None)
        self.REF = row_dict.get('REF', None)
        self.ALT = row_dict.get('ALT', None)
        self.CALL = row_dict.get('CALL', None)
        self.TAG = row_dict.get('TAG', None)
        self.ASSESSMENT = row_dict.get('ASSESSMENT', None)
        self.NOTES = row_dict.get('NOTES', None)
        self.GENE = row_dict.get('GENE', None)
        self.IMPACT = row_dict.get('IMPACT', None)
        self.CONSEQUENCE = row_dict.get('CONSEQUENCE', None)
        self.HGVSp = row_dict.get('HGVSp', None)
        self.T_DP = row_dict.get('T_DP', None)
        self.T_ALT = row_dict.get('T_ALT', None)
        self.T_VAF = row_dict.get('T_VAF', None)
        self.N_DP = row_dict.get('N_DP', None)
        self.N_ALT = row_dict.get('N_ALT', None)
        self.N_VAF = row_dict.get('N_VAF', None)
        self.CLIN_SIG = row_dict.get('CLIN_SIG', None)
        self.gnomAD = row_dict.get('gnomAD', None)
        self.BRCAEx = row_dict.get('BRCAEx', None)
        self.OncoKB = row_dict.get('OncoKB', None)
        self.purecn_probability = row_dict.get('purecn_probability', None)
        self.purecn_status = row_dict.get('purecn_status', None)
        self.purecn_tot_copies = row_dict.get('purecn_tot_copies', None)
        self.include_variant_report_pdf = row_dict.get('include_variant_report_pdf', None)
        self.user_name = row_dict.get('user_name', None)

    def __repr__(self):
        return "<TableIgvSomatic (id='%s', projectid='%s', sdid='%s', captureid='%s')>" % (self.id,self.PROJECT_ID,
                                                                                            self.SDID, self.CAPTURE_ID)

class TableSVS(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "table_svs"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    PROJECT_ID = pssql.Column(pssql.String, nullable=False)
    SDID = pssql.Column(pssql.String, nullable=False)
    CAPTURE_ID = pssql.Column(pssql.String, nullable=False)
    CHROM_A = pssql.Column(pssql.String)
    START_A = pssql.Column(pssql.String)
    END_A = pssql.Column(pssql.String)
    CHROM_B = pssql.Column(pssql.String)
    START_B = pssql.Column(pssql.String)
    END_B = pssql.Column(pssql.String)
    SVTYPE = pssql.Column(pssql.String)
    SV_LENGTH = pssql.Column(pssql.String)
    SUPPORT_READS = pssql.Column(pssql.String)
    TOOL = pssql.Column(pssql.String)
    SAMPLE = pssql.Column(pssql.String)
    GENE_A = pssql.Column(pssql.String)
    IN_DESIGN_A = pssql.Column(pssql.String)
    GENE_B = pssql.Column(pssql.String)
    IN_DESIGN_B = pssql.Column(pssql.String)
    GENE_A_GENE_B_sorted = pssql.Column('GENE_A-GENE_B-sorted', pssql.String)
    CALL = pssql.Column(pssql.String)
    TYPE = pssql.Column(pssql.String)
    SECONDHIT = pssql.Column(pssql.String)
    COMMENT  = pssql.Column(pssql.String)
    ASSESSMENT = pssql.Column(pssql.String)
    CLONALITY  = pssql.Column(pssql.String)
    CONSEQUENCE  = pssql.Column(pssql.String)
    FUNCTIONAL_TYPE  = pssql.Column(pssql.String)
    VARIANT_STRING  = pssql.Column(pssql.String)
    include_variant_report_pdf = pssql.Column(pssql.String)
    user_name = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.PROJECT_ID = row_dict.get('PROJECT_ID', None)
        self.SDID = row_dict.get('SDID', None)
        self.CAPTURE_ID = row_dict.get('CAPTURE_ID', None)
        self.CHROM_A = row_dict.get('CHROM_A', None)
        self.START_A = row_dict.get('START_A', None)
        self.END_A = row_dict.get('END_A', None)
        self.CHROM_B = row_dict.get('CHROM_B', None)
        self.START_B = row_dict.get('START_B', None)
        self.END_B = row_dict.get('END_B', None)
        self.SVTYPE = row_dict.get('SVTYPE', None)
        self.SV_LENGTH = row_dict.get('SV_LENGTH', None)
        self.SUPPORT_READS = row_dict.get('SUPPORT_READS', None)
        self.TOOL = row_dict.get('TOOL', None)
        self.SAMPLE = row_dict.get('SAMPLE', None)
        self.GENE_A = row_dict.get('GENE_A', None)
        self.IN_DESIGN_A = row_dict.get('IN_DESIGN_A', None)
        self.GENE_B = row_dict.get('GENE_B', None)
        self.IN_DESIGN_B = row_dict.get('IN_DESIGN_B', None)
        self.GENE_A_GENE_B_sorted = row_dict.get('GENE_A-GENE_B-sorted', None)
        self.CALL = row_dict.get('CALL', None)
        self.TYPE = row_dict.get('TYPE', None)
        self.SECONDHIT = row_dict.get('SECONDHIT', None)
        self.COMMENT = row_dict.get('COMMENT', None)
        self.ASSESSMENT = row_dict.get('ASSESSMENT', None)
        self.CLONALITY = row_dict.get('CLONALITY', None)
        self.CONSEQUENCE = row_dict.get('CONSEQUENCE', None)
        self.FUNCTIONAL_TYPE = row_dict.get('FUNCTIONAL_TYPE', None)
        self.VARIANT_STRING = row_dict.get('VARIANT_STRING', None)
        self.include_variant_report_pdf = row_dict.get('include_variant_report_pdf', None)
        self.user_name = row_dict.get('user_name', None)


    def __repr__(self):
        return "<TableSVS (id='%s', projectid='%s', sdid='%s', captureid='%s')>" % (self.id,self.PROJECT_ID,
                                                                                            self.SDID, self.CAPTURE_ID)

class TableIgvHotspot(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "hotspot_table"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    gene = pssql.Column(pssql.String)
    aapos = pssql.Column(pssql.String)
    nmut = pssql.Column(pssql.String)
    protmut = pssql.Column(pssql.String)
    prot2mut = pssql.Column(pssql.String)
    dnamut = pssql.Column(pssql.String)
    canmut = pssql.Column(pssql.String)
    conseqmut = pssql.Column(pssql.String)
    transcript = pssql.Column(pssql.String)
    dn_ds = pssql.Column(pssql.String) 
    community_notes = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.id = row_dict.get('id', None)
        self.gene = row_dict.get('gene', None)
        self.aapos = row_dict.get('aapos', None)
        self.nmut = row_dict.get('nmut', None)
        self.protmut = row_dict.get('protmut', None)
        self.prot2mut = row_dict.get('prot2mut', None)
        self.dnamut = row_dict.get('dnamut', None)
        self.canmut = row_dict.get('canmut', None)
        self.conseqmut = row_dict.get('conseqmut', None)
        self.transcript = row_dict.get('transcript', None)
        self.dn_ds = row_dict.get('dn_ds', None) 
        self.community_notes = row_dict.get('community_notes', None)

    def __repr__(self):
        return "<TableIgvHotspot (id='%s', gene='%s', protmut='%s', prot2mut='%s')>" % (self.id,self.gene,
                                                                                            self.protmut,self.prot2mut)

class TableIgvCancerHotspot(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "cancer_hotspot_summary"

    h_id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    gene = pssql.Column(pssql.String)
    hgvsp = pssql.Column(pssql.String)
    amino_acid_position = pssql.Column(pssql.String)
    start_aa = pssql.Column(pssql.Integer)
    end_aa = pssql.Column(pssql.Integer)

    def __init__(self, row_dict):

        self.h_id = row_dict.get('h_id', None)
        self.gene = row_dict.get('gene', None)
        self.hgvsp = row_dict.get('hgvsp', None)
        self.amino_acid_position = row_dict.get('amino_acid_position', None)
        self.start_aa = row_dict.get('start_aa', None)
        self.end_aa = row_dict.get('end_aa', None)

    def __repr__(self):
        return "<TableIgvCancerHotspot (h_id='%s', gene='%s', hgvsp='%s')>" % (self.h_id,self.gene, self.hgvsp)

class TableIgvHotspotUpdate(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "hotspot_summary"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    gene = pssql.Column(pssql.String)
    residue = pssql.Column(pssql.String)
    res_type = pssql.Column(pssql.String)
    variants = pssql.Column(pssql.String)
    variant_arr = pssql.Column(pssql.String)

    def __init__(self, row_dict):
        
        self.residue = row_dict.get('residue', None)
        self.res_type = row_dict.get('res_type', None)
        self.variants = row_dict.get('variants', None)
        self.variant_arr = row_dict.get('variant_arr', None)

    def __repr__(self):
        return "<TableIgvHotspot (id='%s', gene='%s', variant_arr='%s')>" % (self.id,self.gene, self.variant_arr)

class TableIgvWarmspot(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "warmspots_table"

    id = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    gene = pssql.Column(pssql.String)
    aapos = pssql.Column(pssql.String)
    nmut = pssql.Column(pssql.String)
    protmut = pssql.Column(pssql.String)
    prot2mut = pssql.Column(pssql.String)
    dnamut = pssql.Column(pssql.String)
    canmut = pssql.Column(pssql.String)
    conseqmut = pssql.Column(pssql.String)
    transcript = pssql.Column(pssql.String)
    dn_ds = pssql.Column(pssql.String) 
    community_notes = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.id = row_dict.get('id', None)
        self.gene = row_dict.get('gene', None)
        self.aapos = row_dict.get('aapos', None)
        self.nmut = row_dict.get('nmut', None)
        self.protmut = row_dict.get('protmut', None)
        self.prot2mut = row_dict.get('prot2mut', None)
        self.dnamut = row_dict.get('dnamut', None)
        self.canmut = row_dict.get('canmut', None)
        self.conseqmut = row_dict.get('conseqmut', None)
        self.transcript = row_dict.get('transcript', None)
        self.dn_ds = row_dict.get('dn_ds', None) 
        self.community_notes = row_dict.get('community_notes', None)

    def __repr__(self):
        return "<TableIgvWarmspot (id='%s', gene='%s', protmut='%s', prot2mut='%s')>" % (self.id,self.gene,
                                                                                            self.protmut,self.prot2mut)

class TablePsffSummary(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "psff_summary"

    id = pssql.Column(pssql.Integer, primary_key=True, autoincrement=True)
    sample_id = pssql.Column(pssql.String)
    capture_id = pssql.Column(pssql.String)
    basic_qc = pssql.Column(pssql.String)
    franken_plot = pssql.Column(pssql.String)
    ploidy = pssql.Column(pssql.String)
    ctdna_fraction = pssql.Column(pssql.String)
    ctdna_param = pssql.Column(pssql.String)
    ctdna_method = pssql.Column(pssql.String)
    ctdna_category = pssql.Column(pssql.String)
    study_code = pssql.Column(pssql.String)
    study_site = pssql.Column(pssql.String)
    somatic_mutations = pssql.Column(pssql.String) 
    germline_alterations = pssql.Column(pssql.String)
    structural_variants = pssql.Column(pssql.String)
    cnvs = pssql.Column(pssql.String)
    comment_info = pssql.Column(pssql.String)


    def __init__(self, row_dict):

        self.id = row_dict.get('id', None)
        self.sample_id = row_dict.get('sample_id', None)
        self.capture_id = row_dict.get('capture_id', None)
        self.basic_qc = row_dict.get('basic_qc', None)
        self.franken_plot = row_dict.get('franken_plot', None)
        self.ploidy = row_dict.get('ploidy', None)
        self.ctdna_fraction = row_dict.get('ctdna_fraction', None)
        self.ctdna_param = row_dict.get('ctdna_param', None)
        self.ctdna_method = row_dict.get('ctdna_method', None)
        self.ctdna_category = row_dict.get('ctdna_category', None)
        self.study_code = row_dict.get('study_code', None)
        self.study_site = row_dict.get('study_site', None)
        self.somatic_mutations = row_dict.get('somatic_mutations', None)
        self.germline_alterations = row_dict.get('germline_alterations', None)
        self.structural_variants = row_dict.get('structural_variants', None)
        self.cnvs = row_dict.get('cnvs', None)
        self.comment_info = row_dict.get('comment_info', None)



    def __repr__(self):
        return "<TablePsffSummary (id='%s', sample_id='%s', capture_id='%s')>" % (self.id,self.sample_id,self.capture_id)


class TableProbioSummary(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "probio_summary"

    id = pssql.Column(pssql.Integer, primary_key=True, autoincrement=True)
    sample_id = pssql.Column(pssql.String)
    capture_id = pssql.Column(pssql.String)
    basic_qc = pssql.Column(pssql.String)
    franken_plot = pssql.Column(pssql.String)
    ploidy = pssql.Column(pssql.String)
    ctdna_fraction = pssql.Column(pssql.String)
    ctdna_param = pssql.Column(pssql.String)
    ctdna_method = pssql.Column(pssql.String)
    ctdna_category = pssql.Column(pssql.String)
    study_code = pssql.Column(pssql.String)
    study_site = pssql.Column(pssql.String)
    somatic_mutations = pssql.Column(pssql.String) 
    germline_alterations = pssql.Column(pssql.String)
    structural_variants = pssql.Column(pssql.String)
    cnvs = pssql.Column(pssql.String)
    comment_info = pssql.Column(pssql.String)


    def __init__(self, row_dict):

        self.id = row_dict.get('id', None)
        self.sample_id = row_dict.get('sample_id', None)
        self.capture_id = row_dict.get('capture_id', None)
        self.basic_qc = row_dict.get('basic_qc', None)
        self.franken_plot = row_dict.get('franken_plot', None)
        self.ploidy = row_dict.get('ploidy', None)
        self.ctdna_fraction = row_dict.get('ctdna_fraction', None)
        self.ctdna_param = row_dict.get('ctdna_param', None)
        self.ctdna_method = row_dict.get('ctdna_method', None)
        self.ctdna_category = row_dict.get('ctdna_category', None)
        self.study_code = row_dict.get('study_code', None)
        self.study_site = row_dict.get('study_site', None)
        self.somatic_mutations = row_dict.get('somatic_mutations', None)
        self.germline_alterations = row_dict.get('germline_alterations', None)
        self.structural_variants = row_dict.get('structural_variants', None)
        self.cnvs = row_dict.get('cnvs', None)
        self.comment_info = row_dict.get('comment_info', None)



    def __repr__(self):
        return "<TableProbioSummary (id='%s', sample_id='%s', capture_id='%s')>" % (self.id,self.sample_id,self.capture_id)

class TableGenomicProfileSummary(pssql.Model):
    __bind_key__ = 'curation'
    __tablename__ = "genomic_profile_summary"

    id = pssql.Column(pssql.Integer, primary_key=True, autoincrement=True)
    project_name = pssql.Column(pssql.String)
    sample_id = pssql.Column(pssql.String)
    capture_id = pssql.Column(pssql.String)
    study_code = pssql.Column(pssql.String)
    study_site = pssql.Column(pssql.String)
    dob = pssql.Column(pssql.String)
    disease = pssql.Column(pssql.String)
    specimen_assay = pssql.Column(JSON)
    ctdna_param = pssql.Column(pssql.String)
    ctdna_method = pssql.Column(pssql.String)
    genome_wide = pssql.Column(JSON)
    somatic_mutations = pssql.Column(JSON) 
    germline_alterations = pssql.Column(JSON)
    structural_variants = pssql.Column(JSON)
    cnvs = pssql.Column(JSON)
    summary_txt = pssql.Column(pssql.String)

    def __init__(self, row_dict):

        self.id = row_dict.get('id', None)
        self.project_name = row_dict.get('project_name', None)
        self.sample_id = row_dict.get('sample_id', None)
        self.capture_id = row_dict.get('capture_id', None)
        self.study_code = row_dict.get('study_code', None)
        self.study_site = row_dict.get('study_site', None)
        self.dob = row_dict.get('dob', None)
        self.disease = row_dict.get('disease', None)
        self.specimen_assay = row_dict.get('specimen_assay', None)
        self.basic_qc = row_dict.get('basic_qc', None)
        self.franken_plot = row_dict.get('franken_plot', None)
        self.ctdna_param = row_dict.get('ctdna_param', None)
        self.ctdna_method = row_dict.get('ctdna_method', None)
        self.genome_wide = row_dict.get('genome_wide', None)
        self.somatic_mutations = row_dict.get('somatic_mutations', None)
        self.germline_alterations = row_dict.get('germline_alterations', None)
        self.structural_variants = row_dict.get('structural_variants', None)
        self.cnvs = row_dict.get('cnvs', None)
        self.summary_txt = row_dict.get('summary_txt', None)


    def __repr__(self):
        return "<TableGenomicProfileSummary (id='%s', project_name='%s', sample_id='%s', capture_id='%s')>" % (self.id,self.project_name, self.sample_id,self.capture_id)
