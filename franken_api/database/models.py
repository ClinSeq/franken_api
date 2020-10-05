from franken_api.database import db as pssql


class ProbioBloodReferral(pssql.Model):
    __tablename__ = "probio_bloodreferrals"
    crid = pssql.Column(pssql.Integer, primary_key=True, nullable=False)
    pnr  = pssql.Column(pssql.String, nullable=False)
    rid  = pssql.Column(pssql.String, nullable=False)
    datum = pssql.Column(pssql.Date, nullable=False)
    tid  = pssql.Column(pssql.String, nullable=False)
    sign = pssql.Column(pssql.Integer)
    countyletter  = pssql.Column(pssql.String, nullable=False)
    new = pssql.Column(pssql.String, nullable=False)
    progression = pssql.Column(pssql.String, nullable=False)
    follow_up    = pssql.Column(pssql.String, nullable=False)
    cf_dna1 = pssql.Column(pssql.String, nullable=False)
    cf_dna2 = pssql.Column(pssql.String, nullable=False)
    cf_dna3 = pssql.Column(pssql.String, nullable=False)
    kommentar = pssql.Column(pssql.String, nullable=False)
    filnamn = pssql.Column(pssql.String, nullable=False)
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

    def __repr__(self):
        return "<TableSVS (id='%s', projectid='%s', sdid='%s', captureid='%s')>" % (self.id,self.PROJECT_ID,
                                                                                            self.SDID, self.CAPTURE_ID)
