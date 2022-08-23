import logging
from flask import current_app
from flask import request, send_file, make_response, send_from_directory
from flask_restx import Resource
from franken_api.api.franken.parsers import curation_germline_arguments, curation_somatic_arguments, curation_svs_arguments, curation_psff_profile_arguments, curation_psff_profile_id_arguments, curation_probio_profile_id_arguments, curation_probio_profile_arguments, curation_genomic_profile_arguments, curation_genomic_profile_id_arguments, project_list_arguments, cancer_hotspot_arguments
from franken_api.api.restplus import api
from franken_api.api.franken.business import  get_curation_igv_germline, get_curation_igv_somatic, get_curation_svs, post_curation, get_curation_hotspot, get_curation_warmspot, get_curation_psff_profile, curation_update_profile, list_curation_psff_profile, list_curation_probio_profile, get_curation_probio_profile, get_curation_cancer_hotspot, list_curation_genomic_profile, get_curation_genomic_profile, fetch_cancer_hotsport_info
from franken_api.api.franken.serializers import curation_germline, germline_data_list, somatic_data_list, svs_data_list, hotspot_data_list, warmspot_data_list, psff_profile_data_list, probio_profile_data_list, cancer_hotspot_data_list, genomic_profile_data_list

log = logging.getLogger(__name__)
ns3 = api.namespace('curation', description='Curation Database API')

@ns3.route('/')
@api.response(200, 'Check API status')
@api.response(400, '/nfs is not mount locally no data found')
class CurationStatus(Resource):
    def get(self):
        """
        Referral API status check
        ```

        ```
        """
        return 'Working', 200

@ns3.route('/igv/cancer-hotspot')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationCAncerHotspotTable(Resource):
    @api.marshal_with(cancer_hotspot_data_list)
    def get(self):
        """
        Fetch all cancer hotspot information
        ```

        ```
        """
        result, error = get_curation_cancer_hotspot()
        return result, error

    @api.expect(cancer_hotspot_arguments, validate=True)
    @api.marshal_with(cancer_hotspot_data_list)
    def post(self):
        """
           Fetch the patient information
        """
        args = cancer_hotspot_arguments.parse_args()
        gene = args['gene']
        HGVSp = args['HGVSp']
        position = args['position']
        # consequence = args['consequence']
        result, errorcode = fetch_cancer_hotsport_info(gene, HGVSp, position)
        return result, errorcode


@ns3.route('/igv/hotspot')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationHotspotTable(Resource):
    @api.marshal_with(hotspot_data_list)
    def get(self):
        """
        Fetch all hotspot information
        ```

        ```
        """
        result, error = get_curation_hotspot()
        return result, error

@ns3.route('/igv/warmspot')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationWarmspotTable(Resource):
    @api.marshal_with(warmspot_data_list)
    def get(self):
        """
        Fetch all warmspot information
        ```

        ```
        """
        result, error = get_curation_warmspot()
        return result, error

@ns3.route('/igv/germline')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationIgvGermline(Resource):
    @api.marshal_with(germline_data_list)
    @api.expect(project_list_arguments, validate=True)
    def get(self):
        """
        Fetch all Germline based on project access
        ```

        ```
        """
        args = project_list_arguments.parse_args()
        project_ids = args["project_ids"]
        result, error = get_curation_igv_germline(project_ids)
        return result, error

    @api.expect(curation_germline_arguments, validate=True)
    def post(self):
        """
        Fetch all Germline
        ```

        ```
        """
        args = curation_germline_arguments.parse_args()
        result, errorcode = post_curation(dict(args), 'germline')
        return result, 200 #result, errorcode

@ns3.route('/igv/somatic')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationIgvSomatic(Resource):
    @api.marshal_with(somatic_data_list)
    @api.expect(project_list_arguments, validate=True)
    def get(self):
        """
        Fetch all Somatic  based on project access
        ```

        ```
        """
        args = project_list_arguments.parse_args()
        project_ids = args["project_ids"]
        result, error = get_curation_igv_somatic(project_ids)
        return result, error

    @api.expect(curation_somatic_arguments, validate=True)
    def post(self):
        """
        Fetch all Somatic
        ```

        ```
        """
        args = curation_somatic_arguments.parse_args()
        result, errorcode = post_curation(dict(args), 'somatic')
        return result, errorcode

@ns3.route('/igv/svs')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationSVS(Resource):
    @api.marshal_with(svs_data_list)
    @api.expect(project_list_arguments, validate=True)
    def get(self):
        """
        Fetch all SVS  based on project access
        ```

        ```
        """
        args = project_list_arguments.parse_args()
        project_ids = args["project_ids"]
        result, error = get_curation_svs(project_ids)
        return result, error

    @api.expect(curation_svs_arguments, validate=True)
    def post(self):
        """
        Fetch all Germline
        ```

        ```
        """
        args = curation_svs_arguments.parse_args()
        result, errorcode = post_curation(dict(args), 'svs')
        return result, errorcode


@ns3.route('/psff_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationPSFFProfile(Resource):
    @api.marshal_with(psff_profile_data_list)
    @api.expect(curation_psff_profile_id_arguments, validate=True)
    def post(self):
        """
         Get the PSFF Genomic Profiling
        ```

        ```
        """
        args = curation_psff_profile_id_arguments.parse_args()
        result, error = get_curation_psff_profile(dict(args))
        return result, error

@ns3.route('/probio_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationPROBIOProfile(Resource):
    @api.marshal_with(probio_profile_data_list)
    @api.expect(curation_probio_profile_id_arguments, validate=True)
    def post(self):
        """
         Get the PSFF Genomic Profiling
        ```

        ```
        """
        args = curation_probio_profile_id_arguments.parse_args()
        result, error = get_curation_probio_profile(dict(args))
        return result, error

@ns3.route('/genomic_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationGenomicProfile(Resource):
    @api.marshal_with(genomic_profile_data_list)
    @api.expect(curation_genomic_profile_id_arguments, validate=True)
    def post(self):
        """
         Get the Genomic Profiling for all project 
        ```

        ```
        """
        args = curation_genomic_profile_id_arguments.parse_args()
        result, error = get_curation_genomic_profile(dict(args))
        return result, error

@ns3.route('/list_profile/psff')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationPSFFProfile(Resource):
    @api.marshal_with(psff_profile_data_list)
    def get(self):
        """
        Fetch all PSFF Genomic Profiling
        ```

        ```
        """
        result, error = list_curation_psff_profile()
        return result, error

@ns3.route('/list_profile/probio')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationPSFFProfile(Resource):
    @api.marshal_with(probio_profile_data_list)
    def get(self):
        """
        Fetch all PSFF Genomic Profiling
        ```

        ```
        """
        result, error = list_curation_probio_profile()
        return result, error

@ns3.route('/list_profile/genomic')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationListGenomicProfile(Resource):
    @api.marshal_with(genomic_profile_data_list)
    @api.expect(project_list_arguments, validate=True)
    def get(self):
        """
        Fetch all Genomic Profiling Profiling based on projectIds
        ```

        ```
        """
        args = project_list_arguments.parse_args()
        project_ids = args["project_ids"]
        result, error = list_curation_genomic_profile(project_ids)
        return result, error

    # @api.expect(project_list_arguments, validate=True)
    # def post(self):
    #     """
    #      Fetch Genomic Profiling based on projectIds
    #     ```

    #     ```
    #     """
    #     args = project_list_arguments.parse_args()
    #     project_ids = args["project_ids"]
    #     result, errorcode = fetch_genomic_profile(project_ids)
    #     return result, errorcode

@ns3.route('/update_psff_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationUpdatePSFFProfile(Resource):
    @api.expect(curation_psff_profile_arguments, validate=True)
    def post(self):
        """
         Update PSFF Genomic Profiling
        ```

        ```
        """
        args = curation_psff_profile_arguments.parse_args()
        result, errorcode = curation_update_profile(dict(args), 'psff_summary')
        return result, errorcode

@ns3.route('/update_probio_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationUpdatePSFFProfile(Resource):
    @api.expect(curation_probio_profile_arguments, validate=True)
    def post(self):
        """
         Update PSFF Genomic Profiling
        ```

        ```
        """
        args = curation_probio_profile_arguments.parse_args()
        result, errorcode = curation_update_profile(dict(args), 'probio_summary')
        return result, errorcode


@ns3.route('/update_genomic_profile')
@api.response(200, 'Success')
@api.response(400, '/nfs is not mount locally no data found')
class CurationUpdateGenomicProfile(Resource):
    @api.expect(curation_genomic_profile_arguments, validate=True)
    def post(self):
        """
         Update Genomic Profiling for project
        ```

        ```
        """
        args = curation_genomic_profile_arguments.parse_args()
        result, errorcode = curation_update_profile(dict(args), 'genomic_profile_summary')
        return result, errorcode
