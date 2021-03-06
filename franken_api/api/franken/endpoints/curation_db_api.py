import logging
from flask import current_app
from flask import request, send_file, make_response, send_from_directory
from flask_restplus import Resource
from franken_api.api.franken.parsers import curation_germline_arguments, curation_somatic_arguments, curation_svs_arguments, curation_psff_profile_arguments, curation_psff_profile_id_arguments, curation_probio_profile_id_arguments, curation_probio_profile_arguments
from franken_api.api.restplus import api
from franken_api.api.franken.business import  get_curation_igv_germline, get_curation_igv_somatic, get_curation_svs, post_curation, get_curation_hotspot, get_curation_warmspot, get_curation_psff_profile, curation_update_profile, list_curation_psff_profile, list_curation_probio_profile, get_curation_probio_profile
from franken_api.api.franken.serializers import curation_germline, germline_data_list, somatic_data_list, svs_data_list, hotspot_data_list, warmspot_data_list, psff_profile_data_list, probio_profile_data_list

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
    def get(self):
        """
        Fetch all Germline
        ```

        ```
        """
        result, error = get_curation_igv_germline()
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
    def get(self):
        """
        Fetch all Somatic
        ```

        ```
        """
        result, error = get_curation_igv_somatic()
        return result, error

    @api.expect(curation_somatic_arguments, validate=True)
    def post(self):
        """
        Fetch all Germline
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
    def get(self):
        """
        Fetch all Somatic
        ```

        ```
        """
        result, error = get_curation_svs()
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
