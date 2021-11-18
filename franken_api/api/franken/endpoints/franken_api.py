import logging
from flask import current_app
from flask import request, send_file, make_response, send_from_directory
from flask_restplus import Resource
#from franken_api.api.franken.serializers import search_result
from franken_api.api.franken.parsers import pdf_arguments, table_cnv_arguments, search_arguments, capture_arguments, ploturls_arguments, staticplot_arguments, igv_arguments, table_svs_arguments, project_arguments, table_igvnav_arguments, igv_save_file_arguments, table_qc_arguments, purecn_arguments, purecn_max_val_arguments
from franken_api.api.restplus import api
from flask import jsonify
from franken_api.api.franken.business import pdfs_files, get_table_cnv_header, check_nfs_mount, get_sample_ids, get_sample_design_ids, get_static_frankenplot, get_static_image, get_interactive_plot, get_table_svs_header, get_table_igv, save_igvnav_input_file, get_table_qc_header, get_purecn_ctdna, update_pureCN_somatic_germline, get_xml_image
from franken_api.api.franken.serializers import status_result, dropdownlist, dropdownlist_capture, ploturl_list
import io
#import  franken_api.database.models
log = logging.getLogger(__name__)
ns = api.namespace('franken', description='Interactive franken Plots')

@ns.route('/')
@api.response(200, 'Status of API')
@api.response(400, '/nfs is not mount locally')
class frankenStatus(Resource):
    @api.expect(project_arguments, validate=True)
    @api.marshal_with(status_result)
    def get(self):
        """
        Returns status of the endpoint.
        ```
        {
            "server_status": true
        }

        ```
        """
        args = project_arguments.parse_args()
        proj_name = args['project_name']
        status, error_code = check_nfs_mount(current_app.config[proj_name])
        result = {'server_status': True}

        if not status:
            return {'server_status': False}, error_code


        return result, error_code

@ns.route('/samples')
@api.response(200, 'All Samples for dropdown')
@api.response(400, '/nfs is not mount locally no data found')
class DropdownListSample(Resource):
    @api.marshal_with(dropdownlist)
    @api.expect(project_arguments, validate=True)
    def get(self):
        """
        Returns List of sample ids for dropdown in UI.
        ```
        [sdid1, sdid2,......]
        ```
        """
        args = project_arguments.parse_args()
        proj_name = args['project_name']
        result, errorcode = get_sample_ids(current_app.config[proj_name])
        return result, errorcode

@ns.route('/capture')
@api.response(200, 'All Samples for dropdown')
@api.response(400, '/nfs is not mount locally no data found')
class DropdownListCapture(Resource):
    @api.expect(capture_arguments, validate=True)
    @api.marshal_with(dropdownlist_capture)
    def get(self):
        """
        Returns List of capture ids for given sample ids which are used for dropdown in UI.
        ```
        [cpture_id1, capture2,......]
        ```
        """
        args = capture_arguments.parse_args()
        result, errorcode = get_sample_design_ids(current_app.config[args['project_name']], args['sdid'])
        return result, errorcode


@ns.route('/ploturls')
@api.response(200, 'Franken plot urls')
@api.response(400, 'No franken plot images found in qc folder')
class FrankenUrls(Resource):
    @api.expect(ploturls_arguments, validate=True)
    @api.marshal_with(ploturl_list)
    def get(self):
        """
        Returns List of franken plot url.
        ```
        [url1, url2,......]
        ```
        """
        args = ploturls_arguments.parse_args()
        result, errorcode = get_static_frankenplot(current_app.config[args['project_name']], args['project_name'], args['sdid'], args['capture_id'])
        return result, errorcode


@ns.route('/staticimage')
@api.response(200, 'Franken Static plot')
@api.response(400, 'No Static plots found')
class FrankenStaticImages(Resource):
    @api.representation('image/png')
    @api.expect(staticplot_arguments, validate=True)
    def get(self):
        """
        Returns static franken plot.
        ```
        base64 of png
        ```
        """
        args = staticplot_arguments.parse_args()
        result, errorcode = get_static_image(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['imagename'])
        return send_file(result,
                      attachment_filename='frankenplot.png',
                      mimetype='image/png')

@ns.route('/igvXML')
@api.response(200, 'IGV Session XML')
@api.response(400, 'No XML found')
class FrankenStaticImages(Resource):
    @api.expect(staticplot_arguments, validate=True)
    def get(self):
        """
        Returns static franken plot.
        ```
        session XML
        ```
        """
        args = staticplot_arguments.parse_args()
        result, errorcode = get_xml_image(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['imagename'])
        return send_file(result, mimetype='application/xml')


@ns.route('/plot')
@api.response(200, 'Json file to plot')
@api.response(400, 'sample or json file not found')
class frankenPlot(Resource):
    @api.expect(search_arguments, validate=True)
    def get(self):
        """
        Returns Json file to plot.
        ```
        Json data to plot the franken plots
        ```
        """
        args = search_arguments.parse_args()
        result , errocode = get_interactive_plot(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['pname'])
        return result, errocode

@ns.route('/igvsession')
@api.response(200, 'Files to load in igv tracks')
@api.response(400, 'No file found')
class IGVTracksFiles(Resource):
    @api.expect(igv_arguments, validate=True)
    def get(self):
        """
        Returns Files Required To Load in IGV Tracks
        ```
        bam, vcf etc.
        ```
        """
        args = igv_arguments.parse_args()
        if args['filename'] == '1':
            return send_from_directory('/nfs/PROBIO/autoseq-output/P-00356971/PB-P-00356971-CFDNA-03589573-KH20190515-C220190515_PB-P-00356971-N-03589575-KH20190515-C220190515/bams/C2',
                             'PB-P-00356971-CFDNA-03589573-KH-C2-nodups.bam' )

        return send_from_directory('/nfs/PROBIO/autoseq-output/P-00356971/PB-P-00356971-CFDNA-03589573-KH20190515-C220190515_PB-P-00356971-N-03589575-KH20190515-C220190515/bams/C2',
                         'PB-P-00356971-CFDNA-03589573-KH-C2-nodups.bam.bai')
        #return send_from_directory('/nfs/PROBIO/autoseq-output/P-00356971/PB-P-00356971-CFDNA-03589573-KH20190515-C220190515_PB-P-00356971-N-03589575-KH20190515-C220190515/bams/C2/', 'PB-P-00356971-CFDNA-03589573-KH-C2-nodups.bam.bai')


@ns.route('/table/svs')
@api.response(200, 'All Structural Variants')
@api.response(400, '/nfs is not mount locally no data found')
class TableSvs(Resource):
    @api.expect(table_svs_arguments, validate=True)
    def get(self):
        """
        Returns All Structural Variants .
        ```
        { 'header': {
                    columnTitle1:{ title: 'ID', type: 'number', editable:false},
                    columnTitle2:{ title: 'ID', type: 'string', editable:false}
                    },
          'data' : [
                    { columnTitle1: value1,  columnTitle2: value2 }
          ]
        }
        ```
        """
        args = table_svs_arguments.parse_args()
        result, errorcode = get_table_svs_header(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['header'])
        return result, errorcode

@ns.route('/table/igv/<string:variant>')
@api.response(200, 'All Germline and somatic Variants')
@api.response(400, '/nfs is not mount locally no data found')
class TableIgv(Resource):
    @api.expect(table_igvnav_arguments, validate=True)
    def get(self, variant):
        """
        Returns All Structural Variants .
        ```
        { 'header': {
                    columnTitle1:{ title: 'ID', type: 'number', editable:false},
                    columnTitle2:{ title: 'ID', type: 'string', editable:false}
                    },
          'data' : #[
                    { columnTitle1: value1,  columnTitle2: value2 }
          ]
        }
        ```
        """
        args = table_svs_arguments.parse_args()
        result, errorcode = get_table_igv(variant, current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['header'])
        return result, errorcode

@ns.route('/table/qc')
@api.response(200, 'Sample QC Metrics')
@api.response(400, '/nfs is not mount locally no data found')
class TableQc(Resource):
    @api.expect(table_qc_arguments, validate=True)
    def get(self):
        """
        Returns All QC Metrics For Samples .
        ```
        { 'header': {
                    columnTitle1:{ title: 'ID', type: 'number', editable:false},
                    columnTitle2:{ title: 'ID', type: 'string', editable:false}
                    },
          'data' : [
                    { columnTitle1: value1,  columnTitle2: value2 }
          ]
        }
        ```
        """
        args = table_qc_arguments.parse_args()
        result, errorcode = get_table_qc_header(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['header'])
        return result, errorcode

@ns.route('/table/cnv/<string:variant_type>')
@api.response(200, 'CNV Metrics')
@api.response(400, '/nfs is not mount locally no data found')
class TableCNV(Resource):
    @api.expect(table_cnv_arguments, validate=True)
    def get(self, variant_type):
        """
        Returns All QC Metrics For Samples .
        ```
        { 'header': {
                    columnTitle1:{ title: 'ID', type: 'number', editable:false},
                    columnTitle2:{ title: 'ID', type: 'string', editable:false}
                    },
          'data' : [
                    { columnTitle1: value1,  columnTitle2: value2 }
          ]
        }
        ```
        """
        args = table_cnv_arguments.parse_args()
        result, errorcode = get_table_cnv_header(current_app.config[args['project_name']], args['sdid'], args['capture_id'], variant_type, args['header'])
        return result, errorcode



@ns.route('/save/igvinput')
@api.response(200, 'Susscessfully saving igvnav files')
@api.response(400, '/nfs is not mount locally no data found')
class SaveIGVFile(Resource):
    #@api.expect(igv_save_file_arguments, validate=True)
    def post(self):
        """
        Saves IGVnav-input.txt file and structural variant file .
        ```

        ```
        """
        args = request.json
        result, errorcode = save_igvnav_input_file(args['file_name'], args['data'])
        return result, errorcode

# pdf endpoints
@ns.route('/pdf/<string:variant>')
@api.response(200, 'PDF file')
@api.response(400, '/nfs is not mount locally no data found')
class PDFCalls(Resource):
    @api.representation('application/pdf')
    @api.expect(pdf_arguments, validate=True)
    def get(self, variant):
        """
        Returns PDF files .
        """
        args = pdf_arguments.parse_args()
        result, errorcode = pdfs_files(variant, current_app.config[args['project_name']], args['sdid'], args['capture_id'])
        if variant != 'multiqc':
            return send_file(result, attachment_filename=variant+'.pdf', mimetype='application/pdf')
        else:
            return send_file(result, attachment_filename=variant+'.html', mimetype='text/html')

@ns.route('/purecn')
@api.response(200, 'CSV file to purecn')
@api.response(400, 'CSV file not found')
class GetPurecn(Resource):
    @api.expect(purecn_arguments, validate=True)
    def get(self):
        """
            Fetch ploidy and purity from purecn file
        """
        args = purecn_arguments.parse_args()
        result , errocode = get_purecn_ctdna(current_app.config[args['project_name']], args['sdid'], args['capture_id'])
        return result, errocode

@ns.route('/get_purecn_max_val')
@api.response(200, 'Update purecn max value into somatic and germline')
@api.response(400, 'Purecn variant file not found')
class UpdateMaxPurecnVal(Resource):
    @api.expect(purecn_max_val_arguments, validate=True)
    def get(self):
        """
        Get Max PureCN value for each variant and update into somatic and germline table
        """
        args = purecn_max_val_arguments.parse_args()
        result , errocode = update_pureCN_somatic_germline(current_app.config[args['project_name']], args['sdid'], args['capture_id'], args['variant_type'])
        return result, errocode