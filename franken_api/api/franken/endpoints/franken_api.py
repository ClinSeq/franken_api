"""
franken_api.py

"""

import logging

from flask import make_response, request, send_file
from flask_restx import Resource

from franken_api.api.franken.business import (
    check_frankenplot_files,
    check_nfs_mount,
    fetch_curated_pdf,
    fetch_nfs_path,
    fetch_patient_info,
    frankenplot_files,
    generate_curated_json,
    generate_curated_pdf,
    get_cancer_type,
    get_clinical_report_file,
    get_curated_json_file,
    get_igv_report_image,
    get_interactive_plot,
    get_pdf_file,
    get_pdf_file2,
    get_purecn_ctdna,
    get_rna_qc_status,
    get_sample_design_ids,
    get_sample_ids,
    get_static_frankenplot,
    get_static_image,
    get_table_cnv_header,
    get_table_fusion_inspector,
    get_table_fusion_report,
    get_table_igv,
    get_table_qc_header,
    get_table_svs_header,
    pdfs_files,
    read_igv_session_json,
    read_igv_session_xml,
    rna_html_files,
    rna_pdf_files,
    save_igvnav_input_file,
    send_json_mtbp_portal,
    update_curated_info,
    update_pureCN_somatic_germline,
)

# from franken_api.api.franken.serializers import search_result
from franken_api.api.franken.parsers import (
    capture_arguments,
    fetch_patient_info_arguments,
    igv_file_arguments,
    igv_image_report_arguments,
    igv_json_arguments,
    igv_xml_arguments,
    json_urls_arguments,
    pdf_arguments,
    ploturls_arguments,
    project_arguments,
    purecn_arguments,
    purecn_max_val_arguments,
    rna_html_arguments,
    search_arguments,
    send_json_mtbp_arguments,
    staticplot_arguments,
    table_cnv_arguments,
    table_igvnav_arguments,
    table_qc_arguments,
    table_svs_arguments,
    update_curated_info_arguments,
    view_pdf_arguments,
)
from franken_api.api.franken.serializers import (
    dropdownlist,
    dropdownlist_capture,
    ploturl_list,
    status_result,
)
from franken_api.api.restplus import api

# import  franken_api.database.models
log = logging.getLogger(__name__)
ns = api.namespace("franken", description="Interactive franken Plots")


@ns.route("/")
@api.response(200, "Status of API")
@api.response(400, "/nfs is not mount locally")
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
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            status, error_code = check_nfs_mount(nfs_path)
            if not status:
                return {"server_status": False}, error_code
            else:
                return {"server_status": True}, error_code
        else:
            return nfs_path, response_code


@ns.route("/samples")
@api.response(200, "All Samples for dropdown")
@api.response(400, "/nfs is not mount locally no data found")
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
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, error_code = get_sample_ids(nfs_path)
            return result, error_code
        else:
            return nfs_path, response_code


@ns.route("/capture")
@api.response(200, "All Samples for dropdown")
@api.response(400, "/nfs is not mount locally no data found")
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
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_sample_design_ids(nfs_path, args["sdid"])
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/ploturls")
@api.response(200, "Franken plot urls")
@api.response(400, "No franken plot images found in qc folder")
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
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_static_frankenplot(
                nfs_path, proj_name, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/staticimage")
@api.response(200, "Franken Static plot")
@api.response(400, "No Static plots found")
class FrankenStaticImages(Resource):
    @api.representation("image/png")
    @api.expect(staticplot_arguments, validate=True)
    def get(self):
        """
        Returns static franken plot.
        ```
        base64 of png
        ```
        """
        args = staticplot_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_static_image(
                nfs_path, args["sdid"], args["capture_id"], args["imagename"]
            )
            if errorcode == 200:
                return send_file(
                    result, attachment_filename="frankenplot.png", mimetype="image/png"
                )
            else:
                return "", errorcode

        else:
            return nfs_path, response_code


@ns.route("/plot")
@api.response(200, "Json file to plot")
@api.response(400, "sample or json file not found")
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
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_interactive_plot(
                nfs_path, args["sdid"], args["capture_id"], args["pname"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/svs")
@api.response(200, "All Structural Variants")
@api.response(400, "/nfs is not mount locally no data found")
class TableSvsView(Resource):
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
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_svs_header(
                nfs_path, args["sdid"], args["capture_id"], user_id, args["header"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/igv/<string:variant>")
@api.response(200, "All Germline and somatic Variants")
@api.response(400, "/nfs is not mount locally no data found")
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
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_igv(
                variant,
                nfs_path,
                args["sdid"],
                args["capture_id"],
                user_id,
                args["header"],
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/qc")
@api.response(200, "Sample QC Metrics")
@api.response(400, "/nfs is not mount locally no data found")
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
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_qc_header(
                nfs_path, args["sdid"], args["capture_id"], user_id, args["header"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/cnv/<string:variant_type>")
@api.response(200, "CNV Metrics")
@api.response(400, "/nfs is not mount locally no data found")
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
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_cnv_header(
                nfs_path,
                args["sdid"],
                args["capture_id"],
                variant_type,
                user_id,
                args["header"],
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/save/igvinput")
@api.response(200, "Susscessfully saving igvnav files")
@api.response(400, "/nfs is not mount locally no data found")
class SaveIGVFile(Resource):
    # @api.expect(igv_save_file_arguments, validate=True)
    def post(self):
        """
        Saves IGVnav-input.txt file and structural variant file .
        ```

        ```
        """
        args = request.json
        result, errorcode = save_igvnav_input_file(args["file_name"], args["data"])
        return result, errorcode


# pdf endpoints
@ns.route("/pdf/<string:variant>")
@api.response(200, "PDF file")
@api.response(400, "/nfs is not mount locally no data found")
class PDFCalls(Resource):
    @api.representation("application/pdf")
    @api.expect(pdf_arguments, validate=True)
    def get(self, variant):
        """
        Returns PDF files .
        """
        args = pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = pdfs_files(
                variant, nfs_path, args["sdid"], args["capture_id"]
            )
            if variant != "multiqc":
                return send_file(
                    result,
                    attachment_filename=variant + ".pdf",
                    mimetype="application/pdf",
                )
            else:
                return send_file(
                    result, attachment_filename=variant + ".html", mimetype="text/html"
                )
        else:
            return nfs_path, response_code


# RNA PDF & HTML
@ns.route("/RNAReport/<string:file_format>/<string:report_type>")
@api.response(200, "PDF & HTML  file")
@api.response(400, "/nfs is not mount locally no data found")
class FetchRNAReports(Resource):
    # @api.representation('application/pdf')
    @api.expect(pdf_arguments, validate=True)
    def get(self, file_format, report_type):
        """
        Returns PDF & HTML  files .
        """
        args = pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            if file_format == "pdf":
                result, errorcode = rna_pdf_files(
                    file_format, nfs_path, args["sdid"], args["capture_id"]
                )
                return send_file(
                    result,
                    attachment_filename=report_type + ".pdf",
                    mimetype="application/pdf",
                )
            else:
                file_format = (
                    file_format.replace("R", "") + "_fastqc"
                    if report_type == "qc"
                    else file_format
                )
                result, errorcode = rna_html_files(
                    report_type, file_format, nfs_path, args["sdid"], args["capture_id"]
                )
                return send_file(
                    result,
                    attachment_filename=file_format + ".html",
                    mimetype="text/html",
                )
        else:
            return nfs_path, response_code


@ns.route("/RNAReport/html/fusion-report")
@api.response(200, "PDF & HTML  file")
@api.response(400, "/nfs is not mount locally no data found")
class FetchHTMLReports(Resource):
    @api.expect(rna_html_arguments, validate=True)
    def get(self):
        """
        Returns HTML  files .
        """
        args = rna_html_arguments.parse_args()
        proj_name = args["project_name"]
        file_name = args["file_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            file_path = nfs_path + "/" + args["sdid"] + "/fusionreport/" + file_name
            return send_file(
                file_path, attachment_filename=file_name + ".html", mimetype="text/html"
            )
        else:
            return nfs_path, response_code


# @ns.route('/RNAReport/html/<string:file_name>')
# @api.response(200, 'PDF & HTML  file')
# @api.response(400, '/nfs is not mount locally no data found')
# class FetchHTMLReports(Resource):
# 	@api.expect(pdf_arguments, validate=True)
# 	def get(self, file_name):
# 		"""
# 		Returns HTML  files .
# 		"""
# 		args = pdf_arguments.parse_args()
# 		proj_name = args['project_name']
# 		nfs_path,response_code = fetch_nfs_path(proj_name)
# 		if response_code == 200:
# 			file_path = nfs_path + '/' + args['sdid'] + '/' + args['capture_id'] + "fusionreport/" + file_name
# 			return send_file(file_path, attachment_filename=file_name+'.html', mimetype='text/html')
# 		else:
# 			return nfs_path,response_code


## check New Franken Plots
@ns.route("/checkFrankenPlot")
@api.response(200, "check new franken plot generate or not ")
@api.response(400, "/nfs is not mount locally no data found")
class CheckFrankenPlotHtmlFile(Resource):
    @api.expect(pdf_arguments, validate=True)
    def get(self):
        """
        Returns List of sample ids for dropdown in UI.
        ```
        [sdid1, sdid2,......]
        ```
        """
        args = pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = check_frankenplot_files(
                nfs_path, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


# New Franken Plots
@ns.route("/frankenplot")
@api.response(200, "Franken New HTML file")
@api.response(400, "/nfs is not mount locally no data found")
class FrankenPlotHtmlFile(Resource):
    @api.representation("text/html")
    @api.expect(pdf_arguments, validate=True)
    def get(self):
        """
        Returns HTML files .
        """
        args = pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = frankenplot_files(
                nfs_path, args["sdid"], args["capture_id"]
            )
            file_name = args["sdid"] + "_frankenplot.html"
            if errorcode == 200:
                return send_file(
                    result, attachment_filename=file_name, mimetype="text/html"
                )
            else:
                return "", errorcode
        else:
            return nfs_path, response_code


@ns.route("/purecn")
@api.response(200, "CSV file to purecn")
@api.response(400, "CSV file not found")
class GetPurecn(Resource):
    @api.expect(purecn_arguments, validate=True)
    def get(self):
        """
        Fetch ploidy and purity from purecn file
        """
        args = purecn_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_purecn_ctdna(
                nfs_path, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/get_purecn_max_val")
@api.response(200, "Update purecn max value into somatic and germline")
@api.response(400, "Purecn variant file not found")
class UpdateMaxPurecnVal(Resource):
    @api.expect(purecn_max_val_arguments, validate=True)
    def get(self):
        """
        Get Max PureCN value for each variant and update into somatic and germline table
        """
        args = purecn_max_val_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = update_pureCN_somatic_germline(
                nfs_path, args["sdid"], args["capture_id"], args["variant_type"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/get_curated_json")
@api.response(200, "Get the Json format")
@api.response(400, "Json file not found")
class GetJsonFile(Resource):
    @api.expect(json_urls_arguments, validate=True)
    # @api.marshal_with(ploturl_list)
    def get(self):
        """
        Returns Json format

        """
        args = json_urls_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_curated_json_file(
                nfs_path, proj_name, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/get_clinical_report_txt")
@api.response(200, "Get the clinical report txt format")
@api.response(400, "Clinical report txt file not found")
class GetClinicalReprotTxtFile(Resource):
    @api.expect(json_urls_arguments, validate=True)
    # @api.marshal_with(ploturl_list)
    def get(self):
        """
        Returns Json format

        """
        args = json_urls_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_clinical_report_file(
                nfs_path, proj_name, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/generate_json")
@api.response(200, "Generate the Json format")
@api.response(400, "Json file not found")
class GenerateJsonFile(Resource):
    @api.expect(json_urls_arguments, validate=True)
    def get(self):
        """
        Generate the json format
        """
        args = json_urls_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = generate_curated_json(
                nfs_path, proj_name, args["sdid"], args["capture_id"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/patient_info")
@api.response(200, "Fetch the patient information")
@api.response(400, "Json file not found")
class FetchPatientInformation(Resource):
    @api.expect(fetch_patient_info_arguments, validate=True)
    def post(self):
        """
        Fetch the patient information
        """
        args = fetch_patient_info_arguments.parse_args()
        proj_name = args["project_name"]
        result, errorcode = fetch_patient_info(
            proj_name, args["sdid"], args["capture_id"]
        )
        return result, errorcode


@ns.route("/generate_pdf")
@api.response(200, "Generate the PDF format")
@api.response(400, "Pdf file not found")
class GeneratePdfFile(Resource):
    @api.expect(json_urls_arguments, validate=True)
    def get(self):
        """
        Generate the pdf
        """
        args = json_urls_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = generate_curated_pdf(
                nfs_path, proj_name, sample_id, capture_id
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/report_pdf")
@api.response(200, "Fetch the PDF format")
@api.response(400, "Pdf file not found")
class FetchPdfFile(Resource):
    @api.expect(json_urls_arguments, validate=True)
    def get(self):
        """
        Fetch the pdf
        """
        args = json_urls_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = fetch_curated_pdf(
                nfs_path, proj_name, sample_id, capture_id
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/viewPdf")
@api.response(200, "view Pdf")
@api.response(400, "No Static plots found")
class ReportViewPDF(Resource):
    @api.representation("application/pdf")
    @api.expect(view_pdf_arguments, validate=True)
    def get(self):
        """
        View PDF
        """
        args = view_pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_pdf_file(
                nfs_path, args["sdid"], args["capture_id"], args["pdf_name"]
            )
            return send_file(
                result,
                attachment_filename=proj_name + "_" + args["sdid"] + "_report.pdf",
                mimetype="application/pdf",
            )
        else:
            return nfs_path, response_code


@ns.route("/fetch_pdf")
@api.response(200, "view Pdf")
@api.response(400, "No Static plots found")
class ReportFetchPDF(Resource):
    @api.representation("application/pdf")
    @api.expect(pdf_arguments, validate=True)
    def get(self):
        """
        View PDF
        """
        args = pdf_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_pdf_file2(
                nfs_path, args["sdid"], args["capture_id"]
            )
            return send_file(
                result,
                attachment_filename=args["sdid"] + "_report.pdf",
                mimetype="application/pdf",
            )
        else:
            return nfs_path, response_code


@ns.route("/update_curated_info")
@api.response(200, "Update Curation information based on ctdna")
@api.response(400, "No Static plots found")
class UpdateCuratedInfo(Resource):
    @api.expect(update_curated_info_arguments, validate=True)
    def post(self):
        """
        Update Curation information based on ctdna
        """
        args = update_curated_info_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = update_curated_info(
                nfs_path,
                proj_name,
                args["sdid"],
                args["capture_id"],
                args["ctdna_val"],
                args["ctdna_opt"],
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/send_json_mtbp")
@api.response(200, "Send json file to MTBP portal")
@api.response(400, "Json file not found")
class SendJsonMTBP(Resource):
    @api.expect(send_json_mtbp_arguments, validate=True)
    def post(self):
        """
        Send json file to MTBP portal
        """
        args = send_json_mtbp_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        user_name = args["user_name"]
        user_pwd = args["user_pwd"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = send_json_mtbp_portal(
                nfs_path, proj_name, sample_id, capture_id, user_name, user_pwd
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/fusion_inspector")
@api.response(200, "Read fusion inspector tsv")
@api.response(400, "/nfs is not mount locally no data found")
class TableRNAFusionInspector(Resource):
    @api.expect(table_qc_arguments, validate=True)
    def get(self):
        """
        Read fusion inspector
        """
        args = table_qc_arguments.parse_args()
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_fusion_inspector(
                nfs_path, args["sdid"], args["capture_id"], user_id, args["header"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/table/fusion_report")
@api.response(200, "Read fusion Report tsv")
@api.response(400, "/nfs is not mount locally no data found")
class TableRNAFusionReport(Resource):
    @api.expect(table_qc_arguments, validate=True)
    def get(self):
        """
        Read fusion Report
        """
        args = table_qc_arguments.parse_args()
        proj_name = args["project_name"]
        user_id = args["uId"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_table_fusion_report(
                nfs_path, args["sdid"], args["capture_id"], user_id, args["header"]
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/rna_qc_status")
@api.response(200, "check the RNA QC data available or not")
@api.response(400, "No data found")
class CheckRNAQCStatus(Resource):
    @api.expect(capture_arguments, validate=True)
    def get(self):
        """
        check the RNA QC data available or not
        """
        args = capture_arguments.parse_args()
        proj_name = args["project_name"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_rna_qc_status(nfs_path, args["sdid"])
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/cancer_type_list")
@api.response(200, "Fetched all cancer type")
@api.response(400, "No data found")
class GenomicProfileCancerTypeList(Resource):
    def get(self):
        """
        Fetch all cancer records
        ```

        ```
        """
        result, errorcode = get_cancer_type()
        return result, errorcode


@ns.route("/read_igv_json")
@api.response(200, "read the igv session json")
@api.response(400, "No data found")
class ReadIGVSessionJsonFile(Resource):
    @api.expect(igv_json_arguments, validate=True)
    def get(self):
        """
        Read the igv session json
        """
        args = igv_json_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        igv_type = args["igv_type"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = read_igv_session_json(
                nfs_path, sample_id, capture_id, igv_type
            )
            return result, errorcode
        else:
            return nfs_path, response_code


@ns.route("/read_igv_xml")
@api.response(200, "read the igv session XML")
@api.response(400, "No data found")
class ReadIGVSessionXMLFile(Resource):
    @api.expect(igv_xml_arguments, validate=True)
    def get(self):
        """
        Read the igv session xml
        """
        args = igv_xml_arguments.parse_args()
        file_path = args["file_name"]
        result, errorcode = read_igv_session_xml(file_path)
        return result, errorcode


@ns.route("/igvSession/files")
@api.response(200, "Read igv session files")
@api.response(400, "/nfs is not mount locally no data found")
class FetchIGVSessionFiles(Resource):
    @api.expect(igv_file_arguments, validate=True)
    def get(self):
        """
        Read igv session files
        """
        args = igv_file_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        file_name = args["file_name"]

        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            file_path = (
                nfs_path + "/" + sample_id + "/" + capture_id + "/IGVnav/" + file_name
            )
            # Set the mimetype to the appropriate content type
            mimetype = "text/plain"
            # mimetype = 'application/octet-stream'

            response = make_response(
                send_file(
                    file_path,
                    attachment_filename=file_name,
                    mimetype=mimetype,
                    as_attachment=False,
                    conditional=True,
                )
            )
            # response.headers['Content-Encoding'] = 'gzip'
            return response
            # return send_file(file_path, attachment_filename=file_name, mimetype=mimetype, as_attachment=False, conditional=True)
            # return send_file(file_path, as_attachment=True)
        else:
            return nfs_path, response_code


@ns.route("/igv-image/<string:variant_type>")
@api.response(200, "Franken IGV Image")
@api.response(400, "No IGV Image found")
class FrankenIGVStaticImage(Resource):
    @api.representation("image/png")
    @api.expect(igv_image_report_arguments, validate=True)
    def get(self, variant_type):
        """
        Returns staic IGV Image.
        """
        args = igv_image_report_arguments.parse_args()
        proj_name = args["project_name"]
        sample_id = args["sdid"]
        capture_id = args["capture_id"]
        chrom = args["CHROM"]
        start = args["START"]
        end = args["END"]
        nfs_path, response_code = fetch_nfs_path(proj_name)
        if response_code == 200:
            result, errorcode = get_igv_report_image(
                nfs_path, sample_id, capture_id, variant_type, chrom, start, end
            )
            return send_file(
                result, attachment_filename="frankenplot.png", mimetype="image/png"
            )
        else:
            return nfs_path, response_code
