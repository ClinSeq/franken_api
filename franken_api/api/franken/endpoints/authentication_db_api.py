import logging
from flask import current_app
from flask import request, send_file, make_response, send_from_directory
from flask_restplus import Resource
from franken_api.api.franken.parsers import auth_login_arguments, project_list_arguments, auth_register_arguments
from franken_api.api.restplus import api
from franken_api.api.franken.business import  login_validate, fetch_project_list, form_registation, fetch_all_project_list
from franken_api.api.franken.serializers import * 

log = logging.getLogger(__name__)
au = api.namespace('auth', description='Authentication Database API')

@au.route('/')
@api.response(200, 'Check API status')
@api.response(400, '/nfs is not mount locally no data found')
class AuthenticationStatus(Resource):
    def get(self):
        """
        Authentication API status check
        ```

        ```
        """
        return 'Working', 200

@au.route('/login')
@api.response(200, 'login successfully')
@api.response(400, 'login failed')
class AuthenticationLoginValidation(Resource):
    @api.expect(auth_login_arguments, validate=True)
    def post(self):
        """
        Login
        ```

        ```
        """
        args = auth_login_arguments.parse_args()
        email_id = args['email_id']
        passwd = args['passwd']
        result, errorcode = login_validate(email_id, passwd)
        return result, errorcode


@au.route('/project_list')
@api.response(200, 'Fetch the project list')
@api.response(400, 'project list not found')
class AuthenticationProjectList(Resource):
    @api.expect(project_list_arguments, validate=True)
    def post(self):
        """
        Fetch all the project list 
        ```

        ```
        """
        args = project_list_arguments.parse_args()
        project_ids = args['project_ids']
        result, errorcode = fetch_project_list(project_ids)
        return result, errorcode

@au.route('/register')
@api.response(200, 'New user registration successfully')
@api.response(400, 'registration failed')               
class RegisterValidation(Resource):
    @api.expect(auth_register_arguments, validate=True)   
    def post(self):
        """
        New user registration
        ```

        ```
        """
        args = auth_register_arguments.parse_args()
        first_name= args['first_name']
        last_name = args['last_name']
        email_id= args['email_id']
        pwd = args['passwd']
        project_access = args['proj_ids']
        result, errorcode = form_registation(first_name, last_name, email_id, pwd, project_access)
        return result, errorcode 

@au.route('/all_project_list')
@api.response(200, 'Fetch the project list')
@api.response(400, 'project list not found')
class AuthenticationAllProjectList(Resource):
    # @api.expect(project_list_arguments, validate=True)
    def get(self):
        """
        Fetch all the project list 
        ```

        ```
        """
        # args = project_list_arguments.parse_args()
        # project_ids = args['project_ids']
        result, errorcode = fetch_all_project_list()
        return result, errorcode
