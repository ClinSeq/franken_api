CREATE TYPE cur_user_status AS ENUM ('0','1', '2', '-1');
CREATE TYPE cur_user_status_desc AS ENUM ('0 -New user','1 - Active User', '2 - Suspended User', '-1 Delete User');

CREATE TYPE cur_role_status AS ENUM ('0', '1', '-1');
CREATE TYPE cur_role_status_desc AS ENUM('0 - Normal User', '1 - Admin', '-1 Super Admin');

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
    pdf_report proj_status NOT NULL DEFAULT '0';
    sort_order integer NOT NULL DEFAULT 0,
    created_on timestamptz
);

/* INSERT INTO public.cur_projects_t( p_id, project_name, prefix_name, nfs_path, proj_status, created_on) 	VALUES (DEFAULT, 'IPCM', 'iPCM', '/sdata/IPCM/autoseq-output', '1', NOW()); */