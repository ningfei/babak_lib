/**************************************************
   analyze.h     
   Copyright (c) 1996 Babak A. Ardekani
   ALL RIGHTS RESERVED
***************************************************/

/* $Log: analyze.h,v $
 * Revision 3.1.1.1  1996/11/09 13:15:16  babak
 * Akita Univ. version
 *
 * Revision 3.1  1996/11/07 12:47:30  babak
 * first check-in
 * */

/* $Id: analyze.h,v 3.1.1.1 1996/11/09 13:15:16 babak Exp babak $ */

struct header_key 
{
	int sizeof_hdr;
	char data_type[10];
	char db_name[18];
	int extents;
	short session_error;
	char regular;
	char hkey_un0;
};

struct image_dimension
{
	short dim[8];
	short unused8;
	short unused9;
	short unused10;
	short unused11;
	short unused12;
	short unused13;
	short unused14;
	short datatype;
	short bitpix;
	short dim_un0;
	float pixdim[8];
	float funused8;
	float funused9;
	float funused10;
	float funused11;
	float funused12;
	float funused13;
	float compressed;
	float verified;
	int glmax,glmin;
};

struct data_history
{
	char descrip[80];
	char aux_file[24];
	char orient;
	char originator[10];
	char generated[10];
	char scannum[10];
	char patient_id[10];
	char exp_date[10];
	char exp_time[10];
	char hist_un0[3];
	int views;
	int vols_added;
	int start_field;
	int field_skip;
	int omax, omin;
	int smax, smin;
};

struct dsr
{
	struct header_key hk;	
	struct image_dimension dime;
	struct data_history hist;
};
