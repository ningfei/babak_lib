
#define RASPOINT float
#define IMATRIXTYPE short
#define DIMXYTYPE float
#define PIXSIZETYPE float
#define BLOCK char
#define ATOMIC long

struct VARTYP {
	unsigned long	length;
	char		*data;
};

typedef struct VARTYP VARTYPE;

typedef struct {
	BLOCK 		su_id[4];
	short int 	su_uniq;
	char 		su_diskid;
	char 		prodid[13];
	BLOCK 		su_verscre[2];
	BLOCK 		su_verscur[2];
	unsigned long 	su_checksum;
	BLOCK		su_padding[85];
} SUITEDATATYPE;

typedef struct {
	BLOCK		ex_suid[4];
	short 		ex_uniq;	
	char		ex_diskid;
	unsigned short	ex_no;
	char		hospname[33];	
	short		detect;		
	int 		numcells;
	float		zerocell;
	float		cellspace;
	float		srctodet;
	float		srctoiso;
	short		tubetyp;	
	short		dastyp;	
	short		num_dcnk;
	short		dcn_len;
	short		dcn_density;	
	short		dcn_stepsize;
	short 		dcn_shiftcnt;
	int		magstrength;	
	char		patid[13];	// patient ID
	char		patname[25];	// patient name
	short		patage;
	short		patian;
	short		patsex;
	int		patweight;
	short 		trauma;
	char		hist[61];
	char		reqnum[13];
	int		ex_datetime;
	char		refphy[33];
	char		diagrad[33];
	char		op[4];
	char		ex_desc[23];
	char		ex_typ[3];
	short		ex_format;
	double		firstaxtime;
	char		ex_sysid[9];
	int		ex_lastmod;	
	short		protocolflag;
	char		ex_alloc_key[13];
	ATOMIC		ex_delta_cnt;
	BLOCK		ex_verscre[2];
	BLOCK		ex_verscur[2];
	unsigned long	ex_checksum;
	ATOMIC		ex_complete;
	ATOMIC		ex_seriesct;
	ATOMIC		ex_numarch;
	ATOMIC		ex_numseries;
	VARTYPE		ex_series;
	ATOMIC		ex_numunser;
	VARTYPE		ex_unseries;
	ATOMIC		ex_toarchcnt;
	VARTYPE		ex_toarchive;
	ATOMIC		ex_prospcnt;
	VARTYPE		ex_prosp;
	ATOMIC		ex_modelnum;
	ATOMIC		ex_modelcnt;
	VARTYPE		ex_models;	
	short		ex_stat;
	BLOCK		uniq_sys_id[16];
	BLOCK		service_id[16];
	BLOCK		mobile_loc[4];
	BLOCK		study_uid[32];
	short		study_status;
	BLOCK		ex_padding[516];
} EXAMDATATYPE;

typedef struct {
	BLOCK		se_suid[4];
	short		se_uniq;
	char		se_diskid;
	unsigned short	se_exno;
	short		se_no;
	int		se_datetime;
	int		se_actual_dt;
	char		se_desc[30];
	char		pr_sysid[9];
	char		pansysid[9];
	short		se_typ;
	short		se_source;
	short		se_plane;
	short		scan_type;
	int		position;
	int		entry;
	char		anref[3];
	float		lmhor;
	char		prtcl[25];
	short		se_contrast;
	char		start_ras;
	float		start_loc;
	char		end_ras;
	float		end_loc;
	short		se_pseq;
	short		se_sortorder;
	int		se_lndmrkcnt;
	short		se_nacq;
	short		xbasest;
	short		xbaseend;
	short		xenhst;
	short		xenhend;
	int		se_lastmod;
	char		se_alloc_key[13];
	ATOMIC		se_delta_cnt;
	BLOCK		se_verscre[2];
	BLOCK		se_verscur[2];
	float		se_pds_a;
	float		se_pds_c;
	float		se_pds_u;
	unsigned long	se_checksum;
	ATOMIC		se_complete;
	ATOMIC		se_numarch;
	ATOMIC		se_imagect;
	ATOMIC		se_numimages;
	VARTYPE		se_images;
	ATOMIC		se_numunimg;
	VARTYPE		se_unimages;
	ATOMIC		se_toarchcnt;

	VARTYPE		se_toarchive;
	float		echo1_alpha;
	float		echo1_beta;
	unsigned short	echo1_window;
	short		echo1_level;
	float		echo2_alpha;
	float		echo2_beta;
	unsigned short	echo2_window;
	short		echo2_level;
	float		echo3_alpha;
	float		echo3_beta;
	unsigned short	echo3_window;
	short		echo3_level;
	float		echo4_alpha;
	float		echo4_beta;
	unsigned short	echo4_window;
	short		echo4_level;
	float		echo5_alpha;
	float		echo5_beta;
	unsigned short	echo5_window;
	short		echo5_level;
	float		echo6_alpha;
	float		echo6_beta;
	unsigned short	echo6_window;
	short		echo6_level;
	float		echo7_alpha;
	float		echo7_beta;
	unsigned short	echo7_window;
	short		echo7_level;
	float		echo8_alpha;
	float		echo8_beta;
	unsigned short	echo8_window;
	short		echo8_level;
	BLOCK		series_uid[32];
	BLOCK		landmark_uid[32];
	BLOCK		equipmnt_uid[32];
	BLOCK		se_padding[588];
} SERIESDATATYPE;

typedef struct {
	BLOCK		im_suid[4];
	short		im_uniq;
	char		im_diskid;
	unsigned short	im_exno;
	short		im_seno; // series number
	short		im_no;	 // image number
	int		im_datetime;
	int		im_actual_dt;
	float		sctime;
	float		slthick;	// slice thichness (mm)
	IMATRIXTYPE	imatrix_X;	// X image matrix size (nx)
	IMATRIXTYPE	imatrix_Y;	// Y image matrix size (ny)
	float		dfov;		// fov X
	float		dfov_rect;	// fov Y if different
	DIMXYTYPE	dim_X;	// image dimension - X
	DIMXYTYPE	dim_Y;	// image dimension - Y
	PIXSIZETYPE	pixsize_X;	// X pixel size (dx) 
	PIXSIZETYPE	pixsize_Y;	// Y pixel size (dy) 
	BLOCK		pdid[14];
	char		contrastIV[17];
	char		contrastOral[17];
	short		contmode;
	short		serrx;
	short		imgrx;
	short		screenformat;
	short		plane;
	float		scanspacing;	// spacing between scans
	short		im_compress;
	short		im_scouttype;
	char		loc_ras;
	float		loc;
	RASPOINT	ctr_R;	// center R coord of plane image
	RASPOINT	ctr_A;	// center A coord of plane image
	RASPOINT	ctr_S;	// center S coord of plane image
	RASPOINT	norm_R;	// Normal R coord
	RASPOINT	norm_A;	// Normal A coord
	RASPOINT	norm_S;	// Normal S coord
	RASPOINT	tlhc_R;	// R coord of top left hand corner
	RASPOINT	tlhc_A;	// A coord of top left hand corner
	RASPOINT	tlhc_S;	// S coord of top left hand corner
	RASPOINT	trhc_R;	// R coord of top right hand corner
	RASPOINT	trhc_A;	// A coord of top right hand corner
	RASPOINT	trhc_S;	// S coord of top right hand corner
	RASPOINT	brhc_R;	// R coord of bottom right hand corner
	RASPOINT	brhc_A;	// A coord of bottom right hand corner
	RASPOINT	brhc_S;	// S coord of bottom right hand corner
	char		forimgrev[4];
	int		tr;	// TR usec
	int		ti;	// TI usec
	int		te;
	int		te2;
	short		numecho;
	short		echonum;
	float		tbldlta;
	float		nex;
	short		contig;
	short		hrtrate;
	int		tdel;
	float		saravg;
	float		sarpeak;
	short		monsar;
	short		trgwindow;
	float		reptime;	// cardiac repetition time
	short		imgpcyc;	
	short		xmtgain;
	short		rcvgain1;
	short		rcvgain2;
	short		mr_flip;
	int		mindat;
	short		cphase;
	short		swappf;
	short		pauseint;
	float		pausetime;
	int		obplane;
	int		slocfov;
	int		xmtfreq;
	int		autoxmtfreq;
	short		autoxmtgain;
	short		prescan_r1;
	short		prescan_r2;
	int		user_bitmap;
	short		cenfreq;
	short		imode;
	int		iopt;
	short		pseq;
	short		pseqmode;
	char		psdname[33];	// pulse seq. name
	int 		psd_datetime;	
	char		psd_iname[13];
	short		ctyp;
	char		cname[17];
	short		surfctyp;
	short		surfcext;
	int		rawrunnum;
	unsigned long	cal_fldstr;
	short		supp_tech;
	float		vbw;
	short		slquant;
	short		gpre;
	int		intr_del;
	float		user0;
	float		user1;
	float		user2;
	float		user3;
	float		user4;
	float		user5;
	float		user6;
	float		user7;
	float		user8;
	float		user9;
	float		user10;
	float		user11;
	float		user12;
	float		user13;
	float		user14;
	float		user15;
	float		user16;
	float		user17;
	float		user18;
	float		user19;
	float		user20;
	float		user21;
	float		user22;
	float		user23;
	float		user24;
	
	char		im_alloc_key[13];
	int		im_lastmod;
	BLOCK		im_verscre[2];
	BLOCK		im_verscur[2];
	int		im_pds_a;
	int		im_pds_c;
	int		im_pds_u;
	unsigned long	im_checksum;
	ATOMIC		im_archived;
	ATOMIC		im_complete;
	short		satbits;
	short		scic;
	short		satxloc1;
	short		satxloc2;
	short		satyloc1;
	short		satyloc2;
	short		satzloc1;
	short		satzloc2;
	short		satxthick;
	short		satythick;
	short		satzthick;
	short		flax;
	short		venc;
	short		thk_disclmr;	// slice thickness
	short		ps_flag;
	short		ps_status;
	short		image_type;
	short		vas_collapse;
	float		user23n;
	float		user24n;
	short		proj_alg;
	char		proj_name[13];
	float		x_axis_rot;
	float		y_axis_rot;
	float		z_axis_rot;
	int		thresh_min1;	
	int		thresh_max1;	
	int		thresh_min2;	
	int		thresh_max2;	
	short		echo_trn_len;	// echo train length for FSE
	short		frac_echo;
	short		prep_pulse;
	short		cphasenum;
	short		var_echo;
	char		ref_img;
	char		sum_img;
	unsigned short	img_window;
	short		img_level;
	int		slop_int_1;
	int		slop_int_2;
	int		slop_int_3;
	int		slop_int_4;
	int		slop_int_5;
	float		slop_float_1;
	float		slop_float_2;
	float		slop_float_3;
	float		slop_float_4;
	float		slop_float_5;
	char		slop_str_1[16];
	char		slop_str_2[16];
	short		scanactno;
	short		vasflags;
	float		vencscale;
	short		integrity;
	int		fphase;
	short		freq_dir;
	short		vas_mode;
	BLOCK		image_uid[32];
	BLOCK		sop_uid[32];
	short		dont_use_1;
	short		dont_use_2;
	short		dont_use_3;
	short		pscopts;
	short		asoffsetx;
	short		asoffsety;
	short		asoffsetz;
	short		unoriginal;
	short		interleaves;	// number of EPI shots
	short		effechospace;
	short		viewsperseg;
	short		rbpm;
	short		rtpoint;
	short		rcvrtype;
	float		dbdt;
	float		dbdtper;
	float		estdbdtper;
	float		estdbdtts;
	float		saravghead;
	float		neg_scanspacing;
	int		offsetfreq;
	unsigned long	user_usage_tag;
	unsigned long	user_fill_mapMSW;
	unsigned long	user_fill_mapLSW;
	float		user25;
	float		user26;
	float		user27;
	float		user28;
	float		user29;
	float		user30;
	float		user31;
	float		user32;
	float		user33;
	float		user34;
	float		user35;
	float		user36;
	float		user37;
	float		user38;
	float		user39;
	float		user40;
	float		user41;
	float		user42;
	float		user43;
	float		user44;
	float		user45;
	float		user46;
	float		user47;
	float		user48;
	int		slop_int_6;
	int		slop_int_7;
	int		slop_int_8;
	int		slop_int_9;
	BLOCK		mr_padding[32];
} MRIMAGEDATATYPE;

typedef struct {
	int		img_magic;
	int		img_hdr_length;	// header length
	int		img_width;
	int		img_height;
	int		img_depth;
	int		img_compress;
	int		img_dwindow;
	int		img_dlevel;
	int		img_bgshade;
	int 		img_overflow;
	int		img_undflow;
	int		img_top_offset;
	int		img_bot_offset;
	short		img_version;
	unsigned short	img_checksum;
	int		img_p_id;
	int		img_l_id;
	int		img_p_unpack;
	int		img_l_unpack;
	int		img_p_compress;
	int		img_l_compress;
	int		img_p_histo;
	int		img_l_histo;
	int		img_p_text;	
	int		img_l_text;	
	int		img_p_graphics;
	int		img_l_graphics;
	int		img_p_dbHdr;
	int		img_l_dbHdr;
	int		img_levelOffset;
	int		img_p_user;
	int		img_l_user;
	
} PixHdr;
