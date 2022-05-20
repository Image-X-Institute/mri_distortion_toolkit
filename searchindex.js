Search.setIndex({docnames:["Harmonic_sanity_checks","Harmonic_stability_checks","Preamble","_static/themes/README","code_docs","developer_notes","examples","field_calculation","fit_spherical_harmonics","index","marker_extraction","marker_matching","marker_position_errors","reporting"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":5,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":3,"sphinx.domains.rst":2,"sphinx.domains.std":2,"sphinx.ext.todo":2,sphinx:56},filenames:["Harmonic_sanity_checks.md","Harmonic_stability_checks.md","Preamble.md","_static/themes/README.md","code_docs.rst","developer_notes.md","examples.rst","field_calculation.md","fit_spherical_harmonics.md","index.rst","marker_extraction.md","marker_matching.md","marker_position_errors.md","reporting.md"],objects:{"MRI_DistortionQA.FieldAnalysis":[[4,1,1,"","SphericalHarmonicFit"]],"MRI_DistortionQA.FieldAnalysis.SphericalHarmonicFit":[[4,2,1,"","plot_cut_planes"],[4,2,1,"","plot_harmonics_pk_pk"],[4,2,1,"","print_key_harmonics"]],"MRI_DistortionQA.FieldCalculation":[[4,1,1,"","ConvertMatchedMarkersToBz"]],"MRI_DistortionQA.MarkerAnalysis":[[4,1,1,"","MarkerVolume"],[4,1,1,"","MatchedMarkerVolumes"]],"MRI_DistortionQA.MarkerAnalysis.MarkerVolume":[[4,2,1,"","export_to_slicer"],[4,2,1,"","perturb_marker_positions"],[4,2,1,"","plot_3D_markers"],[4,2,1,"","save_dicom_data"]],"MRI_DistortionQA.MarkerAnalysis.MatchedMarkerVolumes":[[4,2,1,"","plot_3D_markers"],[4,2,1,"","plot_compressed_markers"]],"MRI_DistortionQA.Reports":[[4,1,1,"","DefaultTestSuite"],[4,1,1,"","Elekta_Distortion_tests"],[4,1,1,"","MRI_QA_Reporter"],[4,1,1,"","TG284_Tests"]],"MRI_DistortionQA.Reports.DefaultTestSuite":[[4,2,1,"","test_distortion_less_than_2mm_within_r_100"],[4,2,1,"","test_distortion_less_than_4mm_within_r_150"]],"MRI_DistortionQA.Reports.Elekta_Distortion_tests":[[4,2,1,"","test_distortion_less_than_1_mm_within_r_100_mm"],[4,2,1,"","test_distortion_less_than_2_mm_within_r_150_mm"],[4,2,1,"","test_distortion_less_than_4_mm_within_r_200_mm"]],"MRI_DistortionQA.Reports.MRI_QA_Reporter":[[4,2,1,"","write_html_report"]],"MRI_DistortionQA.Reports.TG284_Tests":[[4,2,1,"","test_distortion_less_than_1mm_within_r_100"],[4,2,1,"","test_distortion_less_than_2mm_within_r_200"]],"MRI_DistortionQA.calculate_harmonics":[[4,3,1,"","calculate_harmonics"]],"MRI_DistortionQA.utilities":[[4,1,1,"","bcolors"],[4,3,1,"","build_dicom_affine"],[4,3,1,"","compare_recon_report_with_ground_truth_report"],[4,3,1,"","convert_cartesian_to_spherical"],[4,3,1,"","convert_spherical_harmonics"],[4,3,1,"","convert_spherical_to_cartesian"],[4,3,1,"","dicom_to_numpy"],[4,3,1,"","generate_legendre_basis"],[4,3,1,"","get_all_files"],[4,3,1,"","get_dicom_data"],[4,3,1,"","get_gradient_spherical_harmonics"],[4,3,1,"","plot_MarkerVolume_overlay"],[4,3,1,"","plot_compressed_MarkerVolumes"],[4,3,1,"","reconstruct_Bz"]],MRI_DistortionQA:[[4,0,0,"-","FieldAnalysis"],[4,0,0,"-","FieldCalculation"],[4,0,0,"-","MarkerAnalysis"],[4,0,0,"-","Reports"],[4,0,0,"-","calculate_harmonics"],[4,0,0,"-","utilities"]]},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","function","Python function"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:function"},terms:{"0":[0,1,4,7,8,12,13],"004159":12,"01":[0,1,8],"02":[0,1],"023":8,"03":[0,1],"032":0,"04":[0,1,7,8,10,11,13],"048565":12,"05":[0,1,10,11],"06":[0,1],"07":[0,1],"08":[0,1],"087395":12,"09":[0,1],"09255":12,"1":[1,4,6,8,10,11,12],"10":[0,1,4],"100":[0,1,10,12,13],"1000":4,"1002":4,"104":0,"11":[0,1,11,12],"111853":12,"112788":12,"12":[0,1],"13":[0,1],"134":0,"138":0,"139105":12,"14":[0,1],"14695":4,"14764":4,"15":[0,1],"150":[1,4,8,10,13],"16":[0,1],"165":[0,1],"17":[0,1],"18":[0,1],"19":[0,1],"1mm":1,"2":[1,4,6,8,12],"20":[0,1,4],"200722":12,"20220428":[0,1],"203":4,"20data":[0,1],"20experiment":[0,1],"20linac":[0,1],"20mr":[0,1],"21":0,"22":0,"23":0,"24":0,"25":0,"256":4,"258656":12,"27":0,"275514":12,"288829":12,"2_projectdata":[0,1],"2d":4,"2research":[0,1],"3":[0,1,10,12],"30":0,"300":[8,11],"309114":12,"336":[0,1],"367208":12,"368377":12,"393177":12,"397762":12,"3d":4,"4":[0,1,8,12],"40":0,"400":4,"401":8,"416738":12,"42":4,"43":0,"440425":12,"444":0,"448163":12,"465":0,"4x4":4,"5":[0,1,4,12],"50":[0,1],"508799":12,"537472":12,"559205":12,"586933":12,"587143":12,"5emr":[0,1],"5etest":[0,1],"6":[0,1,4,12],"62996":12,"63493":12,"6603039901":[0,1],"689":0,"693":8,"7":[0,1],"709":0,"720":0,"721":0,"795":8,"796672":12,"8":[0,1,4,8],"806967":12,"86":8,"881676":12,"9":[0,1,4],"900554":12,"921969":12,"957465":12,"96696":12,"\u03bct":[0,8],"abstract":4,"boolean":4,"case":[0,2,4,6,8,10,11],"class":[4,11,12,13],"default":[3,4],"do":[2,4,6,7,8,12,13],"export":[4,11],"final":[4,13],"float":4,"function":[1,4,8],"import":[0,1,4,7,8,10,11,13],"int":4,"long":1,"new":[7,8,10,11,13],"return":[1,4,13],"short":[10,11],"static":12,"switch":0,"throw":4,"true":[0,1,4,10,12,13],"try":[1,4,11,13],"while":11,A:[4,10,11],As:10,At:[4,10],But:[2,11],By:4,For:[4,9,10,11,12],If:[0,4,8,10,11,13],In:[1,2,4,8,10,11],Is:6,It:[1,4],NOT:4,No:10,Not:11,Such:4,The:[0,1,2,4,6,9,10,11,12,13],There:[11,13],These:[2,11],To:[10,11,12,13],With:0,_:4,__file__:0,__resourc:5,_extract_data_from_matchedmarkervolum:13,_residual_pk_pk:1,a00:4,a10:0,a11:[0,8],a20:0,a30:0,a31:0,a50:0,a51:0,a60:0,a80:0,a_1_0:1,a_1_1:1,a_3_0:1,a_3_1:1,a_5_0:1,a_5_1:1,aapm:4,aarnet:[0,1,7,8,10,11],abl:10,about:[4,10,11,12],abov:[1,10,11],abs_di:13,absolut:12,accept:4,access:[4,13],accident:11,accord:4,achiev:11,acquisit:[4,12],actual:[12,13],ad:[4,6],add:[4,10,13],add_arrow:4,addcolorbar:4,addit:[4,11],address:13,affin:4,again:[0,4,11],algorith:4,algorithm:11,align:[4,10,11],all:[0,2,4,8,9,11,12],all_scan:[0,1],allow:[4,7,8,11],allowdoublematch:4,alon:11,alreadi:[4,12],also:[0,1,4,10,11,12,13],although:[4,10],alwai:[4,11],am:13,ammount:1,an:[0,1,4,6,7,10,12,13],analys:11,analysi:[0,4,12],ani:[3,4,9,11],answer:[0,10,11],anterior:0,anywai:11,ap:[0,11,12],ap_hf:12,ap_pa:0,ap_rl:12,ap_volum:12,apart:10,app:[0,1],appear:[4,8,11],append:1,appli:10,applic:4,approach:[4,11],appropri:[8,13],ar:[0,1,2,4,8,10,11,12,13],arbitrari:4,aren:[0,4,10],aris:2,around:[1,12],arrai:[1,4],arrang:1,array_to_normalis:1,arrow:4,assess:[1,5],assessharmonicpk_pk:4,associ:3,assumpt:4,attempt:[4,11],au:[0,1,4,7,8,10,11],autom:8,automat:[4,10,11,13],automatchmark:4,avail:[4,11],ax:1,axi:10,azimuth:4,b0:[1,4,7,10,11,12,13],b0_data:1,b0_harmon:[1,4,8,13],b0_harmonics_ap:0,b0_harmonics_gt:1,b0_harmonics_hf:0,b0_harmonics_rl:0,b11:0,b31:0,b51:0,b:[1,4],b_1_1:1,b_3_1:1,b_5_1:1,b_field:1,b_gx:[1,8],b_gy:[1,8],b_gz:[1,8],back:[4,10,11],back_volum:[0,1],bad:[4,12,13],badg:5,bandwidth:4,bar:4,barplot:4,base:[4,10,11],basic:[4,6,9],bcolor:4,becaus:[2,4,10,11,12,13],been:13,befor:4,being:[0,4],below:[0,1,4,7,8,9,10,11,12,13],best:[4,10,12],better:[0,4,13],between:[0,4,11,12],bfield:[7,8],biggest:1,bit:[4,10,12],block:10,blur:4,bool:[4,13],both:[4,11,13],bottom:0,bound:4,brendan:[7,8,10,11,13],browser:[4,13],build:[4,10],build_dicom_affin:4,built:[10,11],bz:[1,4,8],bz_field:7,bz_recon:4,c:[4,7,8,10,11,13],cad:11,calc:4,calcualt:4,calcul:[0,1,4,6,8,9,11],calculate_harmon:[0,1,8,9],call:[1,4,7,8,10,11,13],can:[0,1,2,4,7,8,9,10,11,12,13],cardin:[4,8],carri:0,cartesian:4,center:[10,11,12],central:12,centroid:[1,8,11],chain:4,challeng:10,chang:[4,10],character:[1,9,11],check:[4,8,10,11,12],chemic:10,circular:4,clean:4,cleaner:4,clearli:12,click:9,close:[8,10,12],cloudstor:[0,1,7,8,10,11],code:[2,6,7,9,10,11,13],coil:[7,8,11],collect:4,color:4,column:[1,4,7,8,10,13],com:[3,4],combin:10,come:11,comfort:1,commonli:4,compar:[1,4,10,11],compare_recon_report_with_ground_truth_report:[1,4],complet:[0,2,13],compon:12,componet:4,compress:4,concern:[0,12],confid:[4,8,10],consist:11,constant:4,constitut:2,conta:13,contain:[4,8,11],content:11,contrast:10,contribut:4,control:[4,8,10],conveni:4,convent:[4,13],convert:[4,7,10],convert_cartesian_to_spher:4,convert_spherical_harmon:4,convert_spherical_to_cartesian:4,convertmatchedmarkerstobz:[0,1,4,7],coord:4,coordiant:4,coordin:[1,4],coordinatematrix:4,copi:[1,4,7,8,10,11,13],coron:0,correct:[4,10],correct_fat_water_shift:[0,1,4,10],correct_fw:[0,1],correspond:[4,8,11],cover:8,creat:[1,2,4,6,7,8,10,13],creation:13,criterion:4,csv:[1,4,7,8,11,13],ct:[0,1,8,11],ct_volum:[0,1],custom:6,customtestsuit:13,cut_off:[4,8],cutoff:4,cutoff_point:[0,1,4],d42ker:3,d:4,dark:4,data:[0,1,2,4,6,7,8,10,12],data_loc:[0,1,7,8,10,11,13],datafram:[4,8,10],datapoint:13,dataset:[4,10,12],dcm:4,deal:12,debug:4,decomposit:4,deepcopi:1,def:[1,13],defaulttestsuit:[4,13],defin:[4,9,11],delet:[4,10],demonstr:[2,10,11,13],dens:1,depend:11,deriv:4,describ:4,design:[2,11,13],detail:[2,11],determin:4,develop:10,diagram:9,dicom:[2,4,6,10],dicom_affin:4,dicom_data:[1,4,7,8,11,13],dicom_data_loc:[7,8,13],dicom_to_numpi:4,dicomfil:4,dict:4,dictionari:4,did:8,differ:[4,10,11,12,13],difficult:10,dir:[0,1],direct:[4,10,11,12,13],directli:[6,11],directori:[5,10,13],discount:[4,11],disort:11,displai:[4,13],distort:[2,4,8,10,11,13],distorted_volum:[4,8,11],distorted_volume_rev:[4,8,11],distorteddata:4,distortoin:11,distribut:[4,12],doc:11,docsrc:5,document:[2,9,11,13],doe:[4,7],doesn:4,doi:4,domin:[0,1,4,8],don:[11,12,13],doubl:4,download:[7,8,10,11,13],drawn:4,drop:12,dsv:10,due:10,dure:4,e:[0,1,4],each:[1,2,4,7,8,9,10,11,12],earlier:4,easi:[6,10],easier:4,easili:11,edg:2,edit:10,edu:[0,1,7,8,10,11],effect:[10,11],either:4,elekta:4,elekta_distortion_test:4,element:4,elev:4,els:13,elsewher:4,enabl:[4,9],encod:[0,10,11,12],encompass:4,end:[2,9],enough:12,ensur:4,enter:4,entir:12,enumer:1,epdf:4,error:[0,1,4,13],especi:1,essenti:[4,11],estim:[0,1,7,10,11,12,13],etc:10,even:[0,1],ever:0,everi:[2,4,8,11],everyth:4,exactli:11,exampl:[2,4,7,9,13],excel:10,except:[0,1,4],excit:0,exhibit:0,exist:4,expect:[0,1,4,8,11],explain:[2,6],explait:11,export_to_slic:[0,4,10,11],express:8,extens:[4,11],extra:11,extract:[1,4,6,9,11,12],extrem:11,eyebal:12,f:[5,11],fail:[6,10,13],fals:[0,1,4,8,10,11,13],fancier:12,far:[0,1,8],fat:[4,6],fat_shift_direct:[0,1,4,10],featur:[4,10],feet:0,few:[1,11],fh:0,field:[0,1,4,6,8,9,11],field_calcul:7,fieldanalysi:[0,1,8,9],fieldcalcul:[0,1,7,9],fielddata:8,fig:1,figsiz:1,figur:[4,8,13],file:[0,1,4,7,8,10,11,13],file_extens:4,filedataset:4,fileid:[0,1],filenam:4,filestoreadin:4,fill:10,filter:4,find:[1,4,6],fine:[4,8],firefox:13,first:[1,2,10],firstli:11,fit:[0,1,4,7,8,13],fit_harmon:8,fly:4,folder:4,follow:[1,5,8,11],formal:4,format:10,former:[2,4,11],forward:[4,10],forward_volum:[0,1],forward_volume_perturb:1,found:4,fov:4,frame:[4,10,13],free:10,freq_encode_direct:0,frequenc:[10,11],from:[0,1,2,3,4,5,6,7,8,9,10,12,13],full:[4,8,11],fulli:4,further:[0,10,11],futur:2,g:[0,1,4],g_x:8,g_x_harmon:[1,4,8,13],g_x_harmonics_ap:0,g_x_harmonics_gt:1,g_x_harmonics_hf:0,g_x_harmonics_rl:0,g_y:8,g_y_harmon:[1,4,8,13],g_y_harmonics_ap:0,g_y_harmonics_gt:1,g_y_harmonics_hf:0,g_y_harmonics_rl:0,g_z:[1,8],g_z_harmon:[1,4,8,13],g_z_harmonics_ap:0,g_z_harmonics_gt:1,g_z_harmonics_hf:0,g_z_harmonics_rl:0,gama:4,gaussian_image_filter_sd:[0,1,4],gener:[1,2,4,5,11,12,13],generate_legendre_basi:4,geometr:[4,9,11],get:[0,1,4,10,11,12,13],get_all_fil:4,get_dicom_data:4,get_gradient_spherical_harmon:4,github:3,give:[4,8,12],given:[4,10],gnl:11,go:[0,1,10,12],goam2:[0,1],good:[1,4,8,10,11,12,13],got:0,gradient:[0,1,4,6,7,8,10],gradient_harmon:[1,4,13],gradxdata:[1,8],gradydata:[1,8],gradzdata:[1,8],gratifi:8,gre_cor_lr:[0,1],gre_cor_rl:[0,1],gre_cor_rl_330:[0,1],gre_sag_fh:[0,1],gre_sag_fh_330:[0,1],gre_sag_hf:[0,1],gre_sag_hf_330:[0,1],gre_tran_ap_large_bw:[0,1],gre_tran_pa_large_bw:[0,1],gre_trans_ap_330:[0,1,7,8,10,11,13],gre_trans_ap_330_f_reset:[0,1],gre_trans_ap_reshim_refreq:[0,1],gre_trans_pa:[0,1],gre_trans_pa_330:[0,1,10,11],gre_trans_pa_reshim_refreq:[0,1],greater:4,grip:12,ground:[4,8,11,13],ground_truth_report:4,ground_truth_volum:[1,4,8,11],groundtruthdata:4,guarante:6,guess:11,gui:10,gx:[1,4],gx_harmon:4,gy:[1,4],gy_harmon:4,gz:4,gz_harmon:4,ha:[4,11],had:3,handi:4,handl:[2,4,6],hard:11,harmon:[4,6,7,9,11],have:[0,1,4,7,8,9,10,11,12,13],haven:13,head:0,help:[4,12],here:[0,1,2,3,4,8,10,11,12,13],hf:[0,12],hf_fh:0,hf_rl:12,high:11,higher:4,highest:4,highli:11,hinder:1,homogen:[4,10,13],honest:13,how:[2,4,6,7,10,13],howev:[1,4,8,11,12],html:[4,13],http:[0,1,3,4,7,8,10,11],hypothes:1,i:[0,1,3,4,6,12,13],idea:4,ideal:11,identifi:4,ignor:[4,5,8,11],ima:4,imag:[1,2,4,9,10,11,12],image_s:4,imagearrai:4,imextens:4,imshow:4,includ:[2,4,11,13],incorpor:6,inde:12,independ:[0,4,9,11],index:9,index_col:[7,8,13],indic:[0,4],influenc:12,info:4,inform:[4,11],inher:1,inhomogen:11,init:5,initi:[0,8],inner:[11,12],inner_most:4,input:[4,8,9],input_bz_data:4,input_data:4,input_data_path:4,input_format:4,inputcoord:4,inputdata:4,insid:4,instanc:[4,10],instanti:10,instruct:11,intend:2,intepret:6,interact:13,interc:4,interest:[1,4],interfac:10,intermidi:7,interpret:4,introduc:[1,12],involv:8,isn:7,issu:[6,11,12],item:11,itself:[1,7],j:4,journei:13,json:[0,1,4,6,7,8,10,13],json_fil:10,json_volum:10,jump:2,just:[4,8,11],k:4,keep:1,kei:11,keyerror:1,know:[4,8,10,11],knowledg:[10,11],larg:1,larger:[1,12],largest:12,later:[7,10,11],latter:[4,11],least:4,left:[0,4],legend:[1,4],legendr:4,legendre_basi:4,leonwong0609:3,less:[1,4],let:[1,10,12,13],librari:10,licens:3,life:4,light:4,like:[0,8,10,12],limit:[4,11,13],linac:[0,1],line:8,linear:11,linspac:1,linux:4,list:[1,4,11],littl:10,load:[4,7,8],localiser_gr:[0,1],locat:[10,11],log:[10,11],look:[0,2,10,11,12,13],lot:[0,4,8,13],low:[1,4,8,10],lower:4,lr:0,ly:4,m:5,magnet:[0,7,9],magneticfield:[1,7],mai:[2,4,10],main:0,major:[1,6],make:[1,4],mani:[4,10],manual:4,map:4,marker:[0,1,4,6,7,8,9,13],marker_match:11,marker_s:4,marker_size_lower_tol:4,marker_size_upper_tol:4,marker_volum:[4,10],marker_volume_back:10,marker_volume_forward:10,markeranalysi:[0,1,8,9,10,11],markercentroid:[10,12],markerextractionexampl:10,markermarker_size_tol:4,markervolum:[0,1,4,6,7,8,10,13],markervolumelist:4,match:[0,1,4,6,7,8,9,10,13],matched_mark:[1,7,11,13],matched_volum:[7,11],matched_volume_no_ref:11,matched_volume_with_rev_data:11,matchedcentroid:[1,4,11],matchedmarkervolum:[0,1,4,11,12,13],mathwork:4,matlab:4,matplotlib:[1,4],matric:4,matrix:4,matter:[4,10,11,12],max:[4,13],maximum:4,mayb:12,me:4,mean:[4,10],measur:[11,13],median:4,median_marker_s:4,mess:11,method:[1,4,5,10,11],middl:10,might:11,min:4,minimum:4,minut:11,mismatch:4,mm:[1,4,8,12],model:4,modul:[9,11],modular:2,moment:4,more:[0,1,2,4,10,11,12,13],most:[1,4,8,10,11,12,13],mostli:4,motion:4,move:[2,4,8,10,11,13],mp:4,mr:[0,1,4,7,8,10,11,13],mr_distortionqa:5,mr_qa_report:4,mri:[0,1,4,11],mri_distortion_qa_sample_data:[7,8,10,11,13],mri_distortionqa:[0,1,4,7,8,10,11,13],mri_qa_recon_report:4,mri_qa_report:[1,4,13],mri_qa_tutori:10,mrk:[0,1,8,10,11],much:[0,12],multipl:10,must:4,n:4,n_markers_expect:[0,1,4],n_order:[1,4,8],name:[4,8],ncol:1,ncoord:4,ndarrai:4,nearest:4,necessari:11,need:[2,4,8,10,11,12],neither:3,nerd:8,never:11,next:[2,6,7,10],nibabel:4,nois:[1,4,10],non:[0,1,11],none:[1,4,8],normal:4,normalis:4,normalise_array_columnwis:1,note:[2,4,8,10,11],noth:4,notic:13,now:[1,8,13],np:[1,4,10],nrow:1,number:[4,12],numer:4,numpi:[1,4,10],nyquist:4,o:5,object:[4,10,11],obtain:11,occasion:13,occur:10,off:[0,11],offset:[4,10,11],oil:10,ok:[0,1,13],onc:[4,7,10],one:[1,4,8,10,11,13],ones:0,onli:[2,4,7,11],onlinelibrari:4,onto:8,oper:4,option:[4,10,11],order:[2,4,8],origin:[0,1,4],other:[1,4,8,9,10],otherwis:4,otsu:4,our:[1,11,13],out:[0,2,4,10,11,13],outlier:[1,11],output:[6,9],output_fold:4,output_format:4,outsid:4,over:[4,8,11,12],overal:0,overlai:[4,12],overlaid:4,overwritten:4,own:13,pa:[0,11],pa_ap:0,page:9,pair:[0,10,11],panda:[4,7,8,10,13],pandas:4,pandas_volum:10,paramet:[4,11],parent:0,parmet:11,part:[4,8,9,10,13],particularli:[7,11],pass:[4,6],path:[0,1,4,7,8,10,11,13],path_to_dicom:4,pathlib:[0,1,4,7,8,10,11,13],pathtodata:4,patient:11,pd:[7,8,10,13],peak:8,percentag:1,perfect:[0,8],perform:[4,11],perturb:[0,4,8],perturb_marker_posit:[1,4],perturbed_report:1,perurb:1,phantom:[1,10,11,12,13],phase:[0,11],phase_encode_direct:0,physic:[0,1],pi:4,pixel:[4,12],pk:[0,1,4,8],place:[4,11],plane:[4,8],pleas:[2,10],plot:[1,4,8,10,11,13],plot_3d_mark:[4,10,11],plot_compressed_mark:4,plot_compressed_markervolum:4,plot_cut_plan:[4,8],plot_harmonics_pk_pk:[4,8],plot_markervolume_overlai:[4,10],plotli:4,plt:1,plu:[0,1,7,8,10,11],point:[2,4,8,10,13],port:4,posit:[1,4,7,10,11,13],possibl:2,posterior:0,preambl:[6,9],predict:[4,10],present:[4,11],pretti:[4,8,13],previosli:[7,8,13],previou:[7,11,13],previous:[4,8,11],principl:4,print:[4,8,11],print_key_harmon:[4,8],prior:[4,11],prj:[0,1],probabl:[1,8,12],problem:11,process:[4,6,10],produc:[4,10],promin:11,provid:[2,8,10,11],pull:[4,11],purpos:[0,1,4,11,13],put:[4,8,10,11],py:[7,8,11,13],pydicom:4,pyplot:[1,4],pytest:5,python:[7,10,11],qa:2,quantifi:[11,12],quantifyfit:4,quantiti:4,quick:[0,4,10],quit:[0,1,4,10,12],r:[0,1,4,7,8,10,11,12,13],r_max:[0,1,4,8,11,13],r_min:4,r_outer:[1,4,8,10,13],radial:[1,4],radii:4,radiu:4,rand:10,random:[1,4,10,12],random_perturb:4,raw:4,reaction:4,read:[4,10,11],read_csv:[7,8,13],read_fil:4,readi:8,readout:11,realli:[7,13],reason:[0,1,4],receiev:[1,8],recomend:4,recommend:[2,11],recon_coord:4,recon_coords_cartesian:4,recon_report:4,reconstruct:[0,1,4,6,8],reconstruct_bz:4,reconstuct:4,recontruct:13,record:1,recreat:11,ref:4,refer:4,referencemark:[1,4,11],rel:[1,4,8],reli:[1,2],reliabl:10,remain:6,remark:0,remov:4,renam:[1,8],repeat:1,report:[0,1,2,6,8,9,11],repres:1,request:[4,11],requir:[4,11],reset_index:12,residu:[0,8],resolut:[4,11],reson:9,result:[0,1,4],return_xyz:4,reusabl:4,revers:[4,6,7,10],reversegradientdata:[1,4,11],rid:11,right:0,rigid:11,rl:[0,12],rl_lr:0,root:5,routin:13,row:4,rpl:[0,1],run:[4,5,13],s:[1,4,7,8,10,11,12,13],sagitt:0,sai:[0,10],same:[1,4,8],sampl:[1,8,10,11],saniti:[4,11],satisfi:[4,10],save:[4,7,8,10,11,13],save_dicom_data:[0,4,7,10,11],save_path:[0,4],scale:4,scan:10,scanner:[10,11],scatter:4,screen:[4,8],screenshot:10,script:4,search:[4,9],second:[1,2,8],section:[2,11,13],see:[1,4,8,10,11,12],seen:[4,13],segment:[1,4],self:13,send:11,sensibl:0,sensit:[1,4],seper:11,sequenc:11,seri:[4,10,12],seriou:12,serv:4,set:[0,4,10,11,12],set_titl:1,set_xlabel:1,set_ylabel:1,setup:4,sever:12,shape:4,share:[0,1],shift:[4,6],should:[1,2,4,6,10,12,13],show:[1,7,8,10,11,12],show_plot:4,shown:4,siemen:11,sign:[0,10],signal:[4,10],similar:0,similarli:4,simpl:[2,6],simplest:10,simpli:[4,11],sinc:[1,10,11],singl:4,sit:13,situat:[2,8,10,11],size:[1,4,10,12],slice:[0,4,10,11,12],slice_direct:0,slicer:[4,10],slicer_centroid:[4,8],slower:4,small:4,snr:[1,10],so:[0,4,8,9,10,11,12,13],softwar:[10,11],solv:11,some:[2,4,8,10,11,12,13],someth:[10,12],sometim:4,somewher:10,sorting_method:[1,4],sourc:[2,11],space:4,specif:11,specifi:4,speed:11,sphere:[1,4,8,13],spheric:[1,4,6,7,9,11,13],sphericalharmonicfit:[0,1,4,8],split:2,spot:4,spuriou:4,squeez:[7,8,13],stabil:[0,1,4],stabl:[0,1],stai:1,standard:11,start:[1,2,4,11],step:[0,2,4,6,7,9,10,13],still:[0,4,10],store:[4,10],str:4,string:[4,13],strong:8,strongli:8,style:4,subject:10,subplot:1,subset:4,substanti:4,succes:11,success:4,suffici:8,suggest:[0,1],suit:13,summar:0,sure:4,surfac:[1,4,8,13],suspect:4,svg:5,system:[4,11,12],systemic_perturb:4,t:[1,4,7,10,11,13],tabl:12,take:[4,10,11],taken:[0,4,12],task:1,technic:11,techniqu:[7,11],tell:[0,8,10,11],tend:12,term:[1,4,11],termin:4,test:[0,1,4,5,6],test_case_1:13,test_case_2:13,test_case_3:13,test_data:[10,13],test_distortion_less_than_1_mm_within_r_100_mm:4,test_distortion_less_than_1mm_within_r_100:4,test_distortion_less_than_2_mm_within_r_150_mm:4,test_distortion_less_than_2mm_within_r_100:4,test_distortion_less_than_2mm_within_r_200:4,test_distortion_less_than_4_mm_within_r_200_mm:4,test_distortion_less_than_4mm_within_r_150:4,tests_to_run:[4,13],tg284_test:4,than:[0,1,4,12],thei:[4,10,11],them:[1,2,4,10,11],theme:3,therefor:[1,8,11,12,13],thi:[0,1,2,4,7,8,9,10,11,13],thick:12,thing:[1,2,4,6,8],think:[0,1,12],third:1,three:[0,4,8,11,12],threshold:4,through:[4,10],thrown:4,tidier:11,tight_layout:1,time:[3,4],titl:4,to_csv:[1,7,8,11],todo:[1,8],toler:4,too:[4,12],took:[3,11],tool:12,total:8,transform:[1,4],transvers:0,trick:2,trimdataby_r_out:4,trivial:1,truth:[4,8,11,13],turn:[0,11],tutori:[9,11,13],two:[0,1,2,3,4,8,11,12,13],type:[4,10,13],typora:3,unabl:[1,11],under:[4,10],underli:10,understood:13,uniform:4,unperturb:1,unperturbed_report:1,untest:4,unzip:[7,8,10,11],up:[0,1,4,8,11,12],updat:[4,5,8,10,11],upper:4,us:[0,1,2,4,7,8,9,10,11,12,13],user:[4,7,8,10,11,13],ut:[0,1,4],util:[1,9,10],utilis:4,valu:[1,4,11],variabl:10,variant:8,variat:1,variou:11,verbos:[0,1,4,8,10,11],veri:[0,4,10,12,13],version:4,versu:6,via:[4,10],visibl:11,visualis:4,vmax:4,vmin:4,volum:[0,1,4,7,8,11,12,13],voxel:[1,4,10,12],vv:5,wa:[0,4,10,11,13],wai:[4,6,10,11,13],want:[0,4,8,10,11],warn:[4,8,10,11,13],warp:0,warpsearchdata:[1,4],wasn:13,water:[4,6],we:[0,1,2,4,6,7,8,11,12,13],welcom:4,well:[1,4,6,9,10],were:4,what:[1,6,10,12,13],whatev:[2,4],when:[0,4,6,12],whenev:11,where:[4,8,10,11],whether:4,which:[0,1,2,4,7,10,11,12,13],why:[0,11],wilei:4,willing:8,window:4,within:[1,4,8,10],wm9vndv47u941ju:[7,8,10,11],won:1,work:[0,4,9,10,13],workflow:[1,2,9,10],world:4,would:[0,1,8,10,11,12],wrapper:[4,8],write:[0,4],write_html_report:[4,13],written:[4,8],x:[0,1,4,8,10,11,12],x_b0:4,x_gnl:4,x_gt:4,y:[0,1,4,8,10,11,12],y_b0:4,y_gnl:4,y_gt:4,yet:[4,10],you:[2,4,7,8,9,10,11,12,13],your:[4,6,11,13],your_home_directori:13,yourself:11,z:[0,1,4,8,10,11,12],z_b0:4,z_gnl:4,z_gt:4,z_max:4,z_min:4,zy:8},titles:["Harmonic sanity checks","Harmonics perturbation analysis","Preamble","&lt;no title&gt;","Code Documentation","Developer notes","Worked Examples","Field Calculation","Spherical Harmonics","MRI Distortion QA","Marker Extraction","Marker Matching","Marker position stability","Reporting"],titleterms:{"1":[0,13],"2":[0,13],"case":13,"do":[10,11],If:12,Is:10,The:8,ad:13,an:11,analysi:1,aren:12,averag:12,b0:0,basic:10,calcul:7,calculate_harmon:4,check:[0,5],code:[1,4,8],comparison:1,conclus:0,content:[6,9],coverag:5,creat:11,custom:13,data:[11,13],develop:5,dicom:11,differ:[0,1],direct:0,directli:13,distort:[1,9],docstr:5,document:4,easi:8,effect:0,error:12,exampl:[6,8,10,11],explain:8,extract:10,fail:11,fat:10,field:7,fieldanalysi:4,fieldcalcul:4,find:10,forward:0,from:11,give:0,gradient:11,guarante:10,gx:0,gy:0,gz:0,handl:10,harmon:[0,1,8,12,13],how:11,i:11,imag:0,incorpor:11,indic:9,intepret:11,interrog:5,intro:0,issu:10,json:11,major:10,marker:[10,11,12],markeranalysi:4,markervolum:11,match:11,mri:9,next:[8,11],note:5,order:0,our:12,output:8,pass:13,perturb:1,plane:0,posit:12,preambl:2,predict:1,process:11,python:0,qa:9,reconstruct:13,remain:8,report:[4,13],residu:1,revers:[0,11],same:0,saniti:0,script:0,set:1,shift:10,should:[0,11],shouldn:0,simpl:[8,11],spheric:8,stabil:12,step:[8,11],t:[0,12],tabl:9,test:13,thi:12,thing:11,util:4,versu:11,wai:8,water:10,we:10,well:11,what:11,when:11,why:12,work:[6,11],wors:12,your:10}})