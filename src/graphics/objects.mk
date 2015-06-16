general_OBJECTS = \
	general/alloc_fac.lo \
	general/atom_handling.lo \
	general/atom_param.lo \
	general/avstruc.lo \
	general/calc_quatfit.lo \
	general/callback.lo \
	general/cluster_manager.lo \
	general/colouring.lo \
	general/coord_man.lo \
	general/correl.lo \
	general/data_structure.lo \
	general/defcell.lo \
	general/diagonalize.lo \
	general/edit_bond.lo \
	general/findhbonds.lo \
	general/findssbonds.lo \
	general/format.lo \
	general/input_gen.lo \
	general/iso_file_io.lo \
	general/ldp_control.lo \
	general/listener.lo \
	general/math_oper.lo \
	general/measure.lo \
	general/monitor.lo \
	general/msd.lo \
	general/nreg.lo \
	general/own_strtok.lo \
	general/plumber.lo \
	general/quatfit.lo \
	general/rdf.lo \
	general/sec_struct.lo \
	general/selections.lo \
	general/signal_handling.lo \
	general/spline.lo \
	general/stdafx.lo \
	general/trace_manager.lo \
	general/trajectory.lo \
	general/utility1.lo \
	general/utility2.lo \
	general/utility3.lo \
	general/utility4.lo \
	general/xasprintf.lo \
	general/xmath.lo \
	general/xprintf.lo \
	general/xscanf.lo \
	general/xtmpfile.lo \
	$(end_of_list)
graphics_OBJECTS = \
	graphics/g_Mmain.lo \
	graphics/g_Mmake.lo \
	graphics/g_clusterdriver.lo \
	graphics/g_coloring.lo \
	graphics/g_contdriver.lo \
	graphics/g_contour.lo \
	graphics/g_contplane.lo \
	graphics/g_font_manager.lo \
	graphics/g_hardcopy.lo \
	graphics/g_isocurve.lo \
	graphics/g_ldpdriver.lo \
	graphics/g_light_models.lo \
	graphics/g_main.lo \
	graphics/g_manipulate.lo \
	graphics/g_opengl_util.lo \
	graphics/g_output.lo \
	graphics/g_plot.lo \
	graphics/g_plot_cell.lo \
	graphics/g_plot_cluster.lo \
	graphics/g_plot_cpk.lo \
	graphics/g_plot_label.lo \
	graphics/g_plot_ldp.lo \
	graphics/g_plot_licorice.lo \
	graphics/g_plot_prop.lo \
	graphics/g_plot_stick.lo \
	graphics/g_plot_trace.lo \
	graphics/g_plotobj.lo \
	graphics/g_plumber.lo \
	graphics/g_quadstereo.lo \
	graphics/g_setup.lo \
	graphics/g_stereo.lo \
	graphics/g_surface_jpr.lo \
	graphics/g_tracedriver.lo \
	graphics/g_trajectory.lo \
	graphics/g_util.lo \
	$(end_of_list)
parser_OBJECTS = \
	parser/p_atom.lo \
	parser/p_calculatec.lo \
	parser/p_center.lo \
	parser/p_contourc.lo \
	parser/p_copy.lo \
	parser/p_definec.lo \
	parser/p_diagonalize.lo \
	parser/p_display.lo \
	parser/p_editc.lo \
	parser/p_exportc.lo \
	parser/p_fillc.lo \
	parser/p_findc.lo \
	parser/p_hardcopyc.lo \
	parser/p_helpc.lo \
	parser/p_manipulatec.lo \
	parser/p_monitorc.lo \
	parser/p_openc.lo \
	parser/p_parser.lo \
	parser/p_plotc.lo \
	parser/p_plumberc.lo \
	parser/p_pythonc.lo \
	parser/p_read_fac.lo \
	parser/p_resetc.lo \
	parser/p_rotate.lo \
	parser/p_s_atomc.lo \
	parser/p_s_cellc.lo \
	parser/p_s_contourc.lo \
	parser/p_s_cutplanec.lo \
	parser/p_s_displaylistsc.lo \
	parser/p_s_gbasisc.lo \
	parser/p_s_lightc.lo \
	parser/p_s_materialc.lo \
	parser/p_s_monitorc.lo \
	parser/p_s_plumberc.lo \
	parser/p_s_windowc.lo \
	parser/p_savec.lo \
	parser/p_scalec.lo \
	parser/p_selectc.lo \
	parser/p_showc.lo \
	parser/p_tcl.lo \
	parser/p_tracec.lo \
	parser/p_trajparse.lo \
	parser/p_translatec.lo \
	parser/p_window.lo \
	$(end_of_list)
readwrite_OBJECTS = \
	readwrite/rw_get_frame_amber.lo \
	readwrite/rw_get_frame_ambera.lo \
	readwrite/rw_get_frame_cerius2.lo \
	readwrite/rw_get_frame_charmm.lo \
	readwrite/rw_get_frame_discover.lo \
	readwrite/rw_get_frame_dl_poly.lo \
	readwrite/rw_get_frame_gromacs.lo \
	readwrite/rw_get_frame_gromos.lo \
	readwrite/rw_get_frame_gromos96a.lo \
	readwrite/rw_get_frame_hyper.lo \
	readwrite/rw_get_frame_mumod.lo \
	readwrite/rw_get_frame_tinker.lo \
	readwrite/rw_get_frame_xmol.lo \
	readwrite/rw_get_frame_xplor.lo \
	readwrite/rw_get_frame_yasp.lo \
	readwrite/rw_gotobmp.lo \
	readwrite/rw_gotosgi.lo \
	readwrite/rw_gototga.lo \
	readwrite/rw_gotoxwd.lo \
	readwrite/rw_ramber.lo \
	readwrite/rw_read_gbl_input.lo \
	readwrite/rw_read_model.lo \
	readwrite/rw_readdict.lo \
	readwrite/rw_rforce.lo \
	readwrite/rw_rgaussian.lo \
	readwrite/rw_rhyper.lo \
	readwrite/rw_rinsight.lo \
	readwrite/rw_rkarp.lo \
	readwrite/rw_rmol2.lo \
	readwrite/rw_rmopac.lo \
	readwrite/rw_rmumod.lo \
	readwrite/rw_ropenmol.lo \
	readwrite/rw_rpdb.lo \
	readwrite/rw_rxmol.lo \
	readwrite/rw_rxyz.lo \
	readwrite/rw_ryasp.lo \
	readwrite/rw_tocolps.lo \
	readwrite/rw_write_coord.lo \
	readwrite/rw_write_model.lo \
	$(end_of_list)
OBJECTS = \
$(general_OBJECTS) $(graphics_OBJECTS) $(parser_OBJECTS) $(readwrite_OBJECTS)
OBJLISTS = \
general/all_lo graphics/all_lo parser/all_lo readwrite/all_lo
