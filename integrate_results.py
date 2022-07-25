import os, ast
import pandas as pd
import cairosvg
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def clear_summary_plots(celltype_repres_dir):
	todel = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)==False]
	for t in todel:
		os.system("rm {}".format(celltype_repres_dir+'/'+t))

def ct_progress_plot(celltype_repres_dir):
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	for l in ls0:
		
		agrfile = celltype_repres_dir+"/"+l+"/post_clustering/agr_Progress.txt"
		lines = open(agrfile, 'r').readlines()
		xlabel = lines[0].replace("\n","")
		x = ast.literal_eval(lines[2])
		heights = ast.literal_eval(lines[3])
		plt.plot(x, heights, label="Overall Agreement", color="yellowgreen", linewidth=3, alpha=0.9)
		plt.fill_between(x, heights, color="yellowgreen", alpha=0.4)
		plt.xlabel(xlabel)
		plt.ylabel("Overall Agreement")
		plt.tight_layout()
		plt.savefig(celltype_repres_dir+"/{}_agr_progress.png".format(l), format="png", dpi=500)
		plt.savefig(celltype_repres_dir+"/{}_agr_progress.svg".format(l), format="svg")
		plt.clf()
		

		oefile = celltype_repres_dir+"/"+l+"/post_clustering/oe_agr_Progress.txt"
		lines = open(oefile, 'r').readlines()
		xlabel = lines[0].replace("\n","")
		x = ast.literal_eval(lines[2])
		heights = ast.literal_eval(lines[3])
		plt.plot(x, heights, label="log(O/E) Agreement", color="palevioletred", linewidth=3, alpha=0.9)
		plt.fill_between(x, heights, color="palevioletred", alpha=0.4)
		plt.xlabel(xlabel)
		plt.ylabel("log(O/E) Agreement")
		plt.tight_layout()
		plt.savefig(celltype_repres_dir+"/{}_OEagr_progress.png".format(l), format="png", dpi=500)
		plt.savefig(celltype_repres_dir+"/{}_OEagr_progress.svg".format(l), format="svg")
		plt.clf()


		ckfile = celltype_repres_dir+"/"+l+"/post_clustering/ck_Progress.txt"
		lines = open(ckfile, 'r').readlines()
		xlabel = lines[0].replace("\n","")
		x = ast.literal_eval(lines[2])
		heights = ast.literal_eval(lines[3])
		plt.plot(x, heights, label="Cohen's Kappa", color="mediumpurple", linewidth=3, alpha=0.9)
		plt.fill_between(x, heights, color="mediumpurple", alpha=0.4)
		plt.xlabel(xlabel)
		plt.ylabel("Cohen's Kappa")
		plt.tight_layout()
		plt.savefig(celltype_repres_dir+"/{}_ck_progress.png".format(l), format="png", dpi=500)
		plt.savefig(celltype_repres_dir+"/{}_ck_progress.svg".format(l), format="svg")
		plt.clf()


	
		

def ct_agr(celltype_repres_dir, ck=True):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		if ck:
			svgfile = celltype_repres_dir+"/"+l+"/16_labels/agr/cohenskappa.svg"
		else:
			svgfile = celltype_repres_dir+"/"+l+"/16_labels/agr/agreement.svg"

		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	if ck:
		# new_im.save(celltype_repres_dir+"/integ_ck_{}.eps".format("_".join(name)))
		new_im.save(celltype_repres_dir+"/integ_ck_{}.png".format("_".join(name)))

	else:
		# new_im.save(celltype_repres_dir+"/integ_agr_{}.eps".format("_".join(name)))
		new_im.save(celltype_repres_dir+"/integ_agr_{}.png".format("_".join(name)))


def ct_cc(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/cc/cc_subplot.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_cc_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_cc_{}.eps".format("_".join(name)))

def ct_clb(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/clb_1/clb_subplot.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_clb_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_clb_{}.eps".format("_".join(name)))

def ct_reprtss(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/tss_rep1/rep_vs_TSS_enrichment_bars.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=600)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_reprtss_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_reprtss_{}.eps".format("_".join(name)))

def ct_enrtss(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/16_labels/tss_rep1/general_TSS_enrichment.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/integ_enrtss_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_enrtss_{}.eps".format("_".join(name)))

def ct_clstrmap(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	agrfiles = {}
	for l in ls0:
		svgfile = celltype_repres_dir+"/"+l+"/post_clustering/clustermap.svg"
		pngfile = svgfile.replace(".svg", ".png")

		if os.path.exists(pngfile)==False:
			cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

		agrfiles[var_setting_dict[l]] = Image.open(pngfile)

	agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
	widths, heights = zip(*(i.size for i in agrfiles.values()))
	total_width = sum(widths)	
	max_height = max(heights)	
	new_im = Image.new('RGB', (total_width, max_height))

	x_offset = 0
	name = []
	for k, im in agrfiles.items():
		name.append(k.split(" ")[-1]) 
		new_im.paste(im, (x_offset,0))
		x_offset += im.size[0]

	new_im.save(celltype_repres_dir+"/clustermap_{}.png".format("_".join(name)))
	# new_im.save(celltype_repres_dir+"/integ_enrtss_{}.eps".format("_".join(name)))

def ct_short_report(celltype_repres_dir):
	var_setting_dict = {
		"concatenated":"Diff. data, Shared training", "rep1_vs_rep2":"Diff. data, Separate training", 
		"rep1_pseudoreps":"Similar data, Shared Training", "rep1_paraminit":"Same data, Separate training"}
	ls0 = [l for l in os.listdir(celltype_repres_dir) if os.path.isdir(celltype_repres_dir+"/"+l)]

	if "chmm" in celltype_repres_dir:
		ls0.remove("rep1_paraminit")
	elif "segway" in celltype_repres_dir:
		ls0.remove("rep1_pseudoreps")

	file_tomerge = [
		"bidir_Cohens_Kappa.svg", "bidir_logOE_Agreement.svg", "bidir_Raw_Agreement.svg"]

	for f in file_tomerge:
		try:
			agrfiles = {}
			for l in ls0:
				svgfile = celltype_repres_dir+"/"+l+"/"+f
				pngfile = svgfile.replace(".svg", ".png")

				if os.path.exists(pngfile)==False:
					cairosvg.svg2png(url=svgfile, write_to=pngfile, dpi=500)

				agrfiles[var_setting_dict[l]] = Image.open(pngfile)

			agrfiles = dict(sorted(agrfiles.items(), key=lambda item: item[0]))
			widths, heights = zip(*(i.size for i in agrfiles.values()))
			total_width = sum(widths)	
			max_height = max(heights)	
			new_im = Image.new('RGB', (total_width, max_height))

			x_offset = 0
			name = []
			for k, im in agrfiles.items():
				name.append(k.split(" ")[-1]) 
				new_im.paste(im, (x_offset,0))
				x_offset += im.size[0]

			new_im.save(celltype_repres_dir+"/integrated_{}_{}.png".format(f.replace(".svg", ""), "_".join(name)))
		except:
			pass
def compare_overalls(res_dir, target_metric="ck"):
	
	var_setting_dict = {
		"concatenated":"S2: Diff. data, Shared train", "rep1_vs_rep2":"S1: Diff. data, Separate train", 
		"rep1_pseudoreps":"S4: Similar data, Shared Train", "rep1_paraminit":"S3s: Same data, Separate train"}

	navig = []
	ct_list = [ct for ct in os.listdir(res_dir) if os.path.isdir(res_dir+"/"+ct)]
	for ct in ct_list:
		settings = [s for s in os.listdir(res_dir+"/"+ct) if os.path.isdir(res_dir+"/"+ct+"/"+s)]

		# if "chmm" in res_dir:
		# 	settings.remove("rep1_paraminit")
		# elif "segway" in res_dir:
		# 	settings.remove("rep1_pseudoreps")

		for s in settings:
			file = res_dir+"/"+ct+"/"+s+"/16_labels/agr/"
			if target_metric=="ck":
				file = file + "cohenskappa.txt"

			elif target_metric=="agr":
				file = file + "agreement.txt"
			
			elif target_metric=="oe":
				file = file + "oe_agreement.txt"

			lines = open(file, 'r').readlines()
			title = lines[0].replace("\n","")
			xlabel = lines[1].replace("\n","")
			ylabel = lines[2].replace("\n","")
			x = ast.literal_eval(lines[3])
			heights = ast.literal_eval(lines[4].replace("-inf", "0"))
			navig.append([ct, var_setting_dict[s], heights[-1]])
	
	navig = pd.DataFrame(navig, columns=["Celltype", "Setting", target_metric]).sort_values(by="Setting")
	# navig.sort_values(by="Setting")

	sns.set_theme(style="whitegrid")
	sns.reset_orig
	plt.style.use('default')

	sns.set(rc={'figure.figsize':(10,8)})
	sns.set_palette(sns.color_palette("deep"))
	sns.barplot(x="Celltype", y=target_metric, hue="Setting", data=navig)

	plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)

	if target_metric=="agr" or target_metric=="ck":
		plt.yticks(np.arange(0, 1.1, step=0.1))
	elif target_metric == "oe":
		plt.yticks(np.arange(0, 4.5, step=0.5))

	plt.tight_layout()
	plt.savefig(res_dir+"/overall_comparison_{}.png".format(target_metric), format="png", dpi=500)
	plt.savefig(res_dir+"/overall_comparison_{}.svg".format(target_metric), format="svg")
	plt.clf()
	

def INTEGRATE_ALL(ct_dir):
	# clear_summary_plots(ct_dir)
	# print("summarizing", ct_dir)
	ct_progress_plot(ct_dir)
	# ct_agr(ct_dir, ck=True)
	# ct_cc(ct_dir)
	# ct_clb(ct_dir)
	# ct_enrtss(ct_dir)
	# ct_reprtss(ct_dir)
	# ct_clstrmap(ct_dir)


if __name__=="__main__":
	ct_list = ['CD14-positive_monocyte', "GM12878", "K562", "HeLa-S3", "MCF-7"]
	segres_dir = "tests/repres_subset/segway/"
	chmmres_dir = "tests/reprod_results/chmm/"
	seg_short = "tests/short_reports/segway/"
	chmm_short = "tests/short_reports/chmm/"

	for ct in ct_list:
		INTEGRATE_ALL("{}/{}".format(segres_dir, ct))
		INTEGRATE_ALL("{}/{}".format(chmmres_dir, ct))

		# ct_short_report(seg_short+ct)
		# ct_short_report(chmm_short+ct)

		
	
	exit()
	compare_overalls(chmmres_dir)
	compare_overalls(segres_dir)
	compare_overalls(chmmres_dir, target_metric="oe")
	compare_overalls(segres_dir, target_metric="oe")
	compare_overalls(chmmres_dir, target_metric="agr")
	compare_overalls(segres_dir, target_metric="agr")
	

