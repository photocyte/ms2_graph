#! /usr/bin/python

##A tool that takes CSV files in the MZMine 2.29 format, and converts it to a graphML with the MS2 similarity in said graph

##TODO, right now the script uses similarity in the RT and m/z to match the MS2 similarity "edge", between the right nodes 
##If a single file/peaklist was searched against itself, then we can simply lookup the feature ID, and save a lot of trouble and potential error?
##If multiple peaklists were searched (e.g. A against B), then we can't rely on the feature ID for lookup, as there may be collisions, 
##and situations like where feature X has an id to match against, but feature id(from A) is present, when we intended feature id(from B), 
##
##Ideally, we could figure out from the CSV which case happened (search against own peaklist, versus seach against 2 peaklists), and have a switch to use the two modes
##But, I believe there isn't any information in the CSV (currently) which could tell you if the same versus different peaklists were searched.
##And currently the similarity search module doesn't export the ID of the matched scan.

import argparse
import time
import math
import csv
print "Basic imports done."

import networkx
print "NetworkX import done."

parser = argparse.ArgumentParser(description='Convert MZMine2 MS2 similarity CSV to graphML. Use the MZMine2 Identification->MS2 similarity search module to to annotate similar ions, then export all the information to CSV with Export/Import->Export CSV, with *all* the exportable options picked, including "Export all IDs for peaks".')

parser.add_argument('-f',metavar="CSVFILE", nargs='+',required=True,help="MZMine2 produced CSV files")
#parser.add_argument('--mz_diff_limit',metavar="0.00015",help="MZMine2 produced CSV files")
#parser.add_argument('--rt_diff_limit',metavar="0.00015",help="MZMine2 produced CSV files")

args = parser.parse_args()

dataframes = dict()
graphs = dict()

class Feature:
	def __init__(self,id,mz,rt,presplit_edges):
		self.id = id
		self.mz = float(mz)
		self.rt = float(rt)
		self.edges = presplit_edges.split(';')
		##convert edges to a list of dictionaries
		for i in range(0,len(self.edges)):
			items = self.edges[i].split(" ")
			
			##Figure out how to parse the edges...
			if len(items) <= 5 or "MS2similarity" not in self.edges[i]:
				print "WARNING! not enough items in identity:",
				print items,len(items)
				self.edges[i] = None
				continue

			elif len(items) > 7:
				print "WARNING! too many items in identity:",
				print items,len(items)
				self.edges[i] = None
				continue

			elif len(items) == 6 and "MS2similarity" in self.edges[i]:
				##print "Have expected number of items."
	                        edge_mz = float(items[1][4:])
        	                s = items[2].find("RT:") + len("RT:")
                	        edge_rt = float(items[2][s:])
                       		edge_score = float(items[3][7:])
                        	edge_num_ions_matched = int(items[4][15:])

                        	s = items[5].find("chedIons:") + len("chedIons:")
                        	edge_matched_ions = items[5][s:].split("_")
                        	self.edges[i] = {'mz':edge_mz,'rt':edge_rt,'score':edge_score,'num_ions_matched':edge_num_ions_matched,'matched_ions':edge_matched_ions}
			elif len(items) == 5 and "MS2" in self.edges[i]:
				##Possibly working with an older version of MZmine that doesn't include the matched ions.
				self.edges[i] = None
				continue

				#edge_mz = float(items[2][4:])
				#s = items[3].find("RT:") + len("RT:")
				#edge_rt = float(items[3][s:])
				#edge_score = float(items[4][7:])
				#edge_num_ions_matched = int(items[5][15:])
				#self.edges[i] = {'mz':edge_mz,'rt':edge_rt,'score':edge_score,'num_ions_matched':edge_num_ions_matched}
	def get_label(self):
		return "ID:"+str(self.id)+" mz:"+str(round(self.mz,4))+" rt:"+str(round(self.rt,1))
		
for file in args.f:
	print "Opening file",file,"..."
	
	handle = open(file,"rU")	 
	csv_rows = csv.reader(handle,delimiter=',')

	graphs[file] = networkx.Graph()
	feature_ids = set()
	feature_ids_with_collisions = set()

	i=0
	for row in csv_rows:
		if i == 0:
			##Header row
			print "CSV header is as follows:"
			print row
			assert 'row ID' in row
			assert 'row m/z' in row
			assert 'row retention time' in row
			assert 'row identity (all IDs)' in row
			row_id_index = row.index('row ID')
			row_mz_index = row.index('row m/z')
			row_rt_index = row.index('row retention time')
			row_ident_index = row.index('row identity (all IDs)')
		elif i > 0:
			id = row[row_id_index]
        		mz = row[row_mz_index]
        		rt = row[row_rt_index]
			print id,mz,rt
        		edges = row[row_ident_index]
			if id in feature_ids:
				feature_ids_with_collisions.add(id)
			else:
				feature_ids.add(id)
        		feature = Feature(id,mz,rt,edges)
        		if feature.edges != None:
				graphs[file].add_node(feature,name=feature.get_label(),mz=feature.mz,rt=feature.rt)	
		i+=1

	#print graphs[file].nodes()
	##Nodes are now setup.
	print "Finished parsing nodes for file",file
	print len(graphs[file]),"nodes parsed."
	print "There were",len(feature_ids_with_collisions),"feature id collisions."
	print "Converting CSV 'edges' to actual NetworkX graph edges"
	print "This could take awhile..."

	##Setup_edges
	node_number = len(graphs[file])
	i_count = 0
	j_count = 0
	last_thing_printed = 0
	for node_a in graphs[file]:
		i_count+=1
		for edge in node_a.edges:
			##Iterate over all the nodes and figure out which edge is equal to which node.
			##Add an edge to node_b if true
			##Ideally, should be just looking up the ID, and adding an edge.
			for node_b in graphs[file]:
				if edge == None:
					break
				mz_diff = math.fabs(node_b.mz - edge['mz'])
				##mz_diff_limit = 3e-6*node_b.mz ##3 PPM
				mz_diff_limit = 0.00015

				rt_diff = math.fabs(node_b.rt - edge['rt'])
				rt_diff_limit = 0.015 ##0.001 mins.
				if mz_diff < mz_diff_limit and rt_diff < rt_diff_limit:
					#print "YES","mz_a:",edge['mz'],"mz_b",node_b.mz,"rt_a",edge['rt'],"rt_b",node_b.rt			
					edge_name = sorted([node_a.get_label(),node_b.get_label()])
					edge_name = edge_name[0]+" and "+edge_name[1]
					if 'matched_ions' in edge.keys():
						graphs[file].add_edge(node_a,node_b,name=edge_name,edge_mz=edge['mz'],edge_rt=edge['rt'],score=edge['score'],num_ions_matched=edge['num_ions_matched'],matched_ions=" ".join(edge['matched_ions']))
					else:
						graphs[file].add_edge(node_a,node_b,name=edge_name,edge_mz=edge['mz'],edge_rt=edge['rt'],score=edge['score'],num_ions_matched=edge['num_ions_matched'])	
					j_count+=1
		percentage_finished = int((float(i_count)/float(node_number))*100)
		if percentage_finished % 10 == 0 and percentage_finished != last_thing_printed:
			print str(percentage_finished)+"% complete"
			last_thing_printed = percentage_finished
					
	print j_count,"CSV edges matched to their putative nodes"
	graph = graphs[file]
	#networkx.write_gml(graph,file+'.gml')
	filename = file+'.graphml'
	print "Writing results to",filename
	networkx.write_graphml(graph,filename)
	print 'done with',file
