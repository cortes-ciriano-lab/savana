"""
Class definitions for SAVANA: ConsensusBreakpoint, PotentialBreakpoint, and Cluster
Created: 13/04/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import json
import uuid

from copy import copy
from statistics import mean, median, pstdev

def generate_uuid():
	""" hex representation of a multiprocessing-safe unique id """
	return uuid.uuid4().hex

class ConsensusBreakpoint():
	""" class for a second-round called breakpoint (stores originating cluster information """
	def __init__(self, locations, source, originating_cluster, end_cluster, labels, breakpoint_notation=None, insert=None):
		self.uid = generate_uuid()
		self.start_chr = locations[0]['chr']
		self.start_loc = int(locations[0]['loc'])
		self.end_chr = locations[1]['chr']
		self.end_loc = int(locations[1]['loc'])
		self.source = source
		if source == 'INS' and not insert:
			raise AttributeError("Must provide an insert for breakpoint type of 'INS")
		self.inserted_sequence = insert
		self.originating_cluster = originating_cluster
		self.end_cluster = end_cluster if end_cluster else originating_cluster
		self.labels = labels
		self.count = None # used later for standardising across output files
		self.local_depths = {} # add later
		self.breakpoint_notation = breakpoint_notation
		# use the labels to calculate support by counting reads
		self.support = {'normal': 0, 'tumour': 0}
		for label, reads in self.labels.items():
			self.support[label]+=len(reads)
		# calculate the length
		self.sv_length = None
		if self.breakpoint_notation == "<INS>":
			self.sv_length = str(len(self.inserted_sequence))
		elif self.start_chr != self.end_chr:
			self.sv_length = 0
		else:
			self.sv_length = abs(int(self.end_loc)-int(self.start_loc))

	def as_dict(self):
		""" return dict representation of breakpoint"""
		self_dict = {
			"start_chr": self.start_chr,
			"start_loc": self.start_loc,
			"end_chr": self.end_chr,
			"end_loc": self.end_loc,
			"inserted_sequence": self.inserted_sequence,
			"source": self.source,
			"originating_cluster": self.originating_cluster.uid,
			"end_cluster": self.end_cluster.uid,
			"labels": "/".join([f'{label}_{str(len(reads))}' for label, reads in self.labels.items()]),
			#"local_depths": self.local_depths, TODO: put this into a string format
			"breakpoint_notation": self.breakpoint_notation
		}
		return self_dict

	def as_bed(self, contig_lengths):
		""" return bed line(s) plus buffer of breakpoint """
		# min is 0 max is length of contig
		bed_lines = [[
			self.start_chr,
			str(min(max(self.start_loc, 0),contig_lengths[self.start_chr])),
			str(min(max(self.start_loc+1, 0),contig_lengths[self.start_chr])),
			str(self.uid),
			'0' # 0th edge
		]]
		if self.breakpoint_notation != "<INS>":
			bed_lines.append([
				self.end_chr,
				str(min(max(self.end_loc, 0), contig_lengths[self.end_chr])),
				str(min(max(self.end_loc+1, 0), contig_lengths[self.end_chr])),
				str(self.uid),
				'1' # 1st edge
			])

		return "\n".join("\t".join(l) for l in bed_lines)+"\n"

	def as_bedpe(self, count):
		""" return bedpe line(s) representation of breakpoint """
		if not self.count:
			self.count = count
		bedpe_line = [
			self.start_chr,
			str(self.start_loc),
			str(self.start_loc),
			self.end_chr,
			str(self.end_loc),
			str(self.end_loc)
		]
		if self.source == "INS":
			# add 1bp to the bedpe location for igv rendering
			bedpe_line[4] = str(int(bedpe_line[4])+1)
			bedpe_line[5] = str(int(bedpe_line[5])+1)
		label_string = "/".join([f'{label.upper()}_{str(len(reads))}' for label, reads in self.labels.items()])
		bp_length = str(abs(int(self.end_loc)-int(self.start_loc)))
		start_cluster = str(self.originating_cluster.uid)
		end_cluster = str(self.end_cluster.uid)
		bedpe_line.append(f'ID_{count}|{bp_length}bp|{start_cluster}/{end_cluster}|{label_string}|{self.breakpoint_notation}')

		return "\t".join(bedpe_line)+"\n"

	def as_read_support(self, count):
		""" return read support line representation of a breakpoint """
		if ',' in ("".join(self.labels.get('tumour',[])+self.labels.get('normal',[]))):
			# quote everything since at least one read name contains a comma
			read_support_line = [
				f'ID_{count}',
				'"'+'","'.join(self.labels.get('tumour',[]))+'"',
				'"'+'","'.join(self.labels.get('normal',[]))+'"'
			]
		else:
			read_support_line = [
				f'ID_{count}',
				'"'+",".join(self.labels.get('tumour',[]))+'"',
				'"'+",".join(self.labels.get('normal',[]))+'"'
			]
		return "\t".join(read_support_line)+"\n"

	def as_variant_stats(self, count, stats_column_order):
		""" return variant line with its stats """
		variant_stats_lines = []
		start_cluster_stats, end_cluster_stats = [],[]
		for key in stats_column_order:
			start_cluster_stats.append(self.originating_cluster.get_stats()[key])
			end_cluster_stats.append(self.end_cluster.get_stats()[key])
		variant_stats_lines.append([
			f'{self.start_chr}:{self.start_loc}',
			f'ID_{count}_1',
			f'{self.breakpoint_notation}',
			f'{self.sv_length}',
			f'{self.originating_cluster.uid}',
			f'{self.end_cluster.uid}',
			f'{self.support["tumour"]}',
			f'{self.support["normal"]}'
		]+[str(s) for s in start_cluster_stats]+[str(s) for s in end_cluster_stats])

		if self.breakpoint_notation != "<INS>":
			variant_stats_lines.append([
				f'{self.end_chr}:{self.end_loc}',
				f'ID_{count}_2',
				f'{self.breakpoint_notation}',
				f'{self.sv_length}',
				f'{self.originating_cluster.uid}',
				f'{self.end_cluster.uid}',
				f'{self.support["tumour"]}',
				f'{self.support["normal"]}'
			]+[str(s) for s in start_cluster_stats]+[str(s) for s in end_cluster_stats])

		variant_stats_str = ''
		for line in variant_stats_lines:
			variant_stats_str+="\t".join(line)+"\n"
		return variant_stats_str

	def get_alts(self, start_base, end_base):
		""" return the vcf ALT notation for a breakpoint (both lines) """
		if self.breakpoint_notation == "<INS>":
			# only one line for insertions
			return ["<INS>"]
		alts = ['', '']
		if self.breakpoint_notation.startswith("+"):
			alts[0]+=f'{start_base}' # start: +
			if self.breakpoint_notation.endswith("+"):
				alts[0]+=f']{self.end_chr}:{self.end_loc}]' # first: ++
				alts[1]+=f'{end_base}]{self.start_chr}:{self.start_loc}]' # second: ++
			else:
				alts[0]+=f'[{self.end_chr}:{self.end_loc}[' # first: +-
				alts[1]+=f']{self.start_chr}:{self.start_loc}]{end_base}' # second: -+
		else:
			if self.breakpoint_notation.endswith("+"):
				alts[0]+=f']{self.end_chr}:{self.end_loc}]' # first: -+
				alts[1]+=f'{end_base}[{self.start_chr}:{self.start_loc}[' # second: +-
			else:
				alts[0]+=f'[{self.end_chr}:{self.end_loc}[' # --
				alts[1]+=f'[{self.start_chr}:{self.start_loc}[{end_base}' # --
			alts[0]+=f'{start_base}' # start: -

		return alts

	def get_stats_str(self):
		""" return the stats of the originating cluster """
		stats_str = f'ORIGINATING_CLUSTER={self.originating_cluster.uid};END_CLUSTER={self.end_cluster.uid};'
		stats_originating = self.originating_cluster.get_stats()
		for key, value in stats_originating.items():
			stats_str+=f'ORIGIN_{key.upper()}={value};'
		stats_end = self.end_cluster.get_stats()
		for key, value in stats_end.items():
			stats_str+=f'END_{key.upper()}={value};'
		return stats_str

	def as_vcf(self, ref_fasta):
		""" return vcf line(s) representation of the breakpoint """
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
		try:
			start_base = ref_fasta.fetch(self.start_chr, self.start_loc - 1, self.start_loc)
		except Exception as e:
			start_base = 'N'
		try:
			end_base = ref_fasta.fetch(self.end_chr, self.end_loc - 1, self.end_loc)
		except Exception as e:
			end_base = 'N'
		alts = self.get_alts(start_base, end_base)
		# construct info column
		gt_tag = ''
		if self.support['normal'] >= 1:
			gt_tag = '0/0'
		else:
			gt_tag = '0/1'
		support_str = ';'.join([f'{label.upper()}_SUPPORT={label_count}' for label, label_count in self.support.items()])
		stats_str = self.get_stats_str()
		info = [f'{support_str};SVLEN={self.sv_length};BP_NOTATION={self.breakpoint_notation};{stats_str}']
		if self.breakpoint_notation == "<INS>":
			info[0] = 'SVTYPE=INS;' + info[0]
			for label, depths in self.local_depths.items():
				for i, bin_label  in enumerate(["BEFORE","AT","AFTER"]):
					info[0]+=f'{label.upper()}_DP_{bin_label}={str(depths[i][0])};'
		else:
			info.append(info[0]) # duplicate info
			# add edge-specific info
			info[0] = f'SVTYPE=BND;MATEID=ID_{self.count}_2;' + info[0]
			for label, depths in self.local_depths.items():
				for i, bin_label  in enumerate(["BEFORE","AT","AFTER"]):
					info[0]+=f'{label.upper()}_DP_{bin_label}={",".join([str(d) for d in depths[i]])};'
					info[1]+=f'{label.upper()}_DP_{bin_label}={",".join([str(d) for d in reversed(depths[i])])};'
			info[1] = f'SVTYPE=BND;MATEID=ID_{self.count}_1;' + info[1]
		# put together vcf line(s)
		vcf_lines = [[
			self.start_chr,
			str(self.start_loc),
			f'ID_{self.count}_1',
			start_base,
			alts[0],
			'.',
			'PASS',
			info[0],
			'GT',
			gt_tag
		]]
		if self.breakpoint_notation != "<INS>":
			vcf_lines.append([
				self.end_chr,
				str(self.end_loc),
				f'ID_{self.count}_2',
				end_base,
				alts[1],
				'.',
				'PASS',
				info[1],
				'GT',
				gt_tag
			])
		vcf_string = ''
		for line in vcf_lines:
			vcf_string+="\t".join(line)+"\n"
		return vcf_string

	# string representation
	def __str__(self):
		return json.dumps(self.as_dict())
	def __repr__(self):
		return self.__str__()

class PotentialBreakpoint():
	""" class for a potential breakpoint identified from a CIGAR string or split-read """
	def __init__(self, locations, source, read_name, read_quality, label, breakpoint_notation, insert=None):
		self.uid = generate_uuid()
		self.start_chr = locations[0]['chr']
		self.start_loc = int(locations[0]['loc'])
		self.end_chr = locations[1]['chr']
		self.end_loc = int(locations[1]['loc'])
		self.inserted_sequence = insert
		self.source = source # SUPP, INS, or DEL
		if source == 'INS' and not insert:
			raise AttributeError("Must provide an insert for breakpoint type of 'INS")
		self.read_name = read_name
		self.mapq = read_quality
		self.label = label
		self.breakpoint_notation = breakpoint_notation

	def as_dict(self):
		""" return dict representation of breakpoint"""
		self_dict = {
			"start_chr": self.start_chr,
			"start_loc": self.start_loc,
			"end_chr": self.end_chr,
			"end_loc": self.end_loc,
			"inserted_sequence": self.inserted_sequence,
			"source": self.source,
			"read_name": self.read_name,
			"mapq": self.mapq,
			"label": self.label,
			"breakpoint_notation": self.breakpoint_notation
		}
		return self_dict

	# string representation
	def __str__(self):
		return json.dumps(self.as_dict())
	def __repr__(self):
		return self.__str__()

	# implemented for interval tree sorting by start
	def __lt__(self, other):
		if self.start_chr == other.start_chr:
			if self.start_loc < other.start_loc:
				return True
		return False
	def __eq__(self, other):
		if self.start_chr == other.start_chr:
			if self.start_loc == other.start_loc:
				return True
		return False
	def __hash__(self):
		return hash(self.uid)
	def __reversed__(self):
		reversed_breakpoint = copy(self)
		setattr(reversed_breakpoint, 'start_chr', self.end_chr)
		setattr(reversed_breakpoint, 'start_loc', self.end_loc)
		setattr(reversed_breakpoint, 'end_chr', self.start_chr)
		setattr(reversed_breakpoint, 'end_loc', self.start_loc)
		return reversed_breakpoint

class Cluster():
	""" class for a cluster containing breakpoint objects within a buffer & sharing an SV type """
	def __init__(self, initial_breakpoint):
		self.uid = generate_uuid()
		self.chr = initial_breakpoint.start_chr
		self.start = initial_breakpoint.start_loc
		if initial_breakpoint.source == "SUPP":
			# only cluster supp breakpoints on start
			self.end = self.start
		else:
			self.end = initial_breakpoint.end_loc
		self.source = initial_breakpoint.source
		self.breakpoints = [initial_breakpoint]
		self.supporting_reads = {initial_breakpoint.read_name}
		self.stats = None

	def overlaps(self, other, buffer):
		""" determine if a breakpoint should be merged to a cluster """
		if self.breakpoints[0].breakpoint_notation != other.breakpoint_notation:
			print(f'Tried to add {other.breakpoint_notaion} breakpoint to {self.breakpoints[0].breakpoint_notation} cluster')
			return False
		self_start_buff = self.start - buffer
		self_end_buff = self.end + buffer
		other_start_buff = other.start_loc - buffer
		other_end_buff = (other.start_loc + buffer) if other.source == "SUPP" else (other.end_loc + buffer)
		if other_start_buff <= self_end_buff and other_start_buff >= self_start_buff:
			return True
		if other_end_buff >= self_start_buff and other_end_buff <= self_end_buff:
			return True
		return False

	def add(self, new_breakpoint):
		""" add a breakpoint to the cluster, updating relevant values """
		# consider whether the cluster is forward or reverse
		if self.start <= self.end:
			self.start = new_breakpoint.start_loc if (new_breakpoint.start_loc < self.start) else self.start
		else:
			# it's reverse
			self.start = new_breakpoint.start_loc if (new_breakpoint.start_loc > self.start) else self.start
		if new_breakpoint.source == "SUPP":
			if self.start <= self.end:
				self.end = 	new_breakpoint.start_loc if (new_breakpoint.start_loc > self.end) else self.end
			else:
				self.end = 	new_breakpoint.start_loc if (new_breakpoint.start_loc < self.end) else self.end
		else:
			if self.start <= self.end:
				self.end = new_breakpoint.end_loc if (new_breakpoint.end_loc > self.end) else self.end
			else:
				self.end = new_breakpoint.end_loc if (new_breakpoint.end_loc < self.end) else self.end
		self.breakpoints.append(new_breakpoint)
		self.supporting_reads.add(new_breakpoint.read_name)
		self.stats = None # reset stats when new breakpoint added

	def as_dict(self):
		""" return dict representation of cluster """
		self_dict = {
			"chr": self.chr,
			"start": self.start,
			"end": self.end,
			"breakpoints": []
		}
		for bp in self.breakpoints:
			self_dict['breakpoints'].append(bp.as_dict())
		return self_dict

	def get_stats(self):
		""" return dict of the stdev of start, event size, and mean size in cluster """
		if not self.stats:
			starts = []
			mapqs = []
			event_sizes = []
			for bp in self.breakpoints:
				starts.append(bp.start_loc)
				mapqs.append(bp.mapq)
				if self.source == "SUPP":
					event_sizes.append(1)
				elif self.source == "INS":
					event_sizes.append(len(bp.inserted_sequence))
				elif self.source == "DEL":
					event_sizes.append(abs(bp.start_loc - bp.end_loc))
			stat_dict = {
				'starts_std_dev': pstdev(starts),
				'mapq_mean': mean(mapqs),
				'event_size_std_dev': pstdev(event_sizes),
				'event_size_median': median(event_sizes),
				'event_size_mean': mean(event_sizes)
			}
			# round to 2 decimal places for stats dict attribute
			self.stats = {}
			for key, value in stat_dict.items():
				try:
					self_value = round(value, 2)
				except Exception as e:
					# skip rounding if error thrown
					self_value = value
				self.stats[key] = self_value

		return self.stats

	# string representation
	def __str__(self):
		return json.dumps(self.as_dict())
	def __repr__(self):
		return self.__str__()

if __name__ == "__main__":
	print("Class definitions for SAVANA")
