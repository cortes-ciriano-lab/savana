"""
Class definitions for SAVANA: ConsensusBreakpoint, PotentialBreakpoint, and Cluster
Created: 13/04/2021
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import sys
import json
import zlib

from copy import copy
from statistics import mean, median, pstdev
#TODO: test using numpy instead of statistics since it's apparently slower
#from numpy import mean, median, std

class ConsensusBreakpoint():
	""" class for a second-round called breakpoint (stores originating cluster information """
	def __init__(self, locations, source, originating_cluster, end_cluster, counts, breakpoint_notation, tumour_only=False, inserts=None):
		self.start_chr = locations[0]['chr']
		self.start_loc = int(locations[0]['loc'])
		self.end_chr = locations[1]['chr']
		self.end_loc = int(locations[1]['loc'])
		if (breakpoint_notation in ["<INS>", "-", "+"]) and not inserts:
			raise AttributeError(f'Must provide an insert for {breakpoint_notation} breakpoints')
		if (breakpoint_notation in ["<INS>", "-", "+"]) and self.start_chr != self.end_chr:
			raise AttributeError(f'{breakpoint_notation} breakpoints must have same start/end location (start_chr={self.start_chr}, end_chr={self.end_chr})')
		if (breakpoint_notation in ["<INS>", "-", "+"]) and self.start_loc != self.end_loc:
			raise AttributeError(f'{breakpoint_notation} breakpoints must have same start/end location (start_loc={self.start_loc}, end_chr={self.end_loc})')
		self.breakpoint_notation = breakpoint_notation
		self.inserted_sequences = inserts
		self.source = source
		self.originating_cluster = originating_cluster
		self.end_cluster = end_cluster if end_cluster else originating_cluster
		self.supporting_reads, self.aln_support_counts, self.phasing, self.total_read_counts = counts
		# use the label counts to calculate read support
		if tumour_only:
			self.read_support_counts = {'tumour': 0}
			for label, reads in self.supporting_reads.items():
				self.read_support_counts[label] += len(reads)
		else:
			self.read_support_counts = {'tumour': 0, 'normal': 0}
			for label, reads in self.supporting_reads.items():
				self.read_support_counts[label] += len(reads)
		# add a count of 0 for aln support if key doesn't exist
		for key in self.read_support_counts.keys():
			if key not in self.aln_support_counts:
				self.aln_support_counts[key] = 0
			if key not in self.total_read_counts:
				self.total_read_counts[key] = 0
		self.count = None # used later for standardising across output files
		# add these all later
		self.phased_local_depths = {}
		self.allele_fractions = {}
		# calculate the length
		self.sv_length = None
		if breakpoint_notation == "<INS>":
			# only use the length of CIGAR insertions
			self.sv_length = str(int(mean([bp.insert_length for bp in self.originating_cluster.breakpoints if bp.breakpoint_notation == "<INS>"])))
		elif breakpoint_notation in ["-", "+"]:
			self.sv_length = 0
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
			"source": self.source,
			"labels": "/".join([f'{label}_{str(len(reads))}' for label, reads in self.labels.items()]),
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
			'0' # 0th edge
		]]
		if (self.breakpoint_notation not in ["<INS>", "-", "+"]):
			bed_lines.append([
				self.end_chr,
				str(min(max(self.end_loc, 0), contig_lengths[self.end_chr])),
				str(min(max(self.end_loc+1, 0), contig_lengths[self.end_chr])),
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
		if self.start_chr == self.end_chr and self.start_loc == self.end_loc:
			# add 1bp to the bedpe location for igv rendering
			bedpe_line[4] = str(int(bedpe_line[4])+1)
			bedpe_line[5] = str(int(bedpe_line[5])+1)
		label_string = "/".join([f'{label.upper()}_{str(count)}' for label, count in self.read_support_counts.items()])
		bp_length = str(abs(int(self.end_loc)-int(self.start_loc)))
		bedpe_line.append(f'ID_{count}|{bp_length}bp|{label_string}|{self.breakpoint_notation}')

		return "\t".join(bedpe_line)+"\n"

	def as_read_support(self, count):
		""" return read support line representation of a breakpoint """
		if ',' in ("".join(self.supporting_reads.get('tumour',[])+self.supporting_reads.get('normal',[]))):
			# quote everything since at least one read name contains a comma
			read_support_line = [
				f'ID_{count}',
				'"'+'","'.join(self.supporting_reads.get('tumour',[]))+'"',
				'"'+'","'.join(self.supporting_reads.get('normal',[]))+'"'
			]
		else:
			read_support_line = [
				f'ID_{count}',
				'"'+",".join(self.supporting_reads.get('tumour',[]))+'"',
				'"'+",".join(self.supporting_reads.get('normal',[]))+'"'
			]
		return "\t".join(read_support_line)+"\n"

	def as_insertion_fasta(self, count):
		""" return fasta line containing inserted sequences for an insertion breakpoint """
		fasta_lines = []
		total_inserts = len(self.inserted_sequences)
		for i, insert in enumerate(self.inserted_sequences, start=1):
			decoded_insert = zlib.decompress(insert).decode()
			fasta_lines.append(f'>ID_{count}-INSSEQ_{str(i).zfill(5)} | {self.breakpoint_notation} | {i}/{total_inserts} | {len(decoded_insert)}bp')
			# split sequence to only have 80 characters per line
			fasta_lines.extend([decoded_insert[i:i+80] for i in range(0, len(decoded_insert), 80)])
		return "\n".join(fasta_lines)+"\n"

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
			f'{self.support["tumour"]}',
			f'{self.support["normal"]}'
		]+[str(s) for s in start_cluster_stats]+[str(s) for s in end_cluster_stats])

		if self.breakpoint_notation not in ["<INS>", "+", "-"]:
			variant_stats_lines.append([
				f'{self.end_chr}:{self.end_loc}',
				f'ID_{count}_2',
				f'{self.breakpoint_notation}',
				f'{self.sv_length}',
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
		elif self.breakpoint_notation in ["-","+"]:
			# only one line for single breakends
			return [f'.{end_base}']
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
		stats_originating = self.originating_cluster.get_stats()
		stats_str = ''
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
			start_base = "N" if start_base == "" else start_base
			if start_base.upper() not in ['N', 'A', 'T', 'C', 'G']:
				if start_base != "":
					print(f'unrecognised base: "{start_base}"')
				start_base = "N"
		except ValueError as _:
			start_base = 'N'
		try:
			end_base = ref_fasta.fetch(self.end_chr, self.end_loc - 1, self.end_loc)
			end_base = "N" if start_base == "" else end_base
			if end_base.upper() not in ['N', 'A', 'T', 'C', 'G']:
				if end_base != "":
					print(f'unrecognised base: "{end_base}"')
					print(f'"{end_base}"')
				end_base = "N"
		except ValueError as _:
			end_base = 'N'

		alts = self.get_alts(start_base, end_base)
		gt_tag = ''
		if 'normal' in self.read_support_counts and self.read_support_counts['normal'] >= 1:
			gt_tag = '0/0'
		else:
			gt_tag = '0/1'
		# construct info column
		support_str = ''
		for key in sorted(self.read_support_counts.keys(), reverse=True):
			support_str += f'{key.upper()}_READ_SUPPORT={self.read_support_counts[key]};'
			support_str += f'{key.upper()}_ALN_SUPPORT={self.aln_support_counts[key]};'
		stats_str = self.get_stats_str()
		info = [f'{support_str};']
		info[0] += f'SVLEN={self.sv_length};'
		info[0] += f'BP_NOTATION={self.breakpoint_notation};'
		info[0] += f'SOURCE={self.source};'
		info[0] += f'CLUSTERED_READS_TUMOUR={self.total_read_counts["tumour"]};'
		if 'normal' in self.total_read_counts:
			info[0] += f'CLUSTERED_READS_NORMAL={self.total_read_counts["normal"]};'
		info[0] += stats_str
		if self.breakpoint_notation == "<INS>":
			info[0] = 'SVTYPE=INS;' + info[0]
			for label, depths in self.phased_local_depths.items():
				for i, bin_label in enumerate(["BEFORE","AT","AFTER"]):
					sum_across_haplotypes = sum([counts[0] for hp, counts in depths[i].items()])
					info[0] += f'{label.upper()}_DP_{bin_label}={sum_across_haplotypes},{sum_across_haplotypes};'
			info[0] = info[0] + f'TUMOUR_AF={",".join([str(af) for af in self.allele_fractions["tumour"]])};'
			if 'normal' in self.allele_fractions:
				info[0] = info[0] + f'NORMAL_AF={",".join([str(af) for af in self.allele_fractions["normal"]])}'
			for label, counts in self.phased_local_depths.items():
				bp_0_phase_at = f'{",".join([str(v[0]) for k,v in counts[1].items()])}'
				info[0] += f';{label.upper()}_TOTAL_HP_AT={bp_0_phase_at};'
		elif self.breakpoint_notation in ["+","-"]:
			info[0] = 'SVTYPE=SBND;' + info[0]
			for label, depths in self.phased_local_depths.items():
				for i, bin_label in enumerate(["BEFORE","AT","AFTER"]):
					sum_across_haplotypes = sum([counts[0] for hp, counts in depths[i].items()])
					info[0] += f'{label.upper()}_DP_{bin_label}={sum_across_haplotypes},{sum_across_haplotypes};'
			info[0] = info[0] + f'TUMOUR_AF={",".join([str(af) for af in self.allele_fractions["tumour"]])};'
			if 'normal' in self.allele_fractions:
				info[0] = info[0] + f'NORMAL_AF={",".join([str(af) for af in self.allele_fractions["normal"]])}'
			for label, counts in self.phased_local_depths.items():
				bp_0_phase_at = f'{",".join([str(v[0]) for k,v in counts[1].items()])}'
				info[0] += f';{label.upper()}_TOTAL_HP_AT={bp_0_phase_at};'
		else:
			info.append(info[0]) # duplicate info
			# add edge-specific info
			info[0] = f'SVTYPE=BND;MATEID=ID_{self.count}_2;' + info[0]
			info[1] = f'SVTYPE=BND;MATEID=ID_{self.count}_1;' + info[1]
			for label, depths in self.phased_local_depths.items():
				for i, bin_label in enumerate(["BEFORE","AT","AFTER"]):
					sum_across_haplotypes_0 = sum([counts[0] for hp, counts in depths[i].items()])
					sum_across_haplotypes_1 = sum([counts[1] for hp, counts in depths[i].items()])
					info[0] += f'{label.upper()}_DP_{bin_label}={sum_across_haplotypes_0},{sum_across_haplotypes_1};'
					info[1] += f'{label.upper()}_DP_{bin_label}={sum_across_haplotypes_1},{sum_across_haplotypes_0};'
			info[0] = info[0] + f'TUMOUR_AF={",".join([str(af) for af in self.allele_fractions["tumour"]])};'
			if 'normal' in self.allele_fractions:
				info[0] = info[0] + f'NORMAL_AF={",".join([str(af) for af in self.allele_fractions["normal"]])}'
			info[1] = info[1] + f'TUMOUR_AF={",".join([str(af) for af in reversed(self.allele_fractions["tumour"])])};'
			if 'normal' in self.allele_fractions:
				info[1] = info[1] + f'NORMAL_AF={",".join([str(af) for af in reversed(self.allele_fractions["normal"])])}'
			for label, counts in self.phased_local_depths.items():
				bp_0_phase_at = f'{",".join([str(v[0]) for k,v in counts[1].items()])}'
				info[0] += f';{label.upper()}_TOTAL_HP_AT={bp_0_phase_at};'
				bp_1_phase_at = f'{",".join([str(v[1]) for k,v in counts[1].items()])}'
				info[1] += f';{label.upper()}_TOTAL_HP_AT={bp_1_phase_at};'
		for label in sorted(self.aln_support_counts.keys(), reverse=True):
			# first line
			if not self.phasing or label not in self.phasing[0]:
				info[0] += f'{label.upper()}_ALT_HP=0,0,{self.aln_support_counts[label]};'
			else:
				bp_0_phase = f'{",".join([str(self.phasing[0][label]["HP"][i]) for i in [1, 2, None]])}'
				info[0] += f';{label.upper()}_ALT_HP={bp_0_phase};'
				if self.phasing[0][label]['PS']:
					info[0] += f'{label.upper()}_PS={",".join([str(ps) for ps in self.phasing[0][label]["PS"]])};'
			# second line (applicable)
			if len(info) > 1:
				if not self.phasing or label not in self.phasing[1]:
					info[1] += f'{label.upper()}_ALT_HP=0,0,{self.aln_support_counts[label]};'
				else:
					bp_1_phase = f'{",".join([str(self.phasing[1][label]["HP"][i]) for i in [1, 2, None]])}'
					info[1] += f';{label.upper()}_ALT_HP={bp_1_phase};'
					if self.phasing[1][label]['PS']:
						info[1] += f'{label.upper()}_PS={",".join([str(ps) for ps in self.phasing[1][label]["PS"]])};'

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
		if self.breakpoint_notation not in ["<INS>", "+", "-"]:
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
	__slots__ = 'start_chr', 'start_loc', 'end_chr', 'end_loc', 'spans_cluster', 'inserted_sequence', 'insert_length', 'source', 'read_name', 'mapq', 'read_quality', 'label', 'breakpoint_notation', 'haplotypes', 'phase_sets', 'insert'
	def __init__(self, locations, source, read_name, read_quality, label, breakpoint_notation, haplotypes, phase_sets, insert=None):
		self.start_chr = locations[0]['chr']
		self.start_loc = int(locations[0]['loc'])
		self.end_chr = locations[1]['chr']
		self.end_loc = int(locations[1]['loc'])
		if (breakpoint_notation in ["<INS>", "-", "+"]) and not insert:
			raise AttributeError(f'Must provide an insert for {breakpoint_notation} breakpoints')
		if (breakpoint_notation in ["<INS>", "-", "+"]) and self.start_chr != self.end_chr:
			raise AttributeError(f'{breakpoint_notation} breakpoints must have same start/end location (start_chr={self.start_chr}, end_chr={self.end_chr})')
		if (breakpoint_notation in ["<INS>", "-", "+"]) and self.start_loc != self.end_loc:
			raise AttributeError(f'{breakpoint_notation} breakpoints must have same start/end location (start_loc={self.start_loc}, end_chr={self.end_loc})')
		self.spans_cluster = True if (self.start_chr == self.end_chr and abs(self.start_loc - self.end_loc) <= 300 and breakpoint_notation not in ["<INS>", "-", "+"]) else False
		self.breakpoint_notation = breakpoint_notation
		self.insert_length = len(insert) if insert else 0
		self.inserted_sequence = zlib.compress(insert.encode()) if insert else None
		if source not in ['SUPPLEMENTARY', 'CIGAR', 'SOFTCLIP']:
			raise AttributeError(f'Invalid PotentialBreakpoint souce: {source}')
		self.source = source
		self.read_name = read_name
		self.mapq = read_quality
		self.label = label
		self.haplotypes = haplotypes if (haplotypes[0] or haplotypes[1]) else None
		self.phase_sets = phase_sets if (phase_sets[0] or phase_sets[1]) else None

	def as_dict(self):
		""" return dict representation of breakpoint"""
		self_dict = {
			"start_chr": self.start_chr,
			"start_loc": self.start_loc,
			"end_chr": self.end_chr,
			"end_loc": self.end_loc,
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
		self.chr = initial_breakpoint.start_chr
		self.start = initial_breakpoint.start_loc
		# use both start and end if initial breakpoint start/end are nearby - separate otherwise
		self.end = initial_breakpoint.end_loc if initial_breakpoint.spans_cluster else self.start
		self.source = initial_breakpoint.source
		self.breakpoints = [initial_breakpoint]
		self.supporting_reads = {initial_breakpoint.read_name}
		self.stats = None

	def overlaps(self, other, buffer):
		""" determine if a breakpoint should be merged to a cluster """
		assert self.start <= other.start_loc,"Overlap check requires breakpoints to be sorted."
		# consider whether the cluster is forward or reverse
		reverse = None
		if self.start <= self.end:
			# it's forward, add the buffer
			#      	 start               end
			# ---------|------------------|-----------
			#      	 start                 end_buff
			# ---------|----------------------|-------
			# then we check if end_buffered >= start
			cluster_end_buffered = self.end + buffer
			reverse = False
		else:
			# it's reverse, subtrack the buffer instead
			#      	  end               start
			# ---------|------------------|-----------
			#   end_buff                start
			# -----|----------------------|-----------
			# then we check if end_buffered <= start
			cluster_end_buffered = self.end - buffer
			reverse = True
		if reverse and cluster_end_buffered <= other.start_loc:
			return True
		elif not reverse and cluster_end_buffered >= other.start_loc:
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
		if new_breakpoint.spans_cluster:
			if self.start <= self.end:
				self.end = new_breakpoint.end_loc if (new_breakpoint.end_loc > self.end) else self.end
			else:
				self.end = new_breakpoint.end_loc if (new_breakpoint.end_loc < self.end) else self.end
		else:
			if self.start <= self.end:
				self.end = 	new_breakpoint.start_loc if (new_breakpoint.start_loc > self.end) else self.end
			else:
				self.end = 	new_breakpoint.start_loc if (new_breakpoint.start_loc < self.end) else self.end
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
				if bp.breakpoint_notation in ["<INS>", "+", "-"]:
					event_sizes.append(bp.insert_length)
				elif (bp.start_chr == bp.end_chr):
					event_sizes.append(abs(bp.start_loc - bp.end_loc))
				else:
					event_sizes.append(0)
			stat_dict = {
				'starts_std_dev': pstdev(starts),
				'mapq_mean': mean(mapqs),
				'event_size_std_dev': pstdev(event_sizes),
				'event_size_median': median(event_sizes),
				'event_size_mean': mean(event_sizes)
			}
			# round to 3 decimal places for stats dict attribute
			self.stats = {}
			for key, value in stat_dict.items():
				try:
					self_value = round(value, 3)
				except Exception as _:
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
