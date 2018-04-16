#!/usr/bin/env python

"""MultiQC module to parse output from OUS variant calling pipeline"""

from __future__ import print_function
from collections import OrderedDict
import logging
import json
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, table
import re

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        super(MultiqcModule, self).__init__(
            name='vcpipe',
            anchor='vcpipe',
            href='https://git.ousamg.io/apps/vcpipe',
            info=' is the variant calling pipeline for Oslo University Hospital'
        )

        self.general_stats_header = OrderedDict()
        self.general_stats_data = dict()
        self.vcpipe_seen = set()
        self.vcpipe_vcfquality_data = dict()
        self.vcpipe_lowcoverage_data = dict()
        self.vcpipe_basequality_data = dict()
        n = dict()

        re_sname = re.compile("\/Diag-(.+?)\/")
        for f in self.find_log_files("vcpipe"):
            sname_match = re.search(pattern=re_sname, string=f["root"])
            if sname_match:
                s_name = sname_match.groups()[0]
                f['s_name'] = self.clean_s_name(s_name, f['root'])
            else:
                raise(ValueError("Couldn't find sample name in path '{}' with re '{}'".format(f["root"], re_sname)))
            n = self.parse_log(f)
            if n > 0:
                self.add_data_source(f)

        log.debug(json.dumps(self.general_stats_header, indent=4))
        self.general_stats_addcols(self.general_stats_data, self.general_stats_header)

        # difference between sum of status and total entries is number of QC failures in that key
        # only show by default if there is a failure in that section
        print_vcf = sum([s["status"] for s in self.vcpipe_vcfquality_data.values()]) < len(self.vcpipe_vcfquality_data)
        if len(self.vcpipe_vcfquality_data) > 0 and print_vcf:
            self.write_data_file(self.vcpipe_vcfquality_data, 'multiqc_vcpipe_vcfquality')
            vq_headers = OrderedDict()
            vq_headers["skewedRatio"] = {
                "title": "Skewed Ratio",
                "min": 0,
                "max": 0.15,
                "scale": "OrRd",
                "format": '{:,.06f}'
            }
            vq_headers["variants"] = {
                "title": "# of Variants",
                "scale": "YlOrBr"
            }
            vq_headers["tiTvRatio"] = {
                "title": "Ti/Tv Ratio",
                "scale": "BrBG",
                "min": 0,
                "max": 5
            }
            vq_headers["contamination"] = {
                "title": "Contamination",
                "scale": "RdYlGn-rev"
            }

            vq_plot = table.plot(self.vcpipe_vcfquality_data, vq_headers)
            self.add_section(
                name="VCF Quality",
                anchor="vcpipe-vcfquality",
                plot=vq_plot
            )

        # print_lowcoverage = sum([s["status"] for s in self.vcpipe_lowcoverage_data.values()]) < len(self.vcpipe_lowcoverage_data)
        # if len(self.vcpipe_lowcoverage_data) > 0 and print_lowcoverage:
        #     self.write_data_file(self.vcpipe_lowcoverage_data, 'multiqc_vcpipe_lowcoverage')

        # if len(self.vcpipe_basequality_data) > 0:
        #     self.write_data_file(self.vcpipe_basequality_data, 'multiqc_vcpipe_basequality')

    def parse_log(self, log_file):
        n = 0
        raw_data = json.loads(log_file["f"])
        log.debug("Processing file {}".format(log_file["root"]))

        if log_file['s_name'] in self.vcpipe_seen:
            log.debug("Duplicate sample name found in {}! Overwriting: {}".format(log_file['root'], log_file['s_name']))
        else:
            self.vcpipe_seen.add(log_file['s_name'])
            self.general_stats_data[log_file['s_name']] = dict()

        if raw_data.get('VcfQuality') is not None:
            self.parse_vcf_quality(log_file['s_name'], *raw_data['VcfQuality'])
            n += 1
        else:
            log.error("Could not find expected key 'VcfQuality' in file '{}'".format(log_file["root"]))

        if raw_data.get('LowCoverageRegions') is not None:
            self.parse_low_coverage(log_file['s_name'], *raw_data['LowCoverageRegions'])
            n += 1
        else:
            log.error("Could not find expected key 'LowCoverageRegions' in file '{}'".format(log_file["root"]))

        if raw_data.get('BaseQuality') is not None:
            self.parse_base_quality(log_file['s_name'], *raw_data['BaseQuality'])
            n += 1
        else:
            log.error("Could not find expected key 'BaseQuality' in file '{}'".format(log_file["root"]))

        # return final_data
        return n

    def parse_vcf_quality(self, s_name, status, data):
        log.debug("Parsing VcfQuality section for {}".format(s_name))

        if 'VcfQuality' not in self.general_stats_header:
            self.general_stats_header["VcfQuality"] = {"title": "VCF Quality"}

        self.general_stats_data[s_name]["VcfQuality"] = status
        self.vcpipe_vcfquality_data[s_name] = data.copy()
        self.vcpipe_vcfquality_data[s_name]["status"] = status

    def parse_low_coverage(self, s_name, status, data):
        log.debug("Parsing LowCoverageRegions section for {}".format(s_name))

        if 'LowCoverageRegions' not in self.general_stats_header:
            self.general_stats_header["LowCoverageRegions"] = {
                "title": "Low Coverage Regions",
                "scale": "RdYlGn-rev"
            }

        if data.get("failed") is not None:
            self.general_stats_data[s_name]['LowCoverageRegions'] = len(data["failed"])
            # self.vcpipe_lowcoverage_data[s_name]
        else:
            self.general_stats_data[s_name]['LowCoverageRegions'] = 0

    def parse_base_quality(self, s_name, status, data):
        log.debug("Parsing BaseQuality section for {}".format(s_name))

        if 'BaseQuality' not in self.general_stats_header:
            self.general_stats_header["BaseQuality"] = {
                "title": "Base Quality (q30)",
                "description": "percentage of bases >= 30x coverage",
                "scale": "OrRd-rev",
                "min": 0,
                "max": 100,
                "suffix": "%"
            }

        self.general_stats_data[s_name]['BaseQuality'] = data["q30_bases_pct"]
