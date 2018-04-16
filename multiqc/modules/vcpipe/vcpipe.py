#!/usr/bin/env python

"""MultiQC module to parse output from OUS variant calling pipeline"""

from __future__ import print_function
from collections import OrderedDict, Counter
import logging
import json
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table
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
        self.vcpipe_lowcoverage_data = dict(depth=dict(), size=dict())
        self.vcpipe_basequality_data = dict()
        n = dict()

        re_sname = re.compile("Diag-(.+?)\/")
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

        self.general_stats_addcols(self.general_stats_data, self.general_stats_header)

        # print_vcf = sum([s["status"] for s in self.vcpipe_vcfquality_data.values()]) < len(self.vcpipe_vcfquality_data)
        if len(self.vcpipe_vcfquality_data) > 0:  # and print_vcf:
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
                "format": "{:.0f}",
                "scale": "RdPu"
            }
            vq_headers["tiTvRatio"] = {
                "title": "Ti/Tv Ratio",
                "scale": "Oranges",
                "format": "{:.03f}",
                "min": 0,
                "max": 5
            }
            vq_headers["contamination"] = {
                "title": "Contamination",
                "format": "{:.0f}",
                "scale": "RdYlGn-rev"
            }

            vq_plot = table.plot(self.vcpipe_vcfquality_data, vq_headers)
            self.add_section(
                name="VCF Quality",
                anchor="vcpipe-vcfquality",
                plot=vq_plot
            )

        # print_lowcoverage = sum([s["status"] for s in self.vcpipe_lowcoverage_data.values()]) < len(self.vcpipe_lowcoverage_data)
        if len(self.vcpipe_lowcoverage_data) > 0:  # and print_lowcoverage:
            log.info(json.dumps(self.vcpipe_lowcoverage_data))
            self.write_data_file(self.vcpipe_lowcoverage_data, 'multiqc_vcpipe_lowcoverage')
            lc_config = {
                "id": "vcpipe-lowcoverage",
                "title": "vcpipe: Regions with Low Coverage by Depth",
                # "cpswitch_c_active": False,
                "ylab": "Read Coverage of Regions",
                "data_labels": ["Regions by Depth", "Regions by Size"]
                # "cpswitch": False
            }
            # lc_cats = [str(x) for x in range(1, 24)] + ["X", "Y", "M"]
            # lc_cats = [str(x) for x in range(1, 30)]
            lc_cats = [
                OrderedDict([(x, {"name": x}) for x in range(1, 21)]),
                ["<10", "10-19", "20-29", "30-49", "50-99", "100+"]
            ]

            lc_plot = bargraph.plot(
                [self.vcpipe_lowcoverage_data["depth"], self.vcpipe_lowcoverage_data["size"]],
                cats=lc_cats,
                pconfig=lc_config
            )
            # lc_plot = bargraph.plot(self.vcpipe_lowcoverage_data, pconfig=lc_config)
            self.add_section(
                name="Low Coverage",
                anchor="vcpipe-lowcoverage",
                plot=lc_plot
            )

        # print_basequality = sum([s["status"] for s in self.vcpipe_basequality_data.values()]) < len(self.vcpipe_basequality_data)
        if len(self.vcpipe_basequality_data) > 0:  # and print_basequality:
            self.write_data_file(self.vcpipe_basequality_data, 'multiqc_vcpipe_basequality')
            bq_headers = OrderedDict()
            bq_headers["q30_bases_pct"] = {
                "title": "Base Quality (q30)",
                "description": "percentage of bases &ge; 30x coverage",
                "scale": "OrRd-rev",
                "suffix": "%",
                "format": "{:.02f}",
                "min": 0,
                "max": 100
            }
            bq_headers["reads"] = {"title": "# of Reads"}
            bq_headers["perfect_index_reads_pct"] = {
                "title": "Perfect Index Reads",
                "format": "{:.02f}",
                "suffix": "%",
                "hidden": True
            }
            bq_headers["mean_qual_score"] = {
                "title": "Mean Quality Score",
                "format": "{:.02f}",
                "hidden": True
            }
            bq_headers["one_mismatch_index_pct"] = {
                "title": "One Mismatch Index",
                "format": "{:.02f}",
                "suffix": "%",
                "hidden": True
            }
            bq_plot = table.plot(self.vcpipe_basequality_data, headers=bq_headers)

            self.add_section(
                name="Base Quality",
                anchor="vcpipe-basequality",
                plot=bq_plot
            )

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
            n += self.parse_vcfquality(log_file['s_name'], *raw_data['VcfQuality'])
        else:
            log.error("Could not find expected key 'VcfQuality' in file '{}'".format(log_file["root"]))

        if raw_data.get('LowCoverageRegions') is not None:
            n += self.parse_lowcoverage(log_file['s_name'], *raw_data['LowCoverageRegions'])
        else:
            log.error("Could not find expected key 'LowCoverageRegions' in file '{}'".format(log_file["root"]))

        if raw_data.get('BaseQuality') is not None:
            n += self.parse_basequality(log_file['s_name'], *raw_data['BaseQuality'])
        else:
            log.error("Could not find expected key 'BaseQuality' in file '{}'".format(log_file["root"]))

        return n

    def parse_vcfquality(self, s_name, status, data):
        log.debug("Parsing VcfQuality section for {}".format(s_name))

        if 'VcfQuality' not in self.general_stats_header:
            self.general_stats_header["VcfQuality"] = {"title": "VCF Quality"}

        self.general_stats_data[s_name]["VcfQuality"] = status
        self.vcpipe_vcfquality_data[s_name] = data.copy()
        self.vcpipe_vcfquality_data[s_name]["status"] = status

        return 1

    def parse_lowcoverage(self, s_name, status, data):
        log.debug("Parsing LowCoverageRegions section for {}".format(s_name))

        if 'LowCoverageRegions' not in self.general_stats_header:
            # self.general_stats_header["LowCoverageRegions"] = {
            #     "title": "Low Coverage Regions",
            #     "scale": "RdYlGn-rev"
            # }
            self.general_stats_header["LowCoverageRegions"] = {"title": "Low Coverage Regions"}

        self.general_stats_data[s_name]['LowCoverageRegions'] = status
        if data.get("failed") is not None:
            # self.general_stats_data[s_name]['LowCoverageRegions'] = len(data["failed"])
            depth_counts = Counter([r["depth"] for r in data["failed"]])
            region_size = Counter([region_cat(r["start"], r["stop"]) for r in data["failed"]])
            self.vcpipe_lowcoverage_data["depth"][s_name] = dict(depth_counts)
            self.vcpipe_lowcoverage_data["size"][s_name] = dict(region_size)
        else:
            log.debug("No depth failures found for {}".format(s_name))
            return 0

        return 1

    def parse_basequality(self, s_name, status, data):
        log.debug("Parsing BaseQuality section for {}".format(s_name))

        if 'BaseQuality' not in self.general_stats_header:
            # self.general_stats_header["BaseQuality"] = {
            #     "title": "Base Quality (q30)",
            #     "description": "percentage of bases >= 30x coverage",
            #     "scale": "OrRd-rev",
            #     "min": 0,
            #     "max": 100,
            #     "suffix": "%"
            # }
            self.general_stats_header["BaseQuality"] = {"title": "Base Quality"}

        # self.general_stats_data[s_name]['BaseQuality'] = data["q30_bases_pct"]
        self.general_stats_data[s_name]['BaseQuality'] = status
        self.vcpipe_basequality_data[s_name] = data.copy()

        return 1


def region_cat(start, end):
    delta = int(end) - int(start) + 1
    if delta < 10:
        return "<10"
    elif delta < 20:
        return "10-19"
    elif delta < 30:
        return "20-29"
    elif delta < 50:
        return "30-49"
    elif delta < 100:
        return "50-99"
    else:
        return "100+"
