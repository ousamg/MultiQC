#!/usr/bin/env python

"""MultiQC module to parse output from OUS variant calling pipeline"""

from __future__ import print_function
from collections import OrderedDict
import logging
import json
from multiqc.modules.base_module import BaseMultiqcModule
import re

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    expected_keys = [
        "CalculateHsMetrics",
        "CollectInsertSizeMetrics",
        "VcfQuality",
        "CollectAlignmentSummaryMetrics",
        "LowCoverageRegions",
        "BaseQuality"
    ]
    process_keys = [
        "VcfQuality",
        "LowCoverageRegions",
        "BaseQuality"
    ]

    def __init__(self):

        super(MultiqcModule, self).__init__(
            name='vcpipe',
            anchor='vcpipe',
            href='https://git.ousamg.io/apps/vcpipe',
            info=' is the variant calling pipeline for Oslo University Hospital'
        )

        self.mod_data = dict()
        re_sname = re.compile("\/Diag-(.+?)\/")
        for f in self.find_log_files("vcpipe"):
            sname_match = re.search(pattern=re_sname, string=f["root"])
            if sname_match:
                s_name = sname_match.groups()[0]
                f['s_name'] = self.clean_s_name(s_name, f['root'])
            else:
                raise(ValueError("Couldn't find sample name in path '{}' with re '{}'".format(f["root"], re_sname)))
            self.mod_data[f['s_name']] = self.parse_log(f)
            self.add_data_source(f)

        gen_stats_header = OrderedDict()
        gen_stats_header["VcfQuality"] = {"title": "VCF Quality"}
        gen_stats_header["LowCoverageRegions"] = {
            "title": "Low Coverage Regions",
            "scale": "RdYlGn-rev"
        }
        gen_stats_header["BaseQuality"] = {
            "title": "Base Quality (q30)",
            "description": "percentage of bases >= 30x coverage",
            "scale": "OrRd-rev",
            "min": 0,
            "max": 100,
            "suffix": "%"
        }
        gen_stats_data = {samp: {col: self.mod_data[samp][col]["gen_stats"] for col in self.mod_data[samp].keys()} for samp in self.mod_data.keys()}
        self.general_stats_addcols(gen_stats_data, gen_stats_header)

    def parse_log(self, log_file):
        raw_data = json.loads(log_file["f"])
        final_data = OrderedDict([(k, {}) for k in self.process_keys])
        for k in self.process_keys:
            if k not in raw_data or raw_data[k] is None:
                log.error("Could not find expected key '{}' in file '{}'".format(k, log_file["root"]))
            elif k == "VcfQuality":
                final_data[k]["gen_stats"] = raw_data[k][0]
                final_data[k]["fields"] = raw_data[k][1]
            elif k == "LowCoverageRegions":
                final_data[k]["gen_stats"] = len(raw_data[k][1]["failed"])
                final_data[k]["fields"] = raw_data[k][1]
            elif k == "BaseQuality":
                final_data[k]["gen_stats"] = raw_data[k][1]["q30_bases_pct"]
                final_data[k]["fields"] = raw_data[k][1]
            else:
                log.error("Key '{}' in process_keys, but don't know what to do with it".format(k))

        return final_data
