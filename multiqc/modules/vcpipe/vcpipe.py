#!/usr/bin/env python

"""MultiQC module to parse output from OUS variant calling pipeline"""

from __future__ import print_function
# from collections import OrderedDict
import logging
import json
from multiqc.modules.base_module import BaseMultiqcModule

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        super(MultiqcModule, self).__init__(
            name='vcpipe',
            anchor='vcpipe',
            href='https://git.ousamg.io/apps/vcpipe',
            info='Variant calling pipeline in Oslo University Hospital'
        )

        for f in self.find_log_files("vcpipe", filehandles=True):
            self.mod_data[f['s_name']] = self.parse_log(f)
            self.add_data_source(f)

    def parse_log(log_fh):
        data = json.load(log_fh)
