#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =========================================================================
#   Program:   S1Processor
#
#   Copyright 2017-2020 (c) CESBIO. All rights reserved.
#
#   This file is part of S1Tiling project
#       https://gitlab.orfeo-toolbox.org/s1-tiling/s1tiling
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# =========================================================================
#
# Authors: Thierry KOLECK (CNES)
#          Luc HERMITTE (CS Group)
#
# =========================================================================

"""
Sub-module that manages decoding of S1Processor options.
"""

import configparser
import copy
import logging
import logging.handlers
# import multiprocessing
import os
import pathlib
import re
import sys
import yaml


def init_logger(mode, paths):
    """
    Initializes logging service.
    """
    # Add the dirname where the current script is
    paths += [pathlib.Path(__file__).parent.parent.absolute()]
    paths = [p / 'logging.conf.yaml' for p in paths]
    cfgpaths = [p for p in paths if p.is_file()]
    # print("from %s, keep %s" % (paths, cfgpaths))

    verbose   = 'debug'   in mode
    log2files = 'logging' in mode
    # print("verbose: ", verbose)
    # print("log2files: ", log2files)
    if cfgpaths:
        with open(cfgpaths[0], 'r') as stream:
            # FullLoader requires yaml 5.1
            # And it SHALL be used, see https://github.com/yaml/pyyaml/wiki/PyYAML-yaml.load(input)-Deprecation
            if hasattr(yaml, 'FullLoader'):
                config = yaml.load(stream, Loader=yaml.FullLoader)
            else:
                print("WARNING - upgrade pyyaml to version 5.1 at least!!")
                config = yaml.load(stream)
        if verbose:
            # Control the maximum global verbosity level
            config["root"]["level"] = "DEBUG"

            # Control the local console verbosity level
            config["handlers"]["console"]["level"] = "DEBUG"
        if log2files:
            if 'file' not in config["root"]["handlers"]:
                config["root"]["handlers"] += ['file']
            if 'important' not in config["root"]["handlers"]:
                config["root"]["handlers"] += ['important']
        main_config = copy.deepcopy(config)
        for _, cfg in main_config['handlers'].items():
            if 'filename' in cfg and '%' in cfg['filename']:
                cfg['filename'] = cfg['filename'] % ('main',)
        logging.config.dictConfig(main_config)
    else:
        # This situation should not happen
        if verbose:
            logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            # os.environ["OTB_LOGGER_LEVEL"]="DEBUG"
        else:
            logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    return config


class Configuration():
    """This class handles the parameters from the cfg file"""
    def __init__(self, configFile):
        config = configparser.ConfigParser(os.environ)
        config.read(configFile)

        # Logs
        self.Mode = config.get('Processing', 'mode')
        self.log_config = init_logger(self.Mode, [pathlib.Path(configFile).parent.absolute()])
        # self.log_queue = multiprocessing.Queue()
        # self.log_queue_listener = logging.handlers.QueueListener(self.log_queue)
        if "debug" in self.Mode and self.log_config['loggers']['s1tiling.OTB']['level'] == 'DEBUG':
            os.environ["OTB_LOGGER_LEVEL"] = "DEBUG"

        # Other options
        self.output_preprocess = config.get('Paths', 'output')
        self.raw_directory     = config.get('Paths', 's1_images')
        self.srtm              = config.get('Paths', 'srtm')
        self.tmpdir            = config.get('Paths', 'tmp')
        if not os.path.isdir(self.tmpdir) and not os.path.isdir(os.path.dirname(self.tmpdir)):
            # Even if tmpdir doesn't exist we should still be able to create it
            logging.critical("ERROR: tmpdir=%s is not a valid path", self.tmpdir)
            sys.exit(1)
        self.GeoidFile         = config.get('Paths', 'geoid_file')
        if config.has_section('PEPS'):
            logging.critical('Since version 2.0, S1Tiling use [DataSource] instead of [PEPS] in config files. Please update your configuration!')
            sys.exit(-1)
        self.eodagConfig       = config.get('DataSource', 'eodagConfig', fallback=None)
        self.download          = config.getboolean('DataSource', 'download')
        self.ROI_by_tiles      = config.get('DataSource', 'roi_by_tiles')
        self.first_date        = config.get('DataSource', 'first_date')
        self.last_date         = config.get('DataSource', 'last_date')
        self.polarisation      = config.get('DataSource', 'polarisation')
        if   self.polarisation == 'VV-VH':
            self.polarisation = 'VV VH'
        elif self.polarisation == 'HH-HV':
            self.polarisation = 'HH HV'
        else:
            logging.critical("Parameter [polarisation] must be HH-HV or VV-VH")
            logging.critical("Please correct it the config file ")
            sys.exit(-1)

        self.type_image         = "GRD"
        self.mask_cond          = config.getboolean('Mask', 'Generate_border_mask')
        self.calibration_type   = config.get('Processing', 'calibration')
        self.removethermalnoise = config.getboolean('Processing', 'remove_thermal_noise')

        self.out_spatial_res    = config.getfloat('Processing', 'output_spatial_resolution')

        self.output_grid        = config.get('Processing', 'tiles_shapefile')
        if not os.path.isfile(self.output_grid):
            logging.critical("ERROR: output_grid=%s is not a valid path", self.output_grid)
            sys.exit(1)

        self.SRTMShapefile = config.get('Processing', 'srtm_shapefile')
        if not os.path.isfile(self.SRTMShapefile):
            logging.critical("ERROR: srtm_shapefile=%s is not a valid path", self.SRTMShapefile)
            sys.exit(-1)
        self.grid_spacing = config.getfloat('Processing', 'orthorectification_gridspacing')
        self.border_threshold = config.getfloat('Processing', 'border_threshold')
        try:
            tiles_file = config.get('Processing', 'tiles_list_in_file')
            self.tile_list = open(tiles_file, 'r').readlines()
            self.tile_list = [s.rstrip() for s in self.tile_list]
            logging.info("The following tiles will be processed: %s", self.tile_list)
        except Exception:  # pylint: disable=broad-except
            tiles = config.get('Processing', 'tiles')
            self.tile_list = [s.strip() for s in re.split(r'\s*,\s*', tiles)]

        self.TileToProductOverlapRatio = config.getfloat('Processing', 'tile_to_product_overlap_ratio')
        self.nb_procs                  = config.getint('Processing', 'nb_parallel_processes')
        self.ram_per_process           = config.getint('Processing', 'ram_per_process')
        self.OTBThreads                = config.getint('Processing', 'nb_otb_threads')
        # self.filtering_activated       = config.getboolean('Filtering', 'filtering_activated')
        # self.Reset_outcore             = config.getboolean('Filtering', 'reset_outcore')
        # self.Window_radius             = config.getint('Filtering', 'window_radius')

    def check_date(self):
        """
        DEPRECATED
        """
        import datetime

        fd = self.first_date
        ld = self.last_date

        try:
            F_Date = datetime.date(int(fd[0:4]), int(fd[5:7]), int(fd[8:10]))
            L_Date = datetime.date(int(ld[0:4]), int(ld[5:7]), int(ld[8:10]))
            return F_Date, L_Date
        except Exception:  # pylint: disable=broad-except
            logging.critical("Invalid date")
            sys.exit()