#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =========================================================================
#   Program:   S1Processor
#
#   Copyright 2017-2022 (c) CNES. All rights reserved.
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
import logging.config
import logging.handlers
import os
from pathlib import Path
import re
import sys
import yaml

from s1tiling.libs import exits
from .otbpipeline import otb_version

RESOURCE_DIR = Path(__file__).parent.parent.absolute() / "resources"


logger = logging.getLogger("s1tiling.configuration")


def init_logger(mode, config_file_path: Path):
    """
    Initializes logging service.
    """
    config_file_path = Path(config_file_path)
    assert config_file_path.is_file(), f"Logging configuration file {config_file_path} does not exist"

    verbose = "debug" in mode
    log2files = "logging" in mode
    with config_file_path.open() as stream:
        if hasattr(yaml, "FullLoader"):
            config = yaml.load(stream, Loader=yaml.FullLoader)
    if verbose:
        # Control the maximum global verbosity level
        config["root"]["level"] = "DEBUG"

        # Control the local console verbosity level
        config["handlers"]["console"]["level"] = "DEBUG"
    if log2files:
        if "file" not in config["root"]["handlers"]:
            config["root"]["handlers"] += ["file"]
        if "important" not in config["root"]["handlers"]:
            config["root"]["handlers"] += ["important"]
    main_config = copy.deepcopy(config)
    for _, cfg in main_config["handlers"].items():
        if "filename" in cfg and "%" in cfg["filename"]:
            cfg["filename"] = cfg["filename"] % ("main",)
    logging.config.dictConfig(main_config)
    return config


class Configuration:
    """This class handles the parameters from the cfg file"""

    def __init__(self, proc_config_file, log_config_file: Path):
        config = configparser.ConfigParser(os.environ)
        config.read(proc_config_file)

        # Logs
        self.logging_mode = config.get("Processing", "mode")
        self.log_config = init_logger(self.logging_mode, log_config_file)
        logger.info("Logging initialized")
        if (
            "debug" in self.logging_mode
            and self.log_config
            and self.log_config["loggers"]["s1tiling.OTB"]["level"] == "DEBUG"
        ):
            os.environ["OTB_LOGGER_LEVEL"] = "DEBUG"

        # Other options
        self.output_preprocess = config.get("Paths", "output")
        self.raw_directory = config.get("Paths", "s1_images")
        self.srtm = config.get("Paths", "srtm")
        self.tmpdir = config.get("Paths", "tmp")
        if not os.path.isdir(self.tmpdir) and not os.path.isdir(
            os.path.dirname(self.tmpdir)
        ):
            # Even if tmpdir doesn't exist we should still be able to create it
            logger.critical("ERROR: tmpdir=%s is not a valid path", self.tmpdir)
            sys.exit(exits.CONFIG_ERROR)
        self.GeoidFile = config.get(
            "Paths", "geoid_file", fallback=str(RESOURCE_DIR / "Geoid/egm96.grd")
        )
        if config.has_section("PEPS"):
            logger.critical(
                "Since version 2.0, S1Tiling use [DataSource] instead of [PEPS] in config files. Please update your configuration!"
            )
            sys.exit(exits.CONFIG_ERROR)
        self.eodagConfig = config.get("DataSource", "eodagConfig", fallback=None)
        self.download = config.getboolean("DataSource", "download")
        self.ROI_by_tiles = config.get("DataSource", "roi_by_tiles")
        self.first_date = config.get("DataSource", "first_date")
        self.last_date = config.get("DataSource", "last_date")
        self.polarisation = config.get("DataSource", "polarisation")
        if self.polarisation == "VV-VH":
            self.polarisation = "VV VH"
        elif self.polarisation == "HH-HV":
            self.polarisation = "HH HV"
        elif self.polarisation not in ["VV", "VH", "HH", "HV"]:
            logger.critical(
                "Parameter [polarisation] must be either HH-HV, VV-VH, HH, HV, VV or VH"
            )
            logger.critical("Please correct the config file ")
            sys.exit(exits.CONFIG_ERROR)
        if self.download:
            self.nb_download_processes = config.getint(
                "DataSource", "nb_parallel_downloads", fallback=1
            )

        self.type_image = "GRD"
        self.mask_cond = config.getboolean("Mask", "generate_border_mask")
        self.cache_srtm_by = config.get(
            "Processing", "cache_srtm_by", fallback="symlink"
        )
        if self.cache_srtm_by not in ["symlink", "copy"]:
            logger.critical(
                "Unexpected value for Processing.cache_srtm_by option: '%s' is neither 'copy' no 'symlink'",
                self.cache_srtm_by,
            )
            sys.exit(exits.CONFIG_ERROR)

        self.calibration_type = config.get("Processing", "calibration")
        self.removethermalnoise = config.getboolean(
            "Processing", "remove_thermal_noise"
        )
        if self.removethermalnoise and otb_version() < "7.4.0":
            logger.critical(
                "ERROR: OTB %s does not support noise removal. Please upgrade OTB to version 7.4.0 or disable 'remove_thermal_noise' in '%s'",
                otb_version(),
                proc_config_file,
            )
            sys.exit(exits.CONFIG_ERROR)

        self.out_spatial_res = config.getfloat(
            "Processing", "output_spatial_resolution"
        )

        self.output_grid = config.get(
            "Processing",
            "tiles_shapefile",
            fallback=str(RESOURCE_DIR / "shapefile/Features.shp"),
        )
        if not os.path.isfile(self.output_grid):
            logger.critical(
                "ERROR: output_grid=%s is not a valid path", self.output_grid
            )
            sys.exit(exits.CONFIG_ERROR)

        self._SRTMShapefile = RESOURCE_DIR / "shapefile" / "srtm_tiles.gpkg"

        self.grid_spacing = config.getfloat(
            "Processing", "orthorectification_gridspacing"
        )
        self.interpolation_method = config.get(
            "Processing", "orthorectification_interpolation_method", fallback="nn"
        )
        try:
            tiles_file = config.get("Processing", "tiles_list_in_file")
            self.tile_list = open(tiles_file, "r").readlines()
            self.tile_list = [s.rstrip() for s in self.tile_list]
            logger.info("The following tiles will be processed: %s", self.tile_list)
        except Exception:  # pylint: disable=broad-except
            tiles = config.get("Processing", "tiles")
            self.tile_list = [s.strip() for s in re.split(r"\s*,\s*", tiles)]

        self.TileToProductOverlapRatio = config.getfloat(
            "Processing", "tile_to_product_overlap_ratio"
        )
        self.nb_procs = config.getint("Processing", "nb_parallel_processes")
        self.ram_per_process = config.getint("Processing", "ram_per_process")
        self.OTBThreads = config.getint("Processing", "nb_otb_threads")
        try:
            self.override_azimuth_cut_threshold_to = config.getboolean(
                "Processing", "override_azimuth_cut_threshold_to"
            )
        except Exception:  # pylint: disable=broad-except
            # We cannot use "fallback=None" to handle ": None" w/ getboolean()
            self.override_azimuth_cut_threshold_to = None

        logger.debug("Running S1Tiling with:")
        logger.debug("[Paths]")
        logger.debug("- geoid_file                     : %s", self.GeoidFile)
        logger.debug("- output                         : %s", self.output_preprocess)
        logger.debug("- s1_images                      : %s", self.raw_directory)
        logger.debug("- srtm                           : %s", self.srtm)
        logger.debug("- tmp                            : %s", self.tmpdir)
        logger.debug("[DataSource]")
        logger.debug("- download                       : %s", self.download)
        logger.debug("- first_date                     : %s", self.first_date)
        logger.debug("- last_date                      : %s", self.last_date)
        logger.debug("- polarisation                   : %s", self.polarisation)
        logger.debug("- roi_by_tiles                   : %s", self.ROI_by_tiles)
        if self.download:
            logger.debug(
                "- nb_parallel_downloads          : %s", self.nb_download_processes
            )
        logger.debug("[Processing]")
        logger.debug("- calibration                    : %s", self.calibration_type)
        logger.debug("- mode                           : %s", self.logging_mode)
        logger.debug("- nb_otb_threads                 : %s", self.OTBThreads)
        logger.debug("- nb_parallel_processes          : %s", self.nb_procs)
        logger.debug("- orthorectification_gridspacing : %s", self.grid_spacing)
        logger.debug("- output_spatial_resolution      : %s", self.out_spatial_res)
        logger.debug("- ram_per_process                : %s", self.ram_per_process)
        logger.debug("- remove_thermal_noise           : %s", self.removethermalnoise)
        logger.debug("- srtm_shapefile                 : %s", self._SRTMShapefile)
        logger.debug(
            "- tile_to_product_overlap_ratio  : %s", self.TileToProductOverlapRatio
        )
        logger.debug("- tiles                          : %s", self.tile_list)
        logger.debug("- tiles_shapefile                : %s", self.output_grid)
        logger.debug("[Mask]")
        logger.debug("- generate_border_mask           : %s", self.mask_cond)

    @property
    def srtm_db_filepath(self):
        """
        Get the SRTMShapefile databe filepath
        """
        return str(self._SRTMShapefile)

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
            logger.critical("Invalid date")
            sys.exit(exits.CONFIG_ERROR)
