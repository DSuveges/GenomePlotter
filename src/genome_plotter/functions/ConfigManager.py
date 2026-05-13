"""Managing configuration settings for the Genome Plotter."""

from __future__ import annotations

import os
import re
from typing import Optional

import yaml
from pydantic import BaseModel, ConfigDict, field_validator


class PlotParameters(BaseModel):
    """Parameters controlling chromosome plot dimensions and shading."""

    model_config = ConfigDict(validate_assignment=True)

    width: int
    pixel_size: int
    dark_start: float
    dark_max: float
    custom_gene_window: int


class BasicParameters(BaseModel):
    """Runtime parameters populated via CLI; mostly None in the bundled config."""

    model_config = ConfigDict(validate_assignment=True)

    chunk_size: int
    missing_tolerance: Optional[float] = None
    plot_folder: Optional[str] = None
    data_folder: Optional[str] = None
    row_length: Optional[int] = None


class ColorSchema(BaseModel):
    """Color definitions for all genomic feature categories."""

    chromosome_colors: dict[str, str]
    cytoband_colors: dict[str, str]
    gwas_point: str
    arrow_colors: dict[str, str]

    @field_validator("gwas_point")
    @classmethod
    def valid_hex(cls, v: str) -> str:
        """Validate that gwas_point is a six-digit hex color."""
        if not re.match(r"^#[0-9A-Fa-f]{6}$", v):
            raise ValueError(f"Invalid hex color: {v}")
        return v


class CytoBandData(BaseModel):
    """Source configuration for cytological band data."""

    url: str
    processed_file: str
    genome_build: Optional[str] = None


class SourcePrototype(BaseModel):
    """Generic source configuration for an FTP-fetched dataset."""

    processed_file: str
    url: Optional[str] = None
    genome_build: Optional[str] = None
    host: Optional[str] = None
    path: Optional[str] = None
    source_file: Optional[str] = None
    arrow_file: Optional[str] = None
    release_date: Optional[str] = None
    release: Optional[int] = None
    version_url: Optional[str] = None
    version: Optional[int] = None


class SourceData(BaseModel):
    """Container for all data source configurations."""

    cytoband_data: CytoBandData
    ensembl_data: SourcePrototype
    gencode_data: SourcePrototype
    gwas_data: SourcePrototype


class Config(BaseModel):
    """Root configuration object for the genome plotter."""

    model_config = ConfigDict(validate_assignment=True)

    plot_parameters: PlotParameters
    basic_parameters: BasicParameters
    color_schema: ColorSchema
    source_data: SourceData

    def save(self, file_path: str) -> None:
        """Save the configuration to a YAML file.

        Args:
            file_path (str): Path to the output file.
        """
        with open(file_path, "w") as f:
            yaml.dump(self.model_dump(), f, default_flow_style=False, sort_keys=False)

    def update_basic_parameters(self, **kwargs) -> None:  # type: ignore[no-untyped-def]
        """Update basic parameters from keyword arguments.

        Args:
            **kwargs: Field names and values to update.
        """
        for key, value in kwargs.items():
            if key in type(self.basic_parameters).model_fields:
                setattr(self.basic_parameters, key, value)

    def get_cytoband_file(self) -> str:
        """Return the path to the cytoband file, raising if it does not exist."""
        file = self.source_data.cytoband_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"
        if not os.path.isfile(full_path):
            raise ValueError(f"Cytological band file ({full_path}) doesn't exists.")
        return full_path

    def get_chromosome_file(self, chromsome: str) -> str:
        """Return the path to a per-chromosome file, raising if it does not exist."""
        file = self.source_data.ensembl_data.processed_file.format(chromsome)
        full_path = f"{self.basic_parameters.data_folder}/{file}"
        if not os.path.isfile(full_path):
            raise ValueError(f"The requested genome file ({full_path}) doesn't exists.")
        return full_path

    def get_gencode_file(self) -> str:
        """Return the path to the GENCODE file, raising if it does not exist."""
        file = self.source_data.gencode_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"
        if not os.path.isfile(full_path):
            raise ValueError(f"Processed GENCODE file ({full_path}) doesn't exists.")
        return full_path

    def get_gwas_file(self) -> str:
        """Return the path to the GWAS file, raising if it does not exist."""
        file = self.source_data.gwas_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"
        if not os.path.isfile(full_path):
            raise ValueError(f"Processed GWAS file ({full_path}) doesn't exists.")
        return full_path

    def get_gencode_arrow_file(self) -> str:
        """Return the path to the GENCODE arrow file, raising if it does not exist."""
        file = self.source_data.gencode_data.arrow_file
        if file is None:
            raise ValueError("GENCODE arrow file not configured.")
        full_path = f"{self.basic_parameters.data_folder}/{file}"
        if not os.path.isfile(full_path):
            raise ValueError(f"GENCODE arrow file ({full_path}) doesn't exists.")
        return full_path
