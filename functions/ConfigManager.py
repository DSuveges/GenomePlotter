"""Managing configuration settings for the Genome Plotter."""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from typing import Optional


# List of dataclasses to describe the configuration file:
@dataclass
class PlotParameters:
    """Dataclass to store plot parameters."""

    width: int
    pixel_size: int
    dark_start: float
    dark_max: float
    custom_gene_window: int


@dataclass
class BasicParameters:
    """Dataclass to store basic parameters."""

    chunk_size: int
    missing_tolerance: Optional[float] = None
    plot_folder: Optional[str] = None
    data_folder: Optional[str] = None
    row_length: Optional[int] = None


@dataclass
class ColorSchema:
    """Dataclass to store color schema."""

    chromosome_colors: dict[str, str]
    cytoband_colors: dict[str, str]
    gwas_point: str
    arrow_colors: dict[str, str]


@dataclass
class CytoBandData:
    """Dataclass to store cytoband data."""

    url: str
    processed_file: str
    genome_build: Optional[str] = None


@dataclass
class SourcePrototype:
    """Dataclass to store source prototype."""

    # Output file is mandatory:
    processed_file: str
    # Optional parameters:
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


@dataclass
class SourceData:
    """Dataclass to store source data."""

    cytoband_data: CytoBandData
    ensembl_data: SourcePrototype
    gencode_data: SourcePrototype
    gwas_data: SourcePrototype

    def __post_init__(self: SourceData) -> None:
        """Convert values to the appropriate data types."""
        for field in self.__dataclass_fields__.keys():
            if isinstance(field_type := self.__dataclass_fields__[field].type, str):
                field_type = globals()[field_type]
            else:
                field_type = self.__dataclass_fields__[field].type

            self.__setattr__(
                field,
                field_type(**self.__getattribute__(field)),
            )


@dataclass
class Config:
    """Dataclass to store the configuration file."""

    plot_parameters: PlotParameters
    basic_parameters: BasicParameters
    color_schema: ColorSchema
    source_data: SourceData

    def __post_init__(self: Config) -> None:
        """Convert values to the appropriate data types."""
        for field in self.__dataclass_fields__.keys():
            if isinstance(field_type := self.__dataclass_fields__[field].type, str):
                field_type = globals()[field_type]
            else:
                field_type = self.__dataclass_fields__[field].type

            self.__setattr__(
                field,
                field_type(**self.__getattribute__(field)),
            )

    # Saving the configuration file:
    def save(self: Config, file_path: str) -> None:
        """Save the configuration file.

        Args:
            file_path (str): Path to the configuration file.
        """
        with open(file_path, "w") as file:
            json.dump(asdict(self), file, indent=3)

    # Updating basic configuration based on command line arguments:
    def update_basic_parameters(self: Config, **kwargs) -> None:  # type: ignore[no-untyped-def]
        """Update basic parameters based on command line arguments.

        Args:
            **kwargs: Command line arguments.
        """
        for key, value in kwargs.items():
            if key in self.basic_parameters.__dataclass_fields__:
                self.basic_parameters.__setattr__(key, value)

    def get_cytoband_file(self: Config) -> str:
        """Get the cytoband file.

        Returns:
            str: Path to the cytoband file if exists.

        Raises:
            ValueError: If the cytoband file does not exist.
        """
        file = self.source_data.cytoband_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"

        if not os.path.isfile(full_path):
            raise ValueError(f"Cytological band file ({full_path}) doesn't exists.")

        return full_path

    def get_chromosome_file(self: Config, chromsome: str) -> str:
        """Get the chromosome file.

        Args:
            chromsome (str): The chromosome number.

        Returns:
            str: Path to the chromosome file if exists.

        Raises:
            ValueError: If the chromosome file does not exist.
        """
        file = self.source_data.ensembl_data.processed_file
        file = file.format(chromsome)
        full_path = f"{self.basic_parameters.data_folder}/{file}"

        if not os.path.isfile(full_path):
            raise ValueError(f"The requested genome file ({full_path}) doesn't exists.")

        return full_path

    def get_gencode_file(self: Config) -> str:
        """Get the GENCODE file.

        Returns:
            str: Path to the GENCODE file if exists.

        Raises:
            ValueError: If the GENCODE file does not exist.
        """
        file = self.source_data.gencode_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"

        if not os.path.isfile(full_path):
            raise ValueError(f"Processed GENCODE file ({full_path}) doesn't exists.")

        return full_path

    def get_gwas_file(self: Config) -> str:
        """Get the GWAS file.

        Returns:
            str: Path to the GWAS file if exists.

        Raises:
            ValueError: If the GWAS file does not exist.
        """
        file = self.source_data.gwas_data.processed_file
        full_path = f"{self.basic_parameters.data_folder}/{file}"

        if not os.path.isfile(full_path):
            raise ValueError(f"Processed GWAS file ({full_path}) doesn't exists.")

        return full_path
