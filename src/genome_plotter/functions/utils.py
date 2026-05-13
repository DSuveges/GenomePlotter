"""Scaling genomics to screen coordinates."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Scaler:
    """Tool to scale genomic to screen coordinates."""

    gene_start: int
    chunk_size: int
    pixel_size: int

    def scale_genomic_coordinate(self: Scaler, coordinate: int) -> int:
        """Scale genomic coordinate.

        Args:
            coordinate (int): coordinate to scale given the init input

        Returns:
            int: scaled coordinate
        """
        return round((coordinate - self.gene_start) / self.chunk_size * self.pixel_size)
