#!/usr/bin/env python3
"""
PipelineOutput.py

This module contains the PipelineOutput class. This class contains information about the state of the pipeline and
the outputs associated.
"""
from dataclasses import dataclass
from typing import *

@dataclass
class PipelineOutput:
    """
    This object will structure outputs from the pipeline.
    """
    error: str = None
    inchi_id: str = None
    product: str = None
    reaction_names: List[str] = None
    reactants: List[Tuple[str]] = None
    hits: List[str] = None
    template_path: str = None
    output_dir: str = None
    route_uuid: str = None
    time_elapsed: float = None
    date_completed: str = None
    additional_info: Dict[str, Any] = None
    source_of_route: str = None
    to_hippo_output_path: str = None
    num_placed: int = 0
    num_successful: int = 0
    previous_output_path: str = None
    new_input_from_errors: bool = False




