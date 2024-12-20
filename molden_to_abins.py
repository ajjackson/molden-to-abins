#! /usr/bin/env python3

from argparse import ArgumentParser
from io import TextIOBase
import json
from pathlib import Path
import re
from typing import Dict, List, Literal, Tuple, TypedDict

from ase.data import atomic_masses, chemical_symbols
from ase.units import Bohr
import numpy as np


def get_parser() -> ArgumentParser:
    parser = ArgumentParser()
    parser.add_argument("filename", type=Path, help="Input file in Molden format")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Path/filename for output JSON file",
    )

    return parser


def _check_first_line(fd: TextIOBase) -> None:
    if not fd.readline().strip() == "[Molden Format]":
        raise ValueError(
            "File is missing [Molden Format] header. Are you sure this is a Molden file?"
        )


def _read_blocks(fd: TextIOBase) -> Dict[str, List[str]]:
    raw_data: Dict[str, List[str]] = {}

    header_re = re.compile(r"\[[\w\-]+\]\s*(AU|Angs)?")

    current_header = None
    current_block = None

    for line in map(str.strip, fd):
        if header_re.match(line):
            if current_header:
                raw_data[current_header] = current_block
            current_header = line
            current_block = []
        else:
            current_block.append(line)

    # Store the last section, as loop is broken at file end
    raw_data[current_header] = current_block

    return raw_data


class AtomData(TypedDict):
    coord: List[float]
    mass: float
    sort: int
    symbol: str


def parse_atoms_data(raw_data: Dict[str, List[str]]) -> Dict[str, AtomData]:
    """Get atomic positions from [FR-COORD] and convert to Angstrom"""
    raw_atom_lines = raw_data["[FR-COORD]"]

    atoms_data = {
        f"atom_{index}": AtomData(
            coord=coord,
            mass=atomic_masses[chemical_symbols.index(symbol)],
            sort=(index),
            symbol=symbol,
        )
        for index, (symbol, *coord) in enumerate(map(_parse_atom_line, raw_atom_lines))
    }

    for atom_data in atoms_data.values():
        atom_data["coord"] = _bohr_to_ang(atom_data["coord"])

    return atoms_data


AtomicDisplacements = List[List[List[List[float]]]]


class KPointsData(TypedDict):
    frequencies: List[List[float]]
    atomic_displacements: AtomicDisplacements
    weights: List[float]
    k_vectors: List[List[float]]
    unit_cell: List[List[float]]


def _parse_displacements(
    raw_data: Dict[str, List[str]]
) -> List[List[List[List[float]]]]:
    """Get displacements into abins-friendly nested array format

    The array indices are (kpt, atom, mode, axis)

    Note that eigenvectors are _usually_ arranged (mode, atom, axis); the Abins
    ordering is optimised for iterating over incoherent atom contributions...

    """

    modes_array = []
    current_mode_data = []

    for line in raw_data["[FR-NORM-COORD]"]:
        if "vibration" in line:
            if current_mode_data:
                modes_array.append(current_mode_data)
                current_mode_data = []
        else:
            # Alternate values with 0. for imaginary component of eigenvector
            x, y, z = line.split()
            current_mode_data.append([float(x), 0.0, float(y), 0.0, float(z), 0.0])

    modes_array.append(current_mode_data)

    # Insert outermost (k-point) axis: we only have one k-point here
    modes_array = np.asarray(modes_array, dtype=float)
    modes_array = np.expand_dims(modes_array, 0)

    # Swap atom and mode indices
    modes_array = modes_array.swapaxes(1, 2)

    return modes_array.tolist()


def normalise_displacements(
    displacements: AtomicDisplacements, atoms_data: dict[str, AtomData]
) -> AtomicDisplacements:
    """Apply normalisation to Abins eigenvector convention"""

    complex_disp = np.asarray(displacements).view(complex).copy()
    masses = np.asarray(
        [atoms_data[f"atom_{i}"]["mass"] for i in range(complex_disp.shape[1])],
        dtype=float,
    )

    renormalised_complex_disp = _normalise_displacements(complex_disp, masses)
    return renormalised_complex_disp.view(float).tolist()


def _normalise_displacements(
    displacements: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    """Renormalise modes, like Abins GAUSSIAN loader"""

    # Sum over squares of cartesian displacement directions
    b_traces = np.einsum("ijkl,ijkl -> ijk", displacements, displacements)

    # Multiply by masses and sum over atoms to obtain mode normalisation
    norm = np.einsum("ijk, j->ik", b_traces, masses)

    # Scale along mode and atom axes
    displacements = np.einsum("ijkl,ik->ijkl", displacements, 1.0 / np.sqrt(norm))
    displacements = np.einsum("ijkl,j->ijkl", displacements, np.sqrt(masses))

    return displacements


def parse_k_points_data(raw_data: Dict[str, List[str]]) -> KPointsData:
    # I think Molden uses recip. cm like Abins, but not obvious from docs
    frequencies = [list(map(float, raw_data["[FREQ]"]))]

    displacements = _parse_displacements(raw_data)

    return KPointsData(
        k_vectors=[[0.0, 0.0, 0.0]],
        unit_cell=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        weights=[1.0],
        frequencies=frequencies,
        atomic_displacements=displacements,
    )


def _bohr_to_ang(coord: List[float]) -> List[float]:
    return [x * Bohr for x in coord]


def _parse_atom_line(line: str) -> Tuple[str, float, float, float]:
    symbol, x, y, z = line.split()
    return symbol, float(x), float(y), float(z)


def main() -> None:
    args = get_parser().parse_args()

    with open(args.filename, "r") as fd:
        _check_first_line(fd)

        raw_data = _read_blocks(fd)

    atoms_data = parse_atoms_data(raw_data)
    k_points_data = parse_k_points_data(raw_data)

    k_points_data["atomic_displacements"] = normalise_displacements(
        displacements=k_points_data["atomic_displacements"], atoms_data=atoms_data
    )

    output_data = {
        "k_points_data": k_points_data,
        "atoms_data": atoms_data,
        "__abins_class__": "AbinsData",
        "__mantid_version__": "6.9",
    }

    if args.output:
        with open(args.output, "wt") as fd:
            json.dump(output_data, fd, indent=4)
    else:
        print(json.dumps(output_data, indent=4))


if __name__ == "__main__":
    main()
