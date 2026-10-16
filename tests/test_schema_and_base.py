import json
import math

import pytest

import chemout
from chemout.parse import parse_text
from chemout.parsers.base import Lines, floats_in, symbol_to_z, to_float
from chemout.schema import Calculation, Termination


def test_symbol_to_z():
    assert symbol_to_z("C") == 6
    assert symbol_to_z("c") == 6
    assert symbol_to_z("Cl") == 17
    assert symbol_to_z("CL") == 17
    assert symbol_to_z("Cl3") == 17
    assert symbol_to_z("C1") == 6
    assert symbol_to_z("He") == 2
    assert symbol_to_z("Xx") is None


def test_to_float_fortran_and_overflow():
    assert to_float("-1.2D-03") == pytest.approx(-0.0012)
    assert to_float("1.5e2") == 150.0
    assert math.isnan(to_float("********"))
    assert floats_in("E = -76.123 au  f=0.0238  nm 0.000") == pytest.approx([-76.123, 0.0238, 0.0])
    assert floats_in("1.0D+02 -3.5e-1 7") == pytest.approx([100.0, -0.35, 7.0])


def test_lines_cursor():
    ls = Lines(["a", "b", "", "c"])
    assert ls.next() == "a"
    assert ls.peek() == ""
    assert ls.find("c")
    assert ls.current == "c"
    ls.next()
    assert ls.done
    assert ls.current == ""  # never raises past EOF


def test_unknown_program_gives_warning_not_exception():
    calc = parse_text("hello\nthis is not a QC output\n")
    assert calc.termination is Termination.UNKNOWN
    assert calc.provenance.warnings[0].section == "sniff"


def test_broken_section_is_a_warning_not_a_crash():
    # An ORCA file whose Mulliken block is corrupted mid-table must still yield the energy.
    text = """
                                 *****************
                                 * O   R   C   A *
                                 *****************
                         Program Version 6.0.1  -   RELEASE   -
================================================================================
                                       INPUT FILE
================================================================================
|  1> ! B3LYP def2-SVP
|  2> * xyz 0 1
|  3> H 0 0 0
|  4> *
|  5>                          ****END OF INPUT****
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  H      0.000000    0.000000    0.000000

----------------
TOTAL SCF ENERGY
----------------

Total Energy       :         -0.50000000000000 Eh             -13.60569 eV

-----------------------
MULLIKEN ATOMIC CHARGES
-----------------------
   0 H :    garbage here
Sum of atomic charges:    0.0000000

FINAL SINGLE POINT ENERGY        -0.500000000000
                             ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 0 minutes 1 seconds 0 msec
"""
    calc = parse_text(text)
    assert calc.program.name.value == "orca"
    assert calc.termination is Termination.NORMAL
    assert calc.energies.total == pytest.approx(-0.5)
    assert calc.molecule.atomic_numbers == [1]
    assert calc.input.method == "B3LYP"
    assert calc.input.basis == "def2-SVP"
    assert calc.populations == []  # corrupted block produced no bogus charges


def test_json_roundtrip(tmp_path):
    text = "ORCA\n * O   R   C   A *\nProgram Version 5.0.4\n****ORCA TERMINATED NORMALLY****\n"
    calc = parse_text(text)
    js = calc.to_json()
    back = Calculation.from_json(js)
    assert back == calc
    assert json.loads(js)["units"]["energy"] == "hartree"


def test_hdf5_roundtrip(tmp_path):
    pytest.importorskip("h5py")
    from chemout.io.hdf5 import read_hdf5

    p = tmp_path / "x.h5"
    calc = parse_text(
        " Entering Gaussian System, Link 0=g16\n Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019\n SCF Done:  E(RHF) =  -1.0  A.U. after    3 cycles\n Normal termination of Gaussian 16 at Tue May  6 22:05:38 2025.\n"
    )
    calc.to_hdf5(str(p))
    back = read_hdf5(str(p))
    assert back.energies.scf == -1.0
    assert back.program.version == "16 RevC.01"


def test_version_exposed():
    assert chemout.__version__
