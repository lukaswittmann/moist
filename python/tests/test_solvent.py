"""Solvent table lookups through the Python API."""

import dataclasses
import io
import math

import pytest

import moist
from moist import SolventDescriptors, get_solvent_descriptors, get_solvent_descriptors_by_id
from moist.library import ffi


def test_fields_mirror_c_struct():
    native = [name for name, _ in ffi.typeof("moist_solvent_descriptors").fields]
    assert native[0] == "struct_size"
    assert [item.name for item in dataclasses.fields(SolventDescriptors)] == native[1:]


def test_water_values():
    water = get_solvent_descriptors("  WATER ")
    assert water.id == 175
    assert water.name == "water"
    # Values the ORCA CPCM(water) references assume
    assert water.epsilon == 78.355
    assert water.refractive_index == 1.3328
    assert water.abraham_alpha == 0.82
    assert water.abraham_beta == 0.35
    assert water.surface_tension_si == pytest.approx(0.07199202032600779, rel=1e-15)
    assert water.mass_density_si == 998.0
    assert water.molar_mass_si == pytest.approx(0.018015, abs=1e-5)
    assert water.mass_density_au == pytest.approx(water.number_density_au * water.molecular_mass_au,
                                                  rel=1e-12)
    assert water.molecular_volume_au * water.number_density_au == pytest.approx(1.0, rel=1e-12)


def test_alias_and_id_agree():
    water = get_solvent_descriptors("water")
    assert get_solvent_descriptors("h2o") == water
    assert get_solvent_descriptors_by_id(water.id) == water
    assert get_solvent_descriptors("DMSO").name == "dimethylsulfoxide"


def test_missing_geometry_gives_nan_derived_values():
    decalin = get_solvent_descriptors("decalin")
    assert decalin.id == 94
    assert math.isfinite(decalin.epsilon) and math.isfinite(decalin.mass_density_au)
    assert math.isnan(decalin.molar_mass_si) and math.isnan(decalin.number_density_au)


def test_descriptors_are_frozen():
    water = get_solvent_descriptors("water")
    with pytest.raises(dataclasses.FrozenInstanceError):
        water.epsilon = 80.0


@pytest.mark.parametrize("lookup, key, message", [
    (get_solvent_descriptors, "FOO", "Unknown solvent: FOO"),
    (get_solvent_descriptors, "   ", "Empty solvent alias"),
    (get_solvent_descriptors_by_id, 0, "Unknown solvent ID: 0"),
    (get_solvent_descriptors_by_id, 9999, "Unknown solvent ID: 9999"),
])
def test_unknown_solvent_raises(lookup, key, message):
    with pytest.raises(RuntimeError, match=message):
        lookup(key)


@pytest.mark.parametrize("lookup, key, error", [
    (get_solvent_descriptors, 175, TypeError),
    (get_solvent_descriptors, "wa\0ter", ValueError),
    (get_solvent_descriptors_by_id, "175", TypeError),
    (get_solvent_descriptors_by_id, True, TypeError),
])
def test_invalid_key_rejected(lookup, key, error):
    with pytest.raises(error):
        lookup(key)


def test_text_and_print():
    water = get_solvent_descriptors("water")
    text = water.text()
    assert text.startswith("Solvent properties:\n")
    assert "water" in text and "78.355" in text
    out = io.StringIO()
    water.print(file=out)
    assert out.getvalue() == text
    # Formats the values passed, without a table lookup
    assert "81.5" in dataclasses.replace(water, epsilon=81.5).text()


def test_exported_from_package():
    assert moist.get_solvent_descriptors is get_solvent_descriptors
    assert "SolventDescriptors" in moist.__all__
