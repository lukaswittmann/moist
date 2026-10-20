"""Bulk descriptors of the native solvent table, by name or table ID."""

from __future__ import annotations

from dataclasses import dataclass, fields

from . import library


@dataclass(frozen=True)
class SolventDescriptors:
    """Bulk properties of one solvent table entry.

    Field names and units follow ``moist_solvent_descriptors``: ``*_si`` in
    SI, ``*_au`` in atomic units. Table values follow the Minnesota Solvent
    Descriptor Database (epsilon and surface tension at 298.15 K, refractive
    index at 293.15 K). A NaN marks a value the table lacks and carries into
    the values derived from it; ``0.0`` is a real value. ``id`` is frozen:
    new solvents are appended.
    """

    id: int
    name: str
    epsilon: float
    refractive_index: float
    abraham_alpha: float
    abraham_beta: float
    surface_tension_si: float
    surface_tension_au: float
    mass_density_si: float
    mass_density_au: float
    molar_mass_si: float
    molecular_mass_au: float
    number_density_si: float
    number_density_au: float
    molecular_volume_si: float
    molecular_volume_au: float

    @classmethod
    def _from_native(cls, native) -> SolventDescriptors:
        values = {item.name: getattr(native, item.name) for item in fields(cls)}
        values["name"] = library.ffi.string(native.name).decode()
        return cls(**values)

    def _as_native(self):
        native = library.ffi.new("moist_solvent_descriptors *")
        native.struct_size = library.ffi.sizeof(native[0])
        for item in fields(self):
            value = getattr(self, item.name)
            setattr(native, item.name, value.encode() if item.name == "name" else value)
        return native

    def text(self) -> str:
        """The "Solvent properties" printout of these values."""
        return library.get_solvent_descriptors_text(self._as_native())

    def print(self, file=None) -> None:
        """Print :meth:`text` to ``file``, standard output by default."""
        print(self.text(), end="", file=file)


def get_solvent_descriptors(name: str) -> SolventDescriptors:
    """Descriptors of a solvent by name or any alias, e.g. ``"H2O"``.

    Case-insensitive, surrounding whitespace ignored; an unknown name raises.
    """
    return SolventDescriptors._from_native(library.get_solvent_descriptors(name))


def get_solvent_descriptors_by_id(solvent_id: int) -> SolventDescriptors:
    """Descriptors of a solvent by its table ID; an unknown ID raises."""
    return SolventDescriptors._from_native(library.get_solvent_descriptors_by_id(solvent_id))
