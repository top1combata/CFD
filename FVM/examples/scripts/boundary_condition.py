from dataclasses import dataclass
from enum import Enum
from typing import TypeVar, Generic


class BoundaryType(Enum):
    fixedValue = 1
    fixedGradient = 2


T = TypeVar("T")
@dataclass
class BoundaryCondition(Generic[T]):
    type: BoundaryType
    value: T


@dataclass
class Boundaries:
    p_boundary: BoundaryCondition[float]
    u_boundary: BoundaryCondition[tuple[float, float]]


def wall_boundaries():
    return Boundaries(
        u_boundary=BoundaryCondition[tuple[float, float]](
            type=BoundaryType.fixedValue,
            value=(0, 0)
        ),
        p_boundary=BoundaryCondition[float](
            type=BoundaryType.fixedGradient,
            value=0
        )
    )


def outlet_boundaries(pressure: float):
    return Boundaries(
        u_boundary=BoundaryCondition[tuple[float, float]](
            type=BoundaryType.fixedGradient,
            value=(0, 0)
        ),
        p_boundary=BoundaryCondition[float](
            type=BoundaryType.fixedValue,
            value=pressure
        )
    )


def inlet_boundaries(velocity: tuple[float, float]):
    return Boundaries(
        u_boundary=BoundaryCondition[tuple[float, float]](
            type=BoundaryType.fixedValue,
            value=velocity
        ),
        p_boundary=BoundaryCondition[float](
            type=BoundaryType.fixedGradient,
            value=0
        )
    )
