
"""atom - Defines classes and function that handle atom and atom groups."""


import numpy


class BaseAtom(object):
    """Base class for an Atom."""

    def __init__(self, index=0, name='XXX', resname='XXX', chain='X',
                 resid=0, charge=0.0, coords=(0, 0, 0)):
        self.name = name
        self.resname = resname
        self.chain = chain
        self.index = index
        self.resid = resid
        self.charge = charge
        self.__set_coords(coords)

    @property
    def coords(self):
        """Get atom cartesian coordinates."""
        return self._coords

    @coords.setter
    def coords(self, pos):
        """Set atom cartesian coordinates."""
        self.__set_coords(pos)

    def __set_coords(self, pos):
        if not isinstance(pos, numpy.ndarray):
            pos = numpy.array(pos, dtype=float)
        if not pos.shape == (3, ):
            err = 'atom coordinates but by a scalar or vector of shape 1 x 3 '\
                  '(found {})'.format(pos.shape)
            raise ValueError(err)
        self._coords = pos


