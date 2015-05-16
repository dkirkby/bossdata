# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

""" Define bit masks used in BOSS data and support symbolic operations on masks.
"""

from __future__ import division,print_function

def define_bitmask(**bits):
    """Define a new type for a bitmask with specified symbolic bit names.

    After defining a bitmask type, its bit names are accessible as class-level
    attributes of the returned type and can be used as integer values, for example::

        >>> Colors = define_bitmask(RED=0,BLUE=1,GREEN=4)
        >>> Colors.BLUE
        2
        >>> '{0:b}'.format(Colors.RED|Colors.GREEN)
        '10001'

    The :func:`decode_bitmask` function is useful for converting an integral value back
    to a list of symbolic bit names.

    Args:
        bits(dict): A dictionary of name,value pairs that define the mapping from
            symbolic bit names to bit offsets, which should be integers >= 0. Although
            this can be passed as a dictionary, the dictionary is usually implicitly
            defined by the argument list, as in the examples above. By convention,
            bit names are all upper case.

    Returns:
        type: A new type with the specified name that has class-level attributes for
            each named bit.

    Raises:
        SyntaxError: bit name is repeated in argument list.
        TypeError: bit name is not mapped to a valid type.
        ValueError: bit name is not mapped to a value value.
    """
    # Check that bit offset values are unique.
    if len(bits) != len(set(bits.values())):
        raise ValueError('Bit offset values must be unique.')
    # Initialize a class dictionary with attributes for each named bit.
    class_dict = dict(zip(bits.keys(),(1<<offset for offset in bits.values())))
    # Add a reverse-lookup mapping from offsets to names.
    class_dict['_reverse_map'] = dict((offset,name) for name,offset in bits.iteritems())
    # Return a new class type for this bitmask.
    return type('BitMask',(),class_dict)

def decode_bitmask(mask,value,strict=True):
    """Decode a integer value into its symbolic bit names.

    Use this function to convert a bitmask value into a list of symbolic bit names, for
    example::

        >>> Colors = define_bitmask(RED=0,BLUE=1,GREEN=4)
        >>> decode_bitmask(Colors,Colors.RED|Colors.BLUE)
        ('RED', 'BLUE')

    For pretty printing, try::

        >>> print('|'.join(decode_bitmask(Colors,Colors.RED|Colors.BLUE)))
        RED|BLUE

    Args:
        mask: A bitmask type, normally created with :func:`create_bitmask`, that defines
            the symbolic bit names to use for the decoding.
        value(int): The integral value to decode.
        strict(bool): If set, then all bits set in value must be defined in the bitmask
            type definition.

    Returns:
        tuple: A tuple of symbolic bit names in order of increasing bit offset. If strict
            is False, then any bits without corresponding symbolic names will appear as
            '1<<n' for offset n.

    Raises:
        AttributeError: mask does not have the attributes necessary to define a bitmask.
        ValueError: value has a bit set that has no symbolic name defined and strict is True.
    """
    names = [ ]
    shifted = int(value)
    offset = 0
    while shifted:
        if shifted & 1:
            name = mask._reverse_map.get(offset)
            if name is None:
                if strict:
                    raise ValueError('Bit {0:d} is set but has no symbolic name defined.'.format(offset))
                name = '1<<{0:d}'.format(offset)
            names.append(name)
        offset += 1
        shifted = shifted >> 1
    return tuple(names)
