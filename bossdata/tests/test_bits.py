# Licensed under a MIT style license - see LICENSE.rst

import pytest

from .. import bits

class TestBitmask(object):
    """Test bitmask definition functionality.
    """
    def test_define_basic(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=1,BIT4=4)
        value = mask.BIT0 | mask.BIT4
        assert value == (1<<0)|(1<<4)

    def test_define_docstring(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=1,BIT4=4)
        assert mask.__doc__.startswith("description")

    def test_define_description(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=(1,"BIT1"),BIT4=4)
        assert mask._description[1] == "BIT1"

    def test_define_large_int(self):
        mask = bits.define_bitmask("MASK","description",BIG=80)
        assert mask.BIG == (1<<80)

    def test_define_missing_bits(self):
        # It is not an error to define a mask with no bits defined.
        mask = bits.define_bitmask("MASK","description")

    def test_define_missing_args(self):
        with pytest.raises(TypeError):
            mask = bits.define_bitmask("MASK",BIT0=0,BIT1=1,BIT4=4)
        with pytest.raises(TypeError):
            mask = bits.define_bitmask(BIT0=0,BIT1=1,BIT4=4)

    def test_define_bad_definition(self):
        with pytest.raises(ValueError):
            mask = bits.define_bitmask("MASK","description",BIT0=0,BAD=-1)
        with pytest.raises(ValueError):
            mask = bits.define_bitmask("MASK","description",BIT0=0,BAD='bad')
        with pytest.raises(ValueError):
            mask = bits.define_bitmask("MASK","description",BIT0=0,BAD=(0,'description','invalid'))

    def test_define_repeated_name(self):
        # This raises a SyntaxError at load/import time.
        # mask = bits.define_bitmask("MASK","description",BIT=0,BIT=1)
        # This will silently overwrite the repeated name. Is there any way to flag this?
        mask = bits.define_bitmask("MASK","description",**{ 'BIT':0,'BIT':1 })

    def test_define_repeated_value(self):
        with pytest.raises(ValueError):
            mask = bits.define_bitmask("MASK","description",A=0,B=0)

class TestDecode(object):
    """Test bitmask decoding functionality.
    """
    def test_decode_basic(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=1,BIT4=4)
        assert bits.decode_bitmask(mask,mask.BIT0|mask.BIT4) == ('BIT0','BIT4')

    def test_decode_attribute_error(self):
        with pytest.raises(AttributeError):
            bits.decode_bitmask('abc',123)

    def test_decode_strict(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=1,BIT4=4)
        with pytest.raises(ValueError):
            bits.decode_bitmask(mask,(1<<2))
        with pytest.raises(ValueError):
            bits.decode_bitmask(mask,mask.BIT0|(1<<2))
        with pytest.raises(ValueError):
            bits.decode_bitmask(mask,mask.BIT4|(1<<2))

    def test_decode_non_strict(self):
        mask = bits.define_bitmask("MASK","description",BIT0=0,BIT1=1,BIT4=4)
        assert bits.decode_bitmask(mask,(1<<2),strict=False) == ('1<<2',)
        assert bits.decode_bitmask(mask,mask.BIT0|(1<<2),strict=False) == ('BIT0','1<<2',)
        assert bits.decode_bitmask(mask,mask.BIT4|(1<<2),strict=False) == ('1<<2','BIT4',)

class TestBOSS(object):
    """Test for expected BOSS bit mask definitions.
    """
    def test_boss(self):
        assert bits.ZWARNING.MANY_OUTLIERS == (1<<4)
        assert bits.BOSS_TARGET1.QSO_CORE == (1<<10)
        assert bits.ANCILLARY_TARGET2.QSO_STD == (1<<20)
