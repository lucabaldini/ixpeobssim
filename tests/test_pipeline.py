#!/usr/bin/env python
#
# Copyright (C) 2016--2020, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from __future__ import print_function, division


import unittest

import ixpeobssim.core.pipeline as pipeline


class TestPipeline(unittest.TestCase):

    """Unit test for the pipeline module.
    """

    def test_setup(self):
        """
        """
        print(pipeline.params())
        test = 'test'
        pipeline.setup(test_field=test)
        self.assertEqual(pipeline.param('test_field'), test)
        print(pipeline.params())
        model = 'my_model'
        pipeline.set_model(model)
        self.assertEqual(pipeline.model(), model)
        print(pipeline.params())

    def test_file_lists(self):
        """
        """
        pipeline.set_model('mymodel')
        for args in ((), ('mcube',), (('sel', 1),), (('sel', 1), 'cmap')):
            print(args, pipeline.file_list(*args, check_files=False))

    def test_pha_start(self):
        pipeline.set_model('mymodel')
        for args in (('pha1*',), ('mcube', 'pha1*')):
            print(args, pipeline.file_list(*args, check_files=False))

    def test_command_line_switches(self):
        """
        """
        args = ['file1', 'file2']
        kwargs = dict(configfile='test.py', duration=100, overwrite=False,
                      ebinning=[2., 4., 8.])
        print(args)
        print(kwargs)
        switches = pipeline._command_line_switches(*args, **kwargs)
        print(switches)

    def test_xpobssim_kwargs(self):
        """
        """
        parser = pipeline.XPOBSSIM_PARSER
        kwargs = dict(configfile='test.py', duration=100, overwrite=False)
        print(kwargs)
        pipeline.setup(overwrite=False)
        switches = pipeline._command_line_switches(**kwargs)
        print(switches)
        kwargs = parser.parse_args(switches).__dict__
        print(kwargs)
        pipeline.setup(overwrite=True)
        switches = pipeline._command_line_switches(**kwargs)
        print(switches)
        kwargs = parser.parse_args(switches).__dict__
        print(kwargs)



if __name__ == '__main__':
    unittest.main()
