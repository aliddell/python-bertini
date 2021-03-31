from os import path as op

from bertini.io.input_file.input_section import BertiniInput
from bertini.io.input_file.config_section import BertiniConfig


class BertiniResult:
    def __init__(self, dirname, **kwargs):
        """The output of a Bertini run.

        Parameters
        ----------
        dirname : str
            Path to directory where run was done.
        kwargs
        """
        self._dirname = dirname

    @property
    def config(self):
        """Configuration needed to reproduce the run."""
        return self._config

    @config.setter
    def config(self, val):
        assert isinstance(val, BertiniConfig)
        self._config = val

    @property
    def dirname(self):
        return self._dirname
    @dirname.setter
    def dirname(self, val):
        assert op.isdir(val)
        self._dirname = val

    @property
    def inputs(self):
        """Inputs needed to reproduce the run."""
        return self._inputs

    @inputs.setter
    def inputs(self, val):
        assert isinstance(val, BertiniInput)
        self._inputs = val