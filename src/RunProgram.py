'''
Created on 03 May 2018

@author: ernesto
'''

import subprocess
import re
import os


class RunProgram(object):
    """
    SuperClass used to run an external program within a Python script
    """

    def __init__(self, program=None, path=None, args=None, arg_sep=None, parameters=None,
                 cmd_line=None, downpipe=None, log_name=None, log_file=None):
        """
        Constructor

        Parameters
        ----------
        program : str, default
                  Program to be run.
        path : str, optional
               Folder containing the 'program'.
        args : list, optional
               List of named tuples tuples formed by an argument (or option) and
               its respective value: i.e. arg=namedtuple('Argument', 'option value').
        parameters : list, optional
                     List of parameters ['c','d'].
        arg_sep : char, default=' '
                  char used as a separator between argument and value.
                  i.e. if '=' then we will get 'a'=1
                  Default is a single whitespace (i.e. ' ').
        cmd_line : str, optional
                   String with command line to run.
        downpipe: list, optional
                  List of RunProgram objects that will be executed in a pipe after
                  self.program has been executed.
        log_name: str, optional
                  Name of the logger.
        log_file: filename, optional
                  Path to the file that will be used by the logging library.
        """

        self.cmd_line = cmd_line
        self.program = program
        if self.program is None and self.cmd_line is None:
            raise ValueError("Parameter 'cmd_line' or 'program' must be provided.")

        self.path = path
        self.args = args
        self.arg_sep = arg_sep if arg_sep is not None else ' '
        self.parameters = parameters
        self.downpipe = downpipe
        self.log_name = log_name
        self.log_file = log_file

        # create the command line if is None
        if self.cmd_line is None:
            self.cmd_line = self.create_command_line()

    def create_command_line(self):
        """

        :return:
        """
        if self.path is not None:
            cmd_line = [os.path.join(self.path, self.program)]
        else:
            cmd_line = [self.program]

        # construct the command line
        if self.args is not None:
            cmd_line.append(' '.join([f"{option}{self.arg_sep}{parameter}"
                                      for option, parameter in self.args]))

        if self.parameters is not None:
            cmd_line.append(' '.join(["{0}".format(param) for param in self.parameters]))

        if self.downpipe is not None:
            cmd_line.append(' '.join(["| {0}".format(runO.cmd_line) for runO in self.downpipe]))

        return ' '.join(cmd_line)


    def run_popen(self, raise_exc=True):
        """
        Run self.program using subprocess Popen method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Parameters
        ----------
        raise_exc: bool, default=True
                   If true, then raise an Exception when error is found.

        Returns
        -------
        tuple
             A tuple containing the STDOUT and STDERR from this program
        """

        log_f = open(self.log_file, 'w') if self.log_file is not None else None

        # execute cmd_line
        p = subprocess.Popen(self.cmd_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, bufsize=256,
                             universal_newlines=True)

        # stderr
        patt = re.compile('#* ERROR|Error')
        is_exception = False
        stderr = ""
        for line in p.stderr:
            line = str(line.rstrip())
            stderr += line + "\n"
            if log_f is not None:
                log_f.write(line + "\n")
            m = patt.match(line)
            if m:
                is_exception = True

        # stdout
        stdout = ""
        for line in p.stdout:
            line = str(line.rstrip())
            stdout += line + "\n"

        if is_exception is True and raise_exc is True:
            raise Exception(stderr)

        if log_f is not None:
            log_f.close()

        return (stdout, stderr, is_exception)

    def run_checkoutput(self):
        """
        Run self.program using subprocess check_output method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Returns
        -------
        Stdout produced after running the program
        """

        try:
            stdout = subprocess.check_output(self.cmd_line, shell=True)
        except subprocess.CalledProcessError as e:
            raise

        return stdout
