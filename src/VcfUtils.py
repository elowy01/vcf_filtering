'''
Created on 03 Aug 2017

@author: ernesto
'''

import os
import subprocess
import gzip
import re
from collections import namedtuple
from Utils.RunProgram import RunProgram

class VcfUtils(object):
    """
    Class to represent a misc of actions that can be done on a single or multiple VCF files
    """
    def __init__(self, vcf=None, vcflist=None, bcftools_folder=None, bgzip_folder=None,
                 gatk_folder=None, java_folder=None, tmp_dir=None):
        """
        Constructor

        Parameters
        ----------
        vcf : str, optional
              Path to gzipped vcf file.
        vcflist : list, optional
                  List of dicts containing setname:vcf_paths (keys:values) pairs.
        bcftools_folder : str, optional
                          Path to folder containing the bcftools binary.
        bgzip_folder : str, optional
                       Path to folder containing the bgzip binary.
        gatk_folder : str, optional
                      Path to folder containing the jar file.
        java_folder : str, optional
                      Path to folder containing the java binary.
        tmp_dir : str, optional
                  Path to java temporary directory. This needs to be
                  set for GATK modules that
                  fail because there is not enough space in the default java tmp dir.

        Imp: Either 'vcf' or 'vcflist' variables should be initialized
        """

        if not vcf and not vcflist:
            raise Exception("Either a vcf file or a list of vcf files should be used\
                             to initialize this class")

        if vcf is not None:
            if os.path.isfile(vcf) is False:
                raise Exception("File does not exist")

        self.vcf = vcf
        self.vcflist = vcflist
        self.bcftools_folder = bcftools_folder
        self.bgzip_folder = bgzip_folder
        self.gatk_folder = gatk_folder
        self.java_folder = java_folder
        self.tmp_dir = tmp_dir

    def reheader(self, newheader, outprefix, samplefile=None, verbose=False):
        """
        Modifiy the VCF's header with the newheader

        Parameters
        ----------
        newheader : str
                    Path to the file containing the new header.
        outprefix : string
                    Prefix for output files
        samplefile : str, optional
                     Path to the file with the sample names that will included
                     in the new header.
        verbose : bool, default=False
                  increase the verbosity.

        Returns
        -------
        outfile : str
                 Path to the VCF with the modified header.
        """
        outfile = outprefix+".reheaded.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-h', newheader), Arg('-o', outfile)]

        if samplefile is not None:
            args.append(Arg('-s', samplefile))

        runner = RunProgram(path=self.bcftools_folder, program='bcftools reheader',
                            args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile

    def add_to_header(self, header_f, outfilename, line_ann):
        """
        Function to add to the header of a VCF the string passed with 'line_ann'

        Parameters
        ----------
        header_f : str
                   Path to file containing the header file that will be modified.
        outfilename : str
                      Path to the new header file that is modified.
        line_ann : str
                   Str with line that will be used to add to the header.

        Returns
        -------
        outfilename : str
                Path to modified header file. The new annotation will be added in the following line
                after the desired annotation.
        """

        of = open(outfilename, 'w')

        # getting the type of line that is being passed
        p = re.compile("^##(\w+)=")

        m1 = p.match(line_ann)
        type_ann = m1.group(1)

        line_seen = False
        with open(header_f) as f:
            for line in f:
                line = line.rstrip("\n")
                m2 = p.match(line)
                if m2 is None:
                    of.write(line+"\n")
                    continue
                type1_ann = m2.group(1)
                if type_ann == type1_ann and line_seen is False:
                    line_seen = True
                    of.write(line+"\n"+line_ann+"\n")
                    continue
                else:
                    of.write(line+"\n")
        of.close()

        return outfilename

    def combine(self, labels, reference, outprefix, compress=False, outdir=None,
                ginterval=None, genotypemergeoption=None, filteredrecordsmergetype=None,
                threads=1, options=None, verbose=False):
        """
        Combine VCFs using GATK's CombineVariants into a single VCF

        Parameters
        ----------
        labels : list
                 List of labels used for each of the VCFs in self.vcflist. The order of the labels
                 should be the same that the VCFs in the list.
        reference : str
                    Path to Fasta file with reference.
        outprefix : str
                    Prefix used for output file.
        compress : bool, default=False
                   Compress the output VCF with bgzip.
        outdir : str, optional
                 Path to folder used to write the results to.
        ginterval : str, optional
                    Genomic interval used to restrict the analysis. i.e. chr20:1000-2000.
        genotypemergeoption : {'UNIQUIFY', 'PRIORITIZE', 'UNSORTED', 'REQUIRE_UNIQUE'}, optional
                              Determines how we should merge genotype records for samples shared
                              across the ROD files.
        filteredrecordsmergetype : {'KEEP_IF_ANY_UNFILTERED', 'KEEP_IF_ANY_UNFILTERED',
        'KEEP_UNCONDITIONAL'}, optional
                                   Determines how we should handle records seen at the
                                   same site in the VCF, but with different FILTER fields.
        threads : int, default=1
                  Number of trades to use.
        options : list, optional
                  List of options. i.e. ['-env','--filteredAreUncalled'].
        verbose : bool, default=False
                  increase the verbosity.

        Returns
        -------
        outfile : str
                Path to the merged VCF.
        """

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-T', 'CombineVariants'), Arg('-R', reference), Arg('-nt', threads)]

        variants_str = ""
        for path, label in zip(self.vcflist, labels):
            if os.path.isfile(path) == False:
                print("Error reading from {0}".format(path))
                raise Exception("File does not exist")
            args.append(Arg('-V:{0}'.format(label), path))

        outfile = ""
        if outdir:
            outfile = "{0}/".format(outdir)
        outfile += "{0}.vcf".format(outprefix)

        if ginterval is not None:
            args.append(Arg('-L', ginterval))

        if genotypemergeoption is not None:
            args.append(Arg('--genotypemergeoption', genotypemergeoption))

        if filteredrecordsmergetype is not None:
            args.append(Arg('--filteredrecordsmergetype', filteredrecordsmergetype))

        params = []
        if options:
            for opt in options:
                params.append(opt)

        pipelist = None
        if compress is True:
            outfile += ".gz"
            compressRunner = RunProgram(path=self.bgzip_folder,
                                        program='bgzip',
                                        parameters=['-c', '>', outfile])
            pipelist = [compressRunner]
        else:
            args.append(Arg('-o', outfile))

        program_str = None
        if self.tmp_dir is not None:
            program_str = "java -Djava.io.tmpdir={0} " \
                          "-jar {1}/GenomeAnalysisTK.jar".format(self.tmp_dir, self.gatk_folder)
        else:
            program_str = "java -jar {0}/GenomeAnalysisTK.jar".format(self.gatk_folder)

        runner = RunProgram(path=self.java_folder, program=program_str,
                            args=args, parameters=params, downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_exc = runner.run_popen()

        return outfile

    def rename_chros(self, chr_types, outfile, compress=True):
        """
        Function to modify the chr names in the VCF file
        For example:
        If file has UCSC-type chr names (i.e. chr1,chr2,...) then this
        function will convert the UCSC-type chr names to Ensembl-type
        chr names (i.e. 1,2,...) or vice-versa

        Parameters
        ----------
        chr_types : {'ucsc','ensembl'}
                    Type of chr names that will be written to the file.
        outfile : str
                  File used for the output VCF.
        compress : bool, default=True

        Returns
        -------
        outfile : str
                 Path to the VCF with the chrosomes renamed.
        """

        command = ""
        if chr_types == 'ensembl':
            if compress is True:
                command += "zcat {0} | awk '{{gsub(/^chr/,\"\"); print}}' - | {1}/bgzip -c > {2}".\
                    format(self.vcf, self.bgzip_folder, outfile)
            else:
                command += "zcat {0} | awk '{{gsub(/^chr/,\"\"); print}}' - > {1}".\
                    format(self.vcf, outfile)
        elif chr_types == 'ucsc':
            if compress is True:
                command += "zcat {0} | awk '{{if($0 !~ /^#/) print \"chr\"$0; " \
                           "else print $0}}' - | {1}/bgzip -c > {2}".format(self.vcf,
                                                                            self.bgzip_folder,
                                                                            outfile)
            else:
                command += "zcat {0} | awk '{{if($0 !~ /^#/) print \"chr\"$0; " \
                           "else print $0}}' - > {1}".format(self.vcf, outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong.\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

        return outfile

    def correct_ambiguity_codes(self, outfile):
        """
        Function to correct the ambiguity bases in the VCF. This ambiguity
        may appear in the REF or ALT columns

        Parameters
        ----------
        outfile : str
                  File where the output VCF will be written.

        Returns
        -------
        outfile : str
                 Path to vcf.gz file compressed with GZIP.
        """
        ref_count = 0
        alt_count = 0

        f = gzip.open(outfile, 'wb')

        with gzip.open(self.vcf, 'r') as fin:
            for line in fin:
                if not line.startswith(b"#"):
                    bits = line.split(b"\t")
                    ref = bits[3].decode("utf-8")
                    alt = bits[4].decode("utf-8")
                    if re.search(r"[^ATGC.,]", ref):
                        ref_count += 1
                        ref = re.sub('[^ACGT.]', 'N', ref)
                    if re.search(r"[^ATGC.,]", alt):
                        alt_count += 1
                        alt = re.sub('[^ACGT.]', 'N', alt)
                    bits[3] = ref.encode('utf-8')
                    bits[4] = alt.encode('utf-8')
                    nline = b'\t'.join(bits)
                    f.write(nline)
                else:
                    f.write(line)
        f.close()

        print("Sites with ambiguous bases in the REF column is:{0}".format(ref_count))
        print("Sites with ambiguous bases in the ALT column is:{0}".format(alt_count))

        return outfile

    def drop_genotypes(self, outfile, verbose=False):
        """
        Function to drop the Genotype information from a VCF.
        This function uses bcftools -G to perform this operation

        Parameters
        ----------
        outfile : str
                  File where the output VCF will be written.
        verbose : bool, default=False
                  increase the verbosity.

        Returns
        -------
        outfile : str
                 Path to the vcf.gz file without the GT information.
        """
        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-o', outfile), Arg('-O', 'z')]

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools view -G', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile

    def drop_info(self, outfile, verbose=False):
        """
        Function to remove the INFO annotation from a VCF.
        This function uses bcftools annotate  to perform this operation

        Parameters
        ----------
        outfile : str
                  File where the output VCF will be written.
        verbose : bool, default=False
                  increase the verbosity.

        Returns
        -------
        outfile : str
                 Path to the vcf.gz file without the INFO annotation.
        """
        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-o', outfile), Arg('-O', 'z')]

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools annotate --remove INFO',
                            args=args,
                            parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile

    def convert_PL2GL(self, outfile, threads=1, verbose=False):
        """
        Function to convert PL fields into GT.
        This function makes use of Bcftools +tag2tag plugin

        Parameters
        ----------
        outfile : str
                  File where the output VCF will be written.
        threads : int, default=1
                  Number of trades to use.
        verbose : bool, default=False
                  increase the verbosity.

        Returns
        -------
        outfile : str
                 Path to the vcf.gz file with the PL fields converted.
        """
        Arg = namedtuple('Argument', 'option value')

        params = [self.vcf, '-Oz', '--', '-r', '--pl-to-gl']

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools +tag2tag',
                            args=[Arg('--threads', threads), Arg('-o', outfile)],
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile
