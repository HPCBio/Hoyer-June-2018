#!/usr/bin/env python3

"""
Run an interative Pilon polishing.  Requires Pilon, bowtie2, samtools

Params:

    trimmed reads
    nanopolished assembly

Note: this is a really whittled down version of Ryan Wick's Pilon script, with
internal copies of his many of his functions which rely on GFA input.  Long-term
it would make sense to somehow add this functionality directly to Unicycler

"""

import argparse
import os
import sys
import subprocess
from collections import defaultdict
from Bio import SeqIO

def main():
    args = get_arguments()
    polish_dir = os.getcwd()

    for i in range(10):
        round_num = i + 1
        input_filename = str(round_num) + '_polish_input.fasta'
        if round_num == 1 and not os.path.islink(input_filename):
            # TODO: needs try/except
            os.symlink(args.input, input_filename)

        if not os.path.exists(input_filename):
            raise CannotPolish("File doesn't exist")

        paired_bam_filename = str(round_num) + '_paired_alignments.bam'
        unpaired_bam_filename = str(round_num) + '_unpaired_alignments.bam'
        change_count = polish_fasta_with_pilon(input_filename,
                                               args,
                                               polish_dir,
                                               round_num,
                                               'snps,indels,gaps')

        # cleanup
        for f in [input_filename + '.1.bt2', input_filename + '.2.bt2',
                  input_filename + '.3.bt2', input_filename + '.4.bt2',
                  input_filename + '.rev.1.bt2', input_filename + '.rev.2.bt2',
                  paired_bam_filename, paired_bam_filename + '.bai',
                  unpaired_bam_filename, unpaired_bam_filename + '.bai',
                  input_filename
                  ]:
            try:
                os.remove(os.path.join(polish_dir, f))
            except FileNotFoundError:
                pass

        if not change_count:
            break


def get_arguments():
    parser = argparse.ArgumentParser(description='Pilon polishing tool for long read assemblies')
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA to be polished')
    parser.add_argument('-o', '--output', required=True,
                        help='Output prefix for FASTA files')

    parser.add_argument('-1', '--short1', required=False,
                        help='FASTQ file of first short reads in each pair')
    parser.add_argument('-2', '--short2', required=False,
                        help='FASTQ file of second short reads in each pair')
    parser.add_argument('-s', '--unpaired', required=False,
                        help='FASTQ file of unpaired short reads')

    parser.add_argument('-a', '--minins', required=False, default=200,
                        help='Avg insert size (go into Bowtie2)')

    parser.add_argument('-d', '--maxins', required=False, default=800,
                        help='Avg insert size stdev (goes into Bowtie2)')

    parser.add_argument('-t', '--threads', type=int, required=False, default=24,
                        help='Number of threads used')

    parser.add_argument('--diploid', action='store_true',
                        help='Is this a diploid assembly')

    parser.add_argument('--min_polish_size', type=int, default=10000,
                        help='Contigs shorter than this value (bp) will not be polished '
                             'using Pilon')

    parser.add_argument('--bowtie2_path', type=str, default='bowtie2',
                        help='Path to the bowtie2 executable')
    parser.add_argument('--bowtie2_build_path', type=str, default='bowtie2-build',
                        help='Path to the bowtie2_build executable')
    parser.add_argument('--pilon_path', type=str, default='pilon',
                        help='Path to a Pilon executable or the Pilon Java archive file')
    parser.add_argument('--java_path', type=str, default='java',
                        help='Path to the java executable')
    parser.add_argument('--samtools_path', type=str, default='samtools',
                        help='Path to the samtools executable')


    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if (args.short1 and not args.short2) or (args.short2 and not args.short1):
        quit_with_error('you must use both --short1 and --short2 or neither')

    if not args.short1 and not args.short2 and not args.unpaired:
        quit_with_error('no input reads provided (--short1, --short2, --unpaired)')

    # Change some arguments to full paths.
    if args.short1:
        args.short1 = os.path.abspath(args.short1)
    if args.short2:
        args.short2 = os.path.abspath(args.short2)
    if args.unpaired:
        args.unpaired = os.path.abspath(args.unpaired)

    pilon_path, _, _ = pilon_path_and_version(args.pilon_path, args.java_path, args)

    args.verbosity = 2
    args.keep = 3

    return args

def polish_fasta_with_pilon(fasta, args, polish_dir, round_num, fix_type):
    """
    Runs Pilon on the graph to hopefully fix up small mistakes.
    """
    #log.log(underline('Pilon polish round ' + str(round_num)))
    print('Pilon polish round ' + str(round_num), file=sys.stderr, flush=True)

    using_paired_reads = bool(args.short1) and bool(args.short2)
    using_unpaired_reads = bool(args.unpaired)
    assert using_paired_reads or using_unpaired_reads

    input_filename = str(round_num) + '_polish_input.fasta'
    output_prefix = str(round_num) + '_pilon'
    pilon_fasta_filename = str(round_num) + '_pilon.fasta'
    pilon_changes_filename = str(round_num) + '_pilon.changes'
    pilon_output_filename = str(round_num) + '_pilon.out'

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [args.bowtie2_build_path, fasta, fasta]

    print(' '.join(bowtie2_build_command), file=sys.stderr, flush=True)

    #log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [args.bowtie2_path, '--local', '--very-sensitive-local',
                       '--threads', str(args.threads), '-I', str(args.minins),
                       '-X', str(args.maxins), '-x', input_filename]
    if using_paired_reads:
        paired_sam_filename = str(round_num) + '_paired_alignments.sam'
        paired_bam_filename = str(round_num) + '_paired_alignments.bam'
        this_bowtie2_command = bowtie2_command + ['-S', paired_sam_filename,
                                                  '-1', args.short1, '-2', args.short2]
        run_bowtie_samtools_commands(args, this_bowtie2_command, paired_sam_filename,
                                     paired_bam_filename)
    else:
        paired_bam_filename = ''

    if using_unpaired_reads:
        unpaired_sam_filename = str(round_num) + '_unpaired_alignments.sam'
        unpaired_bam_filename = str(round_num) + '_unpaired_alignments.bam'
        this_bowtie2_command = bowtie2_command + ['-S', unpaired_sam_filename, '-U', args.unpaired]
        run_bowtie_samtools_commands(args, this_bowtie2_command, unpaired_sam_filename,
                                     unpaired_bam_filename)
    else:
        unpaired_bam_filename = ''

    # Polish with Pilon.
    if args.pilon_path.endswith('.jar'):
        pilon_command = [args.java_path, '-jar', args.pilon_path]
    else:
        pilon_command = [args.pilon_path]
    pilon_command += ['--genome', input_filename, '--changes',
                      '--output', output_prefix, '--outdir', polish_dir, '--fix', fix_type]

    if args.diploid:
        pilon_command += ['--diploid']

    if using_paired_reads:
        pilon_command += ['--frags', paired_bam_filename]

    if using_unpaired_reads:
        pilon_command += ['--unpaired', unpaired_bam_filename]

    print(' '.join(pilon_command), file=sys.stderr, flush=True)
    # log.log(dim('  ' + ' '.join(pilon_command)), 2)
    try:
        pilon_stdout = subprocess.check_output(pilon_command, stderr=subprocess.STDOUT)
        with open(pilon_output_filename, 'wb') as pilon_out:
            pilon_out.write(pilon_stdout)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Pilon encountered an error:\n' + e.output.decode())
    if not os.path.isfile(pilon_fasta_filename):
        raise CannotPolish('Pilon did not produce FASTA file')
    if not os.path.isfile(pilon_changes_filename):
        raise CannotPolish('Pilon did not produce changes file')

    # # Display Pilon changes.
    change_count = defaultdict(int)
    change_lines = defaultdict(list)
    total_count = 0
    pilon_changes = open(pilon_changes_filename, 'rt')
    for line in pilon_changes:
        try:
            seg_name = line.split(':')[0]
            change_count[seg_name] += 1
            total_count += 1
            change_lines[seg_name].append(line.strip())
        except ValueError:
            pass

    # Unicycler polishes the sequences and then stores them internally in  the
    # graph (e.g. mutates the graph each round).  We don't do this; we simply
    # rename the pilon file appropriately
    if total_count == 0:
        print('No Pilon changes, final file generated', file=sys.stderr, flush=True)
        clean_pilon(fasta, 'final.pilon.fasta')
        # return prior input file as final
    elif round_num == 10:
        print('Final round: final file generated', file=sys.stderr, flush=True)
        # return current pilon file as final
        clean_pilon(pilon_fasta_filename, 'final.pilon.fasta', file=sys.stderr, flush=True)
    else:
        print('Total number of changes: ' + int_to_str(total_count), file=sys.stderr, flush=True)
        next_pilon_input = str(round_num + 1) + '_polish_input.fasta'

        clean_pilon(pilon_fasta_filename, next_pilon_input, file=sys.stderr, flush=True)

    if args.keep < 3:
        list_of_files = [input_filename, pilon_fasta_filename, pilon_changes_filename,
                         input_filename + '.1.bt2', input_filename + '.2.bt2',
                         input_filename + '.3.bt2', input_filename + '.4.bt2',
                         input_filename + '.rev.1.bt2', input_filename + '.rev.2.bt2',
                         pilon_output_filename, paired_bam_filename, paired_bam_filename + '.bai',
                         unpaired_bam_filename, unpaired_bam_filename + '.bai']
        for f in list_of_files:
            try:
                os.remove(f)
                print("Removed " + f, file=sys.stderr, flush=True)
            except FileNotFoundError:
                pass

    return total_count

def run_bowtie_samtools_commands(args, bowtie2_command, sam_filename, bam_filename):
    """
    Simple function that wraps the bowtie2 aln + samtools cleanup
    """
    # Run bowtie2 alignment command
    print(' '.join(bowtie2_command), file=sys.stderr, flush=True)
    #log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    # Sort the alignments.
    samtools_sort_command = [args.samtools_path, 'sort', '-@', str(args.threads),
                             '-o', bam_filename, '-O', 'bam', '-T', 'temp', sam_filename]

    print(' '.join(samtools_sort_command), file=sys.stderr, flush=True)
    #log.log(dim('  ' + ' '.join(samtools_sort_command)), 2)
    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # # Always delete the SAM file, regardless of keep level (it's big).
    try:
        os.remove(sam_filename)
    except FileNotFoundError:
        pass

    # Index the alignments.
    samtools_index_command = [args.samtools_path, 'index', bam_filename]
    print(' '.join(samtools_index_command), file=sys.stderr, flush=True)
    #log.log(dim('  ' + ' '.join(samtools_index_command)), 2)

    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

def int_to_str(num, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))

def quit_with_error(message):
    """
    Displays the given message and ends the program's execution.
    """
    #log.log('Error: ' + message, 0, stderr=True)
    print('Error: ' + message, file=sys.stderr, flush=True)
    sys.exit(1)

def clean_pilon(file, name):
    with open(file, "rU") as input_handle, open(name, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if record.id.endswith('_pilon'):
                record.id = record.id[:-6]
            SeqIO.write(record, output_handle, "fasta")

def pilon_path_and_version(pilon_path, java_path, args):
    status = find_pilon(pilon_path, java_path, args)
    if status == 'good':
        pilon_path = args.pilon_path
    else:
        return pilon_path, '', status
    if pilon_path.endswith('.jar'):
        command = [java_path, '-jar', pilon_path, '--version']
    else:
        command = [pilon_path, '--version']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _ = process.communicate()
    version = out.decode().split('Pilon version ')[-1].split()[0]
    try:
        int(version.split('.')[0]), int(version.split('.')[1])
    except (ValueError, IndexError):
        version, status = '?', 'too old'
    return os.path.abspath(pilon_path), version, 'good'

def find_pilon(pilon_path, java_path, args):
    """
    Makes sure the Pilon executable is available. Unlike the other tools, Pilon's target name may
    change based on the version (e.g. pilon-1.18.jar, pilon-1.19.jar, etc.). This function will
    therefore set args.pilon_path to the first matching jar file it finds.
    """
    # If the user specified a Pilon path other than the default, then it must exist.
    if args.pilon_path != 'pilon' and args.pilon_path is not None:
        args.pilon_path = os.path.abspath(args.pilon_path)
        if args.pilon_path.endswith('.jar'):
            if not os.path.isfile(args.pilon_path):
                return 'not found'
        elif shutil.which(args.pilon_path) is None:
            return 'not found'

    # If pilon_path is the default and exists, then that's great!
    elif args.pilon_path == 'pilon' and shutil.which(args.pilon_path) is not None:
        args.pilon_path = shutil.which(args.pilon_path)

    # If the user didn't specify a path and 'pilon' doesn't work, then we need to look for a
    # Pilon jar file.
    else:
        found_pilon_path = get_pilon_jar_path(pilon_path)
        if found_pilon_path:
            args.pilon_path = found_pilon_path
        else:
            return 'not found'

    # Now that we've found Pilon, run the help command to make sure it works.
    if args.pilon_path.endswith('.jar'):
        test_command = [java_path, '-jar', args.pilon_path, '--help']
    else:
        test_command = [args.pilon_path, '--help']
    try:
        pilon_help_out = subprocess.check_output(test_command, stderr=subprocess.STDOUT).decode()
        if 'pilon' not in pilon_help_out.lower():
            raise OSError
    except (FileNotFoundError, OSError, subprocess.CalledProcessError):
        return 'bad'

    return 'good'

def get_pilon_jar_path(pilon_path):
    """
    Returns the path to pilon.jar. If the given path is correct, it just returns that, as an
    absolute path. Otherwise it tries to find it.
    """
    if pilon_path and os.path.isfile(pilon_path):
        return os.path.abspath(pilon_path)
    for directory in os.environ['PATH'].split(':'):
        try:
            path_files = [f for f in os.listdir(directory)
                          if os.path.isfile(os.path.join(directory, f))]
        except FileNotFoundError:
            path_files = []
        pilon_jars = [f for f in path_files if f.startswith('pilon') and f.endswith('.jar')]
        if pilon_jars:
            return os.path.join(directory, sorted(pilon_jars)[-1])  # return the latest version
    return None

class CannotPolish(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

if __name__ == '__main__':
    main()
