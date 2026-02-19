# TODO:
# - [x] replicate index functionality
# - [x] use pyfastx for fasta parsing
# - [x] default to bwa-mem2, always. fallback to bwa. if user specifies override
# - [ ] add rich/tqdm progress indicators at hotspots

import os
import sh
import sys
import shlex
import typer
import pysam
import pyfastx
from xopen import xopen
from pathlib import Path
from operator import attrgetter
from itertools import groupby, chain
from typing import Optional, Annotated, TextIO
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("bwa-meth")
except PackageNotFoundError:
    __version__ = "unknown"

app = typer.Typer(pretty_exceptions_show_locals=True)

BWA_EXTENSIONS     = [".amb", ".ann", ".pac", ".sa", ".bwt"]
BWA_MEM2_EXTENSIONS = [".amb", ".ann", ".pac", ".0123", ".bwt.2bit.64"]

class BWAMethException(Exception):
    """Raise an exception for BWAMeth"""
    pass

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def check_executable(cmd):
    return any([os.access(os.path.join(p, cmd), os.X_OK) for p in os.environ['PATH'].split(":")])

def wrap(text, width=100): # much faster than textwrap
    for s in range(0, len(text), width):
        yield text[s:s+width] + "\n"
 
def bwa_index(orig_fasta: Path, c2t_fasta: Path):
    if check_executable("bwa-mem2"):
        cmd = sh.bwa_mem2.bake("index", "-p") # type: ignore
    elif check_executable("bwa"):
        raise NotImplementedError
    else:
        raise FileNotFoundError("Cannot find bwa or bwa-mem2 executable in PATH")

    index_dir = orig_fasta.with_suffix(".d")
    index_dir.mkdir(exist_ok=True)
    prefix = index_dir / "bwa_meth_index"
    baked_command= cmd.bake(prefix, c2t_fasta)
    eprint("--------------------")
    eprint(f"Running: {str(baked_command)}")

    baked_command(_err=sys.stdout, _out=sys.stderr)

    return 0

def convert_and_write_read(name, seq, qual, read_i, out: TextIO): 
    # Handle non-Illumina formats
    if name.endswith(("_R1", "_R2")):
        name = name[:-3]
    elif name.endswith(("/1", "/2")):
        name = name[:-2]

    seq = seq.upper()

    char_a, char_b = ['CT', 'GA'][read_i]
    # keep original sequence as name.
    name = " ".join(
        (
            name,
            f"YS:Z:{seq}" + f"\tYC:Z:{['CT', 'GA'][read_i]}\n"
        )
    )
    seq = seq.replace(char_a, char_b)
    out.writelines((f"@{name}", f"{seq}\n", "+\n", f"{qual}\n"))

def get_longest_match(cigartuples):
    """Extract longest M/=/X operation from CIGAR."""
    if not cigartuples:
        return 0
    return max((length for op, length in cigartuples if op in (0, 7, 8)), default=0)

def reverse_complement(seq):
    """Complement DNA sequence."""
    comp_table = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(comp_table)

def get_clip_shifts(cigartuples):
    """Calculate hard clip offsets for sequence restoration."""
    if not cigartuples:
        return 0, None

    left = 0
    for op, length in cigartuples:
        if op == 0:  # M - stop at first match
            break
        if op == 5:  # H - count hard clips
            left += length

    right = 0
    for op, length in reversed(cigartuples):
        if op == 0:  # M - stop at last match
            break
        if op == 5:  # H - count hard clips
            right += length

    return left, -right if right else None

def handle_reads_pysam(alignments, set_as_failed, do_not_penalize_chimeras, clean_header):
    def _rebind(aln):
        d = aln.to_dict()
        if d['ref_name'] != '*':
            d['ref_name'] = d['ref_name'][1:]
        if d['next_ref_name'] not in ('=', '*'):
            d['next_ref_name'] = d['next_ref_name'][1:]
        return pysam.AlignedSegment.from_dict(d, clean_header)

    alignments = list(alignments)
    any_failed = False

    new_alignments = []
    for aln in alignments:
        orig_seq = aln.get_tag('YS') if aln.has_tag('YS') else aln.query_sequence
        if aln.has_tag('YS'):
            aln.set_tag('YS', None)

        direction = aln.reference_name[0] if not aln.is_unmapped else None
        aln = _rebind(aln)

        if aln.is_unmapped:
            aln.query_sequence = orig_seq
            new_alignments.append(aln)
            continue

        assert direction == "f" or direction == "r"
        aln.set_tag('YD', direction, value_type='Z')

        if set_as_failed == direction:
            aln.flag |= 0x200
            any_failed = True

        if not do_not_penalize_chimeras:
            if get_longest_match(aln.cigartuples) < (len(orig_seq) * 0.44):
                aln.flag |= 0x200
                aln.flag &= ~0x2
                aln.mapping_quality = min(aln.mapping_quality, 1)
                any_failed = True

        l_shift, r_shift = get_clip_shifts(aln.cigartuples)
        if not aln.is_reverse:
            aln.query_sequence = orig_seq[l_shift:r_shift]
        else:
            aln.query_sequence = reverse_complement(orig_seq[::-1][l_shift:r_shift])

        new_alignments.append(aln)

    if any_failed:
        for aln in new_alignments:
            aln.flag |= 0x200
            aln.flag &= ~0x2

    return new_alignments

@app.command()
def index(reference_fasta: Path):
    c2t_fasta = reference_fasta.with_suffix(".bwameth.c2t.fasta.gz")
    eprint("--------------------")
    eprint(f"Converting C2T in {reference_fasta} to {c2t_fasta}")
    try:
        with xopen(c2t_fasta, "w", format="gz") as fh:
            # for header, seq in fasta_iter(reference_fasta):
            for name, seq in pyfastx.Fasta(str(reference_fasta), build_index=False):
                ########### Reverse ######################
                rheader = f">r{name}\n"
                fh.write(rheader)
                fh.writelines(list(wrap(seq.replace("G", "A"))))

                ########### Forward ######################
                fheader = f">f{name}\n"
                fh.write(fheader)
                fh.writelines(list(wrap(seq.replace("C", "T"))))
    except:
        c2t_fasta.unlink()
        raise

    sys.exit(bwa_index(reference_fasta, c2t_fasta))

@app.command()
def c2t(read_1: Path, read_2: Annotated[Path, typer.Argument()], output: Annotated[Path, typer.Argument()] = Path("/dev/stdout")):
    eprint("--------------------")
    eprint(f"Converting reads in {read_1} and {read_2}")
    fq1_iter = pyfastx.Fastq(str(read_1), build_index=False, full_name=False) # set to false, QNAME need to match between read pairs in interleaved intermidate form
    fq2_iter = pyfastx.Fastq(str(read_2), build_index=False, full_name=False)
    iterleaved_records = chain.from_iterable(zip(fq1_iter, fq2_iter, strict=True))
    
    lt_80bp_count = 0

    with output.open("w") as out:
        for i, (name, seq, qual) in enumerate(iterleaved_records):
            if name == '' or name is None:
                continue
            if len(seq) < 80:
                lt_80bp_count += 1
            convert_and_write_read(name, seq, qual, (i % 2), out)

    if lt_80bp_count > 50:
        eprint(f"WARNING: {lt_80bp_count} reads with length < 80bp")
        eprint( "       : this program is designed for long reads")

    sys.exit(0)

@app.command()
def align(
    index_prefix: Path, 
    threads: Annotated[int, typer.Argument()] = min(4, (os.cpu_count() or 1)),
    # read_group: Annotated[Optional[str], typer.Argument()] = None,
    extra_args: Annotated[Optional[list[str]], typer.Argument()] = None,
    set_as_failed: Annotated[Optional[bool], typer.Argument()] = None,
    do_not_penalize_chimeras: Annotated[Optional[bool], typer.Argument()] = None,
):

    align_cmd = None
    if all(index_prefix.with_suffix(ext).exists() for ext in BWA_EXTENSIONS):
        align_cmd = sh.bwa.bake(shlex.split("mem -T 40 -B 2 -L 10 -CM -U 100 -p")) # type: ignore
        eprint("--------------------")
        eprint("Found bwa index")
    elif all(index_prefix.with_suffix(ext).exists() for ext in BWA_MEM2_EXTENSIONS):
        align_cmd = sh.bwa_mem2.bake(shlex.split("mem -T 40 -B 2 -L 10 -CM -U 100 -p")) # type: ignore
        eprint("--------------------")
        eprint("Found bwa-mem2 index")
    else:
        raise BWAMethException("Index not found, have your ran `bwameth index`?")

    # if not read_group is None and not read_group.startswith('@RG'):
    #     read_group = '@RG\\tID:{rg}\\tSM:{rg}'.format(rg=read_group)
    
    # baked_command = align_cmd.bake("-R", f"'{read_group}'", "-t", threads, *(extra_args or []), index_prefix, "/dev/stdin")
    baked_command = align_cmd.bake("-t", threads, *(extra_args or []), index_prefix, "/dev/stdin")
    eprint(f"Running: {str(baked_command)}")
    eprint("--------------------")

    r_fd, w_fd = os.pipe()
    bg_proc = baked_command(_in=sys.stdin, _out=w_fd, _bg=True)
    os.close(w_fd)    
    
    try:
        with pysam.AlignmentFile(os.fdopen(r_fd, "rb"), "r") as in_samfile:
            header_dict = in_samfile.header.to_dict()
            new_sq = []
            for sq in header_dict.get('SQ', []):
                name = sq['SN']
                if name.startswith('r'):
                    continue
                if name.startswith('f'):
                    sq['SN'] = name[1:]

                new_sq.append(sq)

            header_dict['SQ'] = new_sq

            new_pg = {
                'ID': 'bwa-meth',
                'PN': 'bwa-meth',
                'VN': __version__,
                'CL': " ".join(sys.argv).replace("\t", "\\t")
            }

            header_dict['PG'].append(new_pg)

            clean_header = pysam.AlignmentHeader.from_dict(header_dict)

            with pysam.AlignmentFile("/dev/stdout", "w", header=clean_header) as out_samfile:
                for _, group in groupby(in_samfile, attrgetter('query_name')):
                    processed_alignments = handle_reads_pysam(group, set_as_failed, do_not_penalize_chimeras, clean_header)
                    for aln in processed_alignments:
                        out_samfile.write(aln)
    except Exception:
         bg_proc.kill()
         raise
    finally:
        bg_proc.wait()
    
    sys.exit(0)

def version_callback(value: bool):
    if value:
        typer.echo(f"v{__version__}")
        raise typer.Exit()

@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        help="Show version and exit."
    )
):
    pass

if __name__ == "__main__":
    app()
