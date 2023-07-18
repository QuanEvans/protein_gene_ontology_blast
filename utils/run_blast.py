import subprocess
import time
import pandas as pd
import argparse
from io import StringIO
from Bio import SeqIO
import os

doc = """This script can be imported or run from the command line.
It takes a fasta file and a blast database and runs blastp.
It can either write the output to a file or return a pandas dataframe.

Example usage:
    python run_blast.py query.fasta db.fasta --output_file output.tsv --add_seq --add_self
"""


def get_current_time():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def record_time(func):
    def wrapper(*args, **kwargs):
        func_name = func.__name__
        start_time = time.time()
        print("start %s: %s" % (func_name, get_current_time()))
        result = func(*args, **kwargs)
        end_time = time.time()
        total_time = round(end_time - start_time, 2)
        print("end %s: %s" % (func_name, get_current_time()))
        print("Total time for %s: %s seconds" % (func_name, total_time))
        return result
    return wrapper

@record_time
def run_blast(
        fasta_file,
        database,
        blastp_path:str="blastp",
        evalue:float=10,
        threads:int=8,
        add_seq:bool=False,
        add_self:bool=False,
        self_score_path:str=os.path.join(os.path.dirname(__file__), "SelfScore"),
        output_file:str=None,
        return_df:bool=False,
        processed_file:str=None,
)->pd.DataFrame:
    """
    Run blastp with the given fasta file and database.
    If output file is given, write to it.
    If not, return a pandas dataframe.

    Args:
        fasta_file: fasta file
        database: blast database
        blastp_path: path to blastp, default is blastp
        output_file: output file
        evalue: evalue cutoff, default is 0.1
        threads: number of threads, default is 8
        add_seq: add sequence to the output, default is False
                note, if add seq, a fai file is required on the same dir as the database
                if fai file is not found, it will be created, but samtools is required
        add_self: add self hits to the output, default is False
        self_score_path: path to CPP script SelfScore, default is ./SelfScore
        return_df: return a pandas dataframe, default is False
    """
    print(f"searching {fasta_file} against {database}")
    cmd = [
        blastp_path,
        "-query", fasta_file,
        "-db", database,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-outfmt", "6 qseqid sseqid bitscore pident evalue qlen slen nident",
    ]
    if output_file is not None:
        # if output file is given, write to it
        cmd.extend(["-out", output_file])
    
    P = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = P.communicate()
    if P.returncode != 0:
        raise Exception("Blastp failed with error: %s" % err)
    
    # read in pandas dataframe
    columns = ["query", "target", "bitscore", "pident" , "evalue", "qlen", "slen", "nident"]
    if output_file is not None:
        if processed_file is not None:
            # if unfinished file is given, cat the current output file with the unfinished file
            lines = []
            with open(processed_file, "r") as f:
                lines.extend(f.readlines())
            with open(output_file, "r") as f:
                lines.extend(f.readlines())
            with open(output_file, "w") as f:
                f.writelines(lines)
        df = pd.read_csv(output_file, sep="\t", names=columns)
    else:
        tsv_out = StringIO(out.decode("utf-8"))
        df = pd.read_csv(tsv_out, sep="\t", names=columns)
    
    # if '|'  in df.target.values and len(df.target.values[0].split('|')) > 1:
    #     df.target = df.target.str.split("|").str[1]
    df.target = df.target.apply(lambda x: x.split("|")[1] if "|" in x else x)

    if add_self:
        # add self hits
        run_SelfScore = [
            self_score_path, fasta_file
        ]
        P = subprocess.Popen(run_SelfScore, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = P.communicate()
        if P.returncode != 0:
            raise Exception("SelfScore failed with error: %s" % err)
        str_data = StringIO(out.decode("utf-8"))
        columns = ['target', 'length', 'score', 'bitscore']
        self_df = pd.read_csv(str_data, sep="\t", names=columns, comment="#")
        self_df['qlen'] = self_df['length']
        self_df['slen'] = self_df['length']
        self_df['nident'] = self_df['length']
        self_df = self_df[['target', 'bitscore','qlen', 'slen', 'nident']]
        self_df['query'] = self_df['target'].apply(lambda x: x)
        self_df['pident'] = 100
        self_df['evalue'] = 0
        df = pd.concat([df, self_df], axis=0, ignore_index=True)
    
    if add_seq:
        # add sequence
        fai_file = database + ".fai"
        if not os.path.exists(fai_file):
            # if fai file does not exist, create it
            cmd = ["samtools", "faidx", database]
            P = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = P.communicate()
            if P.returncode != 0:
                raise Exception("Samtools faidx failed with error: %s" % err)
        # read in fasta file
        seq_dict = SeqIO.to_dict(SeqIO.parse(database, "fasta"))
        df['seq'] = df['target'].apply(lambda x: str(seq_dict[x].seq))

    if add_seq:
        # add for hit sequences
        hit_seqs_id = df[df['query'] != df['target']].target.unique()
        fai_file = database + ".fai"
        # cheking if fai file exists
        if not os.path.exists(fai_file):
            print("fai file not found, creating...")
            # using samtools to create fai file
            cmd = [
                "samtools", "faidx", database
            ]
            subprocess.run(cmd, check=True)

        hit_seqs = get_seq(fai_file,database,hit_seqs_id)

        # using biopython to get query sequences
        query_seqs = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            query_seqs[record.id] = str(record.seq)
        # adding sequences to dataframe
        df['target_seq'] = df['target'].apply(lambda x: hit_seqs[x] if x in hit_seqs else query_seqs[x])

    if output_file is not None:
        # if output file is given, write to it
        df.to_csv(output_file, sep="\t", index=False)
    if return_df:
        # if return_df is True, return a pandas dataframe
        return df

    
def get_sequence(
        sequence_id:str,
        database:str,
        blastdbcmd_path:str="blastdbcmd",
        output_file:str=None,
)->str:
    """
    get sequence from blast database using blastdbcmd
    might be slow due to subprocess call and would not be used in the pipeline
    
    Args:
        sequence_id: sequence id
        database: blast database
        blastdbcmd_path: path to blastdbcmd, default is blastdbcmd
        output_file: output file, if not given, return sequence as string
    """
    cmd = [
        blastdbcmd_path,
        "-db", database,
        "-entry", sequence_id,
    ]
    if output_file is not None:
        cmd.extend(["-out", output_file])
        subprocess.run(cmd, check=True)
    else:
        P = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = P.communicate()
        if P.returncode != 0:
            raise Exception("Blastdbcmd failed with error: %s" % err)
        sequence = ''.join(out.decode("utf-8").splitlines()[1:])
        return sequence

def read_fai(fai_file: str) -> pd.DataFrame:
    """
    Read a fai file and return a pandas dataframe
    param fai_file: the fai file path
    return: a pandas dataframe
    """
    df = pd.read_csv(fai_file, sep='\t', header=None)
    df.columns = ['seqid', 'length', 'offset', 'line_bases', 'line_width']
    # because the seqid in the dat file is not the same as the seqid in the fai file
    # we need to convert the seqid in the fai file to the seqid in the dat file
    # example: seqid in the dat file is Q6GZX4, but the seqid in the fai file is sp|Q6GZX4|001R_FRG3G
    # note: this may not work for all the cases; print warning
    if '|' in df['seqid'].values[0]:
        print(f"| present in the seqid in the fai file: {fai_file}, function {read_fai.__name__} with in the module {__name__} may not work properly")
    df['seqid'] = df['seqid'].apply(lambda x: x.split('|')[1] if '|' in x else x)
    # make the seqid the index
    df.set_index('seqid', inplace=True)
    return df

def get_seq(fai_file, fasta_file, seqid_list: list, out_file=None)-> dict:
    """
    Get the sequence by seqid
    param fai_file: the fai file path
    param fasta_file: the fasta file path
    param seqid_list: a list of sequence id
    param out_type: the output type, default is dict, can be csv file
    """
    
    seqid_dict = {} # dict to store the sequence id (key) and sequence (value)
    df = read_fai(fai_file) # read the index file

    with open(fasta_file, 'r') as f:
        for seqid in seqid_list:
            # get the offset and length of the sequence
            if seqid not in df.index:
                seqid_dict[seqid] = None
                continue
            offset = df.loc[seqid, 'offset']
            length = df.loc[seqid, 'length']
            line_bases = df.loc[seqid, 'line_bases']

            read_length = length + (length // line_bases) # the length of the sequence in the fasta file
            # seek to the offset
            f.seek(offset)
            # read the sequence
            seq = f.read(read_length)
            # remove the newline character
            seq = seq.replace('\n', '')
            # store the sequence id and sequence in the dict
            seqid_dict[seqid] = seq

    if out_file is None:
        return seqid_dict
    else:
        assert out_file.endswith('.csv'), 'out_file must be a csv file'
        df = pd.DataFrame.from_dict(seqid_dict, orient='index')
        df.to_csv(out_file, header=False)
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("query", help="fasta file")
    parser.add_argument("db", help="blast database")
    parser.add_argument("-blastp","--blastp_path", help="path to blastp", default="blastp")
    parser.add_argument("-o","--output_file", help="output file")
    parser.add_argument("-e","--evalue", help="evalue cutoff", default=0.1)
    parser.add_argument("-n","--threads", help="number of threads", default=8)
    args = parser.parse_args()
    run_blast(
        args.query,
        args.db,
        blastp_path=args.blastp_path,
        output_file=args.output_file,
        evalue=args.evalue,
        threads=args.threads,
    )