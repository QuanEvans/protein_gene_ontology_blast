import os
from Bio import SeqIO
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm
import argparse
from utils import run_blast, obo_tools

docstring = """
Blast alginment based Protein Gene Ontology Annotation Pipeline (SAGP)

example usage:
    python sagp.py -w /home/username/workdir -f /home/username/seq.fasta -d /home/username/database -g /home/username/goa -t 8

keywords arguments:
    -w, --workdir: working directory for the pipline
    -f, --fasta: path of fasta file
    -d, --database: dir path of database, should include three sub directories: BP, MF, CC, each sub directory should contain a fasta file named SAGP.fasta
    -g, --goa: dir path of goa labels file, should include three files: BP_Term, MF_Term, CC_Term, each file should have multiple lines, each line contains a protein name and its go terms separated by tab, the go terms separated by comma
    -t, --threads: number of threads, default is the number of cores of the machine
"""


class SAGP:

    ### class variables ###
    global oboTools
    oboTools = obo_tools.ObOTools()
    global blast_hits
    blast_hits = dict()
    global aspect_goTerm_labels
    aspect_goTerm_labels = dict()
    ### END ###

    def __init__(self, pipline_wd:str, fasta_file:str, DatabaseDir:str, goa_dir:str, num_threads:int=mp.cpu_count()):
        """init function

        Args:
            pipline_wd (str): the SAGP may called by higher level pipline, this is the working directory of the higher level pipline
            fasta_file (str): path of fasta file
            config_obj (Config): config object that contain the database directory, goa directory, and number of threads information
        """
        
        self.work_dir = pipline_wd
        self.result_file = os.path.join(self.work_dir, "submission.tsv")
        self.fasta_file = fasta_file
        self.DatabaseDir = DatabaseDir
        self.goa_dir = goa_dir
        self.num_threads = num_threads
        self.aspects = ["BP", "MF", "CC"] # biological process, molecular function, cellular component
        self.tmp = os.path.join(self.work_dir, "tmp")
        #shutil.rmtree(self.tmp, ignore_errors=True)
        os.makedirs(self.tmp, exist_ok=True)
    
    def get_seq_dict(self, fast_file:str)->dict:
        """read fasta file and return a dict

        Args:
            fast_file (str): path of fasta file

        Returns:
            dict: dict of fasta file, key is protein name, value is protein sequence
        """
        seq_dict = {}
        for record in SeqIO.parse(fast_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
        return seq_dict

    def read_labels(self, filename:str)->dict:
        """read labels from file
            file should have multiple lines, each line contains a protein name and its go terms separated by tab
            the go terms separated by comma
            example:
                protein1 go1,go2,go3

        Args:
            filename (str): path to the label file

        Returns:
            dict: dict of labels, key is protein name, value is a set of go terms
        """
        with open(filename, 'r') as f:
            go_term_dict = dict()
            for line in f:
                line = line.strip().split()
                if len(line) > 1:
                    name = line[0]
                    go_terms = line[1:]
                    go_terms = [i.split(',') for i in go_terms]
                    go_terms = [j for i in go_terms for j in i]
                    go_term_dict[name] = set(go_terms)
                else:
                    print(line)
        return go_term_dict
    
    def annotate_protein(self, name:str):
        """
        annotate a protein

        the cscore is calculated as follows:
            cscore = sum of bitscore of all hits that have the term / sum of bitscore of all hits

        Args:
            sub_dir (str): protein directory, is a sub directory of SAGP directory
            name (str): protein name
            sequence (str): protein sequence
        
        Returns: None
            Instead of returning, this function will wirte the following results files to the sub_dir:
                seq.fasta: fasta file of the protein
                BP_SAGP_result.tsv: SAGP result of the protein for BP aspect
                MF_SAGP_result.tsv: SAGP result of the protein for MF aspect
                CC_SAGP_result.tsv: SAGP result of the protein for CC aspect
            Note: all the annotation results are sorted by cscore, highest cscore first, and were backpropagated for parent terms
        """
        # annotate protein

        protein_result = dict()

        # for each aspect, get blast hits
        for aspect in self.aspects:
            # get the database labels, key is protein name, value is a set of go terms
            cur_apsect_protein_label = aspect_goTerm_labels[aspect]
            # get current protein blast hits
            cur_blast_hits = blast_hits[aspect].get(name, None)
            if cur_blast_hits is None:
                # no hits for current protein so skip
                continue
            # weighted the bitscore by identity
            cur_blast_hits['ident'] = cur_blast_hits.apply(lambda x: x['nident'] / max(x['qlen'], x['slen']), axis=1)
            #cur_blast_hits['ident'] = cur_blast_hits.apply(lambda x: x['nident'] /  x['qlen'], axis=1)
            cur_blast_hits['bitscore'] = cur_blast_hits['bitscore'] * cur_blast_hits['ident']
            cur_blast_hits = cur_blast_hits[['target','bitscore']]
            # conver to dict
            target_bitscore_dict = dict(zip(cur_blast_hits['target'], cur_blast_hits['bitscore']))

            # annotate protein

            # get all go terms
            term_list = set()
            for target in target_bitscore_dict:
                if target in cur_apsect_protein_label:
                    term_list.update(cur_apsect_protein_label[target])

            result_dict = dict()
            for term in term_list:
                sum_bitscore = sum(target_bitscore_dict.values()) # sum of all bitscore from all targets
                term_sum_bitscore = sum(target_bitscore_dict[target] for target in target_bitscore_dict if term in cur_apsect_protein_label[target]) # sum of all bitscore from all targets that have current term
                result_dict[term] = term_sum_bitscore / sum_bitscore # cscore
            
            # update protein_result
            protein_result[aspect] = result_dict

        protein_result_combined = dict()
        for aspect in self.aspects:
            if aspect not in protein_result:
                continue
            for term, score in protein_result[aspect].items():
                if term in protein_result_combined:
                    print(f"Error: {name} have duplicated term {term}")
                else:
                    protein_result_combined[term] = score
        
        # backpropagate cscore
        protein_result_combined = oboTools.backprop_cscore(protein_result_combined, min_cscore=0.001)
        
        return protein_result_combined, name
                  
    def get_blast_hits(self, workdir:str, fasta_file:str, threads:int=mp.cpu_count())->None:
        """run alignment, use fasta file as query, SAGP as database, extract all hits

        Args:
            workdir (str): working directory for saving results
            fasta_file (str): path of fasta file
            threads (int, optional): number of core to use. Defaults to mp.cpu_count().
        """

        for aspect in self.aspects:
            hits_tsv_name = os.path.join(workdir, f"{aspect}_blast.tsv")
            database = os.path.join(self.DatabaseDir, aspect, 'SAGP.fasta') ## SAGP database
            processed_file = None

            if os.path.exists(hits_tsv_name):
                # check whether it is complete
                first_line = open(hits_tsv_name).readline()
                if first_line.startswith("query"):
                    continue
                
                columns = ["query", "target", "bitscore", "pident" , "evalue", "qlen", "slen", "nident"]
                processed = pd.read_csv(hits_tsv_name, sep='\t', names=columns)
                processed_id = processed['query'].unique().tolist()
                fasta_dict = self.get_seq_dict(fasta_file)
                # remove the last num core query, because it is we don't know whether it is complete
                processed_id = processed_id[:-mp.cpu_count()]
                processed_id_set = set(processed_id)
                fasta_dict = {k: v for k, v in fasta_dict.items() if k not in processed_id_set}
                with open(os.path.join(self.tmp, f"{aspect}_tmp.fasta"), "w") as f:
                    for name, seq in fasta_dict.items():
                        f.write(f">{name}\n{seq}\n")
                fasta_file = os.path.join(self.tmp, f"{aspect}_tmp.fasta")

                processed = processed[processed['query'].isin(processed_id_set)]
                tmp_tsv_name = os.path.join(self.tmp, f"{aspect}_tmp.tsv")
                # check if tmp tsv file exists
                if os.path.exists(tmp_tsv_name):
                    processed_tmp = pd.read_csv(tmp_tsv_name, sep='\t', names=columns)
                    # if exists, then merge
                    processed = pd.concat([processed_tmp, processed])

                processed.to_csv(tmp_tsv_name, sep='\t', index=False, header=False)
                processed_file = tmp_tsv_name

            run_blast.run_blast(fasta_file, database, threads=threads, output_file=hits_tsv_name, processed_file=processed_file)
        
    def main(self):
        """
        sub pipline for sagp
        """
        # read fasta file
        seq_dict = self.get_seq_dict(self.fasta_file)

        # run blast
        self.get_blast_hits(self.work_dir, self.fasta_file, threads=self.num_threads)

        # read blast hits
        # blast_hits is already a global variable
        for aspect in self.aspects:
            hits_tsv_name = os.path.join(self.work_dir, f"{aspect}_blast.tsv")
            blast_hits_aspect_df = pd.read_csv(hits_tsv_name, sep='\t')
            blast_hits_aspect_df_grouped = blast_hits_aspect_df.groupby('query')
            blast_hits_aspect_df_grouped_dict = dict(list(blast_hits_aspect_df_grouped))
            blast_hits[aspect] = blast_hits_aspect_df_grouped_dict

        # read labels
        # aspect_goTerm_labels is already a global variable
        for aspect in self.aspects:
            aspect_goTerm_labels[aspect] = self.read_labels(os.path.join(self.goa_dir, f"{aspect}_Term"))

        # annotate proteins using multiprocessing
        args_list = []
        for name in seq_dict:
            args_list.append(name)

        print("Generating SAGP results...")
        with Pool(self.num_threads) as p:
            final_result = p.map(self.annotate_protein, tqdm(args_list, total=len(args_list), ascii=' >='))
        
        # # write results
        with open(self.result_file, 'w') as f:
            f.write('EntryID\tterm\tscore\n')
            for result, name in final_result:
                for term, score in result.items():
                    f.write(f"{name}\t{term}\t{score}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=docstring, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-w", "--workdir", help="working directory", required=True)
    parser.add_argument("-f", "--fasta", help="path of fasta file", required=True)
    parser.add_argument("-d", "--database", help="dir path of database", required=True)
    parser.add_argument("-g", "--goa", help="dir path of goa labels file", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", default=mp.cpu_count(), type=int)
    args = parser.parse_args()
    sagp = SAGP(args.workdir, args.fasta, DatabaseDir=args.database, goa_dir=args.goa, num_threads=args.threads)
    sagp.main()