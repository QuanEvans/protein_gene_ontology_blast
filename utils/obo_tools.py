import os
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
from collections import defaultdict
from typing import List, Tuple, Dict
import sys
import shutil
import numpy as np

# default paths
file_dir = os.path.dirname(os.path.abspath(__file__))
default_is_a = os.path.join(file_dir, 'obo2csv', 'is_a.csv')
default_names = os.path.join(file_dir, 'obo2csv', 'names.csv')
default_alt_id = os.path.join(file_dir, 'obo2csv', 'alt_id.csv')
default_go_obo = os.path.join(file_dir, 'obo2csv', 'go-basic.obo')
default_obo2csv = os.path.join(file_dir, 'obo2csv', 'obo2csv')

class ObOTools:

    def __init__(self,
        is_a: str = default_is_a, # Path to the is_a.csv file
        names: str = default_names, # Path to the names.txt file
        go_obo: str = default_go_obo, # Path to the go.obo file
        alt_id: str = default_alt_id, # Path to the alt_id.csv file
        obo_url: str = "http://purl.obolibrary.org/obo/go/go-basic.obo", # URL to the go.obo file
        obo2csv: str = default_obo2csv, # Path to the obo2csv cpp file
        add_part_of:bool = True, # whether to add the part_of relationship to the is_a.csv file
        ):
        self.names_path = names
        self.alt_id_path = alt_id
        self.go_obo_path = go_obo
        self.obo2csv_path = obo2csv
        self.obo_url = obo_url
        self.add_part_of = int(add_part_of)

        self.aspects = ['BP', 'CC', 'MF']
        self.check_obo(go_obo, obo_url) # check if the go.obo file exists
        self.check_csv(is_a, names, alt_id, obo2csv, add_part_of, go_obo) # check if the is_a.csv, names.csv, alt_id.csv files exist
        self.is_a_path = is_a
        self.is_a_dict, self.term2aspect, self.alt_id_dict, self.is_dirct, self.is_indirct \
            = self.read_csv(is_a, names, alt_id) # read the is_a.csv, names.csv, alt_id.csv files
        
    def update_go_obo(self, go_obo_file:str=None, obo_url:str=None):
        """
        update the go_obo file and regenerate the is_a.csv, names.csv, alt_id.csv files
        could been updated by the go_obo_file or the obo_url,
        if both are provided, the go_obo_file will be used
        if both are not provided, the default obo_url will be used

        Args:
            go_obo_file: the path to the go_obo file
            obo_url: the url to the go_obo file
        """

        if go_obo_file:
            print(f"Updating the go_obo file using: {go_obo_file}")
            # delete the old go_obo is_a.csv, names.csv, alt_id.csv files
            os.remove(self.go_obo_path)
            os.remove(self.is_a_path)
            os.remove(self.names_path)
            os.remove(self.alt_id_path)
            os.remove(self.part_of_path)
            # cp the go_obo_file to the default_go_obo
            shutil.copyfile(go_obo_file, self.go_obo_path)
            # generate the new is_a.csv, names.csv, alt_id.csv files
            self.check_csv(self.is_a_path, self.names_path, self.alt_id_path, self.obo2csv_path, self.go_obo_path)
            return
        
        if go_obo_file is None and obo_url is None:
            obo_url = self.obo_url

        if obo_url:
            print(f"Updating the go_obo file from url: {obo_url}")
            # delete the old go_obo file is_a.csv, names.csv, alt_id.csv files
            os.remove(self.go_obo_path)
            os.remove(self.is_a_path)
            os.remove(self.names_path)
            os.remove(self.alt_id_path)
            os.remove(self.part_of_path)
            # download the new go_obo file
            self.check_obo(self.go_obo_path, obo_url)
            # generate the new is_a.csv, names.csv, alt_id.csv files
            self.check_csv(self.is_a_path, self.names_path, self.alt_id_path, self.obo2csv_path, self.go_obo_path)
            return

    def check_obo(self, go_obo:str, obo_url:str):
        """
        Check if the go.obo file exists, if not, download it from the obo_url
        """

        if not os.path.exists(go_obo):
            print("Downloading the go.obo file from {}".format(obo_url))
            os.system(f"wget {obo_url} -O {go_obo}")
            print("Download complete!")
    
    def check_csv(self, is_a:str, names:str, alt_id:str, add_part_of:int, obo2csv:str, go_obo):
        """
        Check if the is_a.csv, names.csv, alt_id.csv files exist
        if not exit, generate them
        """
        if not os.path.exists(is_a) or not os.path.exists(names) or not os.path.exists(alt_id):
            print("Generating is_a.csv, names.csv, alt_id.csv files")
            os.system(f"{obo2csv} {go_obo} {is_a} {names} {alt_id} {add_part_of}")
            print("Generation complete!")

    def read_csv(self, is_a:str, names:str, alt_id:str)->Tuple[dict,dict,dict]:
        """
        Read the is_a.csv, names.csv, alt_id.csv files
        """
        # first read the alt_id.csv to a dict
        # this would map the alt_id to the main_id
        alt_id_dict = pd.read_csv(alt_id, header=None, names=['alt_id', 'main_id'], sep='\t'\
                                ).set_index('alt_id').to_dict()['main_id'] # convert to dict

        # read the is_a.csv to a dataframe
        is_a_df = pd.read_csv(is_a, header=None, names=['term', 'aspect', 'direct', 'indicrect'],\
                               sep='\t')
        # split the go term into in dircet and indirect into a set
        is_a_df['direct'] = is_a_df['direct'].apply(lambda x: set(x.split(',')) if type(x) == str else set())
        is_a_df['indicrect'] = is_a_df['indicrect'].apply(lambda x: set(x.split(',')) if type(x) == str else set())
        aspect_map = {'P': 'BP', 'F': 'MF', 'C': 'CC'}
        is_a_df['aspect'] = is_a_df['aspect'].map(aspect_map) # map the aspect to BP, MF, CC
        # create a dict to map the term to the aspect
        term2aspect = is_a_df[['term', 'aspect']].set_index('term').to_dict()['aspect']
        # create a dict that maps the term to its direct parents
        is_dirct = is_a_df[['term', 'direct']].set_index('term').to_dict()['direct']
        is_indirct = is_a_df[['term', 'indicrect']].set_index('term').to_dict()['indicrect']
        # dirct and indirect parents
        is_a_dict = defaultdict(set)
        for d in [is_dirct, is_indirct]:
            for k, v in d.items():
                is_a_dict[k].update(v)
        
        # sort all the dict by the key
        is_a_dict = dict(sorted(is_a_dict.items(), key=lambda x: x[0]))
        term2aspect = dict(sorted(term2aspect.items(), key=lambda x: x[0]))
        alt_id_dict = dict(sorted(alt_id_dict.items(), key=lambda x: x[0]))
        is_dirct = dict(sorted(is_dirct.items(), key=lambda x: x[0]))
        is_indirct = dict(sorted(is_indirct.items(), key=lambda x: x[0]))
        
        return is_a_dict, term2aspect, alt_id_dict, is_dirct, is_indirct

        
    def update_parent(self, protein_name:str, cur_terms:set)->Tuple[str,set]:
        """
        Update the cur_terms to the parents of the cur_terms

        Args:
            protein_name: the name of the protein
            cur_terms: the current terms of the protein
        
        Returns:
            protein_name: the name of the protein
            final_terms: the updated terms of the protein
        """
        final_terms = set()
        for term in cur_terms:
            if term not in self.alt_id_dict:
                # warning: the term is not in the alt_id_dict, the obo file may be outdated
                print(f"Warning: {term} is not in the alt_id_dict, the obo file may be outdated, please check protein: {protein_name}")
                final_terms.add(term)
            else:
                all_parent = self.is_a_dict[self.alt_id_dict[term]].copy()
                all_parent.add(self.alt_id_dict[term])
                final_terms.update(all_parent)
        return protein_name, final_terms
    
    def backprop_cscore(self, go_cscore:dict, min_cscore:float=None, sorting=True)->Dict[str, float]:
        """
        backproprate the child cscore to the parents
        parent cscore should not be smaller than the child cscore
        if the parent cscore is smaller than the child cscore, update the child cscore to the parent cscore

        Args:
            go_cscore: the go_cscore dict where the key is the go term and the value is the cscore
            min_cscore: the minimum cscore, if the cscore is smaller than the min_cscore, the cscore will be set to 0
        return:
            go_cscore: the updated go_cscore dict
        """
        if min_cscore is not None:
            filter_score = {term:cscore for term, cscore in go_cscore.items() if cscore > min_cscore}
            backprop_cscore = filter_score.copy()
        else:
            filter_score = go_cscore
            backprop_cscore = go_cscore.copy()

        for term in filter_score:
            if term not in self.alt_id_dict:
                # warning: the term is not in the alt_id_dict, the obo file may be outdated
                #print(f"Warning: {term} is not in the alt_id_dict, the obo file may be outdated", file=sys.stderr)
                continue
            # replace the alt_id with the main_id, and get all the parents for the term
            all_parent = self.is_a_dict[self.alt_id_dict[term]].copy()
            # loop through all the parents, if the parent is not in the backprop_cscore or the parent cscore is smaller than the child cscore
            # update the child cscore to the parent cscore
            for parent in all_parent:
                if parent not in backprop_cscore or backprop_cscore[parent] < backprop_cscore[term]:
                    backprop_cscore[parent] = backprop_cscore[term]
        
        if sorting:
            # high to low
            backprop_cscore = {k: v for k, v in sorted(backprop_cscore.items(), key=lambda item: item[1], reverse=True)}

        return backprop_cscore
    
    def get_go_terms_by_aspect(self, aspect:str):
        """
        Get all the go terms for the aspect

        Args:
            aspect: the aspect, BP, MF, CC
        
        Returns:
            go_terms: a set of go terms
        """
        go_terms = set()
        for term in self.term2aspect:
            if self.term2aspect[term] == aspect:
                go_terms.add(term)
        return go_terms
    
    def generate_child_matrix(self, term_list:List[str]):
        """
        Generate the child matrix for the aspect

        Args:
            go2vec: the go2vec dict where the key is the go term and the value is the index of the go term in the embedding matrix
        
        Returns:
            child_matrix: the child matrix for the aspect where child_matrix[i][j] = 1 if the jth GO term is a subclass of the ith GO term else 0
        """
        
        training_terms = term_list
        #CM_ij = 1 if the jth GO term is a subclass of the ith GO term
        child_matrix = np.zeros((len(training_terms), len(training_terms)))
        # fill diagonal with 1
        np.fill_diagonal(child_matrix, 1)
        for i, term in enumerate(training_terms):
            for j, child in enumerate(training_terms):
                if i == j:
                    continue
                if term in self.is_a_dict[child]:
                    child_matrix[i][j] = 1
        return child_matrix
    
    def toplogical_child_matrix(self, term_list:List[str]):
        """
        Generate the child matrix for the aspect
        the input term_list should be topologically sorted
        Then we can ignore these i rows on the bottom of the matrix where the sum of the row is 1 (itself)

        Args:
            term_list (List[str]): selected go terms for the aspect
        return:
            child_matrix: the child matrix for the aspect where child_matrix[i][j] = 1 if the jth GO term is a subclass of the ith GO term else 0
        """

        # check if the term_list is topologically sorted
        # if not, raise the error
        sorted_terms,leafs = self.top_sort(term_list, return_leaf=True)
        assert len(sorted_terms) == len(term_list), "The input term_list is not topologically sorted"
        for i, term in enumerate(sorted_terms):
            assert term == term_list[i], f"The {i} idx term:{term_list[i]} is not the same as the sorted term:{term}"

        training_terms = term_list
        #CM_ij = 1 if the jth GO term is a subclass of the ith GO term
        child_matrix = np.zeros((len(training_terms), len(training_terms)))
        # fill diagonal with 1
        np.fill_diagonal(child_matrix, 1)
        for i, term in enumerate(training_terms):
            for j, child in enumerate(training_terms):
                if i == j:
                    continue
                if term in self.is_a_dict[child]:
                    child_matrix[i][j] = 1
        num_leaf = len(leafs)
        child_matrix = child_matrix[:-num_leaf,:]
        return child_matrix
    
    def get_num_leaf(self, term_list:List[str]):
        """
        Get the number of leaf go terms for the aspect

        Args:
            term_list (List[str]): selected go terms for the aspect
        
        Returns:
            num_leaf (int): the number of leaf go terms for the aspect
        """
        term_list = sorted(term_list)
        set_term_list = set(term_list)

        # reverse the is_a_dict to get the child dict
        child_dict = defaultdict(set)
        for child, parent_set in self.is_a_dict.items():
            for parent in parent_set:
                if parent in set_term_list and child in set_term_list:
                    child_dict[child].add(parent)

        # create a dict to store the indegree of each go term
        indegree = {term:0 for term in term_list}
        # sort the indegree dict by the key
        indegree = dict(sorted(indegree.items(), key=lambda item: item[0]))
        
        # loop through the child dict to get the indegree of each go term
        for parent, child_set in child_dict.items():
            for child in child_set:
                indegree[child] += 1

        quene = [ term for term in term_list if indegree[term] == 0 ]

        return len(quene)
    
    def top_sort(self, term_list:List[str], return_leaf:bool=False)->List[str]:
        """
        Topological sort the go terms

        Args:
            term_list (List[str]): selected go terms for the aspect

        Returns:
            sorted_terms (List[str]): sorted go terms
        """
        term_list = sorted(term_list)
        set_term_list = set(term_list)

        # reverse the is_a_dict to get the child dict
        child_dict = defaultdict(set)
        for child, parent_set in self.is_a_dict.items():
            for parent in parent_set:
                if parent in set_term_list and child in set_term_list:
                    child_dict[child].add(parent)

        # create a dict to store the indegree of each go term
        indegree = {term:0 for term in term_list}
        # sort the indegree dict by the key
        indegree = dict(sorted(indegree.items(), key=lambda item: item[0]))
        
        # loop through the child dict to get the indegree of each go term
        for parent, child_set in child_dict.items():
            for child in child_set:
                indegree[child] += 1
        
        # find the go terms with indegree 0 (leaves) and add them to the queue
        indexes = []
        visited = set()

        quene = [ term for term in term_list if indegree[term] == 0 ]
        leafs = quene.copy()

        # for each element of the queue increment visits, add them to the list of ordered nodes
        # and decrease the in-degree of the neighbor nodes
        # and add them to the queue if they reach in-degree == 0
        while quene:
            node = quene.pop(0)
            visited.add(node)
            indexes.append(node)
            parents = self.is_a_dict[node].copy()
            parents = parents.intersection(set_term_list)
            # sort the parents
            parents = sorted(parents)
            if parents:
                for parent in parents:
                    indegree[parent] -= 1
                    if indegree[parent] == 0:
                        quene.append(parent)

        if len(visited) != len(term_list):
            print("Warning: the graph is not a DAG, the topological sort may not be correct", file=sys.stderr)
        else:
            if return_leaf:
                return indexes[::-1], leafs
            return indexes[::-1]

