import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from bml_casp15.tool.foldseek import *
import pickle
import numpy as np
from bml_casp15.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit
from bml_casp15.complex_templates_search.parsers import TemplateHit
from bml_casp15.tertiary_structure_generation.iterative_search_pipeline import build_alignment_indices, PrefilterError
# need to add A if using relaxation in alphafold
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

PDB_CHAIN_IDS_UNRELAX = 'BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


class _FastaChain:
    def __init__(self, sequence, description):
        self.sequence = sequence
        self.description = description


def _combine_pdb(pdbs, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for chain_id, pdb in zip(PDB_CHAIN_IDS_UNRELAX, pdbs):
            for line in open(pdb):
                if not line.startswith('ATOM'):
                    continue
                out.write(line[:21] + chain_id + line[22:])


def _parse_fasta(fasta_file):
    fasta_string = open(fasta_file).read()
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append('')
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line
    return sequences, descriptions


def _combine_a3ms(infiles, outfile):
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
            else:
                seqs += [line]

    with open(outfile, 'w') as fw:
        for (desc, seq) in zip(descriptions, seqs):
            fw.write(f"{desc}\n{seq}\n")


def _convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


def _complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def _cal_tmscore(mmalign_program, inpdb, nativepdb):
    cmd = mmalign_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $2}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    # print(tmscore_contents)
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return tmscore


def _split_pdb(complex_pdb, outdir):
    makedir_if_not_exists(outdir)
    chain_pdbs = {}
    pre_chain = None
    fw = None
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line)
            chain_pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
        elif chain_name == pre_chain:
            fw.write(line)
        else:
            fw.close()
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line)
            chain_pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
            pre_chain = chain_name
    fw.close()
    return chain_pdbs


def _create_template_df(templates):
    row_list = []
    for i in range(len(templates)):
        target = templates.loc[i, 'target']
        qaln = templates.loc[i, 'qaln']
        qstart = int(templates.loc[i, 'qstart'])
        qend = int(templates.loc[i, 'qend'])
        taln = templates.loc[i, 'taln']
        tstart = templates.loc[i, 'tstart']
        tend = templates.loc[i, 'tend']
        evalue = float(templates.loc[i, 'evalue'])
        row_dict = dict(template=target,
                        tpdbcode=target[0:4],
                        aln_temp=taln,
                        tstart=tstart,
                        tend=tend,
                        aln_query=qaln,
                        qstart=qstart,
                        qend=qend,
                        evalue=evalue)
        row_list += [row_dict]
    return pd.DataFrame(row_list)


def _get_sequence(inpdb):
    """Enclosing logic in a function to simplify code"""

    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)

    return ''.join(sequence)


def _get_indices(sequence, start):
    """Returns indices for non-gap/insert residues starting at the given index."""
    indices = []
    counter = start
    for symbol in sequence:
        # Skip gaps but add a placeholder so that the alignment is preserved.
        if symbol == '-':
            indices.append(-1)
        # Skip deleted residues, but increase the counter.
        elif symbol.islower():
            counter += 1
        # Normal aligned residue. Increase the counter and append to indices.
        else:
            indices.append(counter)
            counter += 1
    return indices


def _build_query_to_hit_index_mapping(
        hit_query_sequence,
        hit_sequence,
        indices_hit,
        indices_query,
        original_query_sequence):
    # If the hit is empty (no aligned residues), return empty mapping
    if not hit_query_sequence:
        return {}

    # Remove gaps and find the offset of hit.query relative to original query.
    hhsearch_query_sequence = hit_query_sequence.replace('-', '')
    hit_sequence = hit_sequence.replace('-', '')
    hhsearch_query_offset = original_query_sequence.find(hhsearch_query_sequence)

    # Index of -1 used for gap characters. Subtract the min index ignoring gaps.
    min_idx = min(x for x in indices_hit if x > -1)
    fixed_indices_hit = [
        x - min_idx if x > -1 else -1 for x in indices_hit
    ]

    min_idx = min(x for x in indices_query if x > -1)
    fixed_indices_query = [x - min_idx if x > -1 else -1 for x in indices_query]

    # Zip the corrected indices, ignore case where both seqs have gap characters.
    mapping = {}
    for q_i, q_t in zip(fixed_indices_query, fixed_indices_hit):
        if q_t != -1 and q_i != -1:
            if (q_t >= len(hit_sequence) or
                    q_i + hhsearch_query_offset >= len(original_query_sequence)):
                continue
            mapping[q_i + hhsearch_query_offset] = q_t

    return mapping

class Multimer_iterative_generation_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def search_templates(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = self.params['foldseek_pdb_database']
        foldseek_af_database = self.params['foldseek_af_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program,
                                   databases=[foldseek_pdb_database, foldseek_af_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive=True)

    def check_and_rank_templates(self, template_file, outfile):
        templates = pd.read_csv(template_file, sep='\t')
        sort_indices = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            evalue = float(templates.loc[i, 'evalue'])
            if target.find('.atom.gz') > 0 and evalue < 1e-10:
                sort_indices += [i]
        for i in range(len(templates)):
            if i in sort_indices:
                continue
            sort_indices += [i]
        if len(sort_indices) == 0:
            os.system(f"cp {template_file} {outfile}")
            return False

        templates_sorted = copy.deepcopy(templates.iloc[sort_indices])
        templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
        templates_sorted.reset_index(inplace=True, drop=True)
        templates_sorted.to_csv(outfile, sep='\t')
        return True

    def check_template_overlap_regions(self, template_info1, template_info2, chain_id_map, template_path):
        print(template_info1)
        chain1_id = template_info1['chainid']
        chain1_seq = chain_id_map[chain1_id].sequence
        indices_hit = build_alignment_indices(template_info1['aln_temp'], template_info1['tstart'])
        indices_query = build_alignment_indices(template_info1['aln_query'], template_info1['qstart'])
        mapping = _build_query_to_hit_index_mapping(
            template_info1['aln_query'], template_info1['aln_temp'], indices_hit, indices_query,
            chain1_seq)

        # tstart_1 = np.min(np.array(list(mapping.values())))
        tend_1 = np.max(np.array(list(mapping.values())))

        chain2_id = template_info2['chainid']
        chain2_seq = chain_id_map[chain2_id].sequence
        indices_hit = build_alignment_indices(template_info2['aln_temp'], template_info2['tstart'])
        indices_query = build_alignment_indices(template_info2['aln_query'], template_info2['qstart'])
        mapping = _build_query_to_hit_index_mapping(
            template_info2['aln_query'], template_info2['aln_temp'], indices_hit, indices_query,
            chain2_seq)

        tstart_2 = np.min(np.array(list(mapping.values())))
        # tend_2 = np.max(np.array(list(mapping.values())))

        return tend_1 > tstart_2

    def assess_complex_templates(self, chain_id_map, template_infos, template_path):
        template_names = []
        same_template_infos = {}
        for template_info in template_infos:
            if template_info['template'] not in template_names:
                template_names += [template_info['template']]
                same_template_infos[template_info['template']] = [template_info]
            else:
                same_template_infos[template_info['template']] += [template_info]

        for template in same_template_infos:
            if len(same_template_infos[template]) == 1:
                continue
            else:
                for i in range(len(same_template_infos[template])):
                    for j in range(i + 1, len(same_template_infos[template])):
                        if not self.check_template_overlap_regions(same_template_infos[template][i],
                                                                   same_template_infos[template][j],
                                                                   chain_id_map,
                                                                   template_path):
                            return False
        return True

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      template_files,
                                      start_msa_path,
                                      outpath,
                                      template_path,
                                      iteration,
                                      rank_templates_by_monomers=True):

        prev_df = None
        for i, chain_id in enumerate(chain_id_map):
            templates = pd.read_csv(template_files[i], sep='\t')
            curr_df = _create_template_df(templates)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=('_' + str(i), '_' + str(i + 1)))

        min_evalues = []
        for i in range(len(prev_df)):
            evalues = []
            for j, chainid in enumerate(chain_id_map):
                evalue = float(prev_df.loc[i, f'evalue_{j + 1}'])
                evalues += [evalue]
            min_evalues += [np.min(np.array(evalues))]
        prev_df['min_evalue'] = min_evalues
        prev_df = prev_df.sort_values(by=['min_evalue'])

        keep_indices = []
        chain_template_msas = {}
        for chain_id in chain_id_map:
            chain_template_msas[chain_id] = {'desc': [chain_id_map[chain_id].description],
                                             'seq': [chain_id_map[chain_id].sequence]}

        print(prev_df)
        for i in range(len(prev_df)):
            template_infos = []
            for j, chainid in enumerate(chain_id_map):
                template = prev_df.loc[i, f'template_{j + 1}']
                qaln = prev_df.loc[i, f'aln_query_{j + 1}']
                qstart = int(prev_df.loc[i, f'qstart_{j + 1}'])
                qend = int(prev_df.loc[i, f'qend_{j + 1}'])
                taln = prev_df.loc[i, f'aln_temp_{j + 1}']
                tstart = prev_df.loc[i, f'tstart_{j + 1}']
                tend = prev_df.loc[i, f'tend_{j + 1}']
                evalue = float(prev_df.loc[i, f'evalue_{j + 1}'])
                row_dict = dict(chainid=chainid,
                                template=template,
                                tpdbcode=template[0:4],
                                aln_temp=taln,
                                tstart=tstart,
                                tend=tend,
                                aln_query=qaln,
                                qstart=qstart,
                                qend=qend,
                                evalue=evalue)
                template_infos += [row_dict]

            if not self.assess_complex_templates(chain_id_map, template_infos, template_path):
                continue

            if len(keep_indices) >= self.max_template_count:
                break

            keep_indices += [i]
            for j, chainid in enumerate(chain_id_map):
                query_non_gaps = [res != '-' for res in prev_df.loc[i, f'aln_query_{j + 1}']]
                out_sequence = ''.join(_convert_taln_seq_to_a3m(query_non_gaps, prev_df.loc[i, f'aln_temp_{j + 1}']))
                aln_full = ['-'] * len(chain_id_map[chainid].sequence)

                qstart = int(prev_df.loc[i, f'qstart_{j + 1}'])
                qend = int(prev_df.loc[i, f'qend_{j + 1}'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)
                chain_template_msas[chainid]['desc'] += [prev_df.loc[i, f'template_{j + 1}']]
                chain_template_msas[chainid]['seq'] += [taln_full_seq]

        msa_out_path = outpath  # + '/msas'
        makedir_if_not_exists(msa_out_path)

        out_msas = []
        for chain_id in chain_id_map:
            start_msa = start_msa_path + '/' + chain_id_map[chain_id].description + '.start.a3m'
            fasta_chunks = (f">{chain_template_msas[chain_id]['desc'][i]}\n{chain_template_msas[chain_id]['seq'][i]}"
                            for i in range(len(chain_template_msas[chain_id]['desc'])))
            with open(start_msa + '.temp', 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')
            _combine_a3ms([start_msa, f"{start_msa}.temp"],
                          f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration{iteration}.a3m")
            out_msas += [f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration{iteration}.a3m"]

        interact_dict = {}
        msa_len = -1
        for i in range(0, len(out_msas)):
            current_len = len(open(out_msas[i]).readlines())
            if msa_len == -1:
                msa_len = current_len
            elif current_len != msa_len:
                raise Exception(f"The length of each msas are not equal! {out_msas}")
            interact_dict[f'index_{i + 1}'] = [i + 1 for i in range(int(msa_len / 2) - 1)]
        interact_df = pd.DataFrame(interact_dict)
        interact_csv = outpath + f'/interaction.iteration{iteration}.csv'
        interact_df.to_csv(interact_csv)

        if rank_templates_by_monomers:
            top_template_files = []
            for template_file, chainid in zip(template_files, chain_id_map):
                templates = pd.read_csv(template_file, sep='\t')
                keep_indices = []
                for i in range(len(templates)):
                    hit = TemplateHit(index=i,
                                      name=templates.loc[i, 'target'].split('.')[0],
                                      aligned_cols=int(templates.loc[i, 'alnlen']),
                                      query=templates.loc[i, 'qaln'],
                                      hit_sequence=templates.loc[i, 'taln'],
                                      indices_query=build_alignment_indices(templates.loc[i, 'qaln'],
                                                                             templates.loc[i, 'qstart']),
                                      indices_hit=build_alignment_indices(templates.loc[i, 'taln'],
                                                                           templates.loc[i, 'tstart']),
                                      sum_probs=0.0)
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=chain_id_map[chainid].sequence)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    keep_indices += [i]
                templates_sorted = copy.deepcopy(templates.iloc[keep_indices])
                templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
                templates_sorted.reset_index(inplace=True, drop=True)
                templates_sorted.to_csv(template_file + f'.top{self.max_template_count}', sep='\t')
                top_template_files += [template_file + f'.top{self.max_template_count}']
            return top_template_files, out_msas, interact_csv

        prev_df.iloc[keep_indices].to_csv(outpath + '/complex_templates.csv')
        return [outpath + '/complex_templates.csv'], out_msas, interact_csv

    def copy_atoms_and_unzip(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            if template_pdb.find('.pdb.gz') > 0:
                os.system(f"cp {self.params['foldseek_af_database_dir']}/{template_pdb} {outdir}")
            else:
                os.system(f"cp {self.params['foldseek_pdb_database_dir']}/{template_pdb} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

    def search(self, fasta_file, input_pdb_dir, outdir, native_pdb_dir=""):

        input_pdb_dir = os.path.abspath(input_pdb_dir)

        fasta_file = os.path.abspath(fasta_file)

        targetname = pathlib.Path(fasta_file).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        sequences, descriptions = _parse_fasta(fasta_file)

        chain_id_map = {}
        for chain_id, sequence, description in zip(PDB_CHAIN_IDS_UNRELAX, sequences, descriptions):
            chain_id_map[chain_id] = _FastaChain(sequence=sequence, description=description)

        native_pdb = ""
        if os.path.exists(native_pdb_dir):
            native_pdb = outdir + '/' + '_'.join(
                [chain_id_map[chain_id].description for chain_id in chain_id_map]) + '.pdb'
            _combine_pdb(
                [native_pdb_dir + '/' + chain_id_map[chain_id].description + '.atom' for chain_id in chain_id_map],
                native_pdb)

        iteration_scores = {}

        true_tm_scores = {}

        iteration_result_all = {'targetname': [],
                                'model': [],
                                'start_lddt': [],
                                'end_lddt': [],
                                'start_tmscore': [],
                                'end_tmscore': []}

        iteration_result_avg = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [],
                                'end_tmscore': []}

        cwd = os.getcwd()

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            current_ref_dir = input_pdb_dir
            ref_start_pdb = f"ranked_{i}.pdb"
            ref_start_pkl = f"result_model_{i + 1}_multimer.pkl"

            model_iteration_scores = []
            model_iteration_tmscores = []

            print(f"Start to refine {ref_start_pdb}")

            for num_iteration in range(self.max_iteration):
                os.chdir(cwd)
                current_work_dir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(current_work_dir)

                start_pdb = f"{current_work_dir}/start.pdb"
                start_msa_path = f"{current_work_dir}/start_msas"
                start_pkl = f"{current_work_dir}/start.pkl"

                os.system(f"cp {current_ref_dir}/{ref_start_pdb} {start_pdb}")
                os.system(f"cp {current_ref_dir}/{ref_start_pkl} {start_pkl}")
                if os.path.exists(start_msa_path):
                    os.system(f"rm -rf {start_msa_path}")
                makedir_if_not_exists(start_msa_path)

                for chainid in chain_id_map:
                    os.system(f"cp {current_ref_dir}/msas/{chain_id_map[chainid].description}.paired.a3m "
                              f"{start_msa_path}/{chain_id_map[chainid].description}.start.a3m")

                ref_avg_lddt = 0
                with open(start_pkl, 'rb') as f:
                    prediction_result = pickle.load(f)
                    ref_avg_lddt = np.mean(prediction_result['plddt'])

                ref_tmscore = 0
                if os.path.exists(native_pdb):
                    ref_tmscore = _cal_tmscore(self.params['mmalign_program'], start_pdb, native_pdb)

                model_iteration_scores += [ref_avg_lddt]
                model_iteration_tmscores += [ref_tmscore]

                out_model_dir = f"{current_work_dir}/alphafold"

                if not _complete_result(out_model_dir):

                    chain_pdbs = _split_pdb(start_pdb, current_work_dir)

                    template_files = []

                    out_template_dir = current_work_dir + '/templates'
                    makedir_if_not_exists(out_template_dir)

                    for chain_id in chain_pdbs:

                        if chain_id not in chain_id_map:
                            raise Exception("Multimer fasta file and model doesn't match!")

                        monomer_work_dir = current_work_dir + '/' + chain_id_map[chain_id].description
                        makedir_if_not_exists(monomer_work_dir)
                        os.system(
                            f"mv {chain_pdbs[chain_id]} {monomer_work_dir}/{chain_id_map[chain_id].description}.pdb")
                        foldseek_res = self.search_templates(
                            inpdb=f"{monomer_work_dir}/{chain_id_map[chain_id].description}.pdb",
                            outdir=monomer_work_dir + '/foldseek')

                        if not self.check_and_rank_templates(foldseek_res,
                                                             f"{monomer_work_dir}/structure_templates.csv"):
                            print(
                                f"Cannot find any templates for {chain_id_map[chain_id].description} in iteration {num_iteration + 1}")
                            break

                        template_files += [f"{monomer_work_dir}/structure_templates.csv"]

                        self.copy_atoms_and_unzip(template_csv=f"{monomer_work_dir}/structure_templates.csv",
                                                  outdir=out_template_dir)

                    if len(template_files) != len(chain_id_map):
                        break

                    template_files, msa_files, msa_pair_file = self.concatenate_msa_and_templates(
                        chain_id_map=chain_id_map,
                        template_files=template_files,
                        start_msa_path=start_msa_path,
                        template_path=out_template_dir,
                        outpath=current_work_dir,
                        iteration=num_iteration + 1)

                    makedir_if_not_exists(out_model_dir)

                    if len(template_files) == 1:
                        cmd = f"python run_alphafold_multimer_custom_sim.py " \
                              f"--fasta_path {fasta_file} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--a3ms {','.join(msa_files)} " \
                              f"--msa_pair_file {msa_pair_file} " \
                              f"--temp_struct_csv {template_files[0]} " \
                              f"--struct_atom_dir {out_template_dir} " \
                              f"--output_dir {out_model_dir}"
                    else:
                        cmd = f"python run_alphafold_multimer_custom_sim.py " \
                              f"--fasta_path {fasta_file} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--a3ms {','.join(msa_files)} " \
                              f"--msa_pair_file {msa_pair_file} " \
                              f"--monomer_temp_csvs {','.join(template_files)} " \
                              f"--struct_atom_dir {out_template_dir} " \
                              f"--output_dir {out_model_dir}"

                    try:
                        os.chdir(self.params['alphafold_program_dir'])
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)

                max_lddt_score = 0
                max_index = -1
                for j in range(0, 5):
                    new_pkl = f"{out_model_dir}/result_model_{j + 1}_multimer.pkl"
                    with open(new_pkl, 'rb') as f:
                        new_prediction_result = pickle.load(f)
                        new_avg_lddt = np.mean(new_prediction_result['plddt'])
                        if new_avg_lddt > max_lddt_score:
                            max_lddt_score = new_avg_lddt
                            max_index = j

                print(f'#########Iteration: {num_iteration + 1}#############')
                print(f"plddt before: {ref_avg_lddt}")
                print(f"plddt after: {max_lddt_score}")
                if max_lddt_score > ref_avg_lddt:
                    print("Continue to refine")
                    current_ref_dir = out_model_dir
                    ref_start_pdb = f"ranked_{max_index}.pdb"
                    ref_start_pkl = f"result_model_{max_index + 1}_multimer.pkl"
                    print('##################################################')
                else:
                    # keep the models in iteration 1 even through the plddt score decreases
                    if num_iteration == 0:
                        ref_avg_lddt = 0
                        with open(out_model_dir + '/' + ref_start_pkl, 'rb') as f:
                            prediction_result = pickle.load(f)
                            ref_avg_lddt = np.mean(prediction_result['plddt'])
                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore = _cal_tmscore(self.params['mmalign_program'], out_model_dir + '/' + ref_start_pdb, native_pdb)
                        model_iteration_scores += [ref_avg_lddt]
                        model_iteration_tmscores += [ref_tmscore]
                    break

            # model_iteration_scores += [max_lddt_score]

            if len(model_iteration_scores) > 0:
                iteration_result_all['targetname'] += [targetname]
                iteration_result_all['model'] += [i]
                iteration_result_all['start_lddt'] += [model_iteration_scores[0]]
                iteration_result_all['end_lddt'] += [model_iteration_scores[len(model_iteration_scores) - 1]]
                iteration_result_all['start_tmscore'] += [model_iteration_tmscores[0]]
                iteration_result_all['end_tmscore'] += [model_iteration_tmscores[len(model_iteration_tmscores) - 1]]

            while len(model_iteration_scores) <= self.max_iteration:
                model_iteration_scores += [0]

            while len(model_iteration_tmscores) <= self.max_iteration:
                model_iteration_tmscores += [0]

            iteration_scores[f'model{i + 1}'] = model_iteration_scores
            true_tm_scores[f'model{i + 1}'] = model_iteration_tmscores

        iteration_result_avg['start_lddt'] = [np.mean(np.array(iteration_result_all['start_lddt']))]
        iteration_result_avg['end_lddt'] = [np.mean(np.array(iteration_result_all['end_lddt']))]
        iteration_result_avg['start_tmscore'] = [np.mean(np.array(iteration_result_all['start_tmscore']))]
        iteration_result_avg['end_tmscore'] = [np.mean(np.array(iteration_result_all['end_tmscore']))]

        print(iteration_scores)
        df = pd.DataFrame(iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        df = pd.DataFrame(true_tm_scores)
        df.to_csv(outdir + '/tmscores.csv')

        df = pd.DataFrame(iteration_result_avg)
        df.to_csv(outdir + '/iteration_result_avg.csv')

        df = pd.DataFrame(iteration_result_all)
        df.to_csv(outdir + '/iteration_result_all.csv')
        os.chdir(cwd)

        return iteration_result_all, iteration_result_avg
