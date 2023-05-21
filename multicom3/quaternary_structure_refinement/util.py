import copy
import os
import sys
import time, json
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
import dataclasses
import numpy as np
from multicom3.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit
from multicom3.complex_templates_search.parsers import TemplateHit
from multicom3.monomer_structure_refinement.util import build_alignment_indices, PrefilterError

# need to add A if using relaxation in alphafold
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

PDB_CHAIN_IDS_UNRELAX = 'BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


class FastaChain:
    def __init__(self, sequence, description):
        self.sequence = sequence
        self.description = description


def combine_pdb(pdbs, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for chain_id, pdb in zip(PDB_CHAIN_IDS_UNRELAX, pdbs):
            for line in open(pdb):
                if not line.startswith('ATOM'):
                    continue
                out.write(line[:21] + chain_id + line[22:])


def parse_fasta(fasta_file):
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


def combine_a3ms(infiles, outfile):
    targetname = None
    targetseq = None
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
                if targetname is None:
                    targetname = line
            else:
                seqs += [line]
                if targetseq is None:
                    targetseq = line

    with open(outfile, 'w') as fw:
        fw.write(f"{targetname}\n{targetseq}\n")
        for (desc, seq) in zip(descriptions, seqs):
            if desc == targetname and seq == targetseq:
                continue
            fw.write(f"{desc}\n{seq}\n")


def convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def cal_tmscore(mmalign_program, inpdb, nativepdb):
    cmd = mmalign_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $2}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return tmscore


def cal_tmalign(tmalign_program, inpdb, nativepdb, tmpdir):
    if os.path.exists(tmpdir):
        os.system(f"rm -rf {tmpdir}")

    os.makedirs(tmpdir)

    src_pdb = tmpdir + '/src.pdb'
    native_pdb = tmpdir + '/native.pdb'

    with open(src_pdb, 'w') as fw:
        for line in open(inpdb):
            if not line.startswith('ATOM'):
                continue
            fw.write(line[:21] + 'A' + line[22:])

    with open(native_pdb, 'w') as fw:
        for line in open(nativepdb):
            if not line.startswith('ATOM'):
                continue
            fw.write(line[:21] + 'A' + line[22:])

    cmd = tmalign_program + ' ' + src_pdb + ' ' + native_pdb + " | grep TM-score | awk '{print $2}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    # print(tmscore_contents)
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    return tmscore


def split_pdb(complex_pdb, outdir):
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


def split_pdb_unrelax2relax(complex_pdb, outdir):
    makedir_if_not_exists(outdir)
    chain_pdbs = {}
    pre_chain = None
    chain_idx = 0
    fw = None
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line)
            chain_pdbs[PDB_CHAIN_IDS[chain_idx]] = outdir + '/' + chain_name + '.pdb'
        elif chain_name == pre_chain:
            fw.write(line)
        else:
            fw.close()
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line)
            chain_idx += 1
            chain_pdbs[PDB_CHAIN_IDS[chain_idx]] = outdir + '/' + chain_name + '.pdb'
            pre_chain = chain_name
    fw.close()
    return chain_pdbs


def create_template_df(templates):
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
        aln_len = int(templates.loc[i, 'alnlen'])
        if target.find('.atom.gz') > 0:
            pdbcode = target[0:4]
        else:
            pdbcode = target
        row_dict = dict(template=target.split()[0],
                        tpdbcode=pdbcode,
                        aln_temp=taln,
                        tstart=tstart,
                        tend=tend,
                        aln_query=qaln,
                        qstart=qstart,
                        qend=qend,
                        evalue=evalue,
                        aligned_length=aln_len)
        row_list += [row_dict]
    if len(row_list) == 0:
        return pd.DataFrame(columns=['template', 'tpdbcode', 'aln_temp', 'tstart', 'tend',
                                     'aln_query', 'qstart', 'qend', 'evalue', 'aligned_length'])
    return pd.DataFrame(row_list)


def create_template_df_with_index(templates):
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
        aln_len = int(templates.loc[i, 'alnlen'])
        if target.find('.atom.gz') > 0:
            pdbcode = target[0:4]
        else:
            pdbcode = target
        row_dict = dict(index=i,
                        template=target.split()[0],
                        tpdbcode=pdbcode,
                        aln_temp=taln,
                        tstart=tstart,
                        tend=tend,
                        aln_query=qaln,
                        qstart=qstart,
                        qend=qend,
                        evalue=evalue,
                        aligned_length=aln_len)
        row_list += [row_dict]
    if len(row_list) == 0:
        return pd.DataFrame(columns=['index', 'template', 'tpdbcode', 'aln_temp', 'tstart', 'tend',
                                     'aln_query', 'qstart', 'qend', 'evalue', 'aligned_length'])
    return pd.DataFrame(row_list)

def get_sequence(inpdb):
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


def get_indices(sequence, start):
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


def build_query_to_hit_index_mapping(
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


def check_template_overlap_regions(template_info1, template_info2, chain_id_map, gap):
    chain1_id = template_info1['chainid']
    chain1_seq = chain_id_map[chain1_id].sequence
    indices_hit = build_alignment_indices(template_info1['aln_temp'], template_info1['tstart'])
    indices_query = build_alignment_indices(template_info1['aln_query'], template_info1['qstart'])
    mapping = build_query_to_hit_index_mapping(
        template_info1['aln_query'], template_info1['aln_temp'], indices_hit, indices_query,
        chain1_seq)

    tstart_1 = np.min(np.array(list(mapping.values())))
    tend_1 = np.max(np.array(list(mapping.values())))

    chain2_id = template_info2['chainid']
    chain2_seq = chain_id_map[chain2_id].sequence
    indices_hit = build_alignment_indices(template_info2['aln_temp'], template_info2['tstart'])
    indices_query = build_alignment_indices(template_info2['aln_query'], template_info2['qstart'])
    mapping = build_query_to_hit_index_mapping(
        template_info2['aln_query'], template_info2['aln_temp'], indices_hit, indices_query,
        chain2_seq)

    tstart_2 = np.min(np.array(list(mapping.values())))
    tend_2 = np.max(np.array(list(mapping.values())))
    # print(tend_1)
    # print(gap)
    # print(tstart_2)

    return max(tstart_1, tstart_2) + gap < min(tend_1, tend_2)


def assess_complex_templates(chain_id_map, template_infos, gap=0):
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
                    if check_template_overlap_regions(template_info1=same_template_infos[template][i],
                                                      template_info2=same_template_infos[template][j],
                                                      chain_id_map=chain_id_map,
                                                      gap=gap):
                        return False
    return True


def cal_sequence_identity_mmseq2(mmseq_program, fasta1, fasta2, tmpdir):
    cmd = f"{mmseq_program} easy-search {fasta1} {fasta2} {tmpdir}/alnRes_{chain1}_{chain2} {tmpdir} " \
          + '--format-output "query,target,pident,qstart,qend,qaln,tstart,tend,taln" --alignment-mode 3 >/dev/null 2>&1'
    os.system(cmd)
    contents = open(f"{tmpdir}/alnRes_{chain1}_{chain2}").readlines()
    seqid = 0
    if len(contents) > 0:
        line = contents[0].rstrip('\n')
        aln_name1, aln_name2, seqid, qstart, qend, aln1, tstart, tend, aln2 = line.split('\t')
    return seqid


def assess_complex_templates_homo(chain_id_map, template_infos, mmseq_program, tmpdir, gap=0, seqid=0.95):
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
                    if check_template_overlap_regions(template_info1=same_template_infos[template][i],
                                                      template_info2=same_template_infos[template][j],
                                                      chain_id_map=chain_id_map,
                                                      gap=gap):
                        return False
                with open(f"{tmpdir}/{i}_{template}.fasta", 'w') as fw:
                    fw.write(f">{i}_{template}\n{same_template_infos[template][i]['aln_temp'].replace('-', '')}")

                with open(f"{tmpdir}/{j}_{template}.fasta", 'w') as fw:
                    fw.write(f">{j}_{template}\n{same_template_infos[template][j]['aln_temp'].replace('-', '')}")

                if cal_sequence_identity_mmseq2(mmseq_program, f"{tmpdir}/{i}_{template}.fasta"
                                                f"{tmpdir}/{j}_{template}.fasta", tmpdir) < seqid:
                    return False
    return True


def assess_complex_templates_diff(chain_id_map, template_infos):
    template_names = []
    same_template_infos = {}
    for template_info in template_infos:
        if template_info['template'] not in template_names:
            template_names += [template_info['template']]
            same_template_infos[template_info['template']] = [template_info]
        else:
            same_template_infos[template_info['template']] += [template_info]
    return len(template_names) == len(chain_id_map)


def check_and_rank_monomer_templates_local_or_global(template_result, outfile, query_sequence, max_template_count):
    global_alignment = False

    templates = template_result['local_alignment']
    if len(templates) == 0:
        templates = template_result['global_alignment']
        global_alignment = True

    sort_indices = []
    for i in range(len(templates)):
        target = templates.loc[i, 'target']
        evalue = float(templates.loc[i, 'evalue'])
        if target.find('.atom.gz') > 0:
            if global_alignment and evalue >= 0.8:
                sort_indices += [i]
            if not global_alignment and evalue < 1e-10:
                sort_indices += [i]
    for i in range(len(templates)):
        if i in sort_indices:
            continue
        sort_indices += [i]

    keep_indices = []
    for i in sort_indices:
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
            assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
        except PrefilterError as e:
            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
            print(msg)
            continue
        keep_indices += [i]
        if len(keep_indices) > max_template_count:
            break

    templates_sorted = copy.deepcopy(templates.iloc[keep_indices])
    templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
    templates_sorted.reset_index(inplace=True, drop=True)
    templates_sorted.to_csv(outfile, sep='\t')
    return True


def check_and_rank_monomer_templates_local_and_global(template_result, outfile, query_sequence, max_template_count):
    evalue_keep_indices = []
    for i in range(len(template_result['local_alignment'])):
        hit = TemplateHit(index=i,
                          name=template_result['local_alignment'].loc[i, 'target'].split('.')[0],
                          aligned_cols=int(template_result['local_alignment'].loc[i, 'alnlen']),
                          query=template_result['local_alignment'].loc[i, 'qaln'],
                          hit_sequence=template_result['local_alignment'].loc[i, 'taln'],
                          indices_query=build_alignment_indices(template_result['local_alignment'].loc[i, 'qaln'],
                                                                template_result['local_alignment'].loc[
                                                                    i, 'qstart']),
                          indices_hit=build_alignment_indices(template_result['local_alignment'].loc[i, 'taln'],
                                                              template_result['local_alignment'].loc[i, 'tstart']),
                          sum_probs=0.0)
        try:
            assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
        except PrefilterError as e:
            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
            print(msg)
            continue
        evalue_keep_indices += [i]

    tmscore_keep_indices = []
    for i in range(len(template_result['global_alignment'])):
        hit = TemplateHit(index=i,
                          name=template_result['global_alignment'].loc[i, 'target'].split('.')[0],
                          aligned_cols=int(template_result['global_alignment'].loc[i, 'alnlen']),
                          query=template_result['global_alignment'].loc[i, 'qaln'],
                          hit_sequence=template_result['global_alignment'].loc[i, 'taln'],
                          indices_query=build_alignment_indices(template_result['global_alignment'].loc[i, 'qaln'],
                                                                template_result['global_alignment'].loc[
                                                                    i, 'qstart']),
                          indices_hit=build_alignment_indices(template_result['global_alignment'].loc[i, 'taln'],
                                                              template_result['global_alignment'].loc[i, 'tstart']),
                          sum_probs=0.0)
        try:
            assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
        except PrefilterError as e:
            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
            print(msg)
            continue
        tmscore_keep_indices += [i]

    if len(evalue_keep_indices) == 0 and len(tmscore_keep_indices) == 0:
        return False

    evalue_thresholds = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3]

    evalue_af_indices = []
    evalue_pdb_indices = []
    tmscore_af_indices = []
    tmscore_pdb_indices = []
    for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
        evalue_af_indices = []
        evalue_pdb_indices = []
        for i in evalue_keep_indices:
            target = template_result['local_alignment'].loc[i, 'target']
            evalue = float(template_result['local_alignment'].loc[i, 'evalue'])
            if evalue < evalue_threshold:
                if target.find('.atom.gz') > 0:
                    evalue_pdb_indices += [i]
                else:
                    evalue_af_indices += [i]

        tmscore_af_indices = []
        tmscore_pdb_indices = []
        for i in tmscore_keep_indices:
            target = template_result['global_alignment'].loc[i, 'target']
            evalue = float(template_result['global_alignment'].loc[i, 'evalue'])
            if evalue > tmscore_threshold:
                if target.find('.atom.gz') > 0:
                    tmscore_pdb_indices += [i]
                else:
                    tmscore_af_indices += [i]

        if len(evalue_af_indices) + len(evalue_pdb_indices) \
                + len(tmscore_af_indices) + len(tmscore_pdb_indices) >= max_template_count:
            break

    templates_sorted = copy.deepcopy(template_result['local_alignment'].iloc[evalue_pdb_indices])
    templates_sorted = templates_sorted.append(
        copy.deepcopy(template_result['global_alignment'].iloc[tmscore_pdb_indices]))
    templates_sorted = templates_sorted.append(
        copy.deepcopy(template_result['local_alignment'].iloc[evalue_af_indices]))
    templates_sorted = templates_sorted.append(
        copy.deepcopy(template_result['global_alignment'].iloc[tmscore_af_indices]))

    templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
    templates_sorted.reset_index(inplace=True, drop=True)
    templates_sorted.to_csv(outfile, sep='\t')
    return True


def check_and_rank_complex_templates(chain_id_map, template_results, max_template_count):
    evalue_thresholds = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4]

    complex_templates_df_filtered = None

    for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
        complex_templates_df = None
        for chain_idx, (chain_id, template_result) in enumerate(zip(chain_id_map, template_results)):
            evalue_keep_indices = []
            for i in range(len(template_result['local_alignment'])):
                hit = TemplateHit(index=i,
                                  name=template_result['local_alignment'].loc[i, 'target'].split('.')[0],
                                  aligned_cols=int(template_result['local_alignment'].loc[i, 'alnlen']),
                                  query=template_result['local_alignment'].loc[i, 'qaln'],
                                  hit_sequence=template_result['local_alignment'].loc[i, 'taln'],
                                  indices_query=build_alignment_indices(
                                      template_result['local_alignment'].loc[i, 'qaln'],
                                      template_result['local_alignment'].loc[
                                          i, 'qstart']),
                                  indices_hit=build_alignment_indices(
                                      template_result['local_alignment'].loc[i, 'taln'],
                                      template_result['local_alignment'].loc[i, 'tstart']),
                                  sum_probs=0.0)
                try:
                    assess_hhsearch_hit(hit=hit, query_sequence=chain_id_map[chain_id].sequence)
                except PrefilterError as e:
                    msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                    print(msg)
                    continue
                if template_result['local_alignment'].loc[i, 'evalue'] < evalue_threshold:
                    evalue_keep_indices += [i]

            tmscore_keep_indices = []
            for i in range(len(template_result['global_alignment'])):
                hit = TemplateHit(index=i,
                                  name=template_result['global_alignment'].loc[i, 'target'].split('.')[0],
                                  aligned_cols=int(template_result['global_alignment'].loc[i, 'alnlen']),
                                  query=template_result['global_alignment'].loc[i, 'qaln'],
                                  hit_sequence=template_result['global_alignment'].loc[i, 'taln'],
                                  indices_query=build_alignment_indices(
                                      template_result['global_alignment'].loc[i, 'qaln'],
                                      template_result['global_alignment'].loc[
                                          i, 'qstart']),
                                  indices_hit=build_alignment_indices(
                                      template_result['global_alignment'].loc[i, 'taln'],
                                      template_result['global_alignment'].loc[i, 'tstart']),
                                  sum_probs=0.0)
                try:
                    assess_hhsearch_hit(hit=hit, query_sequence=chain_id_map[chain_id].sequence)
                except PrefilterError as e:
                    msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                    print(msg)
                    continue
                if template_result['global_alignment'].loc[i, 'evalue'] > tmscore_threshold:
                    tmscore_keep_indices += [i]

            templates_filtered = copy.deepcopy(template_result['local_alignment'].iloc[evalue_keep_indices])
            templates_filtered = templates_filtered.append(
                copy.deepcopy(template_result['global_alignment'].iloc[tmscore_keep_indices]))
            templates_filtered.drop(templates_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
            templates_filtered.reset_index(inplace=True, drop=True)

            curr_df = create_template_df(templates_filtered)
            # print(curr_df)
            if complex_templates_df is None:
                complex_templates_df = curr_df
            else:
                complex_templates_df = complex_templates_df.merge(curr_df, how="inner", on='tpdbcode',
                                                                  suffixes=(str(chain_idx), str(chain_idx + 1)))

        keep_indices = []
        for i in range(len(complex_templates_df)):
            template_infos = []
            for j, chain_id in enumerate(chain_id_map):
                template = complex_templates_df.loc[i, f'template{j + 1}']
                qaln = complex_templates_df.loc[i, f'aln_query{j + 1}']
                qstart = int(complex_templates_df.loc[i, f'qstart{j + 1}'])
                qend = int(complex_templates_df.loc[i, f'qend{j + 1}'])
                taln = complex_templates_df.loc[i, f'aln_temp{j + 1}']
                tstart = complex_templates_df.loc[i, f'tstart{j + 1}']
                tend = complex_templates_df.loc[i, f'tend{j + 1}']
                evalue = float(complex_templates_df.loc[i, f'evalue{j + 1}'])
                row_dict = dict(chainid=chain_id,
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

            if not assess_complex_templates(chain_id_map, template_infos):
                continue

            keep_indices += [i]

        complex_templates_df_filtered = copy.deepcopy(complex_templates_df.iloc[keep_indices])
        complex_templates_df_filtered.drop(complex_templates_df_filtered.filter(regex="Unnamed"), axis=1,
                                           inplace=True)
        complex_templates_df_filtered.reset_index(inplace=True, drop=True)

        if len(complex_templates_df_filtered) > max_template_count:
            break

    return complex_templates_df_filtered
