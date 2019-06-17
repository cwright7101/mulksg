#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
from multiprocessing import Pool
import multiprocessing
from subprocess import Popen
import os
import sys
import shutil
import support
import process_cfg
from process_cfg import bool_to_str
from site import addsitedir
from distutils import dir_util
import options_storage

BASE_STAGE = "construction"
READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single", "hq-mate-pairs"]
READS_TYPES_USED_IN_RNA_SEQ = ["paired-end", "single", "trusted-contigs", "untrusted-contigs"]


def prepare_config_mulksg(filename, cfg, log, additional_contigs_fname, K, stage, saves_dir, last_one, execution_home):
    subst_dict = dict()

    subst_dict["K"] = str(K)
    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset)
    subst_dict["output_base"] = process_cfg.process_spaces(cfg.output_dir)
    subst_dict["tmp_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    if additional_contigs_fname:
        subst_dict["additional_contigs"] = process_cfg.process_spaces(additional_contigs_fname)
        subst_dict["use_additional_contigs"] = bool_to_str(True)
    else:
        subst_dict["use_additional_contigs"] = bool_to_str(False)
    subst_dict["main_iteration"] = bool_to_str(last_one)
    subst_dict["entry_point"] = stage
    subst_dict["load_from"] = saves_dir
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["gap_closer_enable"] = bool_to_str(last_one or K >= 55)
    subst_dict["rr_enable"] = bool_to_str(last_one and cfg.rr_enable)
#    subst_dict["topology_simplif_enabled"] = bool_to_str(last_one)
#TODO: If on multi node, max_threads and max_memory need to change!!!
    if K == cfg.iterative_K[-1] :
        subst_dict["max_threads"] = min(multiprocessing.cpu_count(), cfg.max_threads)
        subst_dict["max_memory"] = cfg.max_memory
    else:
        subst_dict["max_threads"] = min(multiprocessing.cpu_count(), cfg.max_threads) #/(len(cfg.iterative_K)-1)
        subst_dict["max_memory"] = cfg.max_memory/(len(cfg.iterative_K)-1)
    subst_dict["save_gp"] = bool_to_str(cfg.save_gp)
    if (not last_one):
        subst_dict["correct_mismatches"] = bool_to_str(False)
    if "resolving_mode" in cfg.__dict__:
        subst_dict["resolving_mode"] = cfg.resolving_mode
    if "pacbio_mode" in cfg.__dict__:
        subst_dict["pacbio_test_on"] = bool_to_str(cfg.pacbio_mode)
        subst_dict["pacbio_reads"] = process_cfg.process_spaces(cfg.pacbio_reads)
    if cfg.cov_cutoff == "off":
        subst_dict["use_coverage_threshold"] = bool_to_str(False)
    else:
        subst_dict["use_coverage_threshold"] = bool_to_str(True)
        if cfg.cov_cutoff == "auto":
            subst_dict["coverage_threshold"] = 0.0
        else:
            subst_dict["coverage_threshold"] = cfg.cov_cutoff
    if cfg.lcer_cutoff is not None:
        subst_dict["lcer_enabled"] = bool_to_str(True)
        subst_dict["lcer_coverage_threshold"] = cfg.lcer_cutoff

    if "series_analysis" in cfg.__dict__:
        subst_dict["series_analysis"] = cfg.series_analysis
    process_cfg.substitute_params(filename, subst_dict, log)


def prepare_config_rnamulksg(filename, log):
    if not options_storage.rna:
        return
    subst_dict = dict()
    subst_dict["ss_enabled"] = bool_to_str(options_storage.strand_specificity is not None)
    subst_dict["antisense"] = bool_to_str(options_storage.strand_specificity == 'rf')
    process_cfg.substitute_params(filename, subst_dict, log)


def prepare_config_construction(filename, log):
    if options_storage.read_cov_threshold is None:
        return
    subst_dict = dict()
    subst_dict["read_cov_threshold"] = options_storage.read_cov_threshold
    process_cfg.substitute_params(filename, subst_dict, log)


def get_read_length(output_dir, K, ext_python_modules_home, log):
    est_params_filename = os.path.join(output_dir, "K%d" % K, "final.lib_data")
    max_read_length = 0
    if os.path.isfile(est_params_filename):
        addsitedir(ext_python_modules_home)
        if sys.version.startswith('2.'):
            import pyyaml2 as pyyaml
        elif sys.version.startswith('3.'):
            import pyyaml3 as pyyaml
        est_params_data = pyyaml.load(open(est_params_filename, 'r'))
        max_read_length = int(est_params_data['nomerge max read length'])
        log.info("Max read length detected as %d" % max_read_length)
    if max_read_length == 0:
        support.error("Failed to estimate maximum read length! File with estimated params: " + est_params_filename, log)
    return max_read_length


def update_k_mers_in_special_cases(cur_k_mers, RL, log, silent=False):
    if options_storage.auto_K_allowed():
        if RL >= 250:
            if not silent:
                log.info("Default k-mer sizes were set to %s because estimated "
                         "read length (%d) is equal to or greater than 250" % (str(options_storage.K_MERS_250), RL))
            return options_storage.K_MERS_250
        if RL >= 150:
            if not silent:
                log.info("Default k-mer sizes were set to %s because estimated "
                         "read length (%d) is equal to or greater than 150" % (str(options_storage.K_MERS_150), RL))
            return options_storage.K_MERS_150
    if RL <= max(cur_k_mers):
        new_k_mers = [k for k in cur_k_mers if k < RL]
        if not silent:
            log.info("K-mer sizes were set to %s because estimated "
                     "read length (%d) is less than %d" % (str(new_k_mers), RL, max(cur_k_mers)))
        return new_k_mers
    return cur_k_mers


def reveal_original_k_mers(RL):
    if options_storage.original_k_mers is None or options_storage.original_k_mers == 'auto':
        cur_k_mers = options_storage.k_mers
        options_storage.k_mers = options_storage.original_k_mers
        original_k_mers = update_k_mers_in_special_cases(options_storage.K_MERS_SHORT, RL, None, silent=True)
        options_storage.k_mers = cur_k_mers
    else:
        original_k_mers = options_storage.original_k_mers
    original_k_mers = [k for k in original_k_mers if k < RL]
    return original_k_mers


def add_configs(command, configs_dir):
    #Order matters here!
    mode_config_mapping = [("single_cell", "mda_mode"), 
                           ("meta", "meta_mode"),
                           ("truseq_mode", "moleculo_mode"),
                           ("rna", "rna_mode"),
                           ("large_genome", "large_genome_mode"),
                           ("plasmid", "plasmid_mode")]
                           #("careful", "careful_mode"),
    for (mode, config) in mode_config_mapping:
        if options_storage.__dict__[mode]:
            if mode == "rna" or mode == "meta":
                command.append(os.path.join(configs_dir, "mda_mode.info"))
            command.append(os.path.join(configs_dir, config + ".info"))
    if options_storage.__dict__["careful"]:
        if options_storage.__dict__["single_cell"]:
            command.append(os.path.join(configs_dir, "careful_mda_mode.info"))
        else:
            command.append(os.path.join(configs_dir, "careful_mode.info"))

    # special case: extra config
    if options_storage.rna and options_storage.fast:
        command.append(os.path.join(configs_dir, "rna_fast_mode.info"))

def create_mulksg_configs(configs_dir, execution_home, cfg, log, K_values_list, additional_contigs_fname, last_one, ext_python_modules_home):
    command_list = []
    contig_files = []
    for K in K_values_list:
        # RL = get_read_length(cfg.output_dir, K, ext_python_modules_home, log)
        # if K >= RL:
        #     support.warning("Value of K (%d) exceeded estimated read length (%d)" %
        #                             (K, RL), log)
        data_dir = os.path.join(cfg.output_dir, "K%d" % K)
        stage = BASE_STAGE
        saves_dir = os.path.join(data_dir, 'saves')
        dst_configs = os.path.join(data_dir, "configs")

        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)
        os.makedirs(data_dir)

        dir_util._path_created = {}  # see http://stackoverflow.com/questions/9160227/dir-util-copy-tree-fails-after-shutil-rmtree
        dir_util.copy_tree(os.path.join(configs_dir, "debruijn"), dst_configs, preserve_times=False)

        log.info("\n== Running assembler: " + ("K%d" % K) + "\n")
        if "read_buffer_size" in cfg.__dict__:
            #FIXME why here???
            process_cfg.substitute_params(os.path.join(dst_configs, "construction.info"), {"read_buffer_size": cfg.read_buffer_size}, log)
        if "scaffolding_mode" in cfg.__dict__:
            #FIXME why here???
            process_cfg.substitute_params(os.path.join(dst_configs, "pe_params.info"), {"scaffolding_mode": cfg.scaffolding_mode}, log)

        prepare_config_rnamulksg(os.path.join(dst_configs, "rna_mode.info"), log)
        prepare_config_construction(os.path.join(dst_configs, "construction.info"), log)
        cfg_fn = os.path.join(dst_configs, "config.info")
        prepare_config_mulksg(cfg_fn, cfg, log, additional_contigs_fname, K, stage, saves_dir, last_one, execution_home)

        command = [os.path.join(execution_home, "mulksg-core"), cfg_fn]
        # print("command is: ", command)
        add_configs(command, dst_configs)
        command_list.append(command)
        contig_files.append(os.path.join(cfg.output_dir, "K%d" % K, "simplified_contigs.fasta"))
    return command_list, contig_files

def run_mulksg_iterations(command, log):
    support.sys_call(command, log)

def prepare_config_scaffold_correction(filename, cfg, log, saves_dir, K):
    subst_dict = dict()

    subst_dict["K"] = str(K)
    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset)
    subst_dict["output_base"] = process_cfg.process_spaces(os.path.join(cfg.output_dir, "SCC"))
    subst_dict["tmp_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    subst_dict["use_additional_contigs"] = bool_to_str(False)
    subst_dict["main_iteration"] = bool_to_str(False)
    subst_dict["entry_point"] = BASE_STAGE
    subst_dict["load_from"] = saves_dir
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory

    #todo
    process_cfg.substitute_params(filename, subst_dict, log)


def run_scaffold_correction(configs_dir, execution_home, cfg, log, latest, K):
    data_dir = os.path.join(cfg.output_dir, "SCC", "K%d" % K)
    saves_dir = os.path.join(data_dir, 'saves')
    dst_configs = os.path.join(data_dir, "configs")
    cfg_file_name = os.path.join(dst_configs, "config.info")

    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)
    os.makedirs(data_dir)

    dir_util.copy_tree(os.path.join(configs_dir, "debruijn"), dst_configs, preserve_times=False)

    log.info("\n== Running scaffold correction \n")
    scaffolds_file = os.path.join(latest, "scaffolds.fasta")
    if not os.path.isfile(scaffolds_file):
        support.error("Scaffodls were not found in " + scaffolds_file, log)
    if "read_buffer_size" in cfg.__dict__:
        construction_cfg_file_name = os.path.join(dst_configs, "construction.info")
        process_cfg.substitute_params(construction_cfg_file_name, {"read_buffer_size": cfg.read_buffer_size}, log)
    process_cfg.substitute_params(os.path.join(dst_configs, "moleculo_mode.info"), {"scaffolds_file": scaffolds_file}, log)
    prepare_config_scaffold_correction(cfg_file_name, cfg, log, saves_dir, K)
    command = [os.path.join(execution_home, "spades-truseq-scfcorrection"), cfg_file_name]
    add_configs(command, dst_configs)
    log.info(str(command))
    support.sys_call(command, log)

def run_mulksg(configs_dir, execution_home, cfg, dataset_data, ext_python_modules_home, log, comm, rank, size):
    #Running this on multiple nodes, need to differentiate the actual work here
    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    while (len(cfg.iterative_K) - 1) % size != 0:
        bgap_idx = 0
        bgap = 0
        for i in range(0,len(cfg.iterative_K) - 1):
            gap = cfg.iterative_K[i+1] - cfg.iterative_K[i]
            if  gap > bgap:
                bgap = gap
                bgap_idx = i
        if bgap == 0:
            break;
        ul = cfg.iterative_K[i+1]
        ll = cfg.iterative_K[i]
        new_k = (int)(ul-ll)/2
        new_k = new_k/2*2 + 1
        cfg.iterative_K.append(new_k)
        cfg.iterative_K = sorted(cfg.iterative_K)
        #add to cfg.iterative_K
        #add a number between the biggest gaps in the current cfg.iterative_K

    if rank == 0: #Setup all the configs, we will put them in a work_list
        # checking and removing conflicting K-mer directories
        if options_storage.restart_from and (options_storage.restart_k_mers != options_storage.original_k_mers):
            processed_K = []
            for k in range(options_storage.MIN_K, options_storage.MAX_K, 2):
                cur_K_dir = os.path.join(cfg.output_dir, "K%d" % k)
                if os.path.isdir(cur_K_dir) and os.path.isfile(os.path.join(cur_K_dir, "final_contigs.fasta")):
                    processed_K.append(k)
            if processed_K:
                RL = get_read_length(cfg.output_dir, processed_K[0], ext_python_modules_home, log)
                needed_K = update_k_mers_in_special_cases(cfg.iterative_K, RL, log, silent=True)
                needed_K = [k for k in needed_K if k < RL]
                original_K = reveal_original_k_mers(RL)

                k_to_delete = []
                for id, k in enumerate(needed_K):
                    if len(processed_K) == id:
                        if processed_K[-1] == original_K[-1]:  # the last K in the original run was processed in "last_one" mode
                            k_to_delete = [original_K[-1]]
                        break
                    if processed_K[id] != k:
                        k_to_delete = processed_K[id:]
                        break
                if not k_to_delete and (len(processed_K) > len(needed_K)):
                    k_to_delete = processed_K[len(needed_K) - 1:]
                if k_to_delete:
                    log.info("Restart mode: removing previously processed directories for K=%s "
                             "to avoid conflicts with K specified with --restart-from" % (str(k_to_delete)))
                    for k in k_to_delete:
                        shutil.rmtree(os.path.join(cfg.output_dir, "K%d" % k))
        bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")
        if os.path.isdir(bin_reads_dir):
            shutil.rmtree(bin_reads_dir)
        cfg.tmp_dir = support.get_tmp_dir(prefix="mulksg_")

        first_kValues = cfg.iterative_K
        last_kValue = []
        last_kValue.append(first_kValues[-1])
        del first_kValues[-1]
        #process all of the files as last_one
        first_commands, contig_files = create_mulksg_configs(configs_dir, execution_home, cfg, log, first_kValues, None, True,ext_python_modules_home)
        work_list = first_commands
        if size > 1 and len(cfg.iterative_K) > 1:
            num_running = 0
            for command in work_list:
                if num_running == (size - 1):
                    support.sys_call(command, log)
                    num_running = 0
                else:
                    node_rank = comm.recv()#source=ANY_SOURCE
                    comm.send(command, node_rank)
                    num_running += 1

            for i in range(1,size):
                node_rank = comm.recv()#source=ANY_SOURCE
                comm.send("Done", node_rank)
                print "Rank", node_rank, "Finished"
        else:
            if len(cfg.iterative_K) > 1:
                pool = Pool(len(first_commands))
                results = pool.map(support.sys_call, first_commands)#(command for command in first_commands))
        #Now we have generated all the contig files, concatenate into one file     
        extra_contig_filename = None
        try:
            extra_contig_filename = cfg.output_dir +"/extracontigs.fa"   #"/K" + last_kValue[0]   
            with open(extra_contig_filename, 'w') as contig_file:
                for fname in contig_files:
                    try:
                        with open(fname) as infile:
                            contig_file.write(infile.read())
                    except:
                        continue
            
        except:
            log.info("Could not create extra contig file!!!!")
            extra_contig_filename = None

        last_command, contig_files = create_mulksg_configs(configs_dir, execution_home, cfg, log, last_kValue, extra_contig_filename, True,ext_python_modules_home)
        last_command = last_command[0]
        support.sys_call(last_command, log)
        K = last_kValue[0]
        latest = os.path.join(cfg.output_dir, "K%d" % K)

        finished_on_stop_after = False

        if cfg.correct_scaffolds and not options_storage.run_completed:
            if options_storage.continue_mode and os.path.isfile(os.path.join(cfg.output_dir, "SCC", "corrected_scaffolds.fasta")) and not options_storage.restart_from == "scc":
                log.info("\n===== Skipping %s (already processed). \n" % "scaffold correction")
            else:
                if options_storage.continue_mode:
                    support.continue_from_here(log)
                run_scaffold_correction(configs_dir, execution_home, cfg, log, latest, 21)
            latest = os.path.join(os.path.join(cfg.output_dir, "SCC"), "K21")
            if options_storage.stop_after == 'scc':
                support.finish_here(log)

        if cfg.correct_scaffolds:
            correct_scaffolds_fpath = os.path.join(latest, "corrected_scaffolds.fasta")
            if os.path.isfile(correct_scaffolds_fpath):
                shutil.copyfile(correct_scaffolds_fpath, cfg.result_scaffolds)
        elif not finished_on_stop_after:  # interupted by --stop-after, so final K is not processed!
            if os.path.isfile(os.path.join(latest, "before_rr.fasta")):
                result_before_rr_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "before_rr.fasta")
                if not os.path.isfile(result_before_rr_contigs) or not options_storage.continue_mode:
                    shutil.copyfile(os.path.join(latest, "before_rr.fasta"), result_before_rr_contigs)
            if options_storage.rna:
                if os.path.isfile(os.path.join(latest, "transcripts.fasta")):
                    if not os.path.isfile(cfg.result_transcripts) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "transcripts.fasta"), cfg.result_transcripts)
                if os.path.isfile(os.path.join(latest, "transcripts.paths")):
                    if not os.path.isfile(cfg.result_transcripts_paths) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "transcripts.paths"), cfg.result_transcripts_paths)
                for filtering_type in options_storage.filtering_types:
                    prefix = filtering_type + "_filtered_"
                    result_filtered_transcripts = os.path.join(cfg.output_dir, prefix + options_storage.transcripts_name)
                    latest_filtered_transcripts = os.path.join(latest, prefix + "final_paths.fasta")
                    if os.path.isfile(latest_filtered_transcripts):
                        if not os.path.isfile(result_filtered_transcripts) or not options_storage.continue_mode:
                            shutil.copyfile(latest_filtered_transcripts, result_filtered_transcripts)
            else:
                if os.path.isfile(os.path.join(latest, "final_contigs.fasta")):
                    if not os.path.isfile(cfg.result_contigs) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "final_contigs.fasta"), cfg.result_contigs)
                if os.path.isfile(os.path.join(latest, "first_pe_contigs.fasta")):
                    result_first_pe_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "first_pe_contigs.fasta")
                    if not os.path.isfile(result_first_pe_contigs) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "first_pe_contigs.fasta"), result_first_pe_contigs)
                if cfg.rr_enable:
                    if os.path.isfile(os.path.join(latest, "scaffolds.fasta")):
                        if not os.path.isfile(cfg.result_scaffolds) or not options_storage.continue_mode:
                            shutil.copyfile(os.path.join(latest, "scaffolds.fasta"), cfg.result_scaffolds)
                    if os.path.isfile(os.path.join(latest, "scaffolds.paths")):
                        if not os.path.isfile(cfg.result_scaffolds_paths) or not options_storage.continue_mode:
                            shutil.copyfile(os.path.join(latest, "scaffolds.paths"), cfg.result_scaffolds_paths)
                if os.path.isfile(os.path.join(latest, "assembly_graph_with_scaffolds.gfa")):
                    if not os.path.isfile(cfg.result_graph_gfa) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "assembly_graph_with_scaffolds.gfa"), cfg.result_graph_gfa)
                if os.path.isfile(os.path.join(latest, "assembly_graph.fastg")):
                    if not os.path.isfile(cfg.result_graph) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "assembly_graph.fastg"), cfg.result_graph)
                if os.path.isfile(os.path.join(latest, "final_contigs.paths")):
                    if not os.path.isfile(cfg.result_contigs_paths) or not options_storage.continue_mode:
                        shutil.copyfile(os.path.join(latest, "final_contigs.paths"), cfg.result_contigs_paths)

        if os.path.isdir(bin_reads_dir):
            shutil.rmtree(bin_reads_dir)
        if os.path.isdir(cfg.tmp_dir):
            shutil.rmtree(cfg.tmp_dir)

        return first_kValues + last_kValue

    elif rank != 0 and len(cfg.iterative_K) > 1:
        comm.send(rank, 0)
        command = comm.recv(source = 0)
        while command != "Done":
            # print rank, command
            support.sys_call(command, log)
            comm.send(rank, 0)
            command = comm.recv(source = 0)
        # comm.send(rank, 0) #tell the master I am done and exiting
        return None