nextflow.enable.dsl = 2

params.test = false
params.resdir = "res"

params.minmaf = 0.01
// assuming homozyogous diploids, eg. min_maf 0.01
params.mincm = 2.0
params.nchroms = 14

params.tpbwt_template_opts = 1
params.tpbwt_Lm = 80
// update on 1/10 from 100 -> 80
params.tpbwt_Lf = params.mincm
params.tpbwt_use_phase_correction = 0

params.hapibd_minoutput = params.mincm
params.hapibd_minseed = params.mincm
params.hapibd_minextend = 1.0
params.hapibd_maxgap = 1000
params.hapibd_minmarkers = 70

params.refinedibd_length = params.mincm
params.refinedibd_lod = 1.6
// 1/9/24, update from 1.1 to 1.6
params.refinedibd_scale = 0
params.refinedibd_window = 40.0
params.refinedibd_trim = 0.15
// in centigmorgans
params.refinedibd_optimize_params_json = ""

params.hmmibd_n = 100
params.hmmibd_m = 5
params.hmmibd_version = ["hmmIBD", "hmmibd2"][0]

params.isorelate_imiss = 0.3
params.isorelate_vmiss = 0.3
params.isorelate_min_snp = 20
// optimized
params.isorelate_minmaf = 0.1
// 0.1, which is different from other callers

// params.tsinferibd_max_tmrca = [1000, 3000]

params.filt_ibd_by_ov = true
// update 1/9/24 false->true
params.peak_validate_meth = 'ihs'
// 'xirs' or 'ihs'
params.ibdne_no_diploid_convertion = "true"
params.ibdne_mincm = params.mincm
params.ibdne_minregion = 10
params.ibdne_flatmeth = ["none", "keep_hap_1_only", "merge"][0]

params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrials = 1000
// params.ifm_mincm = 2.0
params.ifm_mingwcm = 5.0

params.csp_json = ""
// a normal string for single json; a string with wildcard for multiple json files
params.sp_sets_json = ""
// a normal string for single json
params.mp_sets_json = ""
// a normal string for single json


params.meta = ""
// empty for simulation

params.compute_bench_large_size = [0, 1][0]


params.sp_defaults = [
    seqlen: 100 * 15000,
    selpos: Math.round(0.33 * 100 * 15000),
    num_origins: 1,
    N: params.test ? 500 : 10000,
    h: 0.5,
    s: 0.3,
    g_sel_start: 80,
    r: 0.01 / 15000,
    sim_relatedness: 0,
    g_ne_change_start: 200,
    N0: params.test ? 500 : 1000,
    u: 1E-8,
    nsam: params.test ? 100 : 1000,
]

params.mp_defaults = [
    seqlen: 100 * 15000,
    selpos: Math.round(100 * 15000 * 0.33),
    num_origins: 1,
    N: params.test ? 500 : 10000,
    h: 0.5,
    s: 0.3,
    g_sel_start: 80,
    r: 0.01 / 15000,
    sim_relatedness: 0,
    mig: 0.00001,
    sel_mig: 0.01,
    npop: 5,
    nsam: params.test ? 40 : 200,
    Tsplit: 500,
    u: 1E-8,
]


workflow {
    WF_SP()
    WF_MP()
}



workflow WF_SP {
    main:

    def sp_sets = get_sp_sets()

    // *********************** Log Params *************************
    def sp_set_args_keys = sp_sets.collect { _k, v -> v }[0].collect { k, _v -> k }

    ch_sp_sets = Channel.fromList(sp_sets.collect { label, args -> [label, args] })
    if (params.test) {
        ch_sp_sets = ch_sp_sets.take(1)
    }

    Channel
        .value(['label'] + sp_set_args_keys)
        .concat(
            ch_sp_sets.map { k, args ->
                [k] + sp_set_args_keys.collect { i -> args[i] }
            }
        )
        .map { lst -> lst.join("\t") }
        .collectFile(name: "sp_sets.tsv", newLine: true, storeDir: params.resdir, sort: false)

    // ***********************Simulation *************************

    ch_chrs = Channel.fromList(1..params.nchroms)

    ch_sp_input = ch_sp_sets
        .combine(ch_chrs)
        .map { label, args, chrno -> [label, chrno, args] }

    SIM_SP_CHR(ch_sp_input)


    // *********************** Call ibd *************************
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
    // label, chrno, trees, vcf
    // ch_true_ne = SIM_SP_CHR.out.true_ne
    //     .filter { label, chrno,true_ne-> chrno == 1}
    //     .map{label, chrno,true_ne-> [label, true_ne]}

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets, by: 0)
        .map { label, chrno, trees, _vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets, by: 0)
        .map { label, chrno, _trees, vcf, args -> [label, chrno, args, vcf] }

    // ch_in_ibdcall_vcf_with_ne = ch_in_ibdcall_vcf
    //     .combine(ch_true_ne, by: 0)
    //     // label, chrno, args, vcf, true_ne

    // CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)

    // collect IBD and groupby simulation label and ibdcaller
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.ibd
        .map { it + ["hapibd"] }
        .concat(CALL_IBD_TSKIBD.out.ibd.map { it + ["tskibd"] })
        .concat(CALL_IBD_REFINEDIBD.out.ibd.map { it + ["refinedibd"] })
        .concat(CALL_IBD_TPBWT.out.ibd.map { it + ["tpbwt"] })
        .concat(CALL_IBD_HMMIBD.out.ibd.map { it + ["hmmibd"] })
        .concat(CALL_IBD_ISORELATE.out.ibd.map { it + ["isorelate"] })
        .map { label, chrno, ibd, ibdcaller ->
            def key = groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]
        }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect { data -> data[1] }
            return ["${label}", ibdcaller, ibd_lst]
        }
    // ch_out_ibd_grp.map{it-> [it[0], it[1], it[2].size()]}.view()

    ch_true_ibd = CALL_IBD_TSKIBD.out.ibd
        .map { label, chrno, ibd -> [groupKey(label, params.nchroms), [chrno, ibd]] }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { key, data_lst ->
            def label = key.toString()
            def ibd_files_true = data_lst.collect { v -> v[1] }
            [label, ibd_files_true]
        }


    ch_out_ibd_grp = ch_out_ibd_grp
        .combine(ch_true_ibd, by: 0)
        .map { label, ibdcaller, ibd_lst, ibd_list_true ->
            // fix name collision for tskibd as true and inferrre ibd are the same
            [label, ibdcaller, ibd_lst, ibdcaller == "tskibd" ? [] : ibd_list_true]
        }


    ch_out_vcf_grp = SIM_SP_CHR.out.trees_vcf
        .map { label, chrno, _trees, vcf ->
            [groupKey(label, params.nchroms), [chrno, vcf]]
        }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { label, value_lst ->
            def vcf_lst = value_lst.collect { values -> values[1] }
            [label.toString(), vcf_lst]
        }

    ch_grouped_ibd_vcf = ch_out_ibd_grp
        .combine(ch_out_vcf_grp, by: 0)
        .combine(
            ch_sp_sets.map { label, args ->
                def genome_set_id = args.genome_set_id
                [label, genome_set_id]
            },
            by: 0
        )
    // [label, ibdcaller, ibd_lst, ibd_lst_true, vcf_lst, genome_set_id]
    //                                      ^^^^^^^^^^^^^^

    // ********************** Process ibd ***************************
    // NOTE: although the ibd_files_true is passed, but it can still not used 
    // depending the params.filt_ibd_by_ov = True
    PROC_DIST_NE(ch_grouped_ibd_vcf, file("${projectDir}/lib/ibdne.23Apr20.ae9.jar"))
    // out.ne_input_orig
    // out.ne_input_rmpeaks



    // ********************** Run IbdNe ***************************
    ch_in_ibdne = PROC_DIST_NE.out.ne_input_orig
        .map { it + [false] }
        .concat(
            PROC_DIST_NE.out.ne_input_rmpeaks.map { it + [true] }
        )
        .combine(ch_sp_sets, by: 0)
    // add args

    RUN_IBDNE(ch_in_ibdne)


    //  ******************** CMP IBD ********************************

    ch_for_cmp_ibd = ch_out_ibd_grp
        .filter { _label, ibdcaller, _ibd_lst, _ibd_list_true -> ibdcaller != "tskibd" }
        .map { label, ibdcaller, ibd_lst, ibd_list_true ->
            def args = sp_sets[label]
            def csp_id = "no_cspid"
            [label, args, ibdcaller, csp_id, ibd_list_true, ibd_lst]
        }

    CMP_IBD(ch_for_cmp_ibd)

    emit:
    ch_ibdobj_dist = PROC_DIST_NE.out.ibddist_ibd_obj
    ch_ibdobj_ne = PROC_DIST_NE.out.ibdne_ibd_obj
    ch_ibdne = RUN_IBDNE.out
}

workflow WF_MP {
    main:

    def mp_sets = get_mp_sets()
    // *********************** Log params ************************
    def mp_set_args_keys = mp_sets.collect { _k, v -> v }[0].collect { k, _v -> k }

    ch_mp_sets = Channel.fromList(mp_sets.collect { label, args -> [label, args] })
    if (params.test) {
        ch_mp_sets = ch_mp_sets.take(1)
    }

    Channel
        .value(['label'] + mp_set_args_keys)
        .concat(
            ch_mp_sets.map { k, args ->
                [k] + mp_set_args_keys.collect { i -> args[i] }
            }
        )
        .map { lst -> lst.join("\t") }
        .collectFile(name: "mp_sets.tsv", newLine: true, storeDir: params.resdir, sort: false)


    // ***********************Simulation *************************

    ch_chrs = Channel.fromList(1..params.nchroms)

    ch_mp_input = ch_mp_sets
        .combine(ch_chrs)
        .map { label, args, chrno -> [label, chrno, args] }

    SIM_MP_CHR(ch_mp_input)


    // *********************** Call ibd *************************
    ch_trees_vcf = SIM_MP_CHR.out.trees_vcf
    // label, chrno, trees, vcf

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_mp_sets, by: 0)
        .map { label, chrno, trees, _vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_mp_sets, by: 0)
        .map { label, chrno, _trees, vcf, args -> [label, chrno, args, vcf] }

    // ch_in_ibdcall_vcf_with_ne = ch_in_ibdcall_vcf
    //     .map{it + [ [] ] }
    //     // label, chrno, args, vcf, true_ne (empty)
    //     //
    //     //                            ^
    //     //                            |
    //     // Note pass [] to true_ne, this signals multiple-pop simulation,
    //     // a Ne value is caluated from simulation args, see function `mp_ne`
    //     // defined above the process definitions

    // CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)

    // collect IBD and groupby simulation label and ibdcaller
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.ibd
        .map { it + ["hapibd"] }
        .concat(CALL_IBD_TSKIBD.out.ibd.map { it + ["tskibd"] })
        .concat(CALL_IBD_REFINEDIBD.out.ibd.map { it + ["refinedibd"] })
        .concat(CALL_IBD_TPBWT.out.ibd.map { it + ["tpbwt"] })
        .concat(CALL_IBD_HMMIBD.out.ibd.map { it + ["hmmibd"] })
        .concat(CALL_IBD_ISORELATE.out.ibd.map { it + ["isorelate"] })
        .map { label, chrno, ibd, ibdcaller ->
            def key = groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]
        }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect { data -> data[1] }
            return [label, ibdcaller, ibd_lst]
        }

    ch_out_vcf_grp = SIM_MP_CHR.out.trees_vcf
        .map { label, chrno, _trees, vcf ->
            [groupKey(label, params.nchroms), [chrno, vcf]]
        }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { label, value_lst ->
            def vcf_lst = value_lst.collect { values -> values[1] }
            [label.toString(), vcf_lst]
        }

    ch_grouped_ibd_vcf = ch_out_ibd_grp
        .combine(ch_out_vcf_grp, by: 0)
        .combine(
            ch_mp_sets.map { label, args ->
                def genome_set_id = args.genome_set_id
                [label, genome_set_id]
            },
            by: 0
        )
    // [label, ibdcaller, ibd_lst, vcf_lst, genome_set_id]

    // ********************** Process ibd ***************************

    PROC_INFOMAP(ch_grouped_ibd_vcf)
    // out.ifm_orig_ibd_obj
    // out.ifm_rmpeaks_ibd_obj


    // ********************** Run Infomap ***************************
    ch_in_run_infomap = PROC_INFOMAP.out.ifm_orig_ibd_obj
        .map { it + [false] }
        .concat(
            PROC_INFOMAP.out.ifm_rmpeaks_ibd_obj.map { it + [true] }
        )
        .combine(ch_mp_sets, by: 0)

    RUN_INFOMAP(ch_in_run_infomap)


    // ********************** CMP IBD ***************************
    ch_true_ibd = CALL_IBD_TSKIBD.out.ibd
        .map { label, chrno, ibd -> [groupKey(label, params.nchroms), [chrno, ibd]] }
        .groupTuple(by: 0, sort: { a, b -> a[0] <=> b[0] })
        .map { key, data_lst ->
            def label = key.toString()
            def ibd_files_true = data_lst.collect { v -> v[1] }
            [label, ibd_files_true]
        }

    ch_for_cmp_ibd = ch_out_ibd_grp
        .combine(ch_true_ibd, by: 0)
        .filter { _label, ibdcaller, _ibd_lst, _ibd_list_true -> ibdcaller != "tskibd" }
        .map { label, ibdcaller, ibd_lst, ibd_list_true ->
            def args = mp_sets[label]
            def csp_id = "no_cspid"
            [label, args, ibdcaller, csp_id, ibd_list_true, ibd_lst]
        }

    CMP_IBD(ch_for_cmp_ibd)

    emit:
    ch_ifm = RUN_INFOMAP.out
}

workflow PARAM_OPTIMIZATION {
    /////////////////////////////////////////////////////
    // models
    // - read models from json file
    def models = [
        "sp_neu": params.sp_defaults + [s: 0.0, genome_set_id: 10000],
        "mp_neu": params.mp_defaults + [s: 0.0, genome_set_id: 20000],
        "uk_human": [seqlen: 60000000, gc: 0, r: 1E-8, u: 1.00E-8, nsam: 1000, genome_set_id: 60001],
        "uk_human2": [seqlen: 60000000, gc: 0, r: 1E-8, u: 1.38E-8, nsam: 1000, genome_set_id: 70001],
        "uk_pf": [seqlen: 900000, gc: 0, r: 6.6666667E-7, u: 1E-8, nsam: 1000, genome_set_id: 60003],
    ]
    // prepare model input chanel
    ch_input = Channel
        .fromList(models.collect { k, v -> [k, v] })
        .combine(Channel.fromList(1..params.nchroms))
        .map { label, args, chrno -> [label, chrno, args] }
        .branch { label, _chrno, _args ->
            sp: label.startsWith("sp_")
            mp: label.startsWith("mp_")
            uk: label.startsWith("uk_")
        }

    // //////////////////////////////////////////
    // log simulation parameters
    log_sim_params(ch_input.sp, "param_optimization_sp_set.params.tsv")
    log_sim_params(ch_input.mp, "param_optimization_mp_set.params.tsv")
    log_sim_params(ch_input.uk, "param_optimization_uk_set.params.tsv")

    /////////////////////////////////////////////////////
    // simulations and conat vcf/tree results
    SIM_SP_CHR(ch_input.sp)
    SIM_MP_CHR(ch_input.mp)
    SIM_UK_CHR(ch_input.uk)
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf)
        .concat(SIM_UK_CHR.out.trees_vcf)
    // label, chrno, trees, vcf

    ////////////////////////////////////////////////////
    // call IBD with caller-specific parameters
    // - read csp json
    // - verify json data
    // - call ibd using the specified caller and parameters



    // read csp list from json files
    if (params.csp_json != "") {
        ch_csp = Channel
            .fromPath(params.csp_json)
            .flatMap { fn -> read_param_from_json(fn) }
    }
    else {
        // default list for testing purpose
        ch_csp = Channel.from(
            [
                ["hapibd", 3, [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000, minmarkers: 100, minmaf: 0.01]]
            ]
        )
    }


    ch_callibd = ch_trees_vcf
        .map { label, chrno, _trees, vcf -> [label, chrno, models[label], vcf] }
        .combine(ch_csp)
        .branch { label, _chrno, _args, _vcf, caller, _csp_id, _csp_args ->
            hapibd: caller == "hapibd"
            hmmibd: caller == "hmmibd" && (label != "uk_human") && (label != "uk_human2")
            isorelate: caller == "isorelate" && (label != "uk_human") && (label != "uk_human2")
            refinedibd: caller == "refinedibd"
            tpbwt: caller == "tpbwt"
            other: true
        }

    // to avoid issue when resume the pipeline
    println("view one item from the channel to avoid issue when resuming the pipeline")
    ch_callibd.hapibd.take(1).view()
    ch_callibd.hmmibd.take(1).view()
    ch_callibd.isorelate.take(1).view()
    ch_callibd.refinedibd.take(1).view()
    ch_callibd.tpbwt.take(1).view()

    // /////////////////////////////////////////////////
    // Compare IBD inferred vs true
    // inferred ibd
    CALL_IBD_HAPIBD_PARAM(ch_callibd.hapibd)
    CALL_IBD_HMMIBD_PARAM(ch_callibd.hmmibd)
    CALL_IBD_ISORELATE_PARAM(ch_callibd.isorelate)
    CALL_IBD_REFINEDIBD_PARAM(ch_callibd.refinedibd)
    CALL_IBD_TPBWT_PARAM(ch_callibd.tpbwt)

    ///////////////////////////////////////////////////////
    // Log params
    log_ibdcall_csp_params(ch_callibd.hapibd, "param_optimization_hapibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.hmmibd, "param_optimization_hmmibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.isorelate, "param_optimization_isorelate.params.tsv")
    log_ibdcall_csp_params(ch_callibd.refinedibd, "param_optimization_refinedibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.tpbwt, "param_optimization_tpbwt.params.tsv")

    ch_inferred_ibd = CALL_IBD_HAPIBD_PARAM.out
        .mix(CALL_IBD_HMMIBD_PARAM.out)
        .mix(CALL_IBD_ISORELATE_PARAM.out)
        .mix(CALL_IBD_REFINEDIBD_PARAM.out)
        .mix(CALL_IBD_TPBWT_PARAM.out)
        .map { label, chrno, caller, csp_id, ibd -> [groupKey([label, caller, csp_id], params.nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { k, ll -> [k[0], k[1], k[2], ll.collect { it[1] }] }
    // label, caller, csp_id, ibd_lst


    // true ibd
    CALL_IBD_TSKIBD(ch_trees_vcf.map { label, chrno, trees, _vcf -> [label, chrno, models[label], trees] })
    ch_true_ibd = CALL_IBD_TSKIBD.out.ibd
        .map { label, chrno, ibd -> [groupKey(label, params.nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { gkey, ll -> ["${gkey}", ll.collect { it[1] }] }
    // label, ibd_lst. 
    // Note gkey is not a string shoule get a string value of the key to combine with ch_inferred_ibd


    ch_true_inferr = ch_true_ibd
        .combine(ch_inferred_ibd, by: 0)
        .map { label, trueibd_lst, caller, csp_id, inferredibd_lst -> [label, models[label], caller, csp_id, trueibd_lst, inferredibd_lst] }
    //.view()



    // compare ibd
    CMP_IBD(ch_true_inferr)
}

// largely similar to optimize parameter workflow
// with different moddels and parameter parameters
workflow WF_VARY_RECOM_RATE {
    /////////////////////////////////////////////////////
    // models
    // - read models from json file
    def sp_defaults = params.sp_defaults
    def models = [
        "sp_neu_0003": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 3E-9, seqlen: (1.0 / 3E-9).toInteger()],
        "sp_neu_0010": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 1.0E-8, seqlen: (1.0 / 1.0E-8).toInteger()],
        "sp_neu_0030": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 3.0E-8, seqlen: (1.0 / 3.0E-8).toInteger()],
        "sp_neu_0100": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 1.00E-7, seqlen: (1.0 / 1.00E-7).toInteger()],
        "sp_neu_0300": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 3.00E-7, seqlen: (1.0 / 3.00E-7).toInteger()],
        "sp_neu_0667": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 6.67E-7, seqlen: (1.0 / 6.67E-7).toInteger()],
        "sp_neu_1000": sp_defaults + [s: 0.0, genome_set_id: 10000, r: 0.000001000, seqlen: (1.0 / 0.000001000).toInteger()],
    ]

    models.each { }
    // prepare model input chanel
    def nchroms = 4
    // no need to be as large as params.nchrom

    ch_input = Channel
        .fromList(models.collect { k, v -> [k, v] })
        .combine(Channel.fromList(1..nchroms))
        .map { label, args, chrno -> [label, chrno, args] }


    /////////////////////////////////////////////////////
    // log simulation parameters
    log_sim_params(ch_input, "vary_recom_rate_sp_sets.params.tsv")

    /////////////////////////////////////////////////////
    // simulations and conat vcf/tree results
    SIM_SP_CHR(ch_input)

    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf.filter { label, _chrno, _trees, _vcf -> label != "sp_neu_0003" }
    // filter out this one as its genome size in bp is too large for slow ibd callers such as refinedibd

    ch_csp = Channel.from(
        [
            ["hapibd", 1, [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000, minmarkers: 100, minmaf: 0.01]],
            ["hmmibd", 2, [m: 5, n: 100, minmaf: 0.01]],
            ["isorelate", 3, [minmaf: 0.01, min_snp: 20]],
            ["refinedibd", 4, [lod: 4.0, length: 2.0, scale: Math.sqrt(1000 / 100), minmaf: 0.01, window: 40.0, trim: 0.15]],
            ["tpbwt", 5, [Lm: 300, Lf: 2.0, template: 0, use_phase_correction: 0, minmaf: 0.01]],
        ]
    )


    ch_callibd = ch_trees_vcf
        .map { label, chrno, _trees, vcf -> [label, chrno, models[label], vcf] }
        .combine(ch_csp)
        .branch { _label, _chrno, args, _vcf, caller, _csp_id, _csp_args ->
            hapibd: caller == "hapibd"
            hmmibd: caller == "hmmibd"
            isorelate: (caller == "isorelate") && (args.r >= 1.00E-7)
            refinedibd: caller == "refinedibd"
            tpbwt: caller == "tpbwt"
            other: true
        }

    // to avoid issue when resume the pipeline
    println("view one item from the channel to avoid issue when resuming the pipeline")
    ch_callibd.hapibd.take(1).view()
    ch_callibd.hmmibd.take(1).view()
    ch_callibd.isorelate.take(1).view()
    ch_callibd.refinedibd.take(1).view()
    ch_callibd.tpbwt.take(1).view()

    // /////////////////////////////////////////////////
    // Compare IBD inferred vs true
    // inferred ibd
    CALL_IBD_HAPIBD_PARAM(ch_callibd.hapibd)
    CALL_IBD_HMMIBD_PARAM(ch_callibd.hmmibd)
    CALL_IBD_ISORELATE_PARAM(ch_callibd.isorelate)
    CALL_IBD_REFINEDIBD_PARAM(ch_callibd.refinedibd)
    CALL_IBD_TPBWT_PARAM(ch_callibd.tpbwt)

    // Log params
    log_ibdcall_csp_params(ch_callibd.hapibd, "vary_recom_rate_hapibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.hmmibd, "vary_recom_rate_hmmibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.isorelate, "vary_recom_rate_isorelate.params.tsv")
    log_ibdcall_csp_params(ch_callibd.refinedibd, "vary_recom_rate_refinedibd.params.tsv")
    log_ibdcall_csp_params(ch_callibd.tpbwt, "vary_recom_rate_tpbwt.params.tsv")


    ch_inferred_ibd = CALL_IBD_HAPIBD_PARAM.out
        .mix(CALL_IBD_HMMIBD_PARAM.out)
        .mix(CALL_IBD_ISORELATE_PARAM.out)
        .mix(CALL_IBD_REFINEDIBD_PARAM.out)
        .mix(CALL_IBD_TPBWT_PARAM.out)
        .map { label, chrno, caller, csp_id, ibd -> [groupKey([label, caller, csp_id], nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { k, ll -> [k[0], k[1], k[2], ll.collect { it[1] }] }
    // label, caller, csp_id, ibd_lst


    // true ibd
    CALL_IBD_TSKIBD(ch_trees_vcf.map { label, chrno, trees, _vcf -> [label, chrno, models[label], trees] })
    ch_true_ibd = CALL_IBD_TSKIBD.out.ibd
        .map { label, chrno, ibd -> [groupKey(label, nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { gkey, ll -> ["${gkey}", ll.collect { it[1] }] }
    // label, ibd_lst. 
    // Note gkey is not a string shoule get a string value of the key to combine with ch_inferred_ibd


    ch_true_inferr = ch_true_ibd
        .combine(ch_inferred_ibd, by: 0)
        .map { label, trueibd_lst, caller, csp_id, inferredibd_lst -> [label, models[label], caller, csp_id, trueibd_lst, inferredibd_lst] }

    ch_true_ibd.view { it -> it[0] }
    ch_inferred_ibd.view { it -> it[0] }
    ch_true_inferr.view()

    // compare ibd
    CMP_IBD(ch_true_inferr)
}


workflow WF_SP_COMPUTATION_BENCH {

    def sp_sets = get_sp_sets()

    // *********************** Log Params *************************
    def sp_set_args_keys = sp_sets.collect { _k, v -> v }[0].collect { k, _v -> k }
    ch_sp_sets = Channel.fromList(sp_sets.collect { label, args -> [label, args] })

    Channel
        .value(['label'] + sp_set_args_keys)
        .concat(
            ch_sp_sets.map { k, args ->
                [k] + sp_set_args_keys.collect { i -> args[i] }
            }
        )
        .map { lst -> lst.join("\t") }
        .collectFile(name: "sp_sets.tsv", newLine: true, storeDir: params.resdir, sort: false)

    // ***********************Simulation *************************

    ch_chrs = Channel.fromList(1..params.nchroms)

    ch_sp_input = ch_sp_sets
        .combine(ch_chrs)
        .map { label, args, chrno -> [label, chrno, args] }

    SIM_SP_CHR(ch_sp_input)


    // *********************** Call ibd *************************
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
    // label, chrno, trees, vcf

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets, by: 0)
        .map { label, chrno, trees, _vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets, by: 0)
        .map { label, chrno, _trees, vcf, args -> [label, chrno, args, vcf] }

    // CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    if (params.compute_bench_large_size == 0) {
        CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
        CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
        CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
        CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
        CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
        CALL_IBD_HMMIBDRS(ch_in_ibdcall_vcf)
        CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)
    }
    else {
        CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
        // CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
        CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    }
}
workflow THIN_MARKER {

    /////////////////////////////////////////////////////
    // models
    // - read models from json file
    def models = [
        "sp_neu": params.sp_defaults + [s: 0.0, genome_set_id: 10000],
        "mp_neu": params.mp_defaults + [s: 0.0, genome_set_id: 20000],
    ]
    // prepare model input chanel
    ch_input = Channel
        .fromList(models.collect { k, v -> [k, v] })
        .combine(Channel.fromList(1..params.nchroms))
        .map { label, args, chrno -> [label, chrno, args] }
        .branch { label, _chrno, _args ->
            sp: label.startsWith("sp_")
            mp: label.startsWith("mp_")
        }

    /////////////////////////////////////////////////////
    // simulations and conat vcf/tree results
    SIM_SP_CHR(ch_input.sp)
    SIM_MP_CHR(ch_input.mp)
    // SIM_UK_CHR(ch_input.uk)
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf.concat(SIM_MP_CHR.out.trees_vcf)
    // label, chrno, trees, vcf
    // .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    ////////////////////////////////////////////////////
    // call IBD with caller-specific parameters
    // - read csp json
    // - verify json data
    // - call ibd using the specified caller and parameters



    // read csp list from json files
    if (params.csp_json != "") {
        ch_csp = Channel
            .fromPath(params.csp_json)
            .flatMap { fn -> read_param_from_json(fn) }
    }
    else {
        // default list for testing purpose
        ch_csp = Channel.from(
            [
                [
                    "hmmibd",
                    1,
                    [
                        m: 5,
                        n: 100,
                        minmaf: 0.01,
                    ],
                ]
            ]
        )
    }


    ch_callibd = ch_trees_vcf
        .map { label, chrno, _trees, vcf -> [label, chrno, models[label], vcf] }
        .combine(ch_csp)
        .branch { label, _chrno, _args, _vcf, caller, _csp_id, _csp_args ->
            hapibd: caller == "hapibd"
            hmmibd: caller == "hmmibd" && (label != "uk_human") && (label != "uk_human2")
            isorelate: caller == "isorelate" && (label != "uk_human") && (label != "uk_human2")
            refinedibd: caller == "refinedibd"
            tpbwt: caller == "tpbwt"
            other: true
        }

    // to avoid issue when resume the pipeline
    println("view one item from the channel to avoid issue when resuming the pipeline")
    ch_callibd.hapibd.take(1).view()
    ch_callibd.hmmibd.take(1).view()
    ch_callibd.isorelate.take(1).view()
    ch_callibd.refinedibd.take(1).view()
    ch_callibd.tpbwt.take(1).view()

    // /////////////////////////////////////////////////
    // Compare IBD inferred vs true
    // inferred ibd
    CALL_IBD_HAPIBD_PARAM(ch_callibd.hapibd)
    CALL_IBD_HMMIBD_PARAM(ch_callibd.hmmibd)
    CALL_IBD_ISORELATE_PARAM(ch_callibd.isorelate)
    CALL_IBD_REFINEDIBD_PARAM(ch_callibd.refinedibd)
    CALL_IBD_TPBWT_PARAM(ch_callibd.tpbwt)

    ch_inferred_ibd = CALL_IBD_HAPIBD_PARAM.out
        .mix(CALL_IBD_HMMIBD_PARAM.out)
        .mix(CALL_IBD_ISORELATE_PARAM.out)
        .mix(CALL_IBD_REFINEDIBD_PARAM.out)
        .mix(CALL_IBD_TPBWT_PARAM.out)
        .map { label, chrno, caller, csp_id, ibd -> [groupKey([label, caller, csp_id], params.nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { k, ll -> [k[0], k[1], k[2], ll.collect { it[1] }] }
    // label, caller, csp_id, ibd_lst


    // true ibd
    CALL_IBD_TSKIBD(ch_trees_vcf.map { label, chrno, trees, _vcf -> [label, chrno, models[label], trees] })
    ch_true_ibd = CALL_IBD_TSKIBD.out.ibd
        .map { label, chrno, ibd -> [groupKey(label, params.nchroms), [chrno, ibd]] }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { gkey, ll -> ["${gkey}", ll.collect { it[1] }] }
    // label, ibd_lst. 
    // Note gkey is not a string shoule get a string value of the key to combine with ch_inferred_ibd


    ch_true_inferr = ch_true_ibd
        .combine(ch_inferred_ibd, by: 0)
        .map { label, trueibd_lst, caller, csp_id, inferredibd_lst -> [label, models[label], caller, csp_id, trueibd_lst, inferredibd_lst] }
    //.view()



    // compare ibd
    CMP_IBD(ch_true_inferr)
}




process SIM_SP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${params.resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/true_ne/", pattern: "*.true_ne", mode: "symlink"

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*.true_ne"), emit: true_ne

    script:
    def cmd_options = (args + [chrno: chrno])
        .collect { k, v ->
            "--${k} ${v}"
        }
        .join(" ")
    """
    sim_single_pop.py ${cmd_options}

    mkdir tmp; mv tmp_* tmp/
    """

    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

process SIM_MP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${params.resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/demog/", pattern: "*_demog.png", mode: "symlink"

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*_demog.png"), emit: demog

    script:
    def cmd_options = (args + [chrno: chrno])
        .collect { k, v ->
            "--${k} ${v}"
        }
        .join(" ")
    """
    sim_multiple_pop.py ${cmd_options}

    mkdir tmp; mv tmp_* tmp/
    """

    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,_demog.png}
    """
}

process CALL_IBD_HAPIBD {
    tag "${args.genome_set_id}_${chrno}_hapibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/hapibd/", mode: "symlink", pattern: "*_hapibd.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/hapibd/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_hapibd.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        chrno: chrno,
        seqlen: args.seqlen,
        minseed: params.hapibd_minseed,
        minoutput: params.hapibd_minoutput,
        maxgap: params.hapibd_maxgap,
        minextend: params.hapibd_minextend,
        minmarkers: params.hapibd_minmarkers,
        minmaf: params.minmaf,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_hapibd.py ${cmd_options}
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hapibd.ibd time_output.txt"""
}

process CALL_IBD_TSKIBD {
    tag "${args.genome_set_id}_${chrno}_tskibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/tskibd/", mode: "symlink", pattern: "*_tskibd.ibd"

    input:
    tuple val(label), val(chrno), val(args), path(trees)

    output:
    tuple val(label), val(chrno), path("*_tskibd.ibd"), emit: ibd

    script:
    def cmd_options_tskibd = [
        tree: trees,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")

    """
    call_tskibd.py ${cmd_options_tskibd}
    """

    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """touch  ${prefix}_tskibd.ibd time_output.txt"""
}
process CALL_IBD_REFINEDIBD {
    tag "${args.genome_set_id}_${chrno}_refinedibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/refinedibd/", mode: "symlink", pattern: "*_refinedibd.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/refinedibd/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_refinedibd.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        lod: params.refinedibd_lod,
        length: params.refinedibd_length,
        scale: params.refinedibd_scale,
        mem_gb: task.memory.giga,
        minmaf: params.minmaf,
        window: params.refinedibd_window,
        trim: params.refinedibd_trim,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_refinedibd.py ${cmd_options}
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_refinedibd.ibd time_output.txt"""
}

process CALL_IBD_TPBWT {
    tag "${args.genome_set_id}_${chrno}_tpbwt"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/tpbwt/", mode: "symlink", pattern: "*_tpbwt.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/tpbwt/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_tpbwt.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        template: params.tpbwt_template_opts,
        use_phase_correction: params.tpbwt_use_phase_correction,
        minmaf: params.minmaf,
        Lm: params.tpbwt_Lm,
        Lf: params.tpbwt_Lf,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    /usr/bin/time \
    call_tpbwt.py ${cmd_options} 2>time_output.txt
    cat tmp_time_output.txt >> time_output.tx
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_tpbwt.ibd time_output.txt"""
}

process CALL_IBD_HMMIBD {
    tag "${args.genome_set_id}_${chrno}_hmmibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/hmmibd/", mode: "symlink", pattern: "*_hmmibd.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/hmmibd/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_hmmibd.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        n: params.hmmibd_n,
        m: params.hmmibd_m,
        mincm: params.mincm,
        minmaf: params.minmaf,
        genome_set_id: args.genome_set_id,
        num_threads: task.cpus,
        version: params.hmmibd_version,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
        call_hmmibd.py ${cmd_options}
        """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd time_output.txt"""
}

// same as above except that in this version version is always hmmibd2(hmmibdrs) 
process CALL_IBD_HMMIBDRS {
    tag "${args.genome_set_id}_${chrno}_hmmibdrs"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/hmmibdrs/", mode: "symlink", pattern: "*_hmmibdrs.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/hmmibdrs/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_hmmibdrs.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        n: params.hmmibd_n,
        m: params.hmmibd_m,
        mincm: params.mincm,
        minmaf: params.minmaf,
        genome_set_id: args.genome_set_id,
        num_threads: task.cpus,
        version: "hmmibd2",
        optimize_for_large_size: args.optimize_for_large_size ?: 0,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
        call_hmmibd.py ${cmd_options}
        """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibdrs.ibd time_output.txt"""
}

process CALL_IBD_ISORELATE {
    tag "${args.genome_set_id}_${chrno}_isorelate"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/isorelate/", mode: "symlink", pattern: "*_isorelate.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/isorelate/${chrno}/", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf)

    output:
    tuple val(label), val(chrno), path("*_isorelate.ibd"), emit: ibd
    tuple val(label), val(chrno), path("time_output.txt"), emit: time

    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        min_snp: params.isorelate_min_snp,
        min_len_bp: Math.round(params.mincm * (0.01 / args.r)),
        minmaf: params.isorelate_minmaf,
        imiss: params.isorelate_imiss,
        vmiss: params.isorelate_vmiss,
        genome_set_id: args.genome_set_id,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_isorelate.py ${cmd_options} --cpus ${task.cpus}
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_isorelate.ibd time_output.txt"""
}


process PROC_DIST_NE {
    tag "${genome_set_id}"

    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/", pattern: "*.sh", mode: 'symlink'
    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/", pattern: "*.map", mode: 'symlink'
    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/", pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ibddist_ibd/", pattern: "*.ibddist.ibdobj.gz", mode: 'symlink'
    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ibdne_ibd/", pattern: "*.ibdne.ibdobj.gz", mode: 'symlink'

    input:
    tuple val(label), val(ibdcaller), path(ibd_lst), path(ibd_lst_true), path(vcf_lst), val(genome_set_id)
    path ibdne_jar

    output:
    tuple val(label), val(ibdcaller), path("ibdne.jar"), path("*_orig.sh"), path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
    tuple val(label), val(ibdcaller), path("ibdne.jar"), path("*_rmpeaks.sh"), path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
    tuple val(label), val(ibdcaller), path("*.ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
    tuple val(label), val(ibdcaller), path("*.ibdne.ibdobj.gz"), emit: ibdne_ibd_obj

    script:
    // Whether to pass in true ibd list is determined by params.filt_ibd_by_ov 
    def true_ibd_arg = params.filt_ibd_by_ov && (ibdcaller != "tskibd")
        ? [ibd_files_true: "${ibd_lst_true}"]
        : [:]
    def args_local = (true_ibd_arg + [
        ibd_files: "${ibd_lst}",
        vcf_files: "${vcf_lst}",
        genome_set_id: genome_set_id,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_jar: ibdne_jar,
        ibdne_flatmeth: params.ibdne_flatmeth,
        ibdne_flatmeth: params.ibdne_flatmeth,
        ibdne_no_diploid_conversion: params.ibdne_no_diploid_convertion,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    proc_dist_ne.py ${args_local} 
    """

    stub:
    """
    touch ibdne.jar
    touch ${genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${genome_set_id}.ibddist.ibdobj.gz
    touch ${genome_set_id}_orig.ibdne.ibdobj.gz
    touch ${genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
}

process PROC_INFOMAP {
    tag "${genome_set_id}"

    publishDir "${params.resdir}/${genome_set_id}_${label}/${ibdcaller}/ifm_input/", pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
    tuple val(label), val(ibdcaller), path(ibd_lst), path(vcf_lst), val(genome_set_id)

    output:
    tuple val(label), val(ibdcaller), path("*_orig.ifm.ibdobj.gz"), emit: ifm_orig_ibd_obj
    tuple val(label), val(ibdcaller), path("*_rmpeaks.ifm.ibdobj.gz"), emit: ifm_rmpeaks_ibd_obj

    script:
    def args_local = [
        ibd_files: "${ibd_lst}",
        vcf_files: "${vcf_lst}",
        genome_set_id: genome_set_id,
        peak_validate_meth: params.peak_validate_meth,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    proc_infomap.py ${args_local}
    """

    stub:
    """
    touch ${genome_set_id}{_orig.ifm.ibdobj.gz,_rmpeaks.ifm.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag "${args.genome_set_id}_${are_peaks_removed}"

    publishDir "${params.resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ne_output/", mode: 'symlink'

    input:
    tuple val(label), val(ibdcaller), path(ibdne_jar), path(ibdne_sh), path(gmap), path(ibd_gz), val(are_peaks_removed), val(args)

    output:
    tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*.ne")

    script:
    """
    bash ${ibdne_sh}
    """

    stub:
    def src = are_peaks_removed ? "rmpeaks" : "orig"
    """
    touch ${args.genome_set_id}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${args.genome_set_id}_${are_peaks_removed}"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ifm_output/", mode: 'symlink'

    input:
    tuple val(label), val(ibdcaller), path(ibd_obj), val(are_peaks_removed), val(args)

    output:
    tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*_member.pq")

    script:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    def args_local = [
        ibd_obj: ibd_obj,
        npop: args.npop,
        nsam: args.nsam,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
        ntrials: params.ifm_ntrials,
        transform: params.ifm_transform,
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    run_infomap.py ${args_local}
    """

    stub:
    def cut_mode = are_peaks_removed ? 'rmpeaks' : 'orig'
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
}

// Following demographic model used in the hap-ibd paper
process SIM_UK_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${params.resdir}/${args.genome_set_id}_${label}/trees/", pattern: "*.trees", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/vcf/", pattern: "*.vcf.gz", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/daf/", pattern: "*.daf", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/restart_count/", pattern: "*.restart_count", mode: "symlink"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/true_ne/", pattern: "*.true_ne", mode: "symlink"

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*.true_ne"), emit: true_ne

    script:
    // subset arguments as not all keys will be used.
    def cmd_options = (args + [chrno: chrno])
        .findAll { k, _v -> k in ["chrno", "seqlen", "gc", "r", "u", "nsam", "genome_set_id"] }
        .collect { k, v -> "--${k} ${v}" }
        .join(" ")
    """
    sim_uk_pop.py ${cmd_options}
    mkdir tmp; mv tmp_* tmp/
    """

    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}


process CALL_IBD_HAPIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_hapibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/hapibd/${csp_id}", mode: "symlink", pattern: "*_hapibd.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/hapibd/${csp_id}", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(caller), val(csp_id), val(hapibd_args)

    output:
    tuple val(label), val(chrno), val(caller), val(csp_id), path("*_hapibd.ibd")

    script:
    // call args is a dict with keys as follows
    // minmaf: params.minmaf,
    // minseed: params.hapibd_minseed,
    // minoutput: params.hapibd_minoutput,
    // maxgap: params.hapibd_maxgap,
    // minextend: params.hapibd_minextend,
    // minmarkers: params.hapibd_minmarkers,
    def cmd_options = (hapibd_args + [
        vcf: vcf,
        r: args.r,
        chrno: chrno,
        seqlen: args.seqlen,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_hapibd.py ${cmd_options}
    """

    stub:
    def _ibd_dir = String.format("%d", hapibd_args.minseed)
    """touch ${args.genome_set_id}_${chrno}_hapibd.ibd time_output.txt"""
}


process CALL_IBD_TPBWT_PARAM {
    tag "${args.genome_set_id}_${chrno}_tpbwt"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/tpbwt/${csp_id}", mode: "symlink", pattern: "*_tpbwt.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/tpbwt/${csp_id}", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(caller), val(csp_id), val(tpbwt_args)

    output:
    tuple val(label), val(chrno), val(caller), val(csp_id), path("*_tpbwt.ibd")

    script:
    // Lm: params.tpbwt_Lm,
    // Lf: params.tpbwt_Lf,
    // template: params.tpbwt_template_opts,
    // minmaf
    // use_phase_correction
    def cmd_options = (tpbwt_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    /usr/bin/time \
    call_tpbwt.py ${cmd_options} 2>time_output.txt
    cat tmp_time_output.txt >> time_output.txt 
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_tpbwt.ibd time_output.txt"""
}


process CALL_IBD_REFINEDIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_refinedibd"
    publishDir "${params.resdir}//${args.genome_set_id}_${label}/ibd/refinedibd/${csp_id}", mode: "symlink", pattern: "*_refinedibd.ibd"
    publishDir "${params.resdir}//${args.genome_set_id}_${label}/ibd/refinedibd/${csp_id}", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(caller), val(csp_id), val(refinedibd_args)

    output:
    tuple val(label), val(chrno), val(caller), val(csp_id), path("*_refinedibd.ibd")

    script:
    // lod: params.refinedibd_lod,
    // length: params.refinedibd_length,
    // scale: params.refinedibd_scale,
    // minmaf
    // window
    def cmd_options = (refinedibd_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_refinedibd.py ${cmd_options}
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_refinedibd.ibd time_output.txt"""
}

process CALL_IBD_HMMIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_hmmibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/hmmibd/${csp_id}", mode: "symlink", pattern: "*_hmmibd.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/hmmibd/${csp_id}", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(caller), val(csp_id), val(hmmibd_args)

    output:
    tuple val(label), val(chrno), val(caller), val(csp_id), path("*_hmmibd*.ibd")

    script:
    // n: params.hmmibd_n,
    // m: params.hmmibd_m,
    // minmaf:
    def cmd_options = (hmmibd_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
        num_threads: task.cpus,
        version: params.hmmibd_version,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
        call_hmmibd.py ${cmd_options}
        """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd"""
}


process CALL_IBD_ISORELATE_PARAM {
    tag "${args.genome_set_id}_${chrno}_isorelate"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibd/isorelate/${csp_id}", mode: "symlink", pattern: "*_isorelate.ibd"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/ibdcalltime/isorelate/${csp_id}", mode: "symlink", pattern: "time_output.txt"

    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(caller), val(csp_id), val(isorelate_args)

    output:
    tuple val(label), val(chrno), val(caller), val(csp_id), path("*_isorelate.ibd")

    script:
    // maf: params.maf,
    // min_snp: params.isorelate_min_snp,
    def cmd_options = (isorelate_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        min_len_bp: Math.round(params.mincm * (0.01 / args.r)),
        imiss: params.isorelate_imiss,
        vmiss: params.isorelate_vmiss,
        genome_set_id: args.genome_set_id,
    ]).collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    call_isorelate.py ${cmd_options} --cpus ${task.cpus}
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_isorelate.ibd"""
}

process CMP_IBD {
    tag "${args.genome_set_id}"
    publishDir "${params.resdir}/${args.genome_set_id}_${label}/cmpibd/${caller}/${csp_id}", mode: 'symlink'

    input:
    tuple val(label), val(args), val(caller), val(csp_id), path("trueibd_dir/*.ibd"), path("inferredibd_dir/*.ibd")

    output:
    path ("*.ovcsv"), emit: overlap
    path ("*.pairtotibdpq"), emit: pair_totibd
    path ("*.poptotibdcsv"), emit: pop_totibd

    script:
    def cmd_args = [
        id: "${label}_${caller}_${csp_id}",
        nchroms: params.nchroms,
        r: args.r,
        seqlen_bp: args.seqlen,
        trueibd_dir: "trueibd_dir",
        inferredibd_dir: "inferredibd_dir",
    ].collect { k, v -> "--${k} ${v}" }.join(" ")
    """
    cmp_ibd.py ${cmd_args}
    """

    stub:
    """
    touch caller.{ovcsv,pairtotibdpq,poptotibdcsv}
    """
}

// helper function to split string with ','
def to_lst(value) {
    def lst = "${value}".split(",").collect { it }
    return lst
}

def mp_ne(args) {
    def ne = args.N * args.npop + (args.npop - 1) * (args.npop - 1) / 4 / args.npop / args.sel_mig
    return ne
}

def get_maxtmrca(p) {
    def str = p.getName()
    def match = (str =~ /maxtmrca_(\d+)/)
    if (match) {
        def x = match[0][1].toInteger()
        return x
    }
    else {
        return null
    }
}



def read_param_from_json(json_fn) {
    assert json_fn != ""
    def json_slurper = new groovy.json.JsonSlurper()
    def f = file(json_fn)
    def params = json_slurper.parse(f)
    return params
}

def get_sp_sets() {
    def sp_sets = [
        sp_neu: params.sp_defaults + [s: 0.0, genome_set_id: 10000]
    ]
    if (params.sp_sets_json != "") {
        sp_sets = read_param_from_json(params.sp_sets_json)
    }
    return sp_sets
}

def get_mp_sets() {
    def mp_sets = [
        mp_s00: params.mp_defaults + [s: 0.0, genome_set_id: 20000]
    ]

    if (params.mp_sets_json != "") {
        mp_sets = read_param_from_json(params.mp_sets_json)
    }
}


def log_sim_params(ch_sim_input, filename) {
    // get a union set of args keys
    def ch_sim_args_keys = ch_sim_input
        .map { _label, _chrno, args -> args }
        .collect()
        .map { lst ->
            def args_all = [:]
            lst.each { it -> args_all = args_all + it }
            def args_keys = args_all.collect { k, _v -> k }
            [args_keys]
        }
    // preare header line
    def ch_log_header = ch_sim_input
        .take(1)
        .combine(ch_sim_args_keys)
        .map { _label, _chrno, _args, args_keys ->
            ["label", "chrno"] + args_keys
        }
    // prepare value lines
    def ch_log_values = ch_sim_input
        .combine(ch_sim_args_keys)
        .map { label, chrno, args, args_keys ->
            [label, chrno] + args_keys.collect { k -> args[k] }
        }

    // write a file
    def out_file = file(launchDir).resolve(params.resdir) / filename
    ch_log_header
        .concat(ch_log_values)
        .map { it.join("\t") }
        .collectFile(name: out_file, newLine: true, sort: false)
}

def get_global_caller_params(caller) {
    assert caller in ["tpbwt", "hapibd", "refinedibd", "hmmibd", "isorelate"]
    def args = [minmaf: params.minmaf, mincm: params.mincm]
    def caller_global_args = params
        .findAll { k, _v -> k.startsWith(caller + "_") }
        .collectEntries { k, v ->
            [k.replaceFirst(caller + "_", ""), v]
        }
    return args + caller_global_args
}

def log_ibdcall_params(ch_ibdcall_input, caller, filename) {
    def global_caller_args = get_global_caller_params(caller)
    // find all arg keys (union set)
    def ch_sim_args_keys = ch_ibdcall_input
        .map { _label, _chrno, sim_args, _vcf_or_tree -> sim_args }
        .collect()
        .map { lst ->
            def args = [:]
            lst.each { it -> args = args + it }
            def args_keys = args.collect { k, _v -> k }
            [args_keys]
        }
    def ch_caller_args_key = Channel.value(global_caller_args.collect { k, _v -> k })
    def ch_log_header = ch_ibdcall_input
        .take(1)
        .combine(ch_sim_args_keys)
        .combine(ch_caller_args_key)
        .map { _label, _chrno, _sim_args, _vcf_or_tree, sim_args_keys, caller_args_keys ->
            ["label", "chrno"] + sim_args_keys + ["caller"] + caller_args_keys
        }
    def ch_log_values = ch_ibdcall_input
        .combine(ch_sim_args_keys)
        .combine { ch_caller_args_key }
        .map { label, chrno, sim_args, _vcf_or_tree, sim_args_keys, caller_args_keys ->
            def sim_args_values = sim_args_keys.collect { k -> sim_args[k] }
            def caller_args_values = caller_args_keys.collect { k -> global_caller_args[k] }
            [label, chrno] + sim_args_values + [caller] + caller_args_values
        }

    def out_file = file(launchDir).resolve(params.resdir) / filename

    ch_log_header
        .concat(ch_log_values)
        .map { it.join("\t") }
        .collectFile(name: out_file, newLine: true, sort: false)
}

def log_ibdcall_csp_params(ch_ibdcall_csp_input, filename) {
    // find all keys across items
    def ch_sim_args_keys = ch_ibdcall_csp_input
        .map { _label, _chrno, sim_args, _vcf, _caller, _csp_id, _csp_args -> sim_args }
        .collect()
        .map { lst ->
            def args = [:]
            lst.each { it -> args = args + it }
            def args_keys = args.collect { k, _v -> k }
            [args_keys]
        }
    def ch_csp_args_keys = ch_ibdcall_csp_input
        .map { _label, _chrno, _sim_args, _vcf, _caller, _csp_id, csp_args -> csp_args }
        .collect()
        .map { lst ->
            def args = [:]
            lst.each { it -> args = args + it }
            def args_keys = args.collect { k, _v -> k }
            [args_keys]
        }
    def ch_log_header = ch_ibdcall_csp_input
        .take(1)
        .combine(ch_sim_args_keys)
        .combine(ch_csp_args_keys)
        .map { _label, _chrno, _sim_args, _vcf, _caller, _csp_id, _csp_args, sim_args_keys, csp_args_keys ->
            ["label", "chrno"] + sim_args_keys + ["caller", "csp_id"] + csp_args_keys
        }
    def ch_log_values = ch_ibdcall_csp_input
        .combine(ch_sim_args_keys)
        .combine(ch_csp_args_keys)
        .map { label, chrno, sim_args, _vcf, caller, csp_id, csp_args, sim_args_keys, csp_args_keys ->
            def sim_args_values = sim_args_keys.collect { k -> sim_args[k] }
            def csp_args_values = csp_args_keys.collect { k -> csp_args[k] }
            [label, chrno] + sim_args_values + [caller, csp_id] + csp_args_values
        }

    def out_file = file(launchDir).resolve(params.resdir) / filename

    print(out_file)
    ch_log_header
        .concat(ch_log_values)
        .map { it.join("\t") }
        .collectFile(name: out_file, newLine: true, sort: false)
}
