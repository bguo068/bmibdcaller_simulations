nextflow.enable.dsl=2

params.outdir = "results"
params.test = false
params.sim_related = 0
params.simulation_defaults= [
    num_chr : params.test ? 2: 14,
    label : "label",
    chrno : 1,
    s : 0.2,
    h : 0.5,
    selpos_0_1 : 0.33,
    num_origins : 1,
    s_start_g : 50,
    u : 1e-8,  // set to nozero value, used in simulate_chr_mutation
    bp_per_cm : 15000, // set to pf values
    seqlen_in_cm : 100,
    ne_change_start_g : 200,
    N : 10000,
    N0 : 10000,
    nsam : params.test ? 100 : 1000,
    sim_related: params.sim_related,
    test : params.test,
    minregion: params.test ? 10 : 50,
    ibd_cut_mode: "nocut",
    proc_label: "nocut",
    cut_start: 0,
    cut_end: 0
]
params.maf = params.test ? 0.00001: 0.01
params.mincm = 2.0
params.num_neutral_chrs = 0

params.tpbwt_template_opts = 1
params.tpbwt_Lm = params.test ? 152: 300
params.tpbwt_Lf = params.mincm

params.hapibd_minoutput = params.mincm
params.hapibd_minseed = params.mincm
params.hapibd_minextend = 1.0
params.hapibd_maxgap = 1000
params.hapibd_minmarkers = 100

params.refinedibd_length = params.mincm
params.refinedibd_lod = 3.0
params.refinedibd_scale = 0
params.tsktrueibd_mincm= params.mincm
params.ibdne_mincm = params.mincm

params.hmmibd_n = 100
params.hmmibd_m = 5

params.isorelate_imiss = 0.3
params.isorelate_vmiss = 0.3
params.isorelate_min_snp = 100

process SIMULATE_CHR {
    tag "simu_${d.label}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s1_simulations/${d.chrno}", mode: 'copy'
    input:
    val (d)
    output: 
    tuple val(d), path("${d.chrno}.trees"), emit: trees
    tuple val(d), path("tmp_slim_out*.trees"), path("*.daf"), path("*.true_ne"), path("*.restart_count"), emit: others

    script:
    """
    # Note here ignore mu and use 0 for simulating ancestry
    sim_anc.py \
        --chrno=${d.chrno} \
        --s=${d.s} \
        --h=${d.h} \
        --selpos_0_1=${d.selpos_0_1} \
        --num_origins=${d.num_origins} \
        --s_start_g=${d.s_start_g} \
        --u=0 \
        --bp_per_cm=${d.bp_per_cm} \
        --seqlen_in_cm=${d.seqlen_in_cm} \
        --ne_change_start_g=${d.ne_change_start_g} \
        --N=${d.N} \
        --N0=${d.N0} \
        --nsam=${d.nsam} \
        --sim_related=${d.sim_related}
    """
    stub:
    """touch  tmp_slim_out_${d.chrno}.trees ${d.chrno}.daf ${d.chrno}.restart_count ${d.chrno}.trees ${d.chrno}.true_ne"""
}
        

process SIMULATE_CHR_MUTATIONS {
    tag "simu_${d.label}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s1_simulations/${d.chrno}", mode: 'copy'
    input: 
    tuple val(d), path(in_trees_fn)
    output:
    tuple val(d), path("*.vcf.gz")
    script:
    """
    sim_mut.py \
        --in_chrno ${d.chrno} \
        --in_trees ${in_trees_fn} \
        --mu ${d.u} \
        --out_vcf unfilt.vcf.gz
    # filter by minor allele frequency
    bcftools view -q ${params.maf}:minor unfilt.vcf.gz -Oz -o ${d.chrno}.vcf.gz
    rm unfilt.vcf.gz
    """
    stub:
    """touch  ${d.chrno}.vcf.gz"""
}

process TSINFER_TSDATE_PER_CHR {
    tag "simu_${d.label}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s2_tsinfer_tsdate/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf_files), path(true_ne_df)
    output:
    tuple val(d), path("${d.chrno}.trees")
    script: 
    """
    tsinfer_tsdate_gw.py \
        --vcf_files ${vcf_files} \
        --in_bp_per_cm ${d.bp_per_cm} \
        --ne ${true_ne_df} # --out_trees ${d.chrno}.trees
    """
    stub:
    """touch  ${d.chrno}.trees"""
}

process CALC_RAW_IBD {
    tag "simu_${d.label}_${d.src}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s3_calc_raw_ibd/${d.src}/${d.chrno}", mode: 'copy'
    input:
    tuple val (d), path("${d.chrno}.trees")
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log

    script:
    """
    calc_raw_ibd.py --chrno ${d.chrno} --bp_per_cm ${d.bp_per_cm} --mincm ${params.tsktrueibd_mincm} ${d.chrno}.trees
    """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}
process CALL_IBD_HAPIBD {
    tag "simu_${d.label}_${d.chrno}_${d.ibdcaller}"
    publishDir "$params.outdir/${d.label}/s3_seqbased_raw_ibd/${d.ibdcaller}/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf)
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log
    script: 
        """
        call_hapibd.py --vcf ${vcf} --bp_per_cm ${d.bp_per_cm} --seqlen_in_cm ${d.seqlen_in_cm} --chrno ${d.chrno} \\
            --minseed ${params.hapibd_minseed} --minoutput ${params.hapibd_minoutput} --maxgap ${params.hapibd_maxgap} \\
            --minextend ${params.hapibd_minextend} --minmarkers ${params.hapibd_minmarkers} \\
            --mem_gb ${task.memory.giga} --nthreads ${task.cpus}
        """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}
process CALL_IBD_REFINEDIBD {
    tag "simu_${d.label}_${d.chrno}_${d.ibdcaller}"
    publishDir "$params.outdir/${d.label}/s3_seqbased_raw_ibd/${d.ibdcaller}/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf)
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log
    script: 
        """
        call_refinedibd.py --vcf ${vcf} --bp_per_cm ${d.bp_per_cm} --seqlen_in_cm ${d.seqlen_in_cm} --chrno ${d.chrno} \\
            --length ${params.refinedibd_length} --lod ${params.refinedibd_lod} --scale ${params.refinedibd_scale} \\
            --mem_gb ${task.memory.giga} --nthreads ${task.cpus}
        """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}

process CALL_IBD_TPBWT {
    tag "simu_${d.label}_${d.chrno}_${d.ibdcaller}"
    publishDir "$params.outdir/${d.label}/s3_seqbased_raw_ibd/${d.ibdcaller}/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf)
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log
    script: 
        """
        call_tpbwt.py --vcf ${vcf} --bp_per_cm ${d.bp_per_cm} --seqlen_in_cm ${d.seqlen_in_cm} --chrno ${d.chrno} \\
            --template ${params.tpbwt_template_opts} --Lm ${params.tpbwt_Lm} --Lf ${params.tpbwt_Lf} \\
            --mem_gb ${task.memory.giga} --nthreads ${task.cpus}
        """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}

process CALL_IBD_HMMIBD {
    tag "simu_${d.label}_${d.chrno}_${d.ibdcaller}"
    publishDir "$params.outdir/${d.label}/s3_seqbased_raw_ibd/${d.ibdcaller}/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf)
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log
    script: 
        """
        call_hmmibd.py --vcf ${vcf} \\
            --bp_per_cm ${d.bp_per_cm} --seqlen_in_cm ${d.seqlen_in_cm} --chrno ${d.chrno} \\
            --n ${params.hmmibd_n} --m ${params.hmmibd_m} --mincm ${params.mincm}
        """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}

process CALL_IBD_ISORELATE {
    tag "simu_${d.label}_${d.chrno}_${d.ibdcaller}"
    publishDir "$params.outdir/${d.label}/s3_seqbased_raw_ibd/${d.ibdcaller}/${d.chrno}", mode: 'copy'
    input:
    tuple val(d), path(vcf)
    output: 
    tuple val(d), path("${d.chrno}.ibd"), emit: raw_ibd
    tuple val(d), path("${d.chrno}.map"), emit: gmap
    tuple val(d), path("${d.chrno}.log"), emit: log
    script: 
        """
        call_isorelate.py --vcf ${vcf} \\
            --bp_per_cm ${d.bp_per_cm} --seqlen_in_cm ${d.seqlen_in_cm} --chrno ${d.chrno} \\
            --cpus ${task.cpus} --min_snp ${params.isorelate_min_snp} --min_len_bp ${Math.round(params.mincm * d.bp_per_cm)} \\
            --maf ${params.maf} --imiss ${params.isorelate_imiss} --vmiss ${params.isorelate_vmiss}
        """
    stub:
    """touch  ${d.chrno}.ibd ${d.chrno}.map ${d.chrno}.log"""
}

process CMP_TRUE_AND_INFERRED_RAW_IBD {
    tag         "cmp_ibd_${d_infer.label}_${d_infer.chrno}_${d_infer.src}"
    publishDir  "${params.outdir}/${d_infer.label}/S3_cmp_ibd/${d_infer.chrno}/${d_infer.src}", mode: 'copy'
    input:      tuple val(d_infer), path(ibd_infer, stageAs: 'inferred*.ibd'), path(ibd_true, stageAs: 'true*.ibd')
        // use stageAs to avoid name collison
    output: 
                tuple val(d_infer), path('*false_neg.tsv'), path('*false_pos.tsv'), path('*binned_totalibd.tsv'), path('*chr_wide_error.tsv'), emit: tables
                tuple val(d_infer), path('*plot.png'), emit: images

    script: 
    def n_infer_ibd_files = (ibd_infer instanceof List) ? ibd_infer.size() : 1
    def genome_size_cm = d_infer.seqlen_in_cm * n_infer_ibd_files
    def tmrca_opt1 = d_infer.min_tmrca ? "--min_tmrca ${d_infer.min_tmrca}" : ''
    def tmrca_opt2 = d_infer.max_tmrca ? "--max_tmrca ${d_infer.max_tmrca}" : ''
    """
    cmp_ibd.py --raw_true_ibd ${ibd_true} --raw_detected_ibd ${ibd_infer} --bp_per_cm ${d_infer.bp_per_cm} --genome_size_cm ${genome_size_cm} \\
        ${tmrca_opt1} ${tmrca_opt2} --out_prefix ${d_infer.src}_${d_infer.chrno}_
    """
    stub:
    """
    touch  ${d_infer.src}_${d_infer.chrno}_false_neg.tsv
    touch  ${d_infer.src}_${d_infer.chrno}_false_pos.tsv
    touch  ${d_infer.src}_${d_infer.chrno}_binned_totalibd.tsv
    touch  ${d_infer.src}_${d_infer.chrno}_chr_wide_error.tsv
    touch  ${d_infer.src}_${d_infer.chrno}_plot.png
    """
}

process FIND_PEAKS {
    tag "simu_${d.label}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s4_find_peaks_hapibd/${d.chrno}", mode: 'copy'
    input:
    tuple val (d), path(ibd)
    output:
    tuple val(d), path("*peaks.bed"), emit: ibd_peaks
    tuple val(d), path("*peaks.png"), emit: ibd_peak_plot

    script:
    """
    find_ibd_peaks.py --raw_ibd ${ibd} --chr_name ${d.chrno} --chr_end ${d.seqlen_in_cm * d.bp_per_cm} --out_bed ${d.chrno}_peaks.bed
    """
    stub:
    """
    touch  ${d.chrno}_peaks.bed ${d.chrno}_peaks.png
    """
}

process RM_HIGHLY_RELATED {
    tag "simu_${dd[1].label}_${dd[1].src}"
    publishDir "$params.outdir/${dd[1].label}/s4_rm_highly_related/${dd[1].src}", mode: 'copy'
    input:
    tuple val(dd), path(ibds) // should be sorted according chrno 
    output:
    tuple val(dd), path('rm/*.ibd')

    script:
    """
    mkdir rm
    rm_highly_related.py --raw_ibd_files ${ibds} --outdir rm 
    """
    stub:
    """
    mkdir rm
    touch rm/{1..${dd[1].num_chr}}.ibd 
    """
}

process PROC_IBD {
    tag "simu_${d.label}_${d.src}_${d.proc_label}_${d.chrno}"
    publishDir "$params.outdir/${d.label}/s4_proc_ibd/${d.src}/${d.proc_label}/${d.chrno}", mode: 'copy'
    input:
    tuple val (d), path("${d.chrno}.ibd"), path(peaks)
    output:
    tuple val(d), path("${d.chrno}.ibd_proc"),path("${d.chrno}.map_proc"), emit: proc_ibd
    tuple val(d), path("${d.chrno}.targets_cut"), emit: targets

    script:
    def need_to_remove_old_ibd = d.ibd_tmrca_cutoff ? 'true' : 'false'
    def use_peaks = (d.ibd_cut_mode == "autocut_hapibd") && (peaks) 
    def peak_params = use_peaks ? "--peaks_to_rm ${peaks}" : ""
    """
    if ${need_to_remove_old_ibd}; then
        mv ${d.chrno}.ibd all_${d.chrno}.ibd 
        awk 'NR==1 || \$6 <= ${d.ibd_tmrca_cutoff} {print} ' all_${d.chrno}.ibd > ${d.chrno}.ibd
    fi
    if [ \$(head ${d.chrno}.ibd | wc -l) -eq 1 ]; then 
        touch ${d.chrno}.ibd_proc ${d.chrno}.map_proc ${d.chrno}.targets_cut
        exit 0
    fi
    proc_ibd.py --chrno ${d.chrno} --bp_per_cm ${d.bp_per_cm} --ibd_in ${d.chrno}.ibd \
        --ibd_cut_mode ${d.ibd_cut_mode} --cut_start ${d.cut_start} --cut_end ${d.cut_end} \
        ${peak_params}
    """
    stub:
    """
    touch ${d.chrno}.ibd_proc ${d.chrno}.map_proc ${d.chrno}.targets_cut
    """
}


process CALL_IBDNE {
    tag "ne_${d.label}_${d.proc_label}_${d.src}"
    publishDir "${params.outdir}/${d.label}/s5_call_ne/${d.src}/${d.proc_label?: "no_proc"}", mode: 'copy'
    input: 
    tuple val(d), path("ibd_files/*" ), path("map_files/*")
    output:
    tuple val(d), path("ne_res.ne"), emit: ne
    tuple val(d), path("ne_res.region.excl"), emit: region_excl
    tuple val(d), path("ne_res.pair.excl"), emit: pair_excl
    tuple val(d), path("ne_res.boot"),path("ne_res.log"), emit: others

    script:
    """
    set -eux
    cat map_files/{1..${d.num_chr}}.map_proc > all.map
    set +e # stop exit on error
    cat ibd_files/{1..${d.num_chr}}.ibd_proc | java -Xmx${task.memory.giga}G \\
    	-jar ${projectDir}/lib/ibdne.23Apr20.ae9.jar \\
    	map=all.map out=ne_res nthreads=${task.cpus} minregion=${d.minregion} \\
        mincm=${params.ibdne_mincm}  2>captured_err.txt
    #
    # when ibdne fails due to no IBD segment, exit 0 and touch needed files
    ret=\$?
    grep -e 'ERROR: No IBD segments remain after applying filters' captured_err.txt
    match=\$?
    if [ \$ret -ne 0 ] && [ \$match -eq 0 ]; then 
        touch ne_res.ne ne_res.region.excl ne_res.pair.excl ne_res.boot ne_res.log
        >&2 cat captured_err.txt
    fi

    """ 
    stub: 
    """
    touch ne_res.{ne,region.excl,pair.excl,boot,log}
    """
}

process SUMMARIZE_TREES {
    tag "treeseq_summary_${label}"
    publishDir "${params.outdir}/${label}/s6_summarize_treeseq", mode: 'copy'
    publishDir "${params.outdir}/summary/summarize_TreeSeq", mode: 'copy'
    input:
    tuple val(label), stdin
    output:
    tuple path('*.png'), path('*.tsv')
    script:
    """
    cat - > info.tsv
    summarize_treeseq.py --ts_tsv info.tsv
    """
    stub:
    """
    touch out.png out.tsv
    """
}

process SUMMARIZE_IBD {
    tag "ibd_summary_${label}"
    publishDir "${params.outdir}/${label}/s7_summarize_ibd", mode: 'copy'
    publishDir "${params.outdir}/summary/summarize_IBD", mode: 'copy'
    input:
    tuple val(label), stdin
    output:
    tuple path('*.png'), path('*.tsv')
    script:
    """
    cat - > info.tsv
    summarize_ibd.py --ibd_tsv info.tsv
    """
    stub:
    """
    touch out.png out.tsv
    """
}

process SUMMARIZE_NE {
    tag "ne_summary_${label}"
    publishDir "${params.outdir}/${label}/s8_summarize_ne", mode: 'copy'
    publishDir "${params.outdir}/summary/summarize_NE", mode: 'copy'
    input: 
    tuple val(label), stdin 
    file(true_ne_tsv)
    output: 
    tuple path('*.png'), path('*.tsv')
    script:
    """
    cat - > info.tsv
    summarize_ne.py --ne_tsv info.tsv --true_ne_tsv ${true_ne_tsv}
    """
    stub:
    """
    touch out.png out.tsv
    """
}

workflow SELECTION_EFFECTS_USING_TRUEIBD {

    def new_defaults = [s: 0.3, num_origins: 1, s_start_g: 80, N0: 1000]
    def params_dict = params.simulation_defaults + new_defaults

    ch_different_s = channel.fromList([
        [s: 0,   label: 'diff_s_00'],
        [s: 0.1, label: 'diff_s_01'],
        [s: 0.2, label: 'diff_s_02'],
        [s: 0.3, label: 'diff_s_03'],
        [s: 0.4, label: 'diff_s_04'],
    ])

    ch_different_norigs = channel.fromList([
        [ s : 0.0,          label: 'diff_norig_0'],
        [ num_origins: 1,   label: 'diff_norig_1'],
        [ num_origins: 3,   label: 'diff_norig_3'],
        [ num_origins: 9,   label: 'diff_norig_9'],
        [ num_origins: 27,  label: 'diff_norig_27'],
        [ num_origins: 81,  label: 'diff_norig_81'],
        [ num_origins: 243, label: 'diff_norig_243'],
    ])

    ch_different_g = channel.fromList([
        [ s: 0.0,         label: 'diff_g_00'],
        [ s_start_g: 10,  label: 'diff_g_10'],
        [ s_start_g: 40,  label: 'diff_g_40'],
        [ s_start_g: 80,  label: 'diff_g_80'],
        [ s_start_g: 120, label: 'diff_g_120'],
        [ s_start_g: 160, label: 'diff_g_160'],
    ])

    ch_simu_set1 = ch_different_s.map{ d2 -> params_dict + d2 }
    ch_simu_set2 = ch_different_norigs.map{ d2 -> params_dict + d2 }
    ch_simu_set3 = ch_different_g.map{ d2 -> params_dict + d2 }
    ch_simu_sets = ch_simu_set1.concat(ch_simu_set2).concat(ch_simu_set3)

    ch_diff_proc = channel.fromList([
        [ibd_cut_mode: 'nocut',   proc_label: 'nocut'],
        [ibd_cut_mode: 'autocut', proc_label: 'autocut']
    ])

    ch_simu_chr = ch_simu_sets.flatMap{d-> 
        (1..(d.num_chr)).collect{ d + [chrno:it] }
    }

    // allow $params.num_neutral_chrs chromosomes to be neutral while the rest to be under selection
    ch_simu_chr = ch_simu_chr.map {d ->
        if (d.chrno <= params.num_neutral_chrs){
            d.s = 0
            d.num_origins = 0
            d.h = 0
            d.s_start_g = 0
        }
        d
    }

    ch_simu_chr 
    | SIMULATE_CHR

    ch_trees = SIMULATE_CHR.out.trees.map {d, trees -> [ d + ['src': 'true_tree'], trees ]}
    in_ch_for_summarize_trees = ch_trees
        .map {d, trees -> 
            def simlabel = d.label
            def str = "${d.label}\t${d.src}\t${d.chrno}\t${trees.toAbsolutePath()}"
            return [groupKey(simlabel, d.num_chr), str]
        }
        .groupTuple()
        .map {simlabel, str_list-> [simlabel, str_list.sort().join('\n')]}
    
    in_ch_for_summarize_trees 
    | SUMMARIZE_TREES

    ch_trees 
    | CALC_RAW_IBD

    ch_raw_ibd = CALC_RAW_IBD.out.raw_ibd

    in_ch_for_summarize_ibd = ch_raw_ibd
        .map{d, ibd ->
            def simlabel = d.label
            def str = "${d.label}\t${d.src}\t${d.chrno}\t${ibd.toAbsolutePath()}"
            [simlabel, str]
        }
        .groupTuple()
        .map {simlabel, str_list-> [simlabel, str_list.sort().join('\n')]}
    in_ch_for_summarize_ibd
    | SUMMARIZE_IBD

    ///////////////////////////////////////////////
    // group chromosomes to remove highly related, then convert back to normal
    // IBD per chromosome
    if (params.rm_related) {
        ch_raw_ibd | map {d, ibd ->
            def grp = "${d.label}:${d.src}"
            return [groupKey(grp, d.num_chr), [d, ibd]]
        }  | groupTuple(sort: {item->item[0].chrno} ) 
           | map {grp, ll -> 
                def dd = ll.collectEntries{[ it[0].chrno, it[0] ]}
                def ibd_lst = ll.collect{it[1]}
                return [dd, ibd_lst] } 
           | RM_HIGHLY_RELATED
        
        ch_raw_ibd = RM_HIGHLY_RELATED.out
        | map {dd, ibd_unordered ->
            ll = ibd_unordered.collect{
                    [ it.getSimpleName().toInteger(), it] // chrno, ibd
                }
                .sort {a, b-> a[0]<=>b[0]} // sort by chrno
            ll.collect{chrno, ibd -> [dd[chrno], ibd] }
        } | flatMap
    }

    ch_ibd_with_proc = ch_raw_ibd
        .combine(ch_diff_proc) // d, ibd, proc
        .map {d, ibd, proc -> [ d+ proc, ibd] }
        .map {d, ibd -> def peaks = []; [d, ibd, peaks] }

    ch_ibd_with_proc
    | PROC_IBD

    ch_ibdne = PROC_IBD.out.proc_ibd.map{d, ibd, gmap -> 
            [groupKey("${d.label}_${d.proc_label}_${d.src}", d.num_chr), [d, ibd, gmap]]} // use GroupKey to allow groupying without wait for all tasks to be done
        .groupTuple(by:0, sort: {a, b-> a[0].chrno<=>b[0].chrno})//.view() // sort by chrno
        .map{label, items -> 
            def d=items[0][0] + [chrno: 'all', minregion: 10] // only need one param_dict copy
            return [d, items.collect{it[1]}, items.collect{it[2]}];
            }//.dump(tag: 'ne_input')
    CALL_IBDNE(ch_ibdne)

    // -------------------------- Summarize Ne results -----------------------------
    in_ch_summarize_ne = CALL_IBDNE.out.ne
        .map{d, ne -> def str = "${d.label}\t${d.src}\t${d.proc_label}\t${ne.toAbsolutePath()}"
            return [d.label, str]}
        .groupTuple()
        .map {label, str_list -> [label, str_list.sort().join('\n')]}

    in_ch_true_ne = SIMULATE_CHR.out.others
        // each chr given true_ne but they are redudant and only one is needed
        .filter{d, slimtree, daf, truene, restart_count ->d.chrno == 1} 
        .map{d, slimtree, daf, truene, restart_count ->
            def str = "${d.label}\t${truene.toAbsolutePath()}"
        }
        .toList()
        .map {str_list -> str_list.sort().join('\n')}

    SUMMARIZE_NE (in_ch_summarize_ne, in_ch_true_ne)
}

workflow {
    def params_dict = params.simulation_defaults

    // -------------------------------- different selection coefficients ------------------------------
    ch_different_s = channel.fromList([
        [s: 0,                  label: 'Neutral'],
        // [s: 0.2, s_start_g: 80, label: 'selCoeff_0.2_selStartG_80'],
        [s: 0.4, s_start_g: 50, label: 'selCoeff_0.4_selStartG_50'],
    ])

    // -------------------------------- different Ne -------------------- ------------------------------
    ch_desc_ne = channel.fromList([
        [N0: 100000,    label: 'N0_100000'],
        // [N0: 30000,     label: 'N0_30000'],
        [N0: 10000,     label: 'N0_10000'],
        // [N0: 3000,      label: 'N0_3000'],
        [N0: 1000,      label: 'N0_1000'],
    ])

    // -------------------------------- make simulation sets -------------- ------------------------------
    ch_simu_set = ch_different_s
        .combine(ch_desc_ne).map {left, right-> 
            def label_combined  = "${right.label}_${left.label}"
            return params_dict + left + right + [label: label_combined]
        }//.view()

    // -------------------------------- inject different ibd processing------ ------------------------------
    ch_diff_proc = channel.fromList([
        [ibd_cut_mode: 'nocut',     proc_label: "nocut" ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s00e50",    cut_start: "0", cut_end: "50000000", minregion: 20 ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s00e60",    cut_start: "0", cut_end: "60000000", minregion: 20 ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s00e80",    cut_start: "0", cut_end: "80000000", minregion: 10 ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s23e43",   cut_start: "23000000" ,cut_end: "43000000", minregion: 20 ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s13e53",   cut_start: "13000000",cut_end: "53000000", minregion: 20 ],
        // [ibd_cut_mode: 'hardcut',   proc_label: "hardcut_s03e63",    cut_start:  "3000000", cut_end: "63000000", minregion: 20 ],
        [ibd_cut_mode: 'autocut',   proc_label: "autocut_hapibd",                                          minregion:20],
        // [ibd_cut_mode: 'mutcut',    proc_label: "mutcut",                                           minregion:20],
    ])

    // ---------------------------- for testing, also see params.simulation_defaults --------------------
    if (params_dict.test) {
        ch_simu_set = ch_simu_set.take(2)
    }

    // ----------------------------- expand simulation sets into chromosomes --------------------------------
    ch_simu_chr = ch_simu_set.flatMap{d-> 
        (1..(d.num_chr)).collect{ d + [chrno:it] }
    }
    
    // allow $params.num_neutral_chrs chromosomes to be neutral while the rest to be under selection
    ch_simu_chr = ch_simu_chr.map {d ->
        if (d.chrno <= params.num_neutral_chrs){
            d.s = 0
            d.num_origins = 0
            d.h = 0
            d.s_start_g = 0
        }
        d
    }

    // ---------------------------- RUN simulation ancestry ----------------------------------------
    SIMULATE_CHR(ch_simu_chr)

    // ---------------------------- RUN simulation mutation ----------------------------------------
    SIMULATE_CHR_MUTATIONS(SIMULATE_CHR.out.trees)
    ch_sim_mut = SIMULATE_CHR_MUTATIONS.out
    
    // -----------------------------RUN tsinfer tsdate ---------------------------------------------

    // tsinfer tsdate input format: tuple val(d), path(vcf), path(ne_df)
    // combine vcf and ne by simulation label + chrno
    ch_ne = SIMULATE_CHR.out.others.map{d, tree_fn, daf_fn, ne_fn, restart_count -> 
        ["${d.label}:${d.chrno}", ne_fn] }    // add key and extract ne item
    ch_vcf_ne = ch_sim_mut.map{d, vcf -> 
            ["${d.label}:${d.chrno}",d, vcf]} // add key 
        .combine(ch_ne, by: 0)                // combine by key
        .map{tmplabel, d, vcf, ne -> [d, vcf, ne]} // remove tmp key
    TSINFER_TSDATE_PER_CHR(ch_vcf_ne)

    // ------------------------- Gather trees ------------------------------------------------
    ch_inferred_tree = TSINFER_TSDATE_PER_CHR.out.map{d, trees -> def d2 = d+['src': 'tsinfer_tree']; [d2, trees]}
    ch_simulated_tree = SIMULATE_CHR.out.trees.map {d, trees -> def d2 = d + ['src': 'true_tree']; [d2, trees]}
    ch_all_trees = ch_simulated_tree.concat(ch_inferred_tree)

    // ------------------------- Summarize Trees -----------------------------------------------
    in_ch_for_summarize_trees = ch_all_trees
        .map {d, trees -> 
            def simlabel = d.label
            def str = "${d.label}\t${d.src}\t${d.chrno}\t${trees.toAbsolutePath()}"
            return [groupKey(simlabel, d.num_chr), str]
        }
        .groupTuple()
        .map {simlabel, str_list-> [simlabel, str_list.sort().join('\n')]}
    
    SUMMARIZE_TREES(in_ch_for_summarize_trees)

    // ------------------------- -- Run tskibd --------------------------------------------------
    CALC_RAW_IBD(ch_all_trees)

    // ---------------------------- call ibd --------------------------------------------------
    // NOTE: the original general process CALL_IBD was split into multiple
    // So that each process can have separate conda configuration
    ch_in_4_ibd_calling = ch_sim_mut.multiMap{d, vcf-> 
        hapibd: [ ['src': 'hapibd', 'ibdcaller': 'hapibd'] + d, vcf ]
        refinedibd: [ ['src': 'refinedibd', 'ibdcaller': 'refinedibd'] + d, vcf ]
        tpbwt: [ ['src': 'tpbwt', 'ibdcaller': 'tpbwt'] + d, vcf ]
        hmmibd: [ ['src': 'hmmibd', 'ibdcaller': 'hmmibd'] + d, vcf ]
        isorelate: [ ['src': 'isorelate', 'ibdcaller': 'isorelate'] + d, vcf ]
    }

    ch_in_4_ibd_calling.hapibd | CALL_IBD_HAPIBD
    ch_in_4_ibd_calling.refinedibd | CALL_IBD_REFINEDIBD
    ch_in_4_ibd_calling.tpbwt | CALL_IBD_TPBWT
    ch_in_4_ibd_calling.hmmibd | CALL_IBD_HMMIBD
    ch_in_4_ibd_calling.isorelate | CALL_IBD_ISORELATE


    // mix IBD from snp-based ibdcallers
    ch_call_ibd_out = CALL_IBD_HAPIBD.out.raw_ibd.concat(
        CALL_IBD_REFINEDIBD.out.raw_ibd, 
        CALL_IBD_TPBWT.out.raw_ibd,
        CALL_IBD_HMMIBD.out.raw_ibd,
        CALL_IBD_ISORELATE.out.raw_ibd
    )

    // --------------------------- directly compare raw ibd ----------------------------------------
    // Format tuple val(d_infer), path(ibd_infer), path(ibd_true)
    ch_true_raw_ibd = CALC_RAW_IBD.out.raw_ibd.filter {d, ibd-> d.src == 'true_tree'}
    ch_infer_raw_ibd = CALC_RAW_IBD.out.raw_ibd.filter {d, ibd-> d.src != 'true_tree'}
        .concat (ch_call_ibd_out)
    ch_true_infer_raw_ibd = 
            ch_true_raw_ibd .map{d, ibd -> def sim_chr_label = "${d.label}_${d.chrno}"; [sim_chr_label, ibd]}
        .combine(
            ch_infer_raw_ibd.map{d, ibd -> def sim_chr_label = "${d.label}_${d.chrno}"; [sim_chr_label, d, ibd]},
            by: 0
        )
        .map {sim_chr_label,  ibd_true, d_infer, ibd_infer -> 
            def sim_src_label = "${d_infer.label}_${d_infer.src}"
            [sim_src_label, [d_infer.chrno, d_infer, ibd_infer, ibd_true] ]}
        .groupTuple(sort: {a, b -> a[0] <=> b[0]})
        .map {sim_src_label, ll ->
            def d_infer_new = ll[0][1] + ['chrno': 'all']
            def ibd_infer = ll.collect{it[2]}
            def ibd_true = ll.collect{it[3]}
            return [d_infer_new, ibd_infer, ibd_true]
        }
        .flatMap{ d, ibd_infer, ibd_true ->
            if (d.src == 'tsinfer_tree')
                return [
                    [d + [src: 'tsinfer_tree_0_1000', min_tmrca: 0, max_tmrca: 1000], ibd_infer, ibd_true],
                    [d + [src: 'tsinfer_tree_0_3000', min_tmrca: 0, max_tmrca: 3000], ibd_infer, ibd_true],
                    [d + [src: 'tsinfer_tree_0_10000', min_tmrca: 0, max_tmrca: 10000], ibd_infer, ibd_true],
                ]
            else 
                return [
                    [d, ibd_infer, ibd_true]
                ]
        }

    CMP_TRUE_AND_INFERRED_RAW_IBD(ch_true_infer_raw_ibd)


    // --------------------------- Gather all IBD ---------------------------------------
    ch_combined_all_ibd = CALC_RAW_IBD.out.raw_ibd
        .concat(ch_call_ibd_out) // first concat IBD all of all sources: simulated tree, inferred tree, seq-based

    // ---------------------------- Summarize IBD ----------------------------------
    in_ch_for_summarize_ibd = ch_combined_all_ibd
        .map{d, ibd ->
            def simlabel = d.label
            def str = "${d.label}\t${d.src}\t${d.chrno}\t${ibd.toAbsolutePath()}"
            [simlabel, str]
        }
        .groupTuple()
        .map {simlabel, str_list-> [simlabel, str_list.sort().join('\n')]}
    SUMMARIZE_IBD(in_ch_for_summarize_ibd)

    // ---------------------------- Process raw IBD ----------------------------------
    ch_ibd_with_proc = ch_combined_all_ibd
        .combine(ch_diff_proc) // iterate for different proc methods
        .map{d, ibd, proc-> [d + proc, ibd]}
        .map{d, ibd-> [d.src, d, ibd] }
    
    ch_tsinfer_ibd_cutoff = channel.fromList([
        // ['tsinfer_tree', 100],
        // ['tsinfer_tree', 300],
        ['tsinfer_tree', 1000],
        ['tsinfer_tree', 3000],
        ['tsinfer_tree', 10000],
        ['true_tree', 'None' ],
        ['hapibd', 'None' ],
        ['refinedibd', 'None' ],
        ['tpbwt', 'None' ],
        ['hmmibd', 'None' ],
        ['isorelate', 'None' ],
    ])

    CALL_IBD_HAPIBD.out.raw_ibd | FIND_PEAKS
    ch_peaks = FIND_PEAKS.out.ibd_peaks.map{d, peaks -> def grp="${d.label}_${d.chrno}"; [grp,  peaks]}

    in_ch_proc_ibd = ch_ibd_with_proc
        .combine(ch_tsinfer_ibd_cutoff, by: 0)
        .map {src, d, ibd, cutoff ->
            def grp = "${d.label}_${d.chrno}"
            def d2 = d.clone()
            if (cutoff != 'None') {
                d2.src = "${d.src}_upto_n" + String.format("%05d", cutoff)
                d2 = d2+ ['ibd_tmrca_cutoff': cutoff]
            }
            return [grp, d2, ibd]
        }
        .combine(ch_peaks, by: 0)
        .map {grp, d, ibd, peaks -> [d, ibd, peaks]}
        // .view{[it[0].src, it[0].ibd_tmrca_cutoff]}

    PROC_IBD(in_ch_proc_ibd)
    
    //  ---------------------------- Run IBDNe ------------------------------------
    ch_ibdne = PROC_IBD.out.proc_ibd.map{d, ibd, gmap -> 
            [groupKey("${d.label}_${d.proc_label}_${d.src}", d.num_chr), [d, ibd, gmap]]} // use GroupKey to allow groupying without wait for all tasks to be done
        .groupTuple(by:0, sort: {a, b-> a[0].chrno<=>b[0].chrno})//.view() // sort by chrno
        .map{label, items -> 
            def d=items[0][0] + [chrno: 'all', minregion: 10] // only need one param_dict copy
            return [d, items.collect{it[1]}, items.collect{it[2]}];
            }//.view()
    CALL_IBDNE(ch_ibdne)


    // -------------------------- Summarize Ne results -----------------------------
    in_ch_summarize_ne = CALL_IBDNE.out.ne
        .map{d, ne -> 
            def str = "${d.label}\t${d.src}\t${d.proc_label}\t${ne.toAbsolutePath()}"
            return [d.label, str]}
        .groupTuple()
        .map {label, str_list -> [label, str_list.sort().join('\n')]}

    in_ch_true_ne = SIMULATE_CHR.out.others
        // each chr given true_ne but they are redudant and only one is needed
        .filter{d, slimtree, daf, truene, restart_count ->d.chrno == 1} 
        .map{d, slimtree, daf, truene, restart_count ->
            def str = "${d.label}\t${truene.toAbsolutePath()}"
        }
        .toList()
        .map {str_list -> str_list.sort().join('\n')}

    SUMMARIZE_NE (in_ch_summarize_ne, in_ch_true_ne)
}
