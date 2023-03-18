nextflow.enable.dsl=2

include { SIMULATE_CHR }                    from './main.nf'
include { SIMULATE_CHR_MUTATIONS }          from './main.nf'
include { TSINFER_TSDATE_PER_CHR }          from './main.nf'
include { CALC_RAW_IBD }                    from './main.nf'
include { CALL_IBD_HAPIBD }                 from './main.nf'
include { CALL_IBD_REFINEDIBD }             from './main.nf'
include { CALL_IBD_TPBWT }                  from './main.nf'
include { CALL_IBD_HMMIBD }                 from './main.nf'
include { CALL_IBD_ISORELATE }              from './main.nf'
include { CMP_TRUE_AND_INFERRED_RAW_IBD }   from './main.nf'
include { FIND_PEAKS }                      from './main.nf'
include { PROC_IBD }                        from './main.nf'
include { CALL_IBDNE }                      from './main.nf'
include { SUMMARIZE_TREES }                 from './main.nf'
include { SUMMARIZE_IBD }                   from './main.nf'
include { SUMMARIZE_NE }                    from './main.nf'


workflow TEST_SIMULATE_CHR {
    def d =[
        num_chr : 1,
        label : "label",
        chrno : 1,
        s : 0.6,
        h : 0.5,
        selpos_0_1 : 0.33,
        num_origins : 1,
        s_start_g : 50,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10,
        ne_change_start_g : 200,
        N : 10000,
        N0 : 10000,
        nsam : 500
    ]
    channel.value(d)
    | SIMULATE_CHR
}

workflow TEST_SIMULATE_CHR_MUTATIONS {
    def d = [ chrno : 1, u : 1e-8, label: 'label']
    ch_tree= channel.fromPath("${projectDir}/testdata/sim_tree/1.trees")
    channel.value(d).combine(ch_tree) 
    | SIMULATE_CHR_MUTATIONS
}

workflow TEST_TSINFER_TSDATE_PER_CHR {
    def d = [label: 'label', chrno: 1, bp_per_cm: 15000]
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1.vcf.gz")
    ch_true_ne = channel.fromPath("${projectDir}/testdata/sim_tree/1.true_ne")
    channel.value(d).combine(ch_vcf).combine(ch_true_ne)
    | TSINFER_TSDATE_PER_CHR
}

workflow TEST_CALC_RAW_IBD {
    def d = [label: 'label', src: 'sim_tree', chrno: 1, bp_per_cm: 15000]
    ch_tree= channel.fromPath("${projectDir}/testdata/sim_tree/1.trees")
    channel.value(d).combine(ch_tree) 
    | CALC_RAW_IBD
}

workflow TEST_CALL_IBD_HAPIBD {
    def d = [
        label: 'label', ibdcaller: 'hapibd', chrno: 1, bp_per_cm: 15000,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10]
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1.vcf.gz")
    channel.value(d).combine(ch_vcf)
    | CALL_IBD_HAPIBD
}

workflow TEST_CALL_IBD_REFINEDIBD {
    def d = [
        label: 'label', ibdcaller: 'refinedibd', chrno: 1, bp_per_cm: 15000,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10]
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1.vcf.gz")
    channel.value(d).combine(ch_vcf)
    | CALL_IBD_REFINEDIBD
}

workflow TEST_CALL_IBD_TPBWT {
    def d = [
        label: 'label', ibdcaller: 'refinedibd', chrno: 1, bp_per_cm: 15000,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10]
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1.vcf.gz")
    channel.value(d).combine(ch_vcf)
    | CALL_IBD_TPBWT
}

workflow TEST_CALL_IBD_HMMIBD {
    def d = [
        label: 'label', ibdcaller: 'refinedibd', chrno: 1, bp_per_cm: 15000,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10]
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1.vcf.gz")
    channel.value(d).combine(ch_vcf)
    | CALL_IBD_HMMIBD
}

workflow TEST_CALL_IBD_ISORELATE {
    def d = [
        label: 'label', ibdcaller: 'refinedibd', chrno: 1, bp_per_cm: 15000,
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10]
    // isorelate seems slow, using a subset of samples  (20)
    ch_vcf = channel.fromPath("${projectDir}/testdata/sim_tree/1_20samples.vcf.gz")
    channel.value(d).combine(ch_vcf)
    | CALL_IBD_ISORELATE
}

workflow TEST_CMP_TRUE_AND_INFERRED_RAW_IBD {
    def d =[
        num_chr : 1,
        label : "label",
        chrno : 1,
        src: 'infer',
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10,
        min_tmrca: 2
    ]
    ch_true_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/1.ibd")
    ch_infer_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/1.ibd")
    channel.value(d).combine(ch_true_ibd).combine(ch_infer_ibd)
    | CMP_TRUE_AND_INFERRED_RAW_IBD
}

workflow TEST_FIND_PEAKS {
    def d =[
        num_chr : 1,
        label : "label",
        chrno : 1,
        src: 'infer',
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10,
        min_tmrca: 2
    ]
    ch_true_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/1.ibd")
    channel.value(d).combine(ch_true_ibd)
    | FIND_PEAKS
}

workflow TEST_PROC_IBD {
    def d =[
        label : "label",
        proc_label: 'autocut_hapibd',
        chrno : 1,
        src: 'truetree',
        bp_per_cm : 15000, // set to pf values
        seqlen_in_cm : 10,
        ibd_tmrca_cutoff: 100,
        ibd_cut_mode: 'autocut_hapibd', // TODO: change it to 'autocut_hapibd'
        cut_start: 0,
        cut_end: 0,
    ]
    ch_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/1.ibd")
    ch_hapibd_peaks = channel.fromPath("${projectDir}/testdata/sim_tree/1_peaks.bed")
    channel.value(d).combine(ch_ibd).combine(ch_hapibd_peaks)
    | PROC_IBD
}

workflow TEST_CALL_IBDNE {
    def d =[
        num_chr: 2,
        label : "label",
        proc_label: 'autocut_hapibd',
        chrno : 1,
        src: 'truetree',
        minregion: 3
    ]
    ch_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/*ibd_proc").collect()
    ch_map = channel.fromPath("${projectDir}/testdata/sim_tree/*map_proc").collect()
    channel.value(d).combine(ch_ibd.map{[it]}).combine(ch_map.map{[it]})
    | CALL_IBDNE
}

workflow TEST_SUMMARIZE_IBD {
    def sim_label_val = 'label'
    ch_ibd = channel.fromPath("${projectDir}/testdata/sim_tree/1.ibd")
    ch_src = channel.fromList(['src1', 'src2'])
    ch_chrno = channel.fromList(1..3)
    channel.value(sim_label_val).combine(ch_src).combine(ch_chrno).combine(ch_ibd)
    .map {simlabel, src, chrno, ibd->
        def str = "${simlabel}\t${src}\t${chrno}\t${ibd.toAbsolutePath()}"
        return [simlabel, str] }
    .groupTuple()
    .map{simlabel, str_list -> [simlabel, str_list.sort().join('\n')] }
    | SUMMARIZE_IBD
}

workflow TEST_SUMMARIZE_TREES {
    def sim_label_val = 'label'
    ch_tree = channel.fromPath("${projectDir}/testdata/sim_tree/1.trees")
    ch_src = channel.fromList(['src1', 'src2'])
    ch_chrno = channel.fromList(1..3)
    channel.value(sim_label_val).combine(ch_src).combine(ch_chrno).combine(ch_tree)
    .map {simlabel, src, chrno, tree->
        def str = "${simlabel}\t${src}\t${chrno}\t${tree.toAbsolutePath()}"
        return [simlabel, str] }
    .groupTuple()
    .map{simlabel, str_list -> [simlabel, str_list.sort().join('\n')] }
    | SUMMARIZE_TREES
}

workflow TEST_SUMMARIZE_NE {
    def sim_label_val = 'label'
    ch_src = channel.fromList(['src1', 'src2'])
    ch_proc = channel.fromList(['proc1', 'proc2'])
    ch_ne = channel.fromPath("${projectDir}/testdata/sim_tree/ne_res.ne")
    ch_true_ne = channel.fromPath("${projectDir}/testdata/sim_tree/1.true_ne")


    ch_nes = channel.value(sim_label_val).combine(ch_src).combine(ch_proc).combine(ch_ne)
    .map {simlabel, src, proc, ne->
        def str = "${simlabel}\t${src}\t${proc}\t${ne.toAbsolutePath()}"
        return [simlabel, str]} 
    .groupTuple()
    .map {label, str_list -> [label, str_list.sort().join('\n')]}

    ch_true_ne = channel.value(sim_label_val).combine(ch_true_ne)
    .map {simlabel, true_ne -> "${simlabel}\t${true_ne.toAbsolutePath()}"}

    SUMMARIZE_NE(ch_nes, ch_true_ne)
}

workflow {
    // run all test workflow
    TEST_SIMULATE_CHR()
    TEST_SIMULATE_CHR_MUTATIONS()
    TEST_TSINFER_TSDATE_PER_CHR()
    TEST_CALC_RAW_IBD()
    TEST_CALL_IBD_HAPIBD()
    TEST_CALL_IBD_REFINEDIBD()
    TEST_CALL_IBD_TPBWT()
    TEST_CALL_IBD_HMMIBD()
    TEST_CALL_IBD_ISORELATE()
    TEST_CMP_TRUE_AND_INFERRED_RAW_IBD()
    TEST_FIND_PEAKS()
    TEST_PROC_IBD()
    TEST_CALL_IBDNE()
    TEST_SUMMARIZE_IBD()
    TEST_SUMMARIZE_TREES()
    TEST_SUMMARIZE_NE()
}