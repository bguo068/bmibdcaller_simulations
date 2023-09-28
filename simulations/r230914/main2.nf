nextflow.enable.dsl=2

params.test = false
params.resdir = "res"

params.min_mac = 20 // assuming homozyogous diploids, eg. min_maf 20 / 2 / 1000 = 0.01
params.mincm = 2.0
params.nchroms = 14

params.tpbwt_template_opts = 1
params.tpbwt_Lm = params.test ? 152: 100 // optimized
params.tpbwt_Lf = params.mincm
params.tpbwt_use_phase_correction = 0

params.hapibd_minoutput = params.mincm
params.hapibd_minseed = params.mincm
params.hapibd_minextend = 1.0
params.hapibd_maxgap = 1000
params.hapibd_minmarkers = 70

params.refinedibd_length = params.mincm
params.refinedibd_lod = 1.1 // TODO: confirm what does lod mean
params.refinedibd_scale = 0
params.refinedibd_window = 40.0

params.hmmibd_n = 100
params.hmmibd_m = 5

params.isorelate_imiss = 0.3
params.isorelate_vmiss = 0.3
params.isorelate_min_snp = 20 // optimized
params.isorelate_min_mac = 200 // 0.1, which is different from other callers

// params.tsinferibd_max_tmrca = [1000, 3000]

params.filt_ibd_by_ov = false
params.ibdne_mincm = params.mincm
params.ibdne_minregion = 10
params.ibdne_flatmeth = ["none", "keep_hap_1_only", "merge"][0]

params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrials = 1000
params.ifm_mincm = 2.0
params.ifm_mingwcm = 5.0


params.meta = "" // empty for simulation

// helper function to split string with ','
def to_lst(value){
    def lst = "${value}".split(",").collect{it};
    return lst
}

def resdir = params.resdir

def sp_defaults = [
    seqlen : 100 * 15000,
    selpos : Math.round(0.33 * 100 * 15000),
    num_origins : 1,
    N : params.test? 500: 10000,
    h : 0.5,
    s : 0.3,
    g_sel_start : 80,
    r : 0.01 / 15_000,
    sim_relatedness : 0,
    g_ne_change_start : 200,
    N0 : params.test? 500:1000,
    u : 1e-8,
    nsam : params.test?100: 1000, // haploid
]


def mp_ne(args) {
    def ne = args.N * args.npop + (args.npop-1) * (args.npop-1) / 4 / args.npop / args.sel_mig
    return ne
}

def get_maxtmrca(p){
    def str = p.getName()
    def match = (str =~ /maxtmrca_(\d+)/)
    if (match) {
        def x = match[0][1].toInteger()
        return x
    }
    else{
        return null
    }
}

def sp_sets = [

    sp_r0003: sp_defaults + [r: 3e-9,     s:0.0,  seqlen: (1.0/3e-9).toInteger(),     genome_set_id: 40001],
    sp_r0010: sp_defaults + [r: 1e-8,     s:0.0,  seqlen: (1.0/1e-8).toInteger(),     genome_set_id: 40002],
    sp_r0030: sp_defaults + [r: 3e-8,     s:0.0,  seqlen: (1.0/3e-8).toInteger(),     genome_set_id: 40003],
    sp_r0100: sp_defaults + [r: 1e-7,     s:0.0,  seqlen: (1.0/1e-7).toInteger(),     genome_set_id: 40004],
    sp_r0300: sp_defaults + [r: 3e-7,     s:0.0,  seqlen: (1.0/3e-7).toInteger(),     genome_set_id: 40005],
    sp_r0667: sp_defaults + [r: 6.667e-7, s:0.0,  seqlen: (1.0/6.667e-7).toInteger(), genome_set_id: 40006],
    sp_r1000: sp_defaults + [r: 1e-6,     s:0.0,  seqlen: (1.0/1e-6).toInteger(),     genome_set_id: 40007],
]


def sp_set_args_keys = sp_sets.collect{k, v -> v}[0].collect{k, v->k}


process SIM_SP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/trees/", \
        pattern: "*.trees", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/vcf/", \
        pattern: "*.vcf.gz", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/daf/", \
        pattern: "*.daf", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/restart_count/", \
        pattern: "*.restart_count", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/true_ne/", \
        pattern: "*.true_ne", mode: "symlink"

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), \
        emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), \
        emit: restart_count
    tuple val(label), val(chrno), path("*.true_ne"), emit: true_ne

    script:
    def cmd_options = (args + [chrno: chrno]).collect{k, v->
        "--${k} ${v}"}.join(" ")
    """
    sim_single_pop.py $cmd_options

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}

process CALL_IBD_HAPIBD {
    tag "${args.genome_set_id}_${chrno}_hapibd"
    publishDir "${resdir}/${label}/ibd/hapibd/", \
        mode: "symlink", pattern: "*_hapibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_hapibd.ibd")
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
        minmac: params.min_mac,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_hapibd.py $cmd_options
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hapibd.ibd"""
}

process CALL_IBD_TSKIBD {
    tag "${args.genome_set_id}_${chrno}_tskibd"
    publishDir "${resdir}/${label}/ibd/tskibd/", \
        mode: "symlink", pattern: "*_tskibd.ibd"
    input:
    tuple val (label), val(chrno), val(args), path(trees)
    output:
    tuple val(label), val(chrno), path("*_tskibd.ibd")
    script:
    def cmd_options_tskibd = [
        tree: trees,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")

    """
    call_tskibd.py ${cmd_options_tskibd}
    """
    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """touch  ${prefix}_tskibd.ibd"""
}

process CALL_IBD_HMMIBD {
    tag "${args.genome_set_id}_${chrno}_hmmibd"
    publishDir "${resdir}/${label}/ibd/hmmibd/", \
        mode: "symlink", pattern: "*_hmmibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_hmmibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        n: params.hmmibd_n,
        m: params.hmmibd_m,
        mincm: params.mincm,
        minmac: params.min_mac,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    script:
        """
        call_hmmibd.py ${cmd_options}
        """
    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd"""
}


workflow WF_SP {

    main:

    // *********************** Log Params *************************
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    if(params.test) { ch_sp_sets = ch_sp_sets.take(1) }

    Channel.value(['label'] + sp_set_args_keys)
        .concat( ch_sp_sets.map{k, args -> 
            [k] + sp_set_args_keys.collect{i -> args[i]} }
        )
        .map{lst-> lst.join("\t")}
        .collectFile(name: "sp_sets.tsv", newLine: true, storeDir: resdir, sort:false)

    // ***********************Simulation *************************

    ch_chrs = Channel.fromList(1..(params.nchroms))

    ch_sp_input = ch_sp_sets
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    SIM_SP_CHR(ch_sp_input)

    // *********************** Call ibd *************************
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf // label, chrno, trees, vcf
    // ch_true_ne = SIM_SP_CHR.out.true_ne
    //     .filter { label, chrno,true_ne-> chrno == 1}
    //     .map{label, chrno,true_ne-> [label, true_ne]}

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }
        .filter {x -> x[0] != "sp_r0003"} // this r cause very large genome

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
        .filter {x -> x[0] != "sp_r0003"} // this r cause very large genome

    // CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    // CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    // CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    // CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)

}


workflow {

    WF_SP()
}
