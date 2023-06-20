nextflow.enable.dsl=2

params.maf = params.test ? 0.00001: 0.01
params.mincm = 2.0
params.num_neutral_chrs = 0
params.nchroms = 14

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

params.hmmibd_n = 100
params.hmmibd_m = 5

params.isorelate_imiss = 0.3
params.isorelate_vmiss = 0.3
params.isorelate_min_snp = 100

params.tsinferibd_max_tmrca = [1000, 3000]

params.ibdne_mincm = params.mincm
params.ibdne_minregion = 10


params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrails = 1000
params.ifm_mincm = 2.0
params.ifm_mingwcm = 5.0

params.test = false
params.resdir = "res"

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

def mp_defaults = [
    seqlen : 100 * 15000,
    selpos : Math.round(100 * 15000 * 0.33),
    num_origins : 1,
    N : params.test? 500: 10000,
    h : 0.5,
    s : 0.3,
    g_sel_start : 80,
    r : 0.01 / 15_000,
    sim_relatedness : 0,
    mig : 1e-5,
    sel_mig : 0.01,
    npop : 5,
    nsam : params.test? 40: 200, // haploid
    Tsplit : 500,
    u : 1e-8,
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
    sp_neu: sp_defaults + [s: 0.0, genome_set_id: 10000],
    sp_s01: sp_defaults + [s: 0.1, genome_set_id: 10001],
    sp_s02: sp_defaults + [s: 0.2, genome_set_id: 10002],
    sp_s03: sp_defaults + [s: 0.3, genome_set_id: 10003],
    sp_g040:sp_defaults + [g_sel_start: 40, genome_set_id: 10004],
    sp_g080:sp_defaults + [g_sel_start: 80, genome_set_id: 10005],
    sp_g120:sp_defaults + [g_sel_start:120, genome_set_id: 10006],
    sp_o01: sp_defaults + [num_origins: 1, genome_set_id: 10007],
    sp_o03: sp_defaults + [num_origins: 3, genome_set_id: 10008],
    sp_o27: sp_defaults + [num_origins: 27, genome_set_id: 10009],
    sp_rel: sp_defaults + [sim_relatedness: 1, genome_set_id: 30000],
    sp_r0003: sp_defaults + [r: 3e-9, genome_set_id: 40001],
    sp_r0010: sp_defaults + [r: 1e-8, genome_set_id: 40002],
    sp_r0030: sp_defaults + [r: 3e-8, genome_set_id: 40003],
    sp_r0100: sp_defaults + [r: 1e-7, genome_set_id: 40004],
    sp_r0300: sp_defaults + [r: 3e-7, genome_set_id: 40005],
    sp_r0667: sp_defaults + [r: 6.667e-7, genome_set_id: 40006],
    sp_r1000: sp_defaults + [r: 1e-6, genome_set_id: 40007],
]

def mp_sets = [
    mp_s00: mp_defaults + [s:0.0, genome_set_id: 20000],
    mp_s01: mp_defaults + [s:0.1, genome_set_id: 20001],
    mp_s02: mp_defaults + [s:0.2, genome_set_id: 20002],
    mp_s03: mp_defaults + [s:0.3, genome_set_id: 20003],
    mp_rel: mp_defaults + [sim_relatedness: 1, genome_set_id: 30001],
    mp_r0003: mp_defaults + [r: 3e-9, genome_set_id: 50001],
    mp_r0010: mp_defaults + [r: 1e-8, genome_set_id: 50002],
    mp_r0030: mp_defaults + [r: 3e-8, genome_set_id: 50003],
    mp_r0100: mp_defaults + [r: 1e-7, genome_set_id: 50004],
    mp_r0300: mp_defaults + [r: 3e-7, genome_set_id: 50005],
    mp_r0667: mp_defaults + [r: 6.667e-7, genome_set_id: 50006],
    mp_r1000: mp_defaults + [r: 1e-6, genome_set_id: 50007],
]

def sp_set_args_keys = sp_sets.collect{k, v -> v}[0].collect{k, v->k}
def mp_set_args_keys = mp_sets.collect{k, v -> v}[0].collect{k, v->k}


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

process SIM_MP_CHR {
    tag "${args.genome_set_id}_${chrno}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/trees/", \
        pattern: "*.trees", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/vcf/",\
        pattern: "*.vcf.gz", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/restart_count/",\
        pattern: "*.restart_count", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/daf/",\
        pattern: "*.daf", mode: "symlink"
    publishDir "${resdir}/${args.genome_set_id}_${label}/demog/",\
        pattern: "*_demog.png", mode: "symlink"

    input:
    tuple val(label), val(chrno), val(args)

    output:
    tuple val(label), val(chrno), path("*.trees"), path("*.vcf.gz"), \
        emit: trees_vcf
    tuple val(label), val(chrno), path("*.daf"), emit: daf
    tuple val(label), val(chrno), path("*.restart_count"), emit: restart_count
    tuple val(label), val(chrno), path("*_demog.png"), emit: demog

    script:
    def cmd_options = (args + [chrno: chrno]).collect{k, v->
        "--${k} ${v}"}.join(" ")
    """
    sim_multiple_pop.py $cmd_options

    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,_demog.png}
    """
}

process CALL_IBD_TSINFERIBD {
    tag  "${args.genome_set_id}_${chrno}_tsinferibd"
    publishDir "${resdir}/${label}/ibd/tsinferibd/", \
        mode: "symlink", pattern: "*_tsinferibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf), path(true_ne_df)
    output:
    tuple val(label), val(chrno), path("*_tsinferibd.ibd"), emit: ibd
    tuple val(label), val(chrno), path("*_tsinferibd.ibd_maxtmrca_*"), emit: ibd_extra
    tuple val(label), val(chrno), path("*_tsinfer.trees"), emit: trees
    script:
    def cmd_options_tsinfer = [
        vcf_files: vcf,
        r: args.r,
        u: args.u,
        ne: true_ne_df ?: mp_ne(args),
    ].collect{k,v -> "--${k} $v"}.join(" ")
    // Note passing [] to true_ne signals multiple-pop simulation.
    // A Ne value is caluated from simulation args, see function `mp_ne`
    // defined above the process definitions
    def max_tmrmca_str =  params.tsinferibd_max_tmrca.join(" ")
    def cmd_options_tskibd = [
        tree: "${chrno}.trees",
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
        max_tmrca: max_tmrmca_str
    ].collect{k,v -> "--${k} $v"}.join(" ")
    def prefix = "${args.genome_set_id}_${chrno}"
    // first infer trees, then call IBD from inferred trees
    """
    tsinfer_tsdate_gw.py ${cmd_options_tsinfer} # output ${chrno}.trees
    call_tskibd.py ${cmd_options_tskibd}

    mv ${chrno}.trees ${prefix}_tsinfer.trees
    mv ${prefix}_tskibd.ibd ${prefix}_tsinferibd.ibd

    # extra files
    if [ -n "${max_tmrmca_str}" ] ; then
        for tmrca in ${max_tmrmca_str}; do
            mv ${prefix}_tskibd.ibd_maxtmrca_\${tmrca} \
                ${prefix}_tsinferibd.ibd_maxtmrca_\${tmrca}
        done
    fi
    """

    stub:
    def max_tmrmca_str =  params.tsinferibd_max_tmrca.join(" ")
    def prefix = "${args.genome_set_id}_${chrno}"
    """touch  ${prefix}_tsinferibd.ibd ${prefix}_tsinfer.trees

    # extra files
    if [ -n "${max_tmrmca_str}" ] ; then
        for tmrca in ${max_tmrmca_str}; do
            touch  ${prefix}_tsinferibd.ibd_maxtmrca_\${tmrca}
        done
    fi

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
process CALL_IBD_REFINEDIBD {
    tag "${args.genome_set_id}_${chrno}_refinedibd"
    publishDir "${resdir}/${label}/ibd/refinedibd/", \
        mode: "symlink", pattern: "*_refinedibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_refinedibd.ibd")
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
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_refinedibd.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_refinedibd.ibd"""
}

process CALL_IBD_TPBWT {
    tag "${args.genome_set_id}_${chrno}_tpbwtibd"
    publishDir "${resdir}/${label}/ibd/tpbwtibd/", \
        mode: "symlink", pattern: "*_tpbwtibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_tpbwtibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        template: params.tpbwt_template_opts,
        Lm: params.tpbwt_Lm,
        Lf: params.tpbwt_Lf,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_tpbwt.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_tpbwtibd.ibd"""
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
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    script:
        """
        call_hmmibd.py ${cmd_options}
        """
    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd"""
}

process CALL_IBD_ISORELATE {
    tag "${args.genome_set_id}_${chrno}_isorelate"
    publishDir "${resdir}/${label}/ibd/isorelate/", \
        mode: "symlink", pattern: "*_isorelate.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_isorelate.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        min_snp: params.isorelate_min_snp,
        min_len_bp: Math.round(params.mincm * (0.01/args.r)),
        maf: params.maf,
        imiss: params.isorelate_imiss,
        vmiss: params.isorelate_vmiss,
        cpus: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_isorelate.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_isorelate.ibd"""
}

process CMP_TRUE_AND_INFERRED_IBD {
    tag "${args.genome_set_id}_${ibdcaller}_cmpibd"
    publishDir "${resdir}/${label}/ibdcmp/", \
        mode: "symlink", pattern: "*ibdcmpobj.gz"

    input:
    tuple val(label), val(ibdcaller), path(inferred_ibd_lst), \
          path(true_ibd_lst, stageAs: "trueibd??.ibd"), val(args)
    output:
    tuple val(label), val(ibdcaller), path("*.ibdcmpobj.gz")

    script:
    def cmd_options = [
        true_ibd_lst: true_ibd_lst,
        inferred_ibd_lst: inferred_ibd_lst,
        r: args.r,
        seqlen: args.seqlen,
        nchroms: params.nchroms,
        ibdcaller: ibdcaller,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")

    """
    cmp_ibd.py ${cmd_options}
    """
    stub:
    """
    touch ${args.genome_set_id}_${ibdcaller}.ibdcmpobj.gz
    """
}

process PROC_DIST_NE {
    tag "${args.genome_set_id}_${ibdcaller}"

    publishDir "${resdir}/${label}/ne_input/${ibdcaller}", \
            pattern: "*.sh", mode: "symlink"
    publishDir "${resdir}/${label}/ne_input/${ibdcaller}", \
            pattern: "*.map", mode: "symlink"
    publishDir "${resdir}/${label}/ne_input/${ibdcaller}", \
            pattern: "*.ibd.gz", mode: "symlink"
    publishDir "${resdir}/${label}/ibddist_ibd/${ibdcaller}", \
            pattern: "*_ibddist_ibd.pq", mode: "symlink"
    publishDir "${resdir}/${label}/ibdcov/${ibdcaller}",\
            pattern: "*.cov.pq", mode: "symlink"

    input:
        tuple val(label), val(ibdcaller), path(ibd_lst), val(args)
    output:
        tuple val(label), val(ibdcaller), path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label), val(ibdcaller), path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label), val(ibdcaller), path("*_ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
        tuple val(label), val(ibdcaller), path("*.ibdcov.ibdobj.gz"), emit: cov_ibd_obj
        tuple val(label), val(ibdcaller), path("*.ibdne.ibdobj.gz"), emit: ne_ibd_obj
    script:
    def cmd_options = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        genome_set_id: args.genome_set_id,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        r: args.r,
        seqlen: args.seqlen,
        nchroms: params.nchroms,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_dist_ne.py ${cmd_options}
    """
    stub:
    """
    touch ibdne.jar
    touch ${args.genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${args.genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${args.genome_set_id}_ibddist_ibd.pq
    touch ${args.genome_set_id}_ibddist.ibdobj.gz
    touch ${args.genome_set_id}_orig_all.ibdcov.ibdobj.gz
    touch ${args.genome_set_id}_orig_unrel.ibdcov.ibdobj.gz
    touch ${args.genome_set_id}_orig.ibdne.ibdobj.gz
    touch ${args.genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
}

process PROC_INFOMAP {
    tag  "${args.genome_set_id}_${ibdcaller}"

    publishDir "${resdir}/${label}/ifm_input/${ibdcaller}", \
        pattern: "*.ibdobj.gz", mode: "symlink"

    input:
        tuple val(label), val(ibdcaller), path(ibd_lst), val(args)
    output:
        tuple val(label), val(ibdcaller), path("*_ifm_orig.ibdobj.gz"),\
            emit: ifm_orig_ibd_obj
        tuple val(label), val(ibdcaller), path("*_ifm_rmpeaks.ibdobj.gz"),\
            emit: ifm_rmpeaks_ibd_obj
    script:
    def cmd_options = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        genome_set_id: args.genome_set_id,
        r: args.r,
        seqlen: args.seqlen,
        nchroms: params.nchroms,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_infomap.py ${cmd_options}
    """
    stub:
    """
    touch ${args.genome_set_id}{_ifm_orig.ibdobj.gz,_ifm_rmpeaks.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag  "${args.genome_set_id}_${ibdcaller}_${are_peaks_removed}"

    publishDir "${resdir}/${label}/ne_output/${ibdcaller}", mode: "symlink"

    input:
    tuple val(label), val(ibdcaller), \
        path(ibdne_jar), path(ibdne_sh), path(gmap), path(ibd_gz), \
        val(are_peaks_removed), val(args)
    output:
    tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${label}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag  "${args.genome_set_id}_${ibdcaller}_${are_peaks_removed}"

    publishDir "${resdir}/${label}/ifm_output/${ibdcaller}",  mode: "symlink"
    input:
    tuple val(label), val(ibdcaller), path(ibd_obj), val(are_peaks_removed), \
            val(args)
    output:
    tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*_member.pq")
    script:
    def meta = params.meta ? file(params.meta) : ""
    def cut_mode = are_peaks_removed? "rmpeaks": "orig"
    def cmd_options = [
        ibd_obj: ibd_obj,
        meta: meta,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
        ntrails: params.ifm_ntrails,
        transform: params.ifm_transform,
        ifm_mincm: params.ifm_mincm,
        ifm_mingwcm: params.ifm_mingwcm,
    ].collect{k, v-> v ? "--${k} ${v}": " "}.join(" ")
    """
    run_infomap.py ${cmd_options}
    """
    stub:
    def cut_mode = are_peaks_removed? "rmpeaks": "orig"
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
}


/*
process SUMMARIZE_TREES {
    tag "treeseq_summary_${label}"
    publishDir "${params.outdir}/${label}/s6_summarize_treeseq", mode: "copy"
    publishDir "${params.outdir}/summary/summarize_TreeSeq", mode: "copy"
    input:
    tuple val(label), stdin
    output:
    tuple path("*.png"), path("*.tsv")
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
    publishDir "${params.outdir}/${label}/s7_summarize_ibd", mode: "copy"
    publishDir "${params.outdir}/summary/summarize_IBD", mode: "copy"
    input:
    tuple val(label), stdin
    output:
    tuple path("*.png"), path("*.tsv")
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
    publishDir "${params.outdir}/${label}/s8_summarize_ne", mode: "copy"
    publishDir "${params.outdir}/summary/summarize_NE", mode: "copy"
    input:
    tuple val(label), stdin
    file(true_ne_tsv)
    output:
    tuple path("*.png"), path("*.tsv")
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
*/

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
    ch_true_ne = SIM_SP_CHR.out.true_ne
        .filter { label, chrno,true_ne-> chrno == 1}
        .map{label, chrno,true_ne-> [label, true_ne]}

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }

    ch_in_ibdcall_vcf_with_ne = ch_in_ibdcall_vcf
        .combine(ch_true_ne, by: 0)
        // label, chrno, args, vcf, true_ne

    CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)


    // ********************** Compare ibd ***************************


    ch_out_ibd_grp = \
            CALL_IBD_TSINFERIBD.out.ibd.map{it + ["tsinferibd"]}
        .concat(
            // parse extra_files and get maxtrmca and append to ibdcaller name
            CALL_IBD_TSINFERIBD.out.ibd_extra
                .flatMap{ label, chrno, extra_ibd ->
                    // Note: extra_ibd can a path or a list of Path
                    def lst = (extra_ibd instanceof List) ? extra_ibd: [extra_ibd]
                    return lst.collect{ibd ->
                        def maxtmrca = get_maxtmrca(ibd)
                        def ibdcaller = "tskinferibd${maxtmrca}"
                        [label, chrno, ibd, ibdcaller]
                    }
            },
            CALL_IBD_HAPIBD.out.map{it + ["hapibd"]},
            CALL_IBD_TSKIBD.out.map{it + ["tskibd"]},
            CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]},
            CALL_IBD_TPBWT.out.map{it + ["tpbwt"]},
            CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]},
            CALL_IBD_ISORELATE.out.map{it+["isorelate"]}
        )
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]})
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            return [label, ibdcaller, ibd_lst]
        }


    ch_trueibd = ch_out_ibd_grp
        .filter{it[1] == "tskibd"}.map{it-> [it[0], it[2]]} // label, trueibd_list

    ch_in_cmpibd = ch_out_ibd_grp.combine(ch_trueibd, by: 0).combine(ch_sp_sets, by:0)
    // label, ibdcaller, ibd_lst, trueibd_lst, args

    CMP_TRUE_AND_INFERRED_IBD(ch_in_cmpibd)



    // ********************** Process ibd ***************************
    PROC_DIST_NE( ch_out_ibd_grp.combine(ch_sp_sets, by: 0) )



    // ********************** Run IbdNe ***************************
    ch_in_ibdne = PROC_DIST_NE.out.ne_input_orig.map{it + [false]}
        .concat( PROC_DIST_NE.out.ne_input_rmpeaks.map {it + [true]} )
        // label, ibdcaller, jar, sh, map, ibd, are_peaks_removed
        .combine( ch_sp_sets, by: 0) // add args

    RUN_IBDNE(ch_in_ibdne)

    emit: 

    ch_ibdcov = PROC_DIST_NE.out.cov_ibd_obj.combine(ch_sp_sets, by: 0) 
            // label, ibdcaller, ibdobj, args
    ch_ibdcmp = CMP_TRUE_AND_INFERRED_IBD.out.combine(ch_sp_sets, by: 0) 
            // label, ibdcaller, ibdcmpobj, args
    ch_ibdne = RUN_IBDNE.out.combine(ch_sp_sets, by: 0)
            // label, ibdcaller, are_peaks_removed, ne, args,
        

}

workflow WF_MP {

    main:

    // *********************** Log params ************************
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    if(params.test) { ch_mp_sets = ch_mp_sets.take(1) }

    Channel.value(['label'] + mp_set_args_keys)
        .concat( ch_mp_sets.map{k, args -> 
            [k] + mp_set_args_keys.collect{i -> args[i]} }
        )
        .map{lst-> lst.join("\t")}
        .collectFile(name: "mp_sets.tsv", newLine: true, storeDir: resdir, sort:false)
        

    // ***********************Simulation *************************

    ch_chrs = Channel.fromList(1..(params.nchroms))

    ch_mp_input = ch_mp_sets
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    SIM_MP_CHR(ch_mp_input)


    // *********************** Call ibd *************************
    ch_trees_vcf = SIM_MP_CHR.out.trees_vcf // label, chrno, trees, vcf

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_mp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_mp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }

    ch_in_ibdcall_vcf_with_ne = ch_in_ibdcall_vcf
        .map{it + [ [] ] }
        // label, chrno, args, vcf, true_ne (empty)
        //
        //                            ^
        //                            |
        // Note pass [] to true_ne, this signals multiple-pop simulation,
        // a Ne value is caluated from simulation args, see function `mp_ne`
        // defined above the process definitions

    CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)


    // ********************** Compare ibd ***************************
    ch_out_ibd_grp = \
            CALL_IBD_TSINFERIBD.out.ibd.map{it + ["tsinferibd"]}
        .concat(
            // parse extra_files and get maxtrmca and append to ibdcaller name
            CALL_IBD_TSINFERIBD.out.ibd_extra
                .flatMap{label, chrno, extra_ibd ->
                    // Note: extra_ibd can be a path or a list of Path
                    def lst = (extra_ibd instanceof List) ? extra_ibd: [extra_ibd]
                    lst.collect{ibd ->
                        def maxtmrca = get_maxtmrca(ibd)
                        def ibdcaller = "tskinferibd${maxtmrca}"
                        [label, chrno, ibd, ibdcaller]
                    }
            },
            CALL_IBD_HAPIBD.out.map{it + ["hapibd"]},
            CALL_IBD_TSKIBD.out.map{it + ["tskibd"]},
            CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]},
            CALL_IBD_TPBWT.out.map{it + ["tpbwt"]},
            CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]},
            CALL_IBD_ISORELATE.out.map{it+["isorelate"]}
        )
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]})
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            return [label, ibdcaller, ibd_lst]
        }


    ch_trueibd = ch_out_ibd_grp
        .filter{it[1] == "tskibd"}.map{it-> [it[0], it[2]]} // label, trueibd_list

    ch_in_cmpibd = ch_out_ibd_grp.combine(ch_trueibd, by: 0).combine(ch_mp_sets, by:0)
    // label, ibdcaller, ibd_lst, trueibd_lst, args

    CMP_TRUE_AND_INFERRED_IBD(ch_in_cmpibd)



    // ********************** Process ibd ***************************

    PROC_INFOMAP( ch_out_ibd_grp.combine(ch_mp_sets, by: 0) )


    // ********************** Run IbdNe ***************************
    ch_in_run_infomap = PROC_INFOMAP.out.ifm_orig_ibd_obj.map{it + [false]}.concat (
        PROC_INFOMAP.out.ifm_orig_ibd_obj.map{it + [true]}
    ).combine(ch_mp_sets, by: 0)

    RUN_INFOMAP(ch_in_run_infomap)




    emit:

    ch_ibdcmp = CMP_TRUE_AND_INFERRED_IBD.out.combine(ch_mp_sets, by: 0) 
            // label, ibdcaller, ibdcmpobj, args
    ch_ifm = RUN_INFOMAP.out.combine(ch_mp_sets, by: 0)
            // label, ibdcaller, are_peaks_removed, member, args,

}

workflow WF_SUMMARY {
    
}

workflow {

    WF_SP()
    WF_MP()
}
