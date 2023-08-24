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
    sp_const: sp_defaults + [s: 0.0, N0: 10000, genome_set_id: 30005],
    sp_grow: sp_defaults +  [s: 0.0,  N0: 100000, genome_set_id: 30006],
    // sp_r0003: sp_defaults + [r: 3e-9,     s:0.0,  seqlen: (1.0/3e-9).toInteger(),     genome_set_id: 40001],
    // sp_r0010: sp_defaults + [r: 1e-8,     s:0.0,  seqlen: (1.0/1e-8).toInteger(),     genome_set_id: 40002],
    // sp_r0030: sp_defaults + [r: 3e-8,     s:0.0,  seqlen: (1.0/3e-8).toInteger(),     genome_set_id: 40003],
    // sp_r0100: sp_defaults + [r: 1e-7,     s:0.0,  seqlen: (1.0/1e-7).toInteger(),     genome_set_id: 40004],
    // sp_r0300: sp_defaults + [r: 3e-7,     s:0.0,  seqlen: (1.0/3e-7).toInteger(),     genome_set_id: 40005],
    // sp_r0667: sp_defaults + [r: 6.667e-7, s:0.0,  seqlen: (1.0/6.667e-7).toInteger(), genome_set_id: 40006],
    // sp_r1000: sp_defaults + [r: 1e-6,     s:0.0,  seqlen: (1.0/1e-6).toInteger(),     genome_set_id: 40007],
]

def mp_sets = [
    mp_s00: mp_defaults + [s:0.0, genome_set_id: 20000],
    mp_s01: mp_defaults + [s:0.1, genome_set_id: 20001],
    mp_s02: mp_defaults + [s:0.2, genome_set_id: 20002],
    mp_s03: mp_defaults + [s:0.3, genome_set_id: 20003],
    mp_rel: mp_defaults + [sim_relatedness: 1, genome_set_id: 30001],
    // mp_r0003: mp_defaults + [r: 3e-9,     s:0.0, seqlen: (1.0/3e-9).toInteger(),     genome_set_id: 50001],
    // mp_r0010: mp_defaults + [r: 1e-8,     s:0.0, seqlen: (1.0/1e-8).toInteger(),     genome_set_id: 50002],
    // mp_r0030: mp_defaults + [r: 3e-8,     s:0.0, seqlen: (1.0/3e-8).toInteger(),     genome_set_id: 50003],
    // mp_r0100: mp_defaults + [r: 1e-7,     s:0.0, seqlen: (1.0/1e-7).toInteger(),     genome_set_id: 50004],
    // mp_r0300: mp_defaults + [r: 3e-7,     s:0.0, seqlen: (1.0/3e-7).toInteger(),     genome_set_id: 50005],
    // mp_r0667: mp_defaults + [r: 6.667e-7, s:0.0, seqlen: (1.0/6.667e-7).toInteger(), genome_set_id: 50006],
    // mp_r1000: mp_defaults + [r: 1e-6,     s:0.0, seqlen: (1.0/1e-6).toInteger(),     genome_set_id: 50007],
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
        minmac: params.min_mac,
        window: params.refinedibd_window,
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
        use_phase_correction: params.tpbwt_use_phase_correction,
        minmac: params.min_mac,
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
        minmac: params.min_mac,
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
    tag "${genome_set_id}"

    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/", \
        pattern: "*.sh", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/",\
         pattern: "*.map", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/",\
         pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ibddist_ibd/",\
         pattern: "*.ibddist.ibdobj.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ibdne_ibd/",\
         pattern: "*.ibdne.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label),val(ibdcaller), path(ibd_lst), path(ibd_lst_true), path(vcf_lst), val(genome_set_id)
        path(ibdne_jar)
    output:
        tuple val(label),val(ibdcaller),  path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label),val(ibdcaller),  path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label),val(ibdcaller),  path("*.ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
        tuple val(label),val(ibdcaller),  path("*.ibdne.ibdobj.gz"), emit: ibdne_ibd_obj
    script:
    // Whether to pass in true ibd list is determined by params.filt_ibd_by_ov 
    def true_ibd_arg = (params.filt_ibd_by_ov && (ibdcaller!="tskibd")) ? \
         [ibd_files_true: "${ibd_lst_true}"] : [:]
    def args_local = (true_ibd_arg + [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_jar: ibdne_jar,
        ibdne_flatmeth: params.ibdne_flatmeth,
    ]).collect{k, v-> "--${k} ${v}"}.join(" ")
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

    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ifm_input/", \
        pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label),val(ibdcaller),  path(ibd_lst), path(vcf_lst), val(genome_set_id)

    output:
        tuple val(label),val(ibdcaller),  path("*_orig.ifm.ibdobj.gz"), \
                    emit: ifm_orig_ibd_obj
        tuple val(label),val(ibdcaller),  path("*_rmpeaks.ifm.ibdobj.gz"), \
                    emit: ifm_rmpeaks_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
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

    publishDir "${resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ne_output/",  mode: 'symlink'

    input:
        tuple val(label), val(ibdcaller), path(ibdne_jar), path(ibdne_sh), path(gmap),\
            path(ibd_gz), val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${args.genome_set_id}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${args.genome_set_id}_${are_peaks_removed}"
    publishDir "${resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ifm_output/",  mode: 'symlink'
    input:
        tuple val(label), val(ibdcaller), path(ibd_obj), val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*_member.pq")
    script:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    def args_local = [
        ibd_obj: ibd_obj,
        npop: args.npop,
        nsam: args.nsam,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
        ntrials: params.ifm_ntrials,
        transform: params.ifm_transform,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
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

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets, by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }

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
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.map{it + ["hapibd"]}
            .concat( CALL_IBD_TSKIBD.out.map{it + ["tskibd"]})
            .concat(CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]})
            .concat(CALL_IBD_TPBWT.out.map{it + ["tpbwt"]})
            .concat(CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]})
            .concat(CALL_IBD_ISORELATE.out.map{it+["isorelate"]})
            // [label, chrno, ibd, ibdcaller]
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]}) //sort by chrno
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            return ["${label}", ibdcaller, ibd_lst]
        }
    // ch_out_ibd_grp.map{it-> [it[0], it[1], it[2].size()]}.view()

    ch_true_ibd = CALL_IBD_TSKIBD.out // [label, chrno, ibd]
        .map{label, chrno, ibd-> [groupKey(label, params.nchroms), [chrno, ibd]]}
        .groupTuple(by: 0, sort: {a, b-> a[0]<=>b[0]})  // groupby label, sort by chrno
        .map{key, data_lst -> 
            def label = key.toString()
            def ibd_files_true = data_lst.collect{v -> v[1]}
            [label, ibd_files_true]
        }

    ch_out_ibd_grp = ch_out_ibd_grp
        // combine with true ibd list
        .combine(ch_true_ibd, by: 0)
        .map{label, ibdcaller, ibd_lst, ibd_list_true->
            // fix name collision for tskibd as true and inferrre ibd are the same
            [label, ibdcaller, ibd_lst, ibdcaller=="tskibd" ? [] : ibd_list_true]
         }

    ch_out_vcf_grp = SIM_SP_CHR.out.trees_vcf
        .map{label, chrno, _trees, vcf -> 
            [groupKey(label, params.nchroms), [chrno, vcf]]
            } // drop field trees 
            // groupkey                        values
        .groupTuple(by: 0, sort: {a, b -> a[0]<=> b[0]}) // sort by chrno
        .map {label, value_lst -> 
            def vcf_lst = value_lst.collect{values-> values[1]}
            [label.toString(), vcf_lst] // groupKey is not a string
        }

    ch_grouped_ibd_vcf = ch_out_ibd_grp.combine(ch_out_vcf_grp, by: 0)
        // [label, ibdcaller, ibd_lst, ibd_lst_true, vcf_lst]
        .combine(
            // add genome_set_id
            ch_sp_sets.map{label, args -> 
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
    PROC_DIST_NE( ch_grouped_ibd_vcf, file("${projectDir}/lib/ibdne.23Apr20.ae9.jar") )
     // out.ne_input_orig
     // out.ne_input_rmpeaks



    // ********************** Run IbdNe ***************************
    ch_in_ibdne = \
            PROC_DIST_NE.out.ne_input_orig.map{it + [false]}  // orig
        .concat( 
            PROC_DIST_NE.out.ne_input_rmpeaks.map {it + [true]}  // rmpeaks
        )
        // label, ibdcaller, jar, sh, map, ibd, are_peaks_removed
        .combine( ch_sp_sets, by: 0) // add args

    RUN_IBDNE(ch_in_ibdne)

    emit: 

    ch_ibdobj_dist = PROC_DIST_NE.out.ibddist_ibd_obj
    ch_ibdobj_ne   = PROC_DIST_NE.out.ibdne_ibd_obj
    ch_ibdne       = RUN_IBDNE.out
            // label, ibdcaller, are_peaks_removed, ne

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
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.map{it + ["hapibd"]}
            .concat( CALL_IBD_TSKIBD.out.map{it + ["tskibd"]})
            .concat(CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]})
            .concat(CALL_IBD_TPBWT.out.map{it + ["tpbwt"]})
            .concat(CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]})
            .concat(CALL_IBD_ISORELATE.out.map{it+["isorelate"]})
            // [label, chrno, ibd, ibdcaller]
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]}) //sort by chrno
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            return [label, ibdcaller, ibd_lst]
        }

    ch_out_vcf_grp = SIM_MP_CHR.out.trees_vcf
        .map{label, chrno, _trees, vcf -> 
            [groupKey(label, params.nchroms), [chrno, vcf]]} // drop field trees 
            // groupkey                        values
        .groupTuple(by: 0, sort: {a, b -> a[0]<=> b[0]}) // sort by chrno
        .map {label, value_lst -> 
            def vcf_lst = value_lst.collect{values-> values[1]}
            [label.toString(), vcf_lst] // groupKey is not a string
        }

    ch_grouped_ibd_vcf = ch_out_ibd_grp.combine(ch_out_vcf_grp, by: 0)
        // [label, ibdcaller, ibd_lst, vcf_lst]
        .combine(
            // add genome_set_id
            ch_mp_sets.map{label, args -> 
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
    ch_in_run_infomap = \
        PROC_INFOMAP.out.ifm_orig_ibd_obj.map{it + [false]} // orig
    .concat (
        PROC_INFOMAP.out.ifm_rmpeaks_ibd_obj.map{it + [true]} // rmpeaks
    ).combine(ch_mp_sets, by: 0)

    RUN_INFOMAP(ch_in_run_infomap)


    emit:

    ch_ifm = RUN_INFOMAP.out
            // label, ibdcaller, are_peaks_removed, member

}

workflow WF_SUMMARY {
    
}

workflow {

    WF_SP()
    WF_MP()
}

// Following demographic model used in the hap-ibd paper
process SIM_UK_CHR {
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
    // subset arguments as not all keys will be used.
    def cmd_options = (args + [chrno: chrno])\
        .findAll{k, v -> k in ["chrno", "seqlen","gc", "r", "u", "nsam", "genome_set_id"]}\
        .collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    sim_uk_pop.py $cmd_options
    mkdir tmp; mv tmp_* tmp/
    """
    stub:
    def prefix="${args.genome_set_id}_${chrno}"
    """
    touch ${prefix}{.trees,.vcf.gz,.daf,.restart_count,.true_ne}
    """
}


process CALL_IBD_HAPIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_hapibd"
    publishDir "${resdir}/${label}/ibd/hapibd/${arg_sp_dir}", \
        mode: "symlink", pattern: "*_hapibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(hapibd_args), val(arg_sp_dir)
    output:
    tuple val(label), val(chrno), path("*_hapibd.ibd"), val(hapibd_args)
    script:
        // call args is a dict with keys as follows
        // minmac: params.min_mac,
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
    ]).collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_hapibd.py $cmd_options
    """
    stub:
    def ibd_dir = String.format("%d", hapibd_args.minseed)
    """touch ${args.genome_set_id}_${chrno}_hapibd.ibd"""
}

workflow OPTIMIZE_HAPIBD {
    // single population model
    // only simulating neutral situation
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'sp_neu'}

    // multiple population model
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'mp_s00'}

    // UK model with human or Pf recombination rates
    ch_uk_sets = Channel.fromList( [
        // all demographic parameters are hard-coded in the slim script
        ["uk_gc0", [seqlen: 60_000_000, gc:0, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60001]],
        ["uk_gc1", [seqlen: 60_000_000, gc:1, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60002]],
        // the following two lines are models that using UK demograhic patterns but Pf recombination rate
        ["uk_gc0_pf", [seqlen: 900_000, gc:0, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60003]],
        ["uk_gc1_pf", [seqlen: 900_000, gc:1, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60004]],
    ])


    // combine with chrnos 
    ch_chrs = Channel.fromList(1..(params.nchroms))

    // reorder fields
    ch_sp_input = ch_sp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_mp_input = ch_mp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_uk_input = ch_uk_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    // run simulations
    SIM_SP_CHR(ch_sp_input)
    SIM_MP_CHR(ch_mp_input)
    SIM_UK_CHR(ch_uk_input)

    // prepare ibdcaller input
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf) // label, chrno, trees, vcf
        .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    // 
    ch_hapibd_args = Channel.fromList([
        // add minmac 20 (maf>=0.01) default is 2 which would include too many rare variants
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers: 100, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  90, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  80, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  70, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  60, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  50, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  40, minmac:20],
        // [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  30, minmac:20],

        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers: 100, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  30, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  10, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:   3, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers: 100, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:  30, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:  10, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:   3, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers: 100, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:  30, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:  10, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:   3, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers: 100, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:  30, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:  10, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:   3, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers: 100, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:  30, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:  10, minmac: 20],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:   3, minmac: 20],

        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers: 100, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  30, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:  10, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap: 1000,  minmarkers:   3, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers: 100, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:  30, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:  10, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  300,  minmarkers:   3, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers: 100, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:  30, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:  10, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:  100,  minmarkers:   3, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers: 100, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:  30, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:  10, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:   30,  minmarkers:   3, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers: 100, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:  30, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:  10, minmac: 200],
        [minseed: 2, minoutput: 2, minextend: 1, maxgap:    3,  minmarkers:   3, minmac: 200],
    ]).map{args->
        def arg_sp_dir = String.format("ms%d_mo%d_me%d_mg%04d_mm%03d_mac%03d", 
            args.minseed, args.minoutput, args.minextend, args.maxgap, args.minmarkers, args.minmac
        )
        [args, arg_sp_dir]
    }
    if (params.test) {
        ch_hapibd_args = ch_hapibd_args.take(1)
    }

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
    
    ch_in_ibdcall_vcf_with_params = ch_in_ibdcall_vcf 
        .combine(ch_hapibd_args) // label, chrno, args, vcf, hapibd_args, arg_sp_dir

    // call hmmibd
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_HAPIBD_PARAM(ch_in_ibdcall_vcf_with_params)
}


process CALL_IBD_TPBWT_PARAM {
    tag "${args.genome_set_id}_${chrno}_tpbwtibd"
    publishDir "${resdir}/${label}/ibd/tpbwtibd/${arg_sp_dir}", \
        mode: "symlink", pattern: "*_tpbwtibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(tpbwt_args), val(arg_sp_dir)
    output:
    tuple val(label), val(chrno), path("*_tpbwtibd.ibd")
    script:
    // Lm: params.tpbwt_Lm,
    // Lf: params.tpbwt_Lf,
    // template: params.tpbwt_template_opts,
    // minmac
    // use_phase_correction
    def cmd_options = (tpbwt_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_tpbwt.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_tpbwtibd.ibd"""
}

workflow OPTIMIZE_TPBWT {
    // single population model
    // only simulating neutral situation
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'sp_neu'}

    // multiple population model
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'mp_s00'}

    // UK model with human or Pf recombination rates
    ch_uk_sets = Channel.fromList( [
        // all demographic parameters are hard-coded in the slim script
        ["uk_gc0", [seqlen: 60_000_000, gc:0, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60001]],
        // ["uk_gc1", [seqlen: 60_000_000, gc:1, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60002]],
        // the following two lines are models that using UK demograhic patterns but Pf recombination rate
        ["uk_gc0_pf", [seqlen: 900_000, gc:0, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60003]],
        // ["uk_gc1_pf", [seqlen: 900_000, gc:1, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60004]],
    ])


    // combine with chrnos 
    ch_chrs = Channel.fromList(1..(params.nchroms))

    // reorder fields
    ch_sp_input = ch_sp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_mp_input = ch_mp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_uk_input = ch_uk_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    // run simulations
    SIM_SP_CHR(ch_sp_input)
    SIM_MP_CHR(ch_mp_input)
    SIM_UK_CHR(ch_uk_input)

    // prepare ibdcaller input
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf) // label, chrno, trees, vcf
        .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    // 
    ch_tpbwt_args = Channel.fromList([
        // vary Lm
        [Lm: 300, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm: 250, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm: 200, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm: 150, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm: 100, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm:  80, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm:  50, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        [Lm:  30, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        // [Lm:  10, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],
        // [Lm:   3, Lf: 2.0, template:1, use_phase_correction:0, minmac:20],

        // vary Lm, change template
        [Lm: 300, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm: 250, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm: 200, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm: 150, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm: 100, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm:  80, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm:  50, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        [Lm:  30, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        // [Lm:  10, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],
        // [Lm:   3, Lf: 2.0, template:0, use_phase_correction:0, minmac:20],

        // vary Lm, change minmac: 0
        [Lm: 300, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm: 250, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm: 200, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm: 150, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm: 100, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm:  80, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm:  50, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        [Lm:  30, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        // [Lm:  10, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],
        // [Lm:   3, Lf: 2.0, template:1, use_phase_correction:0, minmac: 0],

    ]).map{args->
        def arg_sp_dir = String.format("lm%d_lf%f_tp%d_pc%04d_mac%03d", 
            args.Lm, args.Lf, args.template, args.use_phase_correction, args.minmac
        )
        [args, arg_sp_dir]
    }
    if (params.test) {
        ch_tpbwt_args = ch_tpbwt_args.take(1)
    }

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
    
    ch_in_ibdcall_vcf_with_params = ch_in_ibdcall_vcf 
        .combine(ch_tpbwt_args) // label, chrno, args, vcf, hapibd_args, arg_sp_dir

    // call hmmibd
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_TPBWT_PARAM(ch_in_ibdcall_vcf_with_params)
}


process CALL_IBD_REFINEDIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_refinedibd"
    publishDir "${resdir}/${label}/ibd/refinedibd/${arg_sp_dir}", \
        mode: "symlink", pattern: "*_refinedibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(refinedibd_args), val(arg_sp_dir)

    output:
    tuple val(label), val(chrno), path("*_refinedibd.ibd")
    script:
        // lod: params.refinedibd_lod,
        // length: params.refinedibd_length,
        // scale: params.refinedibd_scale,
        // minmac
        // window
    def cmd_options = (refinedibd_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_refinedibd.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_refinedibd.ibd"""
}

workflow OPTIMIZE_REFINEDIBD {
    // single population model
    // only simulating neutral situation
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'sp_neu'}

    // multiple population model
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'mp_s00'}

    // UK model with human or Pf recombination rates
    ch_uk_sets = Channel.fromList( [
        // all demographic parameters are hard-coded in the slim script
        ["uk_gc0", [seqlen: 60_000_000, gc:0, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60001]],
        // ["uk_gc1", [seqlen: 60_000_000, gc:1, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60002]],
        // the following two lines are models that using UK demograhic patterns but Pf recombination rate
        ["uk_gc0_pf", [seqlen: 900_000, gc:0, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60003]],
        // ["uk_gc1_pf", [seqlen: 900_000, gc:1, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60004]],
    ])


    // combine with chrnos 
    ch_chrs = Channel.fromList(1..(params.nchroms))

    // reorder fields
    ch_sp_input = ch_sp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_mp_input = ch_mp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_uk_input = ch_uk_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    // run simulations
    SIM_SP_CHR(ch_sp_input)
    SIM_MP_CHR(ch_mp_input)
    SIM_UK_CHR(ch_uk_input)

    // prepare ibdcaller input
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf) // label, chrno, trees, vcf
        .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    // 
    ch_refinedibd_args = Channel.fromList([
        [lod: 8.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 40.0],
        [lod: 4.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 40.0],
        [lod: 3.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 40.0],
        [lod: 2.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 40.0],
        [lod: 1.1, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 40.0],
        // change window
        [lod: 8.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 20.0],
        [lod: 4.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 20.0],
        [lod: 3.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 20.0],
        [lod: 2.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 20.0],
        [lod: 1.1, length:2.0, scale: Math.sqrt(1000/100) , minmac: 20 , window: 20.0],
        // change scale
        [lod: 8.0, length:2.0, scale: Math.sqrt(500/100) , minmac: 20 , window: 40.0],
        [lod: 4.0, length:2.0, scale: Math.sqrt(500/100) , minmac: 20 , window: 40.0],
        [lod: 3.0, length:2.0, scale: Math.sqrt(500/100) , minmac: 20 , window: 40.0],
        [lod: 2.0, length:2.0, scale: Math.sqrt(500/100) , minmac: 20 , window: 40.0],
        [lod: 1.1, length:2.0, scale: Math.sqrt(500/100) , minmac: 20 , window: 40.0],
        // change minmac
        [lod: 8.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 2 , window: 40.0],
        [lod: 4.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 2 , window: 40.0],
        [lod: 3.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 2 , window: 40.0],
        [lod: 2.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 2 , window: 40.0],
        [lod: 1.1, length:2.0, scale: Math.sqrt(1000/100) , minmac: 2 , window: 40.0],

        [lod: 8.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 200, window: 40.0],
        [lod: 4.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 200, window: 40.0],
        [lod: 3.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 200, window: 40.0],
        [lod: 2.0, length:2.0, scale: Math.sqrt(1000/100) , minmac: 200, window: 40.0],
        [lod: 1.1, length:2.0, scale: Math.sqrt(1000/100) , minmac: 200, window: 40.0],
    ]).map{args->
        def arg_sp_dir = String.format("lod%f_len%f_scale%f_minmac%d_win%f", 
            args.lod, args.length, args.scale, args.minmac, args.window
        )
        [args, arg_sp_dir]
    }
    if (params.test) {
        ch_refinedibd_args = ch_refinedibd_args.take(1)
    }

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
    
    ch_in_ibdcall_vcf_with_params = ch_in_ibdcall_vcf 
        .combine(ch_refinedibd_args) // label, chrno, args, vcf, hapibd_args, arg_sp_dir

    // call hmmibd
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD_PARAM(ch_in_ibdcall_vcf_with_params)
}



process CALL_IBD_HMMIBD_PARAM {
    tag "${args.genome_set_id}_${chrno}_hmmibd"
    publishDir "${resdir}/${label}/ibd/hmmibd/${arg_sp_dir}", \
        mode: "symlink", pattern: "*_hmmibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf),  val(hmmibd_args), val(arg_sp_dir)
    output:
    tuple val(label), val(chrno), path("*_hmmibd.ibd")
    script:
        // n: params.hmmibd_n,
        // m: params.hmmibd_m,
        // minmac:
    def cmd_options = (hmmibd_args+[
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
    ]).collect{k, v -> "--${k} ${v}"}.join(" ")
    script:
        """
        call_hmmibd.py ${cmd_options}
        """
    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd"""
}


workflow OPTIMIZE_HMMIBD {
    // single population model
    // only simulating neutral situation
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'sp_neu'}

    // multiple population model
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'mp_s00'}

    // UK model with human or Pf recombination rates
    ch_uk_sets = Channel.fromList( [
        // all demographic parameters are hard-coded in the slim script

        // NOTE: comment out the uk-human genome. it might be too big for hmmibd as it cause (memory issue) to run isorelate
        // ["uk_gc0", [seqlen: 60_000_000, gc:0, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60001]],

        // ["uk_gc1", [seqlen: 60_000_000, gc:1, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60002]],
        // the following two lines are models that using UK demograhic patterns but Pf recombination rate
        ["uk_gc0_pf", [seqlen: 900_000, gc:0, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60003]],
        // ["uk_gc1_pf", [seqlen: 900_000, gc:1, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60004]],
    ])


    // combine with chrnos 
    ch_chrs = Channel.fromList(1..(params.nchroms))

    // reorder fields
    ch_sp_input = ch_sp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_mp_input = ch_mp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_uk_input = ch_uk_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    // run simulations
    SIM_SP_CHR(ch_sp_input)
    SIM_MP_CHR(ch_mp_input)
    SIM_UK_CHR(ch_uk_input)

    // prepare ibdcaller input
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf) // label, chrno, trees, vcf
        .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    // 
    ch_hmmibd_args = Channel.fromList([
        [m:  2, n:  10, minmac: 1],
        [m:  2, n:  30, minmac: 1],
        [m:  2, n: 100, minmac: 1],
        [m:  2, n: 300, minmac: 1],
        [m:  2, n: 9999999, minmac: 1],
        [m:  5, n:  10, minmac: 1],
        [m:  5, n:  30, minmac: 1],
        [m:  5, n: 100, minmac: 1],
        [m:  5, n: 300, minmac: 1],
        [m:  5, n: 9999999, minmac: 1],
        [m: 10, n:  10, minmac: 1],
        [m: 10, n:  30, minmac: 1],
        [m: 10, n: 100, minmac: 1],
        [m: 10, n: 300, minmac: 1],
        [m: 10, n: 9999999, minmac: 1],

        [m:  2, n:  10, minmac: 20],
        [m:  2, n:  30, minmac: 20],
        [m:  2, n: 100, minmac: 20],
        [m:  2, n: 300, minmac: 20],
        [m:  2, n: 9999999, minmac: 20],
        [m:  5, n:  10, minmac: 20],
        [m:  5, n:  30, minmac: 20],
        [m:  5, n: 100, minmac: 20],
        [m:  5, n: 300, minmac: 20],
        [m:  5, n: 9999999, minmac: 20],
        [m: 10, n:  10, minmac: 20],
        [m: 10, n:  30, minmac: 20],
        [m: 10, n: 100, minmac: 20],
        [m: 10, n: 300, minmac: 20],
        [m: 10, n: 9999999, minmac: 20],

        [m:  2, n:  10, minmac: 200],
        [m:  2, n:  30, minmac: 200],
        [m:  2, n: 100, minmac: 200],
        [m:  2, n: 300, minmac: 200],
        [m:  2, n: 9999999, minmac: 200],
        [m:  5, n:  10, minmac: 200],
        [m:  5, n:  30, minmac: 200],
        [m:  5, n: 100, minmac: 200],
        [m:  5, n: 300, minmac: 200],
        [m:  5, n: 9999999, minmac: 200],
        [m: 10, n:  10, minmac: 200],
        [m: 10, n:  30, minmac: 200],
        [m: 10, n: 100, minmac: 200],
        [m: 10, n: 300, minmac: 200],
        [m: 10, n: 9999999, minmac: 200],
    ]).map{args->
        def arg_sp_dir = String.format("m%d_n%d_minmac%d", 
            args.m, args.n, args.minmac
        )
        [args, arg_sp_dir]
    }
    if (params.test) {
        ch_hmmibd_args = ch_hmmibd_args.take(1)
    }

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
    
    ch_in_ibdcall_vcf_with_params = ch_in_ibdcall_vcf 
        .combine(ch_hmmibd_args) // label, chrno, args, vcf, hapibd_args, arg_sp_dir

    // call hmmibd
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_HMMIBD_PARAM(ch_in_ibdcall_vcf_with_params)
}

process CALL_IBD_ISORELATE_PARAM {
    tag "${args.genome_set_id}_${chrno}_isorelate"
    publishDir "${resdir}/${label}/ibd/isorelate/${arg_sp_dir}", \
        mode: "symlink", pattern: "*_isorelate.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf), val(isorelate_args), val(arg_sp_dir)
    output:
    tuple val(label), val(chrno), path("*_isorelate.ibd")
    script:
        // maf: params.maf,
        // min_snp: params.isorelate_min_snp,
    def cmd_options = (isorelate_args + [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        min_len_bp: Math.round(params.mincm * (0.01/args.r)),
        imiss: params.isorelate_imiss,
        vmiss: params.isorelate_vmiss,
        cpus: task.cpus,
        genome_set_id: args.genome_set_id,
    ]).collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_isorelate.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_isorelate.ibd"""
}

workflow OPTIMIZE_ISORELATE {
    // single population model
    // only simulating neutral situation
    ch_sp_sets = Channel.fromList(sp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'sp_neu'}

    // multiple population model
    ch_mp_sets = Channel.fromList(mp_sets.collect {label, args->[label, args]})
    .filter{label, args -> label == 'mp_s00'}

    // UK model with human or Pf recombination rates
    ch_uk_sets = Channel.fromList( [
        // all demographic parameters are hard-coded in the slim script

        // NOTE: uk-human model is too big for isorelate (memory error)
        // ["uk_gc0", [seqlen: 60_000_000, gc:0, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60001]],
        // ["uk_gc1", [seqlen: 60_000_000, gc:1, r: 1e-8, u: 1.38e-8, nsam: 1000, genome_set_id: 60002]],
        // the following two lines are models that using UK demograhic patterns but Pf recombination rate
        ["uk_gc0_pf", [seqlen: 900_000, gc:0, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60003]],
        // ["uk_gc1_pf", [seqlen: 900_000, gc:1, r: 6.6666667e-7, u: 1e-8, nsam: 1000, genome_set_id: 60004]],
    ])


    // combine with chrnos 
    ch_chrs = Channel.fromList(1..(params.nchroms))

    // reorder fields
    ch_sp_input = ch_sp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_mp_input = ch_mp_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    ch_uk_input = ch_uk_sets // label, dict (args)
        .combine( ch_chrs )
        .map {label, args, chrno -> [label, chrno, args]}

    // run simulations
    SIM_SP_CHR(ch_sp_input)
    SIM_MP_CHR(ch_mp_input)
    SIM_UK_CHR(ch_uk_input)

    // prepare ibdcaller input
    ch_trees_vcf = SIM_SP_CHR.out.trees_vcf
        .concat(SIM_MP_CHR.out.trees_vcf) // label, chrno, trees, vcf
        .concat(SIM_UK_CHR.out.trees_vcf) // label, chrno, trees, vcf

    // 
    ch_isorelate_args = Channel.fromList([
        [minac: 0, min_snp: 1],
        [minac: 0, min_snp: 3],
        [minac: 0, min_snp: 10],
        [minac: 0, min_snp: 15],
        [minac: 0, min_snp: 20],
        [minac: 0, min_snp: 40],
        [minac: 0, min_snp: 80],
        [minac: 0, min_snp:160],

        [minac: 20, min_snp: 1],
        [minac: 20, min_snp: 3],
        [minac: 20, min_snp: 10],
        [minac: 20, min_snp: 15],
        [minac: 20, min_snp: 20],
        [minac: 20, min_snp: 40],
        [minac: 20, min_snp: 80],
        [minac: 20, min_snp:160],

        [minac: 60, min_snp: 1],
        [minac: 60, min_snp: 3],
        [minac: 60, min_snp: 10],
        [minac: 60, min_snp: 15],
        [minac: 60, min_snp: 20],
        [minac: 60, min_snp: 40],
        [minac: 60, min_snp: 80],
        [minac: 60, min_snp:160],

        [minac: 200, min_snp: 1],
        [minac: 200, min_snp: 3],
        [minac: 200, min_snp: 10],
        [minac: 200, min_snp: 15],
        [minac: 200, min_snp: 20],
        [minac: 200, min_snp: 40],
        [minac: 200, min_snp: 80],
        [minac: 200, min_snp:160],

        [minmac: 600, min_snp: 1],
        [minmac: 600, min_snp: 3],
        [minmac: 600, min_snp: 10],
        [minmac: 600, min_snp: 15],
        [minmac: 600, min_snp: 20],
        [minmac: 600, min_snp: 40],
        [minmac: 600, min_snp: 80],
        [minmac: 600, min_snp:160],
    ]).map{args->
        def nsam = 1000
        def arg_sp_dir = String.format("maf%f_minsnp%d", 
            args.minmac * 2 / nsam , args.min_snp
        )
        [args, arg_sp_dir]
    }
    if (params.test) {
        ch_isorelate_args = ch_isorelate_args.take(1)
    }

    ch_in_ibdcall_trees = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, trees] }

    ch_in_ibdcall_vcf = ch_trees_vcf
        .combine(ch_sp_sets.concat(ch_mp_sets).concat(ch_uk_sets), by:0)
        .map{label, chrno, trees, vcf, args -> [label, chrno, args, vcf] }
    
    ch_in_ibdcall_vcf_with_params = ch_in_ibdcall_vcf 
        .combine(ch_isorelate_args) // label, chrno, args, vcf, hapibd_args, arg_sp_dir

    // call hmmibd
    CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_ISORELATE_PARAM(ch_in_ibdcall_vcf_with_params)
}

