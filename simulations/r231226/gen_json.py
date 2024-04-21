import json
from pathlib import Path


def gen_hapibd_coarse_lst():
    default_dict = dict(
        minseed=2, minoutput=2, minextend=1, maxgap=1000, minmarkers=100, minmaf=0.01
    )
    coarse_lst = [
        dict(maxgap=1000, minmarkers=100),
        dict(maxgap=1000, minmarkers=30),
        dict(maxgap=1000, minmarkers=10),
        dict(maxgap=1000, minmarkers=3),
        dict(maxgap=300, minmarkers=100),
        dict(maxgap=300, minmarkers=30),
        dict(maxgap=300, minmarkers=10),
        dict(maxgap=300, minmarkers=3),
        dict(maxgap=100, minmarkers=100),
        dict(maxgap=100, minmarkers=30),
        dict(maxgap=100, minmarkers=10),
        dict(maxgap=100, minmarkers=3),
        dict(maxgap=30, minmarkers=100),
        dict(maxgap=30, minmarkers=30),
        dict(maxgap=30, minmarkers=10),
        dict(maxgap=30, minmarkers=3),
        dict(maxgap=3, minmarkers=100),
        dict(maxgap=3, minmarkers=30),
        dict(maxgap=3, minmarkers=10),
        dict(maxgap=3, minmarkers=3),
    ]
    # Python 3.9 introduced the merge operator( | ) in the dict class. Using
    # the merge operator is the easiest way to merge dictionaries. The merge
    # operator returns a new dictionary, leaving the original dictionaries
    # unchanged
    full_lst_coarse = [
        ["hapibd", f"coarse_mg{d['maxgap']}_mm{d['minmarkers']}", default_dict | d]
        for d in coarse_lst
    ]
    Path("hapibd_coarse.json").write_text(json.dumps(full_lst_coarse, indent=4))


def gen_hapibd_fine_lst():
    default_dict = dict(
        minseed=2, minoutput=2, minextend=1, maxgap=1000, minmarkers=100, minmaf=0.01
    )
    fine_lst = [
        dict(maxgap=1000, minmarkers=100),
        dict(maxgap=1000, minmarkers=90),
        dict(maxgap=1000, minmarkers=80),
        dict(maxgap=1000, minmarkers=70),
        dict(maxgap=1000, minmarkers=60),
        dict(maxgap=1000, minmarkers=50),
        dict(maxgap=1000, minmarkers=40),
        dict(maxgap=1000, minmarkers=30),
    ]
    # Python 3.9 introduced the merge operator( | ) in the dict class. Using
    # the merge operator is the easiest way to merge dictionaries. The merge
    # operator returns a new dictionary, leaving the original dictionaries
    # unchanged
    full_lst_fine = [
        ["hapibd", f"fine_mg{d['maxgap']}_mm{d['minmarkers']}", default_dict | d]
        for d in fine_lst
    ]
    Path("hapibd_fine.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_hmmibd_origalg_lst():
    default_dict = dict(m=5, n=9999999, minmaf=0.01)
    fine_lst = [
        dict(m=2, n=10, minmaf=0.01),
        dict(m=2, n=30, minmaf=0.01),
        dict(m=2, n=100, minmaf=0.01),
        dict(m=2, n=300, minmaf=0.01),
        dict(m=2, n=9999999, minmaf=0.01),
        dict(m=5, n=10, minmaf=0.01),
        dict(m=5, n=30, minmaf=0.01),
        dict(m=5, n=100, minmaf=0.01),
        dict(m=5, n=300, minmaf=0.01),
        dict(m=5, n=9999999, minmaf=0.01),
        dict(m=10, n=10, minmaf=0.01),
        dict(m=10, n=30, minmaf=0.01),
        dict(m=10, n=100, minmaf=0.01),
        dict(m=10, n=300, minmaf=0.01),
        dict(m=10, n=9999999, minmaf=0.01),
    ]
    # Python 3.9 introduced the merge operator( | ) in the dict class. Using
    # the merge operator is the easiest way to merge dictionaries. The merge
    # operator returns a new dictionary, leaving the original dictionaries
    # unchanged
    full_lst_fine = [
        ["hmmibd", f"origalg_m{d['m']}_n{d['n']}", default_dict | d] for d in fine_lst
    ]
    Path("hmmibd_origalg.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_hmmibd_origalgmaf001_lst():
    default_dict = dict(m=5, n=9999999, minmaf=0.01)
    fine_lst = [
        dict(m=2, n=10, minmaf=0.001),
        dict(m=2, n=30, minmaf=0.001),
        dict(m=2, n=100, minmaf=0.001),
        dict(m=2, n=300, minmaf=0.001),
        dict(m=2, n=9999999, minmaf=0.001),
        dict(m=5, n=10, minmaf=0.001),
        dict(m=5, n=30, minmaf=0.001),
        dict(m=5, n=100, minmaf=0.001),
        dict(m=5, n=300, minmaf=0.001),
        dict(m=5, n=9999999, minmaf=0.001),
        dict(m=10, n=10, minmaf=0.001),
        dict(m=10, n=30, minmaf=0.001),
        dict(m=10, n=100, minmaf=0.001),
        dict(m=10, n=300, minmaf=0.001),
        dict(m=10, n=9999999, minmaf=0.001),
    ]
    # Python 3.9 introduced the merge operator( | ) in the dict class. Using
    # the merge operator is the easiest way to merge dictionaries. The merge
    # operator returns a new dictionary, leaving the original dictionaries
    # unchanged
    full_lst_fine = [
        ["hmmibd", f"origalgmaf001_m{d['m']}_n{d['n']}", default_dict | d]
        for d in fine_lst
    ]
    Path("hmmibd_origalgmaf001.json").write_text(json.dumps(full_lst_fine, indent=4))
    # NOTES:
    # under the sp model, include rare variants does not seems increase much of
    # of the marker count, for example, one chrom has 29k when maf>0, and 25k
    # when maf >0.01. And using rare variants seems to slightly increase error rates


def gen_isorelate_grpa_lst():
    default_dict = dict(minmaf=0.01, min_snp=20)
    fine_lst = [
        dict(minmaf=0.01, min_snp=1),
        dict(minmaf=0.01, min_snp=3),
        dict(minmaf=0.01, min_snp=10),
        dict(minmaf=0.01, min_snp=15),
        dict(minmaf=0.01, min_snp=20),
        dict(minmaf=0.01, min_snp=40),
        dict(minmaf=0.01, min_snp=80),
        dict(minmaf=0.01, min_snp=160),
        dict(minmaf=0.03, min_snp=1),
        dict(minmaf=0.03, min_snp=3),
        dict(minmaf=0.03, min_snp=10),
        dict(minmaf=0.03, min_snp=15),
        dict(minmaf=0.03, min_snp=20),
        dict(minmaf=0.03, min_snp=40),
        dict(minmaf=0.03, min_snp=80),
        dict(minmaf=0.03, min_snp=160),
        dict(minmaf=0.1, min_snp=1),
        dict(minmaf=0.1, min_snp=3),
        dict(minmaf=0.1, min_snp=10),
        dict(minmaf=0.1, min_snp=15),
        dict(minmaf=0.1, min_snp=20),
        dict(minmaf=0.1, min_snp=40),
        dict(minmaf=0.1, min_snp=80),
        dict(minmaf=0.1, min_snp=160),
    ]
    full_lst_fine = [
        ["isorelate", f"grpa_mm{d['minmaf']}_ms{d['min_snp']}", default_dict | d]
        for d in fine_lst
    ]
    Path("isorelate_grpa.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_refinedibd_win40_lst():
    import math

    default_dict = dict(
        lod=4.0, length=2.0, scale=math.sqrt(1000 / 100), minmaf=0.01, window=40.0
    )
    fine_lst = [
        dict(lod=8.0, minmaf=0.01),
        dict(lod=4.0, minmaf=0.01),
        dict(lod=3.0, minmaf=0.01),
        dict(lod=2.0, minmaf=0.01),
        dict(lod=1.8, minmaf=0.01),
        dict(lod=1.6, minmaf=0.01),
        dict(lod=1.4, minmaf=0.01),
        dict(lod=1.2, minmaf=0.01),
        dict(lod=1.1, minmaf=0.01),
        #
        dict(lod=8.0, minmaf=0.1),
        dict(lod=4.0, minmaf=0.1),
        dict(lod=3.0, minmaf=0.1),
        dict(lod=2.0, minmaf=0.1),
        dict(lod=1.8, minmaf=0.1),
        dict(lod=1.6, minmaf=0.1),
        dict(lod=1.4, minmaf=0.1),
        dict(lod=1.2, minmaf=0.1),
        dict(lod=1.1, minmaf=0.1),
    ]
    full_lst_fine = [
        ["refinedibd", f"win40_lod{d['lod']}_mm{d['minmaf']}", default_dict | d]
        for d in fine_lst
    ]
    Path("refinedibd_win40.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_refinedibd_win20_lst():
    import math

    default_dict = dict(
        lod=4.0, length=2.0, scale=math.sqrt(1000 / 100), minmaf=0.01, window=40.0
    )
    fine_lst = [
        dict(lod=8.0, minmaf=0.01, window=20),
        dict(lod=4.0, minmaf=0.01, window=20),
        dict(lod=3.0, minmaf=0.01, window=20),
        dict(lod=2.0, minmaf=0.01, window=20),
        dict(lod=1.8, minmaf=0.01, window=20),
        dict(lod=1.6, minmaf=0.01, window=20),
        dict(lod=1.4, minmaf=0.01, window=20),
        dict(lod=1.2, minmaf=0.01, window=20),
        dict(lod=1.1, minmaf=0.01, window=20),
        #
        dict(lod=8.0, minmaf=0.1, window=20),
        dict(lod=4.0, minmaf=0.1, window=20),
        dict(lod=3.0, minmaf=0.1, window=20),
        dict(lod=2.0, minmaf=0.1, window=20),
        dict(lod=1.8, minmaf=0.1, window=20),
        dict(lod=1.6, minmaf=0.1, window=20),
        dict(lod=1.4, minmaf=0.1, window=20),
        dict(lod=1.2, minmaf=0.1, window=20),
        dict(lod=1.1, minmaf=0.1, window=20),
    ]
    full_lst_fine = [
        ["refinedibd", f"win20_lod{d['lod']}_mm{d['minmaf']}", default_dict | d]
        for d in fine_lst
    ]
    Path("refinedibd_win20.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_refinedibd_trim_minmaf_lst():
    import math

    default_dict = dict(
        lod=4.0,
        length=2.0,
        scale=math.sqrt(1000 / 100),
        minmaf=0.01,
        window=40.0,
        trim=0.15,
    )
    fine_lst = [
        dict(trim=0.15, minmaf=0.01),
        dict(trim=0.12, minmaf=0.01),
        dict(trim=0.10, minmaf=0.01),
        dict(trim=0.08, minmaf=0.01),
        dict(trim=0.05, minmaf=0.01),
        dict(trim=0.02, minmaf=0.01),
        dict(trim=0.01, minmaf=0.01),
        #
        dict(trim=0.15, minmaf=0.1),
        dict(trim=0.12, minmaf=0.1),
        dict(trim=0.10, minmaf=0.1),
        dict(trim=0.08, minmaf=0.1),
        dict(trim=0.05, minmaf=0.1),
        dict(trim=0.02, minmaf=0.1),
        dict(trim=0.01, minmaf=0.1),
    ]
    full_lst_fine = [
        ["refinedibd", f"trimminmaf_trim{d['trim']}_mm{d['minmaf']}", default_dict | d]
        for d in fine_lst
    ]
    Path("refinedibd_trimminmaf.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_tpbwt_tp0_lst():
    default_dict = dict(Lm=300, Lf=2.0, template=0, use_phase_correction=0, minmaf=0.01)
    fine_lst = [
        dict(Lm=300, template=0, minmaf=0.001),
        dict(Lm=250, template=0, minmaf=0.001),
        dict(Lm=200, template=0, minmaf=0.001),
        dict(Lm=150, template=0, minmaf=0.001),
        dict(Lm=130, template=0, minmaf=0.001),
        dict(Lm=110, template=0, minmaf=0.001),
        dict(Lm=100, template=0, minmaf=0.001),
        dict(Lm=90, template=0, minmaf=0.001),
        dict(Lm=80, template=0, minmaf=0.001),
        dict(Lm=50, template=0, minmaf=0.001),
        ##
        dict(Lm=300, template=0, minmaf=0.01),
        dict(Lm=250, template=0, minmaf=0.01),
        dict(Lm=200, template=0, minmaf=0.01),
        dict(Lm=150, template=0, minmaf=0.01),
        dict(Lm=130, template=0, minmaf=0.01),
        dict(Lm=110, template=0, minmaf=0.01),
        dict(Lm=100, template=0, minmaf=0.01),
        dict(Lm=90, template=0, minmaf=0.01),
        dict(Lm=80, template=0, minmaf=0.01),
        dict(Lm=50, template=0, minmaf=0.01),
        ##
        dict(Lm=300, template=0, minmaf=0.1),
        dict(Lm=250, template=0, minmaf=0.1),
        dict(Lm=200, template=0, minmaf=0.1),
        dict(Lm=150, template=0, minmaf=0.1),
        dict(Lm=130, template=0, minmaf=0.1),
        dict(Lm=110, template=0, minmaf=0.1),
        dict(Lm=100, template=0, minmaf=0.1),
        dict(Lm=90, template=0, minmaf=0.1),
        dict(Lm=80, template=0, minmaf=0.1),
        dict(Lm=50, template=0, minmaf=0.1),
    ]
    full_lst_fine = [
        ["tpbwt", f"tp0_lm{d['Lm']}_mm{d['minmaf']}", default_dict | d]
        for d in fine_lst
    ]
    Path("tpbwt_tp0.json").write_text(json.dumps(full_lst_fine, indent=4))


def gen_tpbwt_tp1_lst():
    default_dict = dict(Lm=300, Lf=2.0, template=0, use_phase_correction=0, minmaf=0.01)
    fine_lst = [
        dict(Lm=300, template=1, minmaf=0.001),
        dict(Lm=250, template=1, minmaf=0.001),
        dict(Lm=200, template=1, minmaf=0.001),
        dict(Lm=150, template=1, minmaf=0.001),
        dict(Lm=130, template=1, minmaf=0.001),
        dict(Lm=110, template=1, minmaf=0.001),
        dict(Lm=100, template=1, minmaf=0.001),
        dict(Lm=90, template=1, minmaf=0.001),
        dict(Lm=80, template=1, minmaf=0.001),
        dict(Lm=50, template=1, minmaf=0.001),
        ##
        dict(Lm=300, template=1, minmaf=0.01),
        dict(Lm=250, template=1, minmaf=0.01),
        dict(Lm=200, template=1, minmaf=0.01),
        dict(Lm=150, template=1, minmaf=0.01),
        dict(Lm=130, template=1, minmaf=0.01),
        dict(Lm=110, template=1, minmaf=0.01),
        dict(Lm=100, template=1, minmaf=0.01),
        dict(Lm=90, template=1, minmaf=0.01),
        dict(Lm=80, template=1, minmaf=0.01),
        dict(Lm=50, template=1, minmaf=0.01),
        ##
        dict(Lm=300, template=1, minmaf=0.1),
        dict(Lm=250, template=1, minmaf=0.1),
        dict(Lm=200, template=1, minmaf=0.1),
        dict(Lm=150, template=1, minmaf=0.1),
        dict(Lm=130, template=1, minmaf=0.1),
        dict(Lm=110, template=1, minmaf=0.1),
        dict(Lm=100, template=1, minmaf=0.1),
        dict(Lm=90, template=1, minmaf=0.1),
        dict(Lm=80, template=1, minmaf=0.1),
        dict(Lm=50, template=1, minmaf=0.1),
    ]
    full_lst_fine = [
        ["tpbwt", f"tp1_lm{d['Lm']}_mm{d['minmaf']}", default_dict | d]
        for d in fine_lst
    ]
    Path("tpbwt_tp1.json").write_text(json.dumps(full_lst_fine, indent=4))


if __name__ == "__main__":
    gen_hapibd_coarse_lst()
    gen_hapibd_fine_lst()
    gen_hmmibd_origalg_lst()
    gen_hmmibd_origalgmaf001_lst()
    gen_isorelate_grpa_lst()
    gen_refinedibd_win20_lst()
    gen_refinedibd_win40_lst()
    gen_refinedibd_trim_minmaf_lst()
    gen_tpbwt_tp0_lst()
    gen_tpbwt_tp1_lst()
