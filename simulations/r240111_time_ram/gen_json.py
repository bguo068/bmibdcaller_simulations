import json
from pathlib import Path

SP_DEFAULTS = dict(
    seqlen=100 * 15000,
    selpos=round(0.33 * 100 * 15000),
    num_origins=1,
    N=10000,
    h=0.5,
    s=0.3,
    g_sel_start=80,
    r=0.01 / 15_000,
    sim_relatedness=0,
    g_ne_change_start=200,
    N0=1000,
    u=1e-8,
    nsam=1000,  # // haploid
)


def write_json_for_test():
    d = dict(
        sp_neu1=SP_DEFAULTS | dict(s=0.0, nsam=30, genome_set_id=20000),
    )
    Path("test.json").write_text(json.dumps(d, indent=2))


def write_json_for_demographic_models_small_size():
    d = dict(
        sp_neu1=SP_DEFAULTS | dict(s=0.0, nsam=30, genome_set_id=20000),
        sp_neu2=SP_DEFAULTS | dict(s=0.0, nsam=100, genome_set_id=30000),
        sp_neu3=SP_DEFAULTS | dict(s=0.0, nsam=300, genome_set_id=40000),
        sp_neu4=SP_DEFAULTS | dict(s=0.0, nsam=1000, genome_set_id=50000),
    )
    Path("demog_small_size.json").write_text(json.dumps(d, indent=2))


def write_json_for_demographic_models_small_size2():
    d = dict(
        sp_neu1a=SP_DEFAULTS | dict(s=0.0, nsam=30, genome_set_id=20001),
        sp_neu2a=SP_DEFAULTS | dict(s=0.0, nsam=100, genome_set_id=30001),
        sp_neu3a=SP_DEFAULTS | dict(s=0.0, nsam=300, genome_set_id=40001),
        sp_neu4a=SP_DEFAULTS | dict(s=0.0, nsam=1000, genome_set_id=50001),
    )
    Path("demog_small_size2.json").write_text(json.dumps(d, indent=2))


def write_json_for_demographic_models_large_size():
    # note the different N0 compared to small size model
    d = dict(
        sp_neu5=SP_DEFAULTS | dict(s=0.0, N0=100000, nsam=1000, genome_set_id=60000),
        sp_neu6=SP_DEFAULTS | dict(s=0.0, N0=100000, nsam=3000, genome_set_id=70000),
        sp_neu7=SP_DEFAULTS | dict(s=0.0, N0=100000, nsam=10000, genome_set_id=80000),
        sp_neu8=SP_DEFAULTS | dict(s=0.0, N0=100000, nsam=30000, genome_set_id=90000),
    )
    Path("demog_large_size.json").write_text(json.dumps(d, indent=2))


if __name__ == "__main__":
    write_json_for_test()
    write_json_for_demographic_models_large_size()
    write_json_for_demographic_models_small_size()
    write_json_for_demographic_models_small_size2()
