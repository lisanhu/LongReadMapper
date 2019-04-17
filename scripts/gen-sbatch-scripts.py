#! /usr/bin/python3
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def gen_script(opts, job_lines='', file=sys.stdout):
    content = """\
#!/bin/sh
{header}
{job_header}\

{job_lines}
{job_footer}\
"""
    header = ""
    for k in opts:
        tpl, val = opts[k]
        if val is not None:
            tpl = "#SBATCH " + tpl
            header += tpl.format(val) + "\n"
        else:
            continue

    job_header = 'echo "Job started"'
    job_footer = 'echo "Job ended"'
    content = content.format(header = header, job_header = job_header, job_footer=job_footer, job_lines=job_lines)
    print(content, file=file)


def build_opts(opts=None):
    opts = opts.copy()
    if opts is None:
        opts = {
        "name": "test",
        "nodes": 1,
        "ntasks": 1,
        "cpus": 32,
        "email": "lisanhu@udel.edu",
        "mail_type": "END,FAIL",
        "mem": "200G",
        "time": "20:00:00",
        }

    tpls = {
    "name": "-J {}",
    "nodes": "-N {}",
    "ntasks": "-n {}",
    "cpus": "-c {}",
    "email": "--mail-user={}",
    "mail_type": "--mail-type={}",
    "mem": "--mem={}",
    "time": "--time={}",
    "ofile": "-o {}",
    }

    to_be_del = []

    for k in opts:
        if opts[k] is None:
            to_be_del.append(k)
        else:
            opts[k] = (tpls[k], opts[k])

    for k in to_be_del:
        del opts[k]
    return opts


def gen_jobs():
    import itertools as it

    batch_sizes = [1000000]
    # seed_lens = [20]
    # threses = [100,200]
    seed_lens = list(range(12,33))
    threses = list(range(100,1001,100))

    # print(len(list(it.product(batch_sizes, seed_lens, threses)))) 210 jobs
    lines = """\
export BATCH_SIZE={b}
export SEED_LEN={sl}
export THRES={t}
export ACC_NUM_CORES={c}
export COMMAND={cmd}
export FA_FILE={fa}
export FQ_FILES={fqs}
export SAM_PREFIX=${{SEED_LEN}}_${{THRES}}_${{BATCH_SIZE}}_${{SLURM_JOBID}}

${{COMMAND}} ${{FA_FILE}} ${{FQ_FILES}} ${{BATCH_SIZE}} ${{SEED_LEN}} ${{THRES}} > ${{SAM_PREFIX}}.sam\
"""
    opts = {
    "name": None,
    "nodes": None,
    "ntasks": None,
    "cpus": 36,
    "email": "hcarter@udel.edu",
    "mail_type": "END,FAIL",
    "mem": "200G",
    "time": "20:00:00",
    }
    c = 0
    for b, sl, t in it.product(batch_sizes, seed_lens, threses):
        opts = opts.copy()
        opts["ofile"] = f"{sl}_{t}_{b}_%j.out"
        # print()
        fn = f"{sl}_{t}_{b}.job"
        fp = open(fn, "w")
        l = lines.format(b=b, sl=sl, t=t, c=opts["cpus"], cmd="./accaln", fa="ucsc.hg19.fa", fqs="ucsc.hg19_single.fq")
        gen_script(build_opts(opts), job_lines=l, file=fp)
        fp.close()


if __name__ == '__main__':
    gen_jobs()
