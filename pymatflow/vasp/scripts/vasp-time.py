import argparse
import datetime
import matplotlib.pyplot as plt


def get_info(outcar):
    """
    get the timing information from OUTCAR of vasp run
    :return out a dict -> {
        "job_done": bool,
        "start_time": xxx,
        "total_cpu_time": float,
        "elapsed_time": float,
    }
    """
    with open(outcar, 'r') as fin:
        lines = fin.readlines()
    # check whether calculation is finished
    if len(lines[-1].split()) == 4 and lines[-1].split()[0] == "Voluntary" and lines[-1].split()[1] == "context":
        job_done = True
    else:
        job_done = False

    out = {}
    out["job_done"] = job_done

    for line in lines:
        # if it is an empty line continue to next line
        if len(line.split()) == 0:
            continue
        if line.split()[0] == "executed" and line.split()[1] == "on" and line.split()[3] == "date":
            out["start_time"] = line.split("\n")[0]
        if line.split()[0] == "Total" and line.split()[1] == "CPU" and line.split()[2] == "time":
            out["total_cpu_time"] = float(line.split()[5]) # in unit of second
        if line.split()[0] == "Elapsed" and line.split()[1] == "time":
            out["elapsed_time"] = float(line.split()[3])

    return out




def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--outcar", type=str, default="./OUTCAR",
        help="OUTCAR of vasp run")
    
    # ==========================================================
    args = parser.parse_args()            

    info = get_info(outcar=args.outcar)

    # calculate the running time and print it out
    # Importante: the length of the time string might be different, depending
    # on the value of hours and minutes and seconds. if they are two digits
    # number, they will be divided like: '11: 6: 2', only when they all are
    # two digtis number, they will not be divided '11:16:12'
    # so we have to preprocess it to build the right time string to pass into
    # datetime.datetime.strptime()
    start_str = info["start_time"].split()[4]+"-"+info["start_time"].split()[5]

    start = datetime.datetime.strptime(start_str, "%Y.%m.%d-%H:%M:%S")
    

    print("---------------------------\n")
    print("Time analysis for vasp run\n")
    print("---------------------------\n")    
    if info["job_done"] == False:
        print("Job is not finished yet!\n")
        return
    print("**Time consuming:**\n")
    print("- job starts at %s\n" % start)
    print("- CPU time:\n")
    print("  - in unit of sec  :  %.3f\n" % (info["total_cpu_time"]))
    print("  - in unit of min  :  %.3f\n" % (info["total_cpu_time"]/60))
    print("  - in unit of hour :  %.3f\n" % (info["total_cpu_time"]/3600))
    print("- Elapsed time:\n")
    print("  - in unit of sec  :  %.3f\n" % (info["elapsed_time"]))
    print("  - in unit of min  :  %.3f\n" % (info["elapsed_time"]/60))
    print("  - in unit of hour :  %.3f\n" % (info["elapsed_time"]/3600))


if __name__ == "__main__":
    main()