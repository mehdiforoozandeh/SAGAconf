import os, time

def scan_old_files(subdir, t_cutoff):
    """returns list of files in subdir that are older than t_cutoff"""
    oldfiles = []

    now = time.time()
    l_files = os.listdir(subdir)

    for filename in l_files:
        filestamp = os.stat(os.path.join(subdir, filename)).st_mtime
        filecompare = now - t_cutoff

        if  filestamp < filecompare:
            oldfiles.append("{}/{}".format(subdir, filename))

    return oldfiles

def main():
    t_1hr_cutoff = 3600
    main_dir = "segway112samples/"
    list_of_running_samples = [s for s in os.listdir(main_dir) if os.path.exists("{}/{}/call_segway_train".format(main_dir, s))]

    for s in list_of_running_samples:
        to_remove = []
        for c in os.listdir("{}/{}/call_segway_train/cmdline".format(main_dir, s)):
            to_remove += scan_old_files("{}/{}/call_segway_train/cmdline/{}".format(main_dir, s, c), 3*t_1hr_cutoff)
        
        for c in os.listdir("{}/{}/call_segway_train/output/o".format(main_dir, s)):
            to_remove += scan_old_files("{}/{}/call_segway_train/output/o/{}".format(main_dir, s, c), 3*t_1hr_cutoff)
        
        for c in os.listdir("{}/{}/call_segway_train/output/e".format(main_dir, s)):
            to_remove += scan_old_files("{}/{}/call_segway_train/output/e/{}".format(main_dir, s, c), 3*t_1hr_cutoff)
        
        for f in to_remove:
            os.remove(f)

if __name__ == "__main__":
    main_dir = "segway112samples/"
    
    # while there are running jobs
    while len([s for s in os.listdir(main_dir) if os.path.exists("{}/{}/call_segway_train".format(main_dir, s))])>0:
        main()
        # every 30min
        time.sleep(1800)
