import os
import multiprocessing
import numpy as np

def generate_seeds(seed, num_seeds):
    seed_seq = np.random.SeedSequence(seed)
    seeds = seed_seq.generate_state(num_seeds)
    return seeds

def run_boost2pandas(bstr_idx, seed):
    os.system(f"python boost2pandas.py --bstr_idx {bstr_idx} --seed {seed}")
    return

if __name__ == '__main__':
    

    seed = 260292
    n_bstrp = 10

    seeds = generate_seeds(seed, 400)

    # for i in range(n_bstrp):
    #     print(i,seeds[i])
    #     os.system(f"python boost2pandas.py --bstr_idx {i} --seed {seeds[i]}")
    for j in range(40):
        input_data = [(i, seeds[i]) for i in range(j*n_bstrp,(j+1)*n_bstrp)]
        pool = multiprocessing.Pool()
        pool.starmap(run_boost2pandas, input_data)

        pool.close()
        pool.join()
