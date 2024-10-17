"""
Module containing functions for absolute copy number fitting
Created: 27/05/2024
Python 3.9.7
Carolin Sauer
"""
#!/usr/bin/env python3

#----
# 1. Import modules
#----
import numpy as np
from scipy.stats import gaussian_kde
import statistics
import math
import copy

#----
# 2. Functions for purity estimation
#----
def is_constant(x):
    return np.all(x == x[0])


def Modes(x, min_size=0.1):
    '''
    The Modes function is a simple, deterministic function that differences the kernel density of x and reports a number of modes equal to half the number of changes in direction, 
    although the min.size function can be used to reduce the number of modes returned, and defaults to 0.1, eliminating modes that do not have at least 10% of the distributional area. 
    The Modes function returns a list with three components: modes, modes.dens, and size. 
    The elements in each component are ordered according to the decreasing density of the modes. The modes component is a vector of the values of x associated with the modes. 
    The modes.dens component is a vector of the kernel density estimates at the modes. 
    The size component is a vector of the proportion of area underneath each mode. 
    Function adapted from the Modes() function of the "LaplacesDemon" R package
    '''
    if x is None:
        raise ValueError("The x argument is required.")
    # Convert input to a numpy array of floats and filter out non-finite values
    x = np.asarray(x).astype(float)
    x = x[np.isfinite(x)]
    # Check if the input array is constant
    if is_constant(x):
        return {'modes': [np.nan], 'mode_dens': [np.nan], 'size': [1]}
    # Compute the kernel density estimate
    density = gaussian_kde(x)
    # Create a grid of 1000 points over the range of the data
    x_grid = np.linspace(np.min(x), np.max(x), 1000)
    # Evaluate the density over the grid
    dens_y = density(x_grid)
    # Compute the differences between consecutive density values
    dens_y_diff = np.diff(dens_y)
    # Identify where the density is increasing (1) or not increasing (0)
    incr = np.where(dens_y_diff > 0, 1, 0)
    # Initialize the list of segment boundaries
    begin = [0]
    for i in range(1, len(incr)):
        if incr[i] != incr[i - 1]:
            begin.append(i)
    # Add the end of the array as the final boundary
    begin.append(len(incr))
    #
    size = []
    modes = []
    mode_dens = []
    # Sum of all density values for normalization
    dens_y_sum = np.sum(dens_y)
    #
    j = 0
    while j < len(begin) - 1:
        start_idx = begin[j] # Define the start and end of the current segment
        end_idx = begin[j + 2] if j + 2 < len(begin) else len(dens_y) # Extract the segment of the density
        segment = dens_y[start_idx:end_idx]  # Calculate the size of the segment
        segment_size = np.sum(segment) / dens_y_sum
        if segment_size >= min_size: # Define the grid points and density values for the segment
            kde_x = x_grid[start_idx:end_idx]
            kde_y = segment
            mode_idx = np.argmax(kde_y) # Find the index of the mode (maximum density) within the segment
            modes.append(kde_x[mode_idx]) # Store the mode and its density
            mode_dens.append(kde_y[mode_idx])
            size.append(segment_size) # Store the size of the segment
        j += 2
    # Order the results by density in descending order
    order = np.argsort(mode_dens)[::-1]
    size = np.array(size)[order]
    modes = np.array(modes)[order]
    mode_dens = np.array(mode_dens)[order]
    # Normalize the sizes to sum to 1 if their sum exceeds 1
    if np.sum(size) > 1:
        size = size / np.sum(size)
    # Return dictionary
    return {'modes': modes, 'mode_dens': mode_dens, 'size': size}


def is_unimodal(x, min_size=0.1):
    '''
    Function to test if data is unimodal
    '''
    modes_result = Modes(x, min_size)
    if len(modes_result['modes']) == 1 and not np.any(np.isnan(modes_result['modes'])):
        return True
    else:
        return False


def skew(x, na_rm=True, type=3):
    '''
    Calculate the skewness of the input array x based on the specified type.
    '''
    if na_rm:
        x = x[~np.isnan(x)]
    n = len(x)
    if n == 0:
        return np.nan
    #
    if type == 1:
        skewer = np.sqrt(n) * (np.sum((x - np.mean(x))**3) / (np.sum((x - np.mean(x))**2)**(3/2)))
    elif type == 2:
        skewer = n * np.sqrt(n - 1) * (np.sum((x - np.mean(x))**3) / ((n - 2) * np.sum((x - np.mean(x))**2)**(3/2)))
    elif type == 3:
        skewer = np.sum((x - np.mean(x))**3) / (n * np.std(x)**3)
    else:
        raise ValueError("Type must be 1, 2, or 3")
    #
    return skewer


def kurtosi(x, na_rm=True, type=3):
    '''
    Calculate the kurtosis of the input array x based on the specified type.
    '''
    if na_rm:
        x = x[~np.isnan(x)]
    n = len(x)
    if n == 0:
        return np.nan
    #
    if type == 1:
        kurt = (np.sum((x - np.mean(x))**4) * n) / (np.sum((x - np.mean(x))**2)**2) - 3
    elif type == 2:
        kurt = (n * (n + 1) * np.sum((x - np.mean(x))**4)) / ((n - 1) * (n - 2) * (n - 3) * (np.sum((x - np.mean(x))**2) / (n - 1))**2) - 3 * (n - 1)**2 / ((n - 2) * (n - 3))
    elif type == 3:
        kurt = (np.sum((x - np.mean(x))**4) / (n * np.std(x)**4)) - 3
    else:
        raise ValueError("Type must be 1, 2, or 3")
    #
    return kurt


def bimodality_coefficient(x, na_rm=False):
    '''
    Calculate the bimodality coefficient using the skew and kurtosi functions.
    '''
    if na_rm:
        x = x[~np.isnan(x)]
    n = len(x)
    if n == 0:
        return np.nan
    #
    m3 = skew(x, na_rm=na_rm, type=2)
    # m3 = scipy_skew(x)
    m4 = kurtosi(x, na_rm=na_rm, type=2)
    # m4 = scipy_kurtosis(x)
    bc = (m3**2 + 1) / (m4 + 3 * ((n - 1)**2 / ((n - 2) * (n - 3))))
    return bc

#----
# 3. Functions for ploidy estimation
#----
def relative_to_absolute_CN(relative_CN, purity, ploidy):
    '''
    Convert relative copy number values to absolute copy number values using purity and ploidy
    '''
    acn = ploidy + (relative_CN - 1)*(ploidy+(2/purity)-2) 
    return acn


def acn_distance(relative_CN, purity, ploidy, distance_function="RMSD", weights=None): # function needs list of relative seg CNs and list of weights
    '''
    Estimate the distance function (goodness of fit) using RMSD (default) or MAD. 
    '''
    acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
    differences = [abs(x - round(x)) for x in acn]
    if weights == None:
        weights = [1]*len(relative_CN)
    if distance_function == "MAD":
        # estimate weighted distances using mean absolute deviation (MAD)
        distance = sum(differences[i] * weights[i] for i in range(len(differences))) / sum(weights)
    if distance_function == "RMSD":
        distance = math.sqrt(sum(differences[i]**2 * weights[i] for i in range(len(differences))) / sum(weights))
    return distance


def define_search_space(min, max, by=0.01):
    '''
    Define search space within which to look for potential purity-ploidy solutions.
    '''
    digs = len(str(by))-2 if isinstance(by,int) != True else 1
    # print(by,digs)
    if by <= 0:
        raise ValueError("by must be >= zero")
    result = []
    current = min
    if (min < max and by > 0):
        while (current <= max):
            result.append(current)
            current += by
            current = round(current,digs)
    return result


def build_search_grid(purity_seq, ploidy_seq):
    '''
    Define grid using purity and ploidy search spaces across which distance functions for purity and ploidy solutions will be estimated.
    '''
    unique_combinations = []
    for i in range(len(purity_seq)):
        for j in range(len(ploidy_seq)):
            unique_combinations.append([purity_seq[i], ploidy_seq[j], i+1, j+1])
    return unique_combinations


def estimate_grid_distances(min_cellularity, max_cellularity, cellularity_step, min_ploidy, max_ploidy, ploidy_step, relative_CN, weights=None, distance_function="RMSD"):
    '''
    Estimate goodness of fits (distance functions) across search grid.
    '''
    purs = define_search_space(min_cellularity,max_cellularity,by=cellularity_step)
    plois = define_search_space(min_ploidy,max_ploidy,by=ploidy_step)
    grid = build_search_grid(purs,plois)
    fits = []
    for pair in grid:
        cur_pur = pair[0]
        cur_ploi = pair[1]
        d = acn_distance(relative_CN, cur_pur, cur_ploi, weights=weights, distance_function=distance_function)
        pair.append(d)
        fits.append(pair)
    return fits


def reduce_grid(fits,distance_filter_scale_factor = 1.25):
    '''
    Reduce grid space by summarising fits (i.e. only keep peaks of hot areas of heatmap). In addition, only return fits that are within distance_filter_scale_factor * min(distance).
    '''
    reduced_grid = copy.deepcopy(fits)
    for xdelta in [-1, 0, 1]:
        for ydelta in [-1, 0, 1]:
            if xdelta != 0 or ydelta != 0:
                keep_idxs = []
                for idx,sol in enumerate(reduced_grid):
                    xc = sol[-2] + xdelta
                    yc = sol[-3] + ydelta
                    # Check if there's a corresponding (xc, yc) in fits
                    dc = next((d[-1] for d in fits if d[-2] == xc and d[-3] == yc),None)
                    if dc is None or sol[-1] <= dc:
                        keep_idxs.append(idx)
                reduced_grid = [reduced_grid[i] for i in keep_idxs]
    if isinstance(distance_filter_scale_factor, (int,float)):
        scaler = min([x[-1] for x in reduced_grid]) * distance_filter_scale_factor
        reduced_grid = [x for x in reduced_grid if x[-1] < scaler]
    return reduced_grid


def is_acceptable_fit(purity, ploidy, relative_CN, weights, max_proportion_zero = 0.1, min_proportion_close_to_whole_number = 0.5, max_distance_from_whole_number = 0.25, main_cn_step_change = 1):
    '''
    Test if given fit is acceptable using parameters looking at the proportion of zero or negative segments, proportion of genome fitted close to an integer, and distance between most common CN state peaks.
    '''
    acn = [relative_to_absolute_CN(x, purity, ploidy) for x in relative_CN]
    acn_int = [round(x) for x in acn]
    differences = [abs(x - round(x)) for x in acn]
    if weights == None:
        weights = [1]*len(relative_CN)
    # Filter fits based on fraction genome fitted to ACN of >= 0
    zeros_idxs =  [i for i in range(len(acn_int)) if acn_int[i] <= 0]
    prop_zero = sum([weights[i] for i in zeros_idxs]) / sum(weights)
    if prop_zero > max_proportion_zero:
        print(f"Fit of purity={purity} and ploidy={ploidy} is NOT an acceptable solution, as proportion of segments fitted to zero = {prop_zero}." )
        return False
    # Filter fits based on fraction genome that is fitted to < min_proportion_close_to_whole_number
    state_idxs = [i for i in range(len(differences)) if differences[i] < max_distance_from_whole_number] 
    prop_state = sum([weights[i] for i in state_idxs]) / sum(weights)
    if prop_state < min_proportion_close_to_whole_number:
        print(f"Fit of purity={purity} and ploidy={ploidy} is NOT an acceptable solution, as proportion of segments close to whole number  = {prop_state}." )
        return False
    # Remove fits which result in overstretchin/overfitting of copy number states (i.e. copy number state skipping) normally caused by oversegmentation
    most_common = statistics.mode(acn_int)
    second_most_common = statistics.mode([x for x in acn_int if x != most_common])
    if abs(most_common - second_most_common) > main_cn_step_change:
        print(f"Fit of purity={purity} and ploidy={ploidy} is NOT an acceptable solution, as main CN change step = {abs(most_common - second_most_common)}." )
        return False
    else:
        return True


def viable_solutions(fits_r, relative_CN, weights, max_proportion_zero = 0.1, min_proportion_close_to_whole_number = 0.5, max_distance_from_whole_number = 0.25, main_cn_step_change = 1):
    '''
    Return acceptable solution candidates.
    '''
    solutions = []
    for sol in fits_r:
        purity,ploidy,d = sol[0],sol[1],sol[-1]
        if is_acceptable_fit(purity, ploidy, relative_CN, weights,
                        max_proportion_zero = max_proportion_zero,
                        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
                        max_distance_from_whole_number = max_distance_from_whole_number, main_cn_step_change = main_cn_step_change) == True:
            solutions.append([sol[0],sol[1],sol[-1]])
    # sort solutions by distance function
    return solutions


def rank_solutions(solutions,distance_precision=3):
    '''
    Rank acceptable solution candidates.
    '''
    ranked = copy.deepcopy(solutions)
    # round distance function
    for x in ranked:
        x[-1] = round(x[-1],distance_precision)
    ranked = sorted(ranked,key=lambda x: (x[-1],x[1])) # sort by rounded distance, and second by ploidy
    # add rank position to output
    for idx,x in enumerate(ranked):
        x.append(idx+1)
    return ranked


def relative_to_absolute_minor_total_CN(chrom, rel_copy_number_segments, allele_counts, fitted_purity, fitted_ploidy):
    print(f"    ... estimating absolute minor and major allele copy number for {chrom} ...")
    acn_minor_major = []
    for x in rel_copy_number_segments:
        s_start, s_end, rcn, sid = x[1],x[2],x[-1],x[3]
        afs = [abs(0.5-float(x[10])) for x in allele_counts if int(x[1]) >= s_start and int(x[2]) <= s_end]
        if len(afs) < 1:
            print(f'        BAF and minor allele copy number cannot be estimated for segment {sid}... Total absolute copy number estimated only.')
            minorCN = ''
            baf_mean = ''
            totalCN = round(relative_to_absolute_CN(rcn, fitted_purity, fitted_ploidy),4)
            x[-1] = totalCN if totalCN > 0 else 0
            # x.append(minorCN)
        elif len(afs) >= 1:
            baf_mean = 0.5 + statistics.mean(afs)
            CN_a = (fitted_purity - 1 + (rcn*(1-baf_mean)*(2*(1-fitted_purity)+fitted_purity*fitted_ploidy))) / fitted_purity
            CN_b = (fitted_purity - 1 + (rcn*(baf_mean)*(2*(1-fitted_purity)+fitted_purity*fitted_ploidy))) / fitted_purity
            minorCN = round(min(CN_a,CN_b),4) if round(min(CN_a,CN_b),4) > 0 else 0
            totalCN = round((CN_a + CN_b),4) if round((CN_a + CN_b),4) > 0 else 0
            x[-1] = totalCN 
        x.extend((minorCN,baf_mean,len(afs)))
        acn_minor_major.append(x)
    return acn_minor_major


if __name__ == "__main__":
	print("Helper functions for Copy Number fitting")