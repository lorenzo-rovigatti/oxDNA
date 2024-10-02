import numpy as np
import argparse
from os import environ, remove, path
from sys import stderr
from typing import List
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.cluster import DBSCAN
from matplotlib import animation
from json import dump, load
from oxDNA_analysis_tools.config import check
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo, TopInfo
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, linear_read, conf_to_str, describe, write_conf

def split_trajectory(traj_info, top_info, labs):
    """
    Splits the trajectory into the clustered trajectories

    Parameters:
        traj_info (TrajInfo): Metadata on the trajectory file
        top_info (TopInfo): Metadata on the topology file
        labs (numpy.array): The cluster each configuration belongs to.
    """

    #How many are in each cluster?
    print ("cluster\tmembers")

    slabs = set(labs)

    for cluster in slabs:
        in_cluster = list(labs).count(cluster)
        print ("{}\t{}".format(cluster, in_cluster))

        #Clear old trajectory files
        try:
            remove("cluster_"+str(cluster)+".dat")
        except: pass

    print ("INFO: splitting trajectory...", file=stderr)
    print ("INFO: Trajectories for each cluster will be written to cluster_<cluster number>.dat", file=stderr)

    fnames = ["cluster_"+str(cluster)+".dat" for cluster in slabs]
    files = [open(f, 'w+') for f in fnames]
    i = 0

    for chunk in linear_read(traj_info, top_info):
        for conf in chunk:
            files[labs[i]].write(conf_to_str(conf, include_vel=traj_info.incl_v))
            i += 1
    [f.close() for f in files]

    print(f"INFO: Wrote trajectory files: {fnames}", file=stderr)

    return

def find_element(n, x, array):
    """
    Finds the id of the nth time element x appears in an array.
    """
    c = 0
    for i, j in enumerate(array):
        if (j == x):
            if (c == n):
                return i
            c += 1
    return -1

def get_centroid(points, metric_name, labs, traj_info, top_info) -> List[int]:
    """
    Takes the output from DBSCAN and produces the trajectory and centroid from each cluster.

    Parameters:
        points (numpy.array): The points fed to the clstering algorithm.
        metric_name (str): The type of data the points represent.
        labs (numpy.array): The cluster each point belongs to.
        traj_info (TrajInfo): Trajectory metadata.
        tpo_file (TopInfo): Topology metadata.
    """

    if metric_name == 'euclidean':
        points = points[np.newaxis,:,:] - points[:,np.newaxis,:]
        points = np.sqrt(np.sum(points**2, axis=2))    
    
    print("INFO: Finding cluster centroid...", file=stderr)
    cids = []
    for cluster in (set(labs)):
        masked = points[labs == cluster]
        in_cluster_id = np.sum(masked, axis = 1).argmin()

        centroid_id = find_element(in_cluster_id, cluster, labs)
        cids.append(centroid_id)

        centroid = get_confs(top_info, traj_info, centroid_id, 1)[0]
        fname = "centroid_"+str(cluster)+".dat"
        write_conf(fname, centroid, include_vel=traj_info.incl_v)
        print(f"INFO: Wrote centroid file {fname}", file=stderr)

    return cids

def make_plot(op, labels, centroid_ids):
    # Prepping a plot of the first 3 dimensions of the provided op
    dimensions = []
    x = []
    y = []
    dimensions.append(x)
    dimensions.append(y)

    # if the op is 1-dimensional add a time dimension
    add_time = False
    if op.shape[1] == 1:
        add_time = True
        op = np.hstack((op, np.arange(op.shape[0]).reshape(op.shape[0], 1)))

    if op.shape[1] > 2:
        z = []
        dimensions.append(z)
    
    for i in op:
        for j, dim in enumerate(dimensions):
            dim.append(i[j])

    dimensions = np.array(dimensions)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    print("INFO: Making cluster plot...")
    if len(dimensions) == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel("OP0")
    ax.set_ylabel("OP1")

    if len(dimensions) == 3:
        ax.set_zlabel("OP2")
        #to show the plot immediatley and interactivley
        '''a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', 7))
        b = fig.colorbar(a, ax=ax)
        plt.show()'''
        
        #to make a video showing a rotating plot
        plot_file = "animated.mp4"
        def init():
            nonlocal labels, dimensions, n_clusters, centroid_ids
            a = ax.scatter(x, y, z, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters+1))
            cen = ax.scatter(dimensions[0][centroid_ids], dimensions[1][centroid_ids], dimensions[2][centroid_ids], s=1.5, c=[0 for _ in centroid_ids], cmap=ListedColormap(['black']))
            fig.colorbar(a, ax=ax)
            handles, _ = cen.legend_elements(prop="colors", num = 1)
            l = ax.legend(handles, ['Centroids'])
            return [fig]

        def animate(i):
            ax.view_init(elev=10., azim=i)
            return [fig]

        try:
            anim = animation.FuncAnimation(fig, animate, init_func=init, frames=range(360), interval=20, blit=True)
            anim.save(plot_file, fps=30, extra_args=['-vcodec', 'libx264'])
        except:
            print("WARNING: ffmpeg not found, cannot make animated plot, opening interactivley instead", file=stderr)
            f = init()
            plt.show()

    else:
        plot_file = "cluster_plot.png"
        if add_time:
            a = ax.scatter(dimensions[1], dimensions[0], s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters+1), vmin=min(labels)-0.5, vmax=max(labels)+0.5)
            cen = ax.scatter(dimensions[1][centroid_ids], dimensions[0][centroid_ids], s=1.5, c=[0 for _ in centroid_ids], cmap=ListedColormap(['black']))
            ax.set_xlabel("conf id")
            ax.set_ylabel("OP0")
        else:
            a = ax.scatter(x, y, s=2, alpha=0.4, c=labels, cmap=plt.get_cmap('tab10', n_clusters+1), vmin=min(labels)-0.5, vmax=max(labels)+0.5)
            cen = ax.scatter(dimensions[0][centroid_ids], dimensions[1][centroid_ids], s=1.5, c=[0 for _ in centroid_ids], cmap=ListedColormap(['black']))

        b = fig.colorbar(a, ax=ax, ticks=list(set(labels)))
        handles, _ = cen.legend_elements(prop="colors", num = 1)
        l = ax.legend(handles, ['Centroids'])
        ax.add_artist(l)
        plt.tight_layout()
        plt.savefig(plot_file)
    print("INFO: Saved cluster plot to {}".format(plot_file), file=stderr)


def perform_DBSCAN(traj_info:TrajInfo, top_info:TopInfo, op:np.ndarray, metric:str, eps:float, min_samples:int):
    """
    Use the DBSCAN algorithm to identify clusters of configurations based on a given order parameter.

    Parameters:
        traj_info (TrajInfo): Information about the trajectory
        top_info (TopInfo): Information about the topology
        op (np.ndarray): The order parameter(s) to use
        metric (str): Either 'euclidean' or 'precomputed' for whether the distance needs to be calculated
        eps (float): The maximum distance between two points to be considered in the same neighborhood
        min_samples (int): The minimum number of points to be considered a neighborhood
    """
    
    check(["python", "sklearn", "matplotlib"])
    
    #dump the input as a json file so you can iterate on eps and min_samples
    dump_file = "cluster_data.json"
    print("INFO: Serializing input data to {}".format(dump_file), file=stderr)
    print("INFO: Run  `oat clustering {} -e<eps> -m<min_samples>`  to adjust clustering parameters".format(dump_file), file=stderr)
    out = {
        "data": op.tolist(), 
        "traj" : traj_info.path,
        "metric" : metric
    }
    dump(out, open(dump_file, 'w+'))


    
    print("INFO: Running DBSCAN...", file=stderr)

    #DBSCAN parameters:
    #eps: the pairwise distance that configurations below are considered neighbors
    #min_samples: The smallest number of neighboring configurations required to start a cluster
    #metric: If the matrix fed in are points in n-dimensional space, then the metric needs to be "euclidean".
    #        If the matrix is already a square distance matrix, the metrix needs to be "precomputed".
    #the eps and min_samples need to be determined for each input based on the values of the input data
    #If you're making your own multidimensional data, you probably want to normalize your data first.
    print("INFO: Current values: eps={}, min_samples={}".format(eps, min_samples))
    db = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(op) 
    labels = db.labels_
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print ("Number of clusters:", n_clusters_)

    split_trajectory(traj_info, top_info, labels)
    centroid_ids =  get_centroid(op, metric, labels, traj_info, top_info)
    make_plot(op, labels, centroid_ids)
    
    print("INFO: Run  `oat clustering {} -e<eps> -m<min_samples>`  to adjust clustering parameters".format(dump_file), file=stderr)

    return labels

def cli_parser(prog="clustering.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Calculates clusters based on provided order parameters.")
    parser.add_argument('serialized_data', type=str, nargs=1, help="The json-formatted input file")
    parser.add_argument('-e', '--eps', type=float, nargs=1, help="The epsilon parameter for DBSCAN (maximum distance to be considered a 'neighbor')")
    parser.add_argument('-m', '--min_samples', type=int, nargs=1, help="The min_samples parameter for DBSCAN (number of neighbors which define a point as a central point in a cluster)")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()
    data_file = args.serialized_data[0]
    if args.eps:
        eps = args.eps[0]
    else:
        eps = 12
    if args.min_samples:
        min_samples = args.min_samples[0]
    else:
        min_samples = 8

    #load a previously serialized dataset
    with open(data_file, 'r') as f:
        data = load(f)
    points = np.array(data["data"])
    top_info, traj_info = describe(None, data["traj"])
    metric = data["metric"]
    labels = perform_DBSCAN(traj_info, top_info, points, metric, eps, min_samples)

if __name__ == '__main__':
    main()