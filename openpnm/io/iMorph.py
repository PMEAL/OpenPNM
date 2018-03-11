import os as os
import scipy as sp
from openpnm.core import logging
from openpnm.io import GenericIO
from openpnm.network import GenericNetwork
from openpnm.topotools import extend, trim
from openpnm.network import GenericNetwork
logger = logging.getLogger(__name__)


class iMorph(GenericIO):
    r"""
    Combines two output files from the iMorph program to build a pore network.
    throats_cellsThroatsGraph_Nodes.txt - stores node shape and type
    information throats_cellsThroatsGraph.txt - stores node connectivity
    """

    @classmethod
    def load(cls, path,
             node_file="throats_cellsThroatsGraph_Nodes.txt",
             graph_file="throats_cellsThroatsGraph.txt",
             network=None, voxel_size=None, return_geometry=False):
        r"""
        Loads network data from an iMorph processed image stack

        Parameters
        ----------
        path : string
            The path of the folder where the subfiles are held

        node_file : string
            The file that describes the pores and throats, the
            default iMorph name is: throats_cellsThroatsGraph_Nodes.txt

        graph_file : string
            The file that describes the connectivity of the network, the
            default iMorph name is: throats_cellsThroatsGraph.txt

        network : OpenPNM Network Object
            The OpenPNM Network onto which the data should be loaded.  If no
            network is supplied then an empty import network is created and
            returned.

        voxel_size : float
            Allows the user to define a voxel size different than what is
            contained in the node_file. The value must be in meters.

        return_geometry : Boolean
            If True, then all geometrical related properties are removed from
            the Network object and added to a GenericGeometry object.  In this
            case the method returns a tuple containing (network, geometry). If
            False (default) then the returned Network will contain all
            properties that were in the original file.  In this case, the user
            can call the ```split_geometry``` method explicitly to perform the
            separation.

        Returns
        -------
        If no Network object is supplied then one will be created and returned.

        If return_geometry is True, then a tuple is returned containing both
        the network and a geometry object.
        """
        #
        node_file = os.path.join(path.resolve(), node_file)
        graph_file = os.path.join(path.resolve(), graph_file)
        # parsing the nodes file
        with open(node_file, 'r') as file:
            Np = sp.fromstring(file.readline().rsplit('=')[1], sep='\t',
                               dtype=int)[0]
            vox_size = sp.fromstring(file.readline().rsplit(')')[1], sep='\t',)[0]

            # network always recreated to prevent errors
            network = GenericNetwork(Np=Np, Nt=0)

            # Define expected properies
            network['pore.volume'] = sp.nan
            scrap_lines = [file.readline() for line in range(4)]
            while True:
                vals = file.readline().split('\t')
                if len(vals) == 1:
                    break
                network['pore.volume'][int(vals[0])] = float(vals[3])
                if 'pore.'+vals[2] not in network.labels():
                    network['pore.'+vals[2]] = False
                network['pore.'+vals[2]][int(vals[0])] = True

        if voxel_size is None:
            voxel_size = vox_size * 1.0E-6  # file stores value in microns

        if voxel_size < 0:
            raise(Exception('Error - Voxel size must be specfied in ' +
                            'the Nodes file or as a keyword argument.'))

        # parsing the graph file
        with open(graph_file, 'r') as file:
            # Define expected properties
            network['pore.coords'] = sp.zeros((Np, 3))*sp.nan
            network['pore.types'] = sp.nan
            network['pore.color'] = sp.nan
            network['pore.radius'] = sp.nan
            network['pore.dmax'] = sp.nan
            network['pore.node_number'] = sp.nan
            # Scan file to get pore coordinate data
            scrap_lines = [file.readline() for line in range(3)]
            line = file.readline()
            xmax = 0.0
            ymax = 0.0
            zmax = 0.0
            node_num = 0
            while line != 'connectivity table\n':
                vals = sp.fromstring(line, sep='\t')
                xmax = vals[1] if vals[1] > xmax else xmax
                ymax = vals[2] if vals[2] > ymax else ymax
                zmax = vals[3] if vals[3] > zmax else zmax
                network['pore.coords'][int(vals[0]), :] = vals[1:4]
                network['pore.types'][int(vals[0])] = vals[4]
                network['pore.color'][int(vals[0])] = vals[5]
                network['pore.radius'][int(vals[0])] = vals[6]
                network['pore.dmax'][int(vals[0])] = vals[7]
                network['pore.node_number'][int(vals[0])] = node_num
                node_num += 1
                line = file.readline()
            # Scan file to get to connectivity data
            scrap_lines.append(file.readline())  # Skip line
            # Create sparse lil array incrementally build adjacency matrix
            lil = sp.sparse.lil_matrix((Np, Np), dtype=int)
            while True:
                vals = sp.fromstring(file.readline(), sep='\t', dtype=int)
                if len(vals) <= 1:
                    break
                lil.rows[vals[0]] = vals[2:]
                lil.data[vals[0]] = sp.ones(vals[1])

        # fixing any negative volumes or distances so they are 1 voxel/micron
        network['pore.volume'][sp.where(network['pore.volume'] < 0)[0]] = 1.0
        network['pore.radius'][sp.where(network['pore.radius'] < 0)[0]] = 1.0
        network['pore.dmax'][sp.where(network['pore.dmax'] < 0)[0]] = 1.0

        # Add adjacency matrix to OpenPNM network
        conns = sp.sparse.triu(lil, k=1, format='coo')
        network.update({'throat.all': sp.ones(len(conns.col), dtype=bool)})
        network['throat.conns'] = sp.vstack([conns.row, conns.col]).T

        network['pore.to_trim'] = False
        network['pore.to_trim'][network.pores('*throat')] = True
        Ts = network.pores('to_trim')
        new_conns = network.find_neighbor_pores(pores=Ts, flatten=False)
        extend(network=network, throat_conns=new_conns, labels='new_conns')
        for item in network.props('pore'):
            item = item.split('.')[1]
            arr = sp.ones_like(network['pore.'+item])[0]
            arr = sp.tile(A=arr, reps=[network.Nt, 1])*sp.nan
            network['throat.'+item] = sp.squeeze(arr)
            network['throat.'+item][network.throats('new_conns')] = \
                network['pore.'+item][Ts]
        trim(network=network, pores=Ts)

        # setting up boundary pores
        x_coord, y_coord, z_coord = sp.hsplit(network['pore.coords'], 3)
        network['pore.front_boundary'] = sp.ravel(x_coord == 0)
        network['pore.back_boundary'] = sp.ravel(x_coord == xmax)
        network['pore.left_boundary'] = sp.ravel(y_coord == 0)
        network['pore.right_boundary'] = sp.ravel(y_coord == ymax)
        network['pore.bottom_boundary'] = sp.ravel(z_coord == 0)
        network['pore.top_boundary'] = sp.ravel(z_coord == zmax)

        # removing any pores that got classified as a boundary pore but
        # weren't labled a border_cell_face
        ps = sp.where(~sp.in1d(network.pores('*_boundary'),
                               network.pores('border_cell_face')))[0]
        ps = network.pores('*_boundary')[ps]
        for side in ['front', 'back', 'left', 'right', 'top', 'bottom']:
            network['pore.'+side+'_boundary'][ps] = False
        # setting internal label
        network['pore.internal'] = False
        network['pore.internal'][network.pores('*_boundary', mode='not')] = True

        # adding props to border cell face throats and from pores
        Ts = sp.where(network['throat.conns'][:, 1] >
                       network.pores('border_cell_face')[0] - 1)[0]
        faces = network['throat.conns'][Ts, 1]
        for item in network.props('pore'):
            item = item.split('.')[1]
            network['throat.'+item][Ts] = network['pore.'+item][faces]
        network['pore.volume'][faces] = 0.0

        # applying unit conversions
        # TODO: Determine if radius and dmax are indeed microns and not voxels
        network['pore.coords'] = network['pore.coords'] * 1e-6
        network['pore.radius'] = network['pore.radius'] * 1e-6
        network['pore.dmax'] = network['pore.dmax'] * 1e-6
        network['pore.volume'] = network['pore.volume'] * voxel_size**3
        network['throat.coords'] = network['throat.coords'] * 1e-6
        network['throat.radius'] = network['throat.radius'] * 1e-6
        network['throat.dmax'] = network['throat.dmax'] * 1e-6
        network['throat.volume'] = network['throat.volume'] * voxel_size**3

        return network.project
