# -*- coding: utf-8 -*-
"""
===============================================================================
MixedInvasionPercolationCoop: IP allowing pores and throats to invade separately
With added cooperative filling algorithms
===============================================================================

"""
import time
import logging
import heapq as hq
import scipy as sp
import numpy as np
from scipy.sparse import coo_matrix, dok_matrix
from openpnm.algorithms import MixedInvasionPercolation
from transforms3d._gohlketransforms import angle_between_vectors

logger = logging.getLogger(__name__)


class MixedInvasionPercolationCoop(MixedInvasionPercolation):
    r"""
    An implemetation of invasion percolation which can invade bonds, sites or a
    mixture of both. Inlets can be treated as individual injection points that
    share a common pressure or have their own and progess independently.
    Inlets can also be single pores or clusters.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which the invasion should occur.

    Notes
    ----
    n/a

    """

    def __init__(self, settings={}, **kwargs):
        def_set = {
            "pore_entry_pressure": "pore.entry_pressure",
            "throat_entry_pressure": "throat.entry_pressure",
            "snap_off": "",
            "invade_isolated_Ts": False,
            "late_pore_filling": "",
            "late_throat_filling": "",
            "gui": {
                "setup": {
                    "pore_entry_pressure": "",
                    "throat_entry_pressure": "",
                    "snap_off": "",
                    "invade_isolated_Ts": "",
                },
                "set_inlets": {"pores": None, "clusters": None},
                "set_outlets": {"pores": None, "overwrite": False},
                "apply_flow": {"flowrate": None},
                "apply_trapping": {"partial": False},
                "set_residual": {"pores": None, "overwrite": False},
            },
        }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)

    def setup(
        self,
        phase=None,
        pore_entry_pressure="pore.entry_pressure",
        throat_entry_pressure="throat.entry_pressure",
        snap_off="",
        invade_isolated_Ts=False,
        late_pore_filling="",
        late_throat_filling="",
        cooperative_pore_filling="",
    ):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        pore_entry_pressure : string
            The dictionary key on the Phase object where the pore entry
            pressure values are stored.  The default is
            'pore.entry_pressure'.

        throat_entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.entry_pressure'.

        snap_off : string
            The dictionary key on the Phase object where the throat snap-off
            pressure values are stored.

        invade_isolated_Ts : boolean
            If True, isolated throats are invaded at the higher invasion
            pressure of their connected pores.

        late_pore_filling : string
            The name of the model used to determine late pore filling as
            a function of applied pressure.

        late_throat_filling : string
            The name of the model used to determine late throat filling as
            a function of applied pressure.

        cooperative_pore_filling : string
            The name of the model used to determine the meniscus properties
            required for assessing cooperative pore filling.
        """
        if phase:
            self.settings["phase"] = phase.name
        if throat_entry_pressure:
            self.settings["throat_entry_pressure"] = throat_entry_pressure
            phase = self.project.find_phase(self)
        self["throat.entry_pressure"] = phase[self.settings["throat_entry_pressure"]]
        if len(np.shape(self["throat.entry_pressure"])) > 1:
            self._bidirectional = True
        else:
            self._bidirectional = False
        if pore_entry_pressure:
            self.settings["pore_entry_pressure"] = pore_entry_pressure
            phase = self.project.find_phase(self)
        self["pore.entry_pressure"] = phase[self.settings["pore_entry_pressure"]]
        if snap_off:
            self.settings["snap_off"] = snap_off
        if invade_isolated_Ts:
            self.settings["invade_isolated_Ts"] = invade_isolated_Ts
        if late_pore_filling:
            self.settings["late_pore_filling"] = late_pore_filling
        if late_throat_filling:
            self.settings["late_throat_filling"] = late_throat_filling
        if cooperative_pore_filling:
            self.settings["cooperative_pore_filling"] = cooperative_pore_filling
        self.reset()

    def _max_pressure(self):
        phase = self.project.find_phase(self)
        if self.settings["throat_entry_pressure"]:
            max_tPc = np.max(phase[self.settings["throat_entry_pressure"]])
        else:
            max_tPc = 0.0
        if self.settings["pore_entry_pressure"]:
            max_pPc = np.max(phase[self.settings["pore_entry_pressure"]])
        else:
            max_pPc = 0.0
        return np.max([max_tPc, max_pPc])

    def _my_dot(self, a, b):
        return a[:, 0] * b[:, 0] + a[:, 1] * b[:, 1] + a[:, 2] * b[:, 2]

    def trilaterate_v(self, P1, P2, P3, r1, r2, r3):
        r"""
        Find whether 3 spheres intersect
        """
        temp1 = P2 - P1
        e_x = temp1 / np.linalg.norm(temp1, axis=1)[:, np.newaxis]
        temp2 = P3 - P1
        i = self._my_dot(e_x, temp2)[:, np.newaxis]
        temp3 = temp2 - i * e_x
        e_y = temp3 / np.linalg.norm(temp3, axis=1)[:, np.newaxis]
        d = np.linalg.norm(P2 - P1, axis=1)[:, np.newaxis]
        j = self._my_dot(e_y, temp2)[:, np.newaxis]
        x = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
        y = (r1 * r1 - r3 * r3 - 2 * i * x + (i * i) + (j * j)) / (2 * j)
        temp4 = r1 * r1 - x * x - y * y
        return temp4 >= 0

    def _get_throat_pairs(self):
        r"""
        Generate an array of pores with all connected throats and pairs of
        throats that connect to the same pore
        """
        network = self.project.network
        # Collect all throat pairs sharing a pore as list of lists
        neighbor_Ts = network.find_neighbor_throats(
            pores=network.pores(), flatten=False
        )
        # Pores associated to throat
        Ps = []
        # Throats connected to each pore
        Ts = []
        # Pairs of throats sharing a pore
        T1 = []
        T2 = []
        start = 0
        # Build lookup pair index arrays for each coordination number up to the
        # Maximum coordination max_c
        max_c = np.amax(network.num_neighbors(pores=network.Ps, flatten=False))
        pair_T1 = []
        pair_T2 = []
        logger.info("Building throat pair matrices")
        for num_t in range(max_c + 1):
            temp1 = []
            temp2 = []
            for t1 in range(num_t)[:-1]:
                for t2 in range(num_t)[t1 + 1 :]:
                    temp1.append(t1)
                    temp2.append(t2)
            pair_T1.append(np.asarray(temp1))
            pair_T2.append(np.asarray(temp2))
        for p, nTs in enumerate(neighbor_Ts):
            num_t = len(nTs)
            for i in range(num_t):
                Ps.append(p)
                Ts.append(nTs[i])
            # Pair indices into Ts
            tempt1 = pair_T1[num_t] + start
            tempt2 = pair_T2[num_t] + start
            for i in range(len(tempt1)):
                T1.append(tempt1[i])
                T2.append(tempt2[i])
            start += num_t
        # Nt * 2 long
        Ps = np.asarray(Ps)
        Ts = np.asarray(Ts)
        # indices into the above arrays based on throat pairs
        T1 = np.asarray(T1)
        T2 = np.asarray(T2)

        return Ps, Ts, T1, T2

    def _apply_cen_to_throats(self, p_cen, t_cen, t_norm, men_cen):
        r"""
        Take the pore center and throat center and work out which way
        the throat normal is pointing relative to the vector between centers.
        Offset the meniscus center along the throat vector in the correct
        direction
        """
        v = p_cen - t_cen
        sign = np.sign(np.sum(v * t_norm, axis=1))
        c3 = np.vstack((men_cen * sign, men_cen * sign, men_cen * sign)).T
        coords = t_cen + c3 * t_norm
        return coords

    def _transform_point_normal(self, point, normal):
        r"""
        Transforms point normal plane definition to parametric form
        Ax + By +Cz + D = 0
        """
        return [
            a
            for a in zip(
                normal[:, 0], normal[:, 1], normal[:, 2], -self._my_dot(point, normal)
            )
        ]

    def _plane_intersect(self, a, b):
        """
        a, b   4-tuples/lists
               Ax + By +Cz + D = 0
               A, B, C, D in order
        output: 2 points on line of intersection, np.arrays, shape (3,)
        https://bit.ly/2LkBEyc
        """
        a_vec, b_vec = np.array(a[:3]), np.array(b[:3])
        aXb_vec = np.cross(a_vec, b_vec)
        A = np.array([a_vec, b_vec, aXb_vec])
        d = np.array([-a[3], -b[3], 0.0]).reshape(3, 1)
        if np.linalg.det(A) == 0:
            return None, None
        else:
            p_inter = np.linalg.solve(A, d).T
            return p_inter[0], (p_inter + aXb_vec)[0]

    def _t(self, p, q, r):
        r"""
        Equation of line passing through points p and q parameterized by t
        https://bit.ly/2EpQ6DD
        """
        x = p - q
        return np.dot(r - q, x) / np.dot(x, x)

    def _distance(self, p, q, r):
        r"""
        Shortest distance between line passing through p and q and point r
        https://bit.ly/2EpQ6DD
        """
        return np.linalg.norm(self._t(p, q, r) * (p - q) + q - r)

    def _pair_permutations(self, n):
        perms = []
        for i in range(n):
            for j in range(n - (i + 1)):
                perms.append([i, j + i + 1])
        return perms

    def _perpendicular_vector(self, v, v_ref=None):
        if v_ref is None:
            np.array([1.0, 0.0, 0.0])
        return np.cross(v, v_ref)

    def _throat_pair_angle(self, t1, t2, pore, network):
        conn1 = network["throat.conns"][t1]
        conn2 = network["throat.conns"][t2]
        lookup = np.vstack((pore, pore)).T
        p1 = conn1[conn1 != lookup]
        p2 = conn2[conn2 != lookup]
        v1 = network["pore.coords"][pore] - network["pore.coords"][p1]
        v2 = network["pore.coords"][pore] - network["pore.coords"][p2]
        v_ref = self._perpendicular_vector(v1, v2)
        v1_perp = self._perpendicular_vector(v1, v_ref)
        v2_perp = self._perpendicular_vector(v2, v_ref)
        return angle_between_vectors(v1_perp, v2_perp, axis=1)

    def setup_coop_filling(self, inv_points=None):
        r"""
        Populate the coop filling throat-throat pair matrix
        """
        self._setup_coop_filling_creep(inv_points)
        self._setup_coop_filling_bulge(inv_points)

    def _setup_coop_filling_creep(self, inv_points=None):
        r"""
        This coop filling model iterates through the invasion pressures set by
        the inv_points variable and calls the meniscus model whose dictionary
        key must be given in the algorithm's setup.
        The meniscus model supplies the position of the meniscus inside each
        throat as if there were a meniscus in every throat for a given pressure
        The contact line of the meniscus traces a circle around the inner
        surface of the throat which is assumed to be toroidal.
        The contact circle lies on a plane that is defined by the throat's
        normal vector which runs along the axis of symmetry of the torus.
        For every pore, every connecting throat is compared with each of it's
        neighboring throats connected to the same pore. If the planes intersect
        then the meniscus contact circles may eventually touch if they can
        advance enough. For highly wetting fluid the contact point may be
        advanced well into the throat whilst still being at negative capillary
        pressure.
        """
        start = time.time()
        net = self.project.network
        phase = self.project.find_phase(self)
        all_phys = self.project.find_physics(phase=phase)
        if inv_points is None:
            inv_points = np.arange(0, 1.01, 0.01) * self._max_pressure()
        # Throat centroids
        try:
            t_centroids = net["throat.centroid"]
        except KeyError:
            t_centroids = np.mean(net["pore.coords"][net["throat.conns"]], axis=1)
        # Throat normal vector
        try:
            t_norms = net["throat.normal"]
        except KeyError:
            t_norms = (
                net["pore.coords"][net["throat.conns"][:, 1]]
                - net["pore.coords"][net["throat.conns"][:, 0]]
            )
        cpf = self.settings["cooperative_pore_filling"]
        model = all_phys[0].models[cpf]
        # Throat Diameter and fiber radius
        try:
            t_rad = net[model["throat_diameter"]] / 2 + model["r_toroid"]
        except KeyError:
            t_rad = net["throat.diameter"] / 2
        # Equations of throat planes at the center of each throat
        planes = self._transform_point_normal(t_centroids, t_norms)
        Nt = net.Nt
        adj_mat = dok_matrix((Nt, Nt), dtype=int)
        # Run through and build the throat-throat pair adjacency matrix
        # If planes of throats intersect then meniscii in throats may also
        # Intersect at a given pressure.
        for pore in net.Ps:
            ts = net.find_neighbor_throats(pores=pore)
            perms = self._pair_permutations(len(ts))
            for perm in perms:
                ta, tb = ts[perm]
                p, q = self._plane_intersect(planes[ta], planes[tb])
                if p is not None:
                    d1 = self._distance(p, q, t_centroids[ta])
                    d2 = self._distance(p, q, t_centroids[tb])
                    if (t_rad[ta] >= d1) and (t_rad[tb] >= d2):
                        adj_mat[ta, tb] = pore + 1
        # Meniscus Filling Angle
        tfill_angle = cpf + ".alpha"
        # Capillary pressure adjacency maxtrix
        pairs = np.asarray([list(key) for key in adj_mat.keys()])
        pores = np.asarray([adj_mat[key] for key in adj_mat.keys()]) - 1
        self.tt_Pc = adj_mat.copy().astype(float).tocsr()
        # Initialize the pressure matrix with nans
        # This is used to check for the first intersection pressure and
        # Prevent overwriting
        self.tt_Pc[pairs[:, 0], pairs[:, 1]] = np.nan
        angles = self._throat_pair_angle(pairs[:, 0], pairs[:, 1], pores, net)
        T1 = pairs[:, 0]
        T2 = pairs[:, 1]
        hits = []
        for Pc in inv_points:
            # Don't use zero as can get strange numbers in menisci data
            if Pc == 0.0:
                Pc = 1e-6
            # regenerate model with new target Pc
            for phys in all_phys:
                phys.models[cpf]["target_Pc"] = Pc
                phys.regenerate_models(propnames=cpf)
            # check whether this throat pair already has a coop value
            check_nans = np.asarray(np.isnan(self.tt_Pc[T1, T2]).tolist()[0])
            fill_angle_sum = np.sum(phase[tfill_angle][pairs], axis=1)
            coalescence = fill_angle_sum >= angles
            mask = check_nans * coalescence
            if np.any(mask):
                self.tt_Pc[T1[mask], T2[mask]] = Pc
                hits.append(Pc)
        # Change to lil for single throat lookups
        #        self.tt_Pc = self.tt_Pc.tolil()
        logger.info(
            "Coop filling finished in " + str(np.around(time.time() - start, 2)) + " s"
        )
        logger.info("Coop Hits", np.unique(np.asarray(hits)))

    def _setup_coop_filling_bulge(self, inv_points=None):
        r"""
        Evaluate the cooperative pore filling condition that the combined
        filling angle in next neighbor throats cannot exceed the geometric
        angle between their throat planes.
        This is used when the invading fluid has access to multiple throats
        connected to a pore

        Parameters
        ----------
        inv_points : array_like
            The invasion pressures at which to assess coopertive pore filling.
        """
        net = self.project.network
        phase = self.project.find_phase(self)
        all_phys = self.project.find_physics(phase=phase)
        if inv_points is None:
            inv_points = np.arange(0, 1.01, 0.01) * self._max_pressure()

        start = time.time()
        cpf = self.settings["cooperative_pore_filling"]
        tfill_angle = cpf + ".alpha"
        tmen_rad = cpf + ".radius"
        tmen_cen = cpf + ".center"
        try:
            # The following properties will all be there for Voronoi
            p_centroids = net["pore.centroid"]
            t_centroids = net["throat.centroid"]
            p_rad = net["pore.indiameter"] / 2
            t_norms = net["throat.normal"]
        except KeyError:
            # Chances are this isn't Voronoi so calculate or replace all
            p_centroids = net["pore.coords"]
            temp = net["pore.coords"][net["throat.conns"]]
            t_centroids = np.mean(temp, axis=1)
            p_rad = net["pore.diameter"] / 2
            t_norms = net["throat.normal"]

        Ps, Ts, T1, T2 = self._get_throat_pairs()
        # Network indices for the pairs
        pps = Ps[T1]
        pt1 = Ts[T1]
        pt2 = Ts[T2]
        # Arrays to get upper and lower triangular for reference to later
        # When checking for phase presence one throat at a time - need to
        # Check all throat pairs
        T_all_t1 = Ts[np.concatenate((T1, T2))]
        # T_all_t2 = Ts[np.concatenate((T2, T1))]
        # Throat-Throat cooperative filling pressure
        data = np.ones(len(T_all_t1), dtype=float)
        data.fill(np.nan)
        # coo --> useful for building matrix
        # coo = coo_matrix(
        #     (data, (T_all_t1, T_all_t2)), shape=(len(T_all_t1), len(T_all_t2))
        # )
        # csr --> for indexing into entire matrix when populating with
        # actual coop filling values
        #        self.tt_Pc = coo.tocsr()
        # Pair pore center and radius
        pp_cen = p_centroids[pps]
        pp_rad = p_rad[pps]
        # Make sure throat normals are unit vector
        unit = np.linalg.norm(t_norms, axis=1)
        t_norms /= np.vstack((unit, unit, unit)).T

        for Pc in inv_points:
            # regenerate model with new target Pc
            for phys in all_phys:
                phys.models[cpf]["target_Pc"] = Pc
                phys.regenerate_models(propnames=cpf)
            men_cen_dist = phase[tmen_cen]
            # Work out meniscii coord for each direction along the throat
            men_cen_coord = self._apply_cen_to_throats(
                p_centroids[Ps], t_centroids[Ts], t_norms[Ts], men_cen_dist[Ts]
            )
            # Pair centers
            pc1 = men_cen_coord[T1]
            pc2 = men_cen_coord[T2]
            # Center to center vector between neighboring meniscii
            c2c = pc1 - pc2
            dist = np.linalg.norm(c2c, axis=1)
            # Pair meniscii radii
            pr1 = phase[tmen_rad][pt1]
            pr2 = phase[tmen_rad][pt2]
            # nans may exist if pressure is outside the range
            # set these to zero to be ignored by next step without
            # causing RuntimeWarning
            pr1[np.isnan(pr1)] = 0
            pr2[np.isnan(pr2)] = 0
            # Negative mensicii radii means positive pressure
            # Assume meniscii only interact when bulging into pore
            check_neg = np.logical_and(pr1 < 0, pr2 < 0)
            # simple initial distance check on sphere rads
            check_rads = (np.abs(pr1 + pr2)) >= dist
            # check whether the filling angle is ok at this Pc
            check_alpha_T1 = ~np.isnan(phase[tfill_angle][pt1])
            check_alpha_T2 = ~np.isnan(phase[tfill_angle][pt2])
            check_alpha = check_alpha_T1 * check_alpha_T2
            # check whether this throat pair already has a coop value
            check_nans = np.asarray(np.isnan(self.tt_Pc[pt1, pt2]).tolist()[0])
            mask = check_neg * check_alpha * check_nans * check_rads
            # if all checks pass
            if np.any(mask):
                # Check if intersecting circle lies within pore
                inter = self.trilaterate_v(
                    P1=pc1[mask],
                    P2=pc2[mask],
                    P3=pp_cen[mask],
                    r1=pr1[mask][:, np.newaxis],
                    r2=pr2[mask][:, np.newaxis],
                    r3=pp_rad[mask][:, np.newaxis],
                )
                inter = inter.flatten()
                if np.any(inter):
                    self.tt_Pc[pt1[mask][inter], pt2[mask][inter]] = Pc
                    self.tt_Pc[pt2[mask][inter], pt1[mask][inter]] = Pc
            # Save meniscii data for each throat into each connecting pore
            men_data = {}
            men_data["Ps"] = Ps
            men_data["Ts"] = Ts
            men_data["cen"] = men_cen_coord
            men_data["rad"] = phase[tmen_rad][Ts]

        # Change to lil for single throat lookups
        self.tt_Pc = self.tt_Pc.tolil()
        logger.info(
            "Coop filling finished in " + str(np.around(time.time() - start, 2)) + " s"
        )

    def _check_coop(self, pore, queue):
        r"""
        Method run in loop after every pore invasion. All connecting throats
        are now given access to the invading phase. Two throats with access to
        the invading phase can cooperatively fill any pores that they are both
        connected to, common pores.
        The invasion of theses throats connected to the common pore is handled
        elsewhere.
        """
        net = self.project.network
        t_inv = "throat.invasion_sequence"
        p_inv = "pore.invasion_sequence"
        for throat in net.find_neighbor_throats(pores=pore):
            # A pore has just been invaded, all it's throats now have
            # An interface residing inside them
            if self[t_inv][throat] == -1:
                # If the throat is not the invading throat that gave access
                # to this pore, get the pores that this throat connects with
                a = set(net["throat.conns"][throat])
                # Get a list of pre-calculated coop filling pressures for all
                # Throats this throat can coop fill with
                ts_Pc = self.tt_Pc.data[throat]
                # Network indices of throats that can act as filling pairs
                ts = self.tt_Pc.rows[throat]
                # If there are any potential coop filling throats
                if np.any(~np.isnan(ts_Pc)):
                    ts_Pc = np.asarray(ts_Pc)
                    ts = np.asarray(ts)
                    ts = ts[~np.isnan(ts_Pc)]
                    ts_Pc = ts_Pc[~np.isnan(ts_Pc)]
                    # For each throat find the common pore and the uncommon
                    # pores
                    for i, t in enumerate(ts):
                        # Find common pore (cP) and uncommon pores (uPs)
                        b = set(net["throat.conns"][t])
                        cP = list(a.intersection(b))
                        uPs = list(a.symmetric_difference(b))
                        # If the common pore is not invaded but the others are
                        # The potential coop filling event can now happen
                        # Add the coop pressure to the queue
                        if (np.all(self[p_inv][uPs] > -1)) and (self[p_inv][cP] == -1):
                            # Coop pore filling fills the common pore
                            # The throats that gave access are not invaded now
                            # However, isolated throats between invaded pores
                            # Are taken care of elsewhere...
                            hq.heappush(queue, [ts_Pc[i], list(cP), "pore"])
