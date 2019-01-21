import scipy as sp


class BoundingBox(list):
    r"""
    Creates a bounding box around the given points, and offers several
    useful functions for querying the bounding box, such as determining the
    distance between a given of points and a given face of the box, or the
    area of a given face.
    """

    norm_xy = sp.array([0, 0, 1])
    norm_xz = sp.array([0, 1, 0])
    norm_yz = sp.array([1, 0, 0])

    def __init__(self, coords):
        super().__init__()
        if coords.shape[0] < 4:
            raise Exception('Must give 4 or more points to define a box')
        if coords.shape[1] != 3:
            raise Exception('Given coordinates must be 3D')
        temp = [sp.amin(coords, axis=0), sp.amax(coords, axis=0)]
        self.extend(temp)
        self.extents = sp.array([self[1][0] - self[0][0],
                                 self[1][1] - self[0][1],
                                 self[1][2] - self[0][2]])
        self.offset = sp.array([self[0][0], self[0][1], self[0][2]])

    @property
    def length_top_to_bottom(self):
        L = self.extents[self.norm_xy.astype(bool)][0]
        return L

    @property
    def length_left_to_right(self):
        L = self.extents[self.norm_xz.astype(bool)][0]
        return L

    @property
    def length_front_to_back(self):
        L = self.extents[self.norm_yz.astype(bool)][0]
        return L

    @property
    def area_of_top(self):
        r"""
        Calculates the area of the top face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_xy.astype(bool)])
        return A

    @property
    def area_of_bottom(self):
        r"""
        Calculates the area of the bottom face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_xy.astype(bool)])
        return A

    @property
    def area_of_left(self):
        r"""
        Calculates the area of the left face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_xz.astype(bool)])
        return A

    @property
    def area_of_right(self):
        r"""
        Calculates the area of the right face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_xz.astype(bool)])
        return A

    @property
    def area_of_front(self):
        r"""
        Calculates the area of the front face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_yz.astype(bool)])
        return A

    @property
    def area_of_back(self):
        r"""
        Calculates the area of the back face of the bounding box
        """
        A = sp.prod(self.extents[~self.norm_yz.astype(bool)])
        return A

    def dist_to_left(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the left face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        return sp.dot(points - self.offset, self.norm_xz)

    def dist_to_right(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the right face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        d = self.dist_to_left(points)
        return sp.absolute(d - self.extents[self.norm_xz.astype(bool)])

    def dist_to_front(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the front face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        return sp.dot(points - self.offset, self.norm_yz)

    def dist_to_back(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the back face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        d = self.dist_to_front(points)
        return sp.absolute(d - self.extents[self.norm_yz.astype(bool)])

    def dist_to_bottom(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the bottom face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        return sp.dot(points - self.offset, self.norm_xy)

    def dist_to_top(self, points):
        r"""
        Given a set of points, calculates the distance between each point
        and the top face of the bounding box

        Parameters
        ----------
        points : array_like

        Returns
        -------
        An ND-array containing the absolute distance between each supplied
        point and the given face
        """
        d = self.dist_to_bottom(points)
        return sp.absolute(d - self.extents[self.norm_xy.astype(bool)])
