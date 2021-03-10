# Generate Random Networks Based on Delaunay Triangulations,  Voronoi Tessellations, or Both

(New in Version 1.6)  A random network offers several advantages over the traditional Cubic arrangement: the topology is more 'natural' looking, and a wider pore size distribution can be achieved since pores are not constrained by the lattice spacing.  Random networks can be tricky to generate however, since the connectivity between pores is difficult to determine or define.  One surprisingly simple option is to use Delaunay triangulation to connect base points (which become pore centers) that are randomly distributed in space.  The Voronoi tessellation is a complementary graph that arises directly from the Delaunay graph which also connects essentially randomly distributed points in space into a network.  OpenPNM offers both of these types of network, plus the ability to create a network containing *both* including interconnections between Delaunay and Voronoi networks via the ```DelaunayVoronoiDual``` class.  In fact, creating the dual, interconnected network is the starting point, and any unwanted elements can be easily trimmed.

## Generate a Square Network with a DelaunayVoronoiDual Topology

``` python
>>> import OpenPNM as op
>>> pn = op.Network.DelaunayVoronoiDual(num_points=100, domain_size=[1, 1, 1])

```

The above line of code is deceptively simple.  The returned network (```pn```) contains a fully connected Delaunay network, its complementary Voronoi network, and interconnecting throats (or bonds) between each Delaunay pore (node) and its neighboring Voronoi pores.  Such a highly complex network would be useful for modeling pore phase transport (i.e. diffusion) on one network (i.e. Delaunay), solid phase transport (i.e. heat transfer) on the other network (i.e. Voronoi), and exchange of a species (i.e. heat) between the solid and void phases via the interconnecting bonds.  Each pore and throat is labelled accordingly (i.e. 'pore.delaunay', 'throat.voronoi'), and the interconnecting throats are labelled 'throat.interconnect'.  Moreover, pores and throats lying on the surface of the network are labelled 'surface'.  

A quick visualization of this network can be accomplished using OpenPNM's built-in graphing tool.  The following shows only the Voronoi connections that lie on the surface of the cube:

``` python
>>> Ts = pn.throats(['voronoi', 'surface'], mode='intersection')
>>> fig = op.Network.tools.plot_connections(network=pn, throats=Ts)

```

![](https://i.imgur.com/YDBtnFK.png)

One central feature of these networks are the flat boundaries, which are essential when performing transport calculations since they provide well-defined control surfaces for calculating flux.  This flat surfaces are accomplished by reflecting the base points across each face prior to performing the tessellations.  

## Obtain a Solo Delaunay (or Voronoi) Network

It is simple to delete one network (or the other) by trimming all of the other network's pores, which also removes all connected throats including the interconnections:

``` python
>>> pn.trim(pores=pn.pores('voronoi'))
>>> fig = op.Network.tools.plot_connections(network=pn)

```

## Create Random Networks of Spherical or Cylindrical Shape

Many porous materials come in spherical or cylindrical shapes, such as catalyst pellets.  The ```DelaunayVoronoiDual``` Network class can produce these geometries by specifying the ```domain_size``` in cylindrical [r, z] or spherical [r] coordinates:

``` python
>>> cyl = op.Network.DelaunayVoronoiDual(num_points=500, domain_size=[1, 5])
>>> fig = op.Network.tools.plot_connections(network = cyl,
...                                         throats=cyl.throats('surface'))

```

Which should yield this:

![](http://i.imgur.com/zzegpym.png)

``` python
>>> sph = op.Network.DelaunayVoronoiDual(num_points=500, domain_size=[2])
>>> fig = op.Network.tools.plot_connections(network = sph,
...                                         throats=sph.throats('surface'))

```

Which should result in this:

![](http://i.imgur.com/Lg9S8q1.png)

Note that the cylindrical and spherical networks don't look very nice when too few points are used, so at least about 200 is recommended.  

## Assign Pore Sizes to the Random Network

With pore centers randomly distributed in space it becomes challenging to know what pore size to assign to each location.  Assigning pores that are too large results in overlaps, which makes it impossible to properly account for porosity and transport lengths.  OpenPNM includes a Geometry model called ```largest_sphere``` that solves this problem.  Let's assign the largest possible pore size to each Voronoi node in the ```sph``` network just created:

``` python
>>> Ps = sph.pores('voronoi')
>>> Ts = sph.throats('voronoi')
>>> geom = op.Geometry.GenericGeometry(network=sph, pores=Ps, throats=Ts)
>>> mod = op.Geometry.models.pore_diameter.largest_sphere
>>> geom.models.add(propname='pore.diameter', model=mod)
>>> mod = op.Geometry.models.throat_length.straight
>>> geom.models.add(propname='throat.length', model=mod)
>>> mod = op.Geometry.models.throat_diameter.minpore
>>> geom.models.add(propname='throat.diameter', model=mod, factor=0.5)

```

The resulting geometrical properties can be viewed with ```geom.plot_histograms()``` (note that each realization will differ slightly):

![](https://i.imgur.com/4XQTWfi.png)
