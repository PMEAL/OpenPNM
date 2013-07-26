r"""
Very much under construction.
Tried to combine examples from:
     * http://docs.enthought.com/mayavi/mayavi/auto/example_protein.html#example-protein
     * http://docs.enthought.com/mayavi/mayavi/auto/example_tvtk_segmentation.html#example-tvtk-segmentation
     * http://vtk.1045678.n5.nabble.com/Fill-closed-surface-with-voxels-td1234737.html
     * http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataToImageData
"""


#import OpenPNM
import scipy as sp
from mayavi import mlab
from tvtk.api import tvtk



class genGDLfromFibres(object):
    def __init__(self, shape=[0.,1.,0.,1.,0.,1.],num_fibres=100, diam=0.05,length=0.4):
        self._shape=sp.asarray(shape)
        self._num_fibres=num_fibres
        self._diam=diam
        self._length = length
        
    def gen_start_points(self):
        return sp.rand(self._num_fibres,3)*(self._shape[1::2]-self._shape[0::2]) \
                +self._shape[0::2]
    
    def gen_fibres(self):
        start_points = self.gen_start_points()
        fibre_directions = self.gen_fibre_directions()
        end_points   = start_points + fibre_directions
        for ii in range(0,3,1):
            ind = end_points[:,ii] < self._shape[ii*2]
            end_points[ind,ii]= start_points[ind,ii]-fibre_directions[ind,ii]
        for ii in range(0,3,1):
            ind = end_points[:,ii] > self._shape[1+ii*2]
            end_points[ind,ii] = start_points[ind,ii]-fibre_directions[ind,ii]
            
                
        self._points = sp.zeros((2*self._num_fibres,3),float)
        self._points[0::2,:] = start_points
        self._points[1::2,:] = end_points
        self._connections = sp.zeros((self._num_fibres,2),dtype=int)
        self._connections[:,0] = sp.arange(0,2*self._num_fibres,2,dtype=int)
        self._connections[:,1] = self._connections[:,0]+1
        
        
    def gen_fibre_directions(self):
        theta = sp.rand(self._num_fibres)*sp.pi
        phi   = sp.rand(self._num_fibres)*2*sp.pi
        directions = sp.zeros((self._num_fibres,3))
        directions[:,0] = sp.sin(theta)*sp.cos(phi)
        directions[:,1] = sp.sin(theta)*sp.sin(phi)
        directions[:,2] = sp.cos(theta)
        directions = directions*self._length
        return directions



if __name__ == '__main__':

    gen=genGDLfromFibres()
    gen.gen_fibres()
    
    from mayavi import mlab
    mlab.figure(1, bgcolor=(0, 0, 0))
    mlab.clf()
    pts = mlab.points3d(gen._points[:,0],gen._points[:,1],gen._points[:,2])
    pts.trait_set(visible=False)
    pts.mlab_source.dataset.lines = gen._connections
    tube = mlab.pipeline.tube(pts, tube_radius=0.025,tube_sides=12)
    tube.filter.capping=True
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))
    
    image = tvtk.ImageData(spacing=(gen._shape[1::2]-gen._shape[0::2])/50, origin=gen._shape[0::2])
    image.dimensions=[50,50,50]
    image.point_data.scalars=sp.zeros((50,50,50),dtype=int).ravel()
    image.point_data.scalars.name = 'scalars'
    
    stencil = tvtk.PolyDataToImageStencil()
    
    
    
    #mlab.view(49, 31.5, 52.8, (4.2, 37.3, 20.6))
    mlab.savefig('test.png')
    mlab.show()
    
    
    
#        // Stencil is created from PolyData and empty ImageData
#        vtkImageData* vi = vtkImageData::New();
#        vi->SetOrigin( 0,0,0 ); // adjust these to your needs
#        vi->SetSpacing( 1,1,1 ); // adjust these to your needs
#        vi->SetDimensions( 181,217,181 ); // adjust these to your needs
#        vi->SetScalarTypeToUnsignedChar ();
#        vi->AllocateScalars();
#        // outputMesh is of vtkPolyData* type and contains your mesh data
#        vtkPolyDataToImageStencil* pti = vtkPolyDataToImageStencil::New();
#        pti->SetInput( outputMesh );
#        pti->Update();
#        vtkImageStencil* is = vtkImageStencil::New();
#        is->SetInput( vi );
#        is->SetStencil( pti->GetOutput() );
#        is->ReverseStencilOff();
#        is->SetBackgroundValue(1);
#        is->Update();
#        // is->GetOutput() returns your image data as vtkImageData*
    
#    def image_data():
#    data = random.random((3, 3, 3))
#    i = tvtk.ImageData(spacing=(1, 1, 1), origin=(0, 0, 0))
#    i.point_data.scalars = data.ravel()
#    i.point_data.scalars.name = 'scalars'
#    i.dimensions = data.shape
#    return i
