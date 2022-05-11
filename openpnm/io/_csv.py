from openpnm.io._pandas import project_to_pandas
from openpnm.io import _parse_filename


def project_to_csv(project, filename=''):
    r"""
    Save all the pore and throat data on the Network and Phase objects to a CSV
    file

    Parameters
    ----------
    project : list
        An openpnm ``project`` object
    filename : str or path object
        The name of the file to store the data

    """
    network = project.network
    phases = project.phases
    df = project_to_pandas(network=network, phases=phases,
                           join=True, delim='.')
    if filename == '':
        filename = project.name
    fname = _parse_filename(filename=filename, ext='csv')
    df.to_csv(fname, index=False)
