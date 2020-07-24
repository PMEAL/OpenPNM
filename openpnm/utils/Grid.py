from terminaltables import AsciiTable


class Grid():

    def __init__(self, rows=1, cols=1, blank='---'):
        super().__init__()
        self.blank = blank
        header = [blank for i in range(cols)]
        self._grid = []
        [self._grid.append(header.copy()) for _ in range(rows)]

    def __getitem__(self, row):
        return self._grid[row]

    def __setitem__(self, row, cols):
        if len(cols) != self.ncols:
            raise Exception('Must write entire row')
        self._grid[row] = cols

    def _set_blank(self, blank):
        if hasattr(self, '_blank'):
            for i, row in enumerate(self._grid):
                for j, col in enumerate(self._grid[i]):
                    if self._grid[i][j] == self._blank:
                        self._grid[i][j] = blank
        self._blank = blank

    def _get_blank(self):
        if hasattr(self, '_blank'):
            return self._blank

    blank = property(fset=_set_blank, fget=_get_blank)

    @property
    def ncols(self):
        header = [self.blank for i in self._grid[0]]
        return len(header)

    @property
    def nrows(self):
        return len(self._grid)

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def size(self):
        return self.nrows * self.ncols

    def add_row(self, row=None):
        header = [self.blank for i in self._grid[0]]
        if row:
            self._grid.insert(row, header)
        else:
            self._grid.append(header)

    def add_col(self, col=None):
        for row in self._grid:
            if col:
                row.inert(col, self.blank)
            else:
                row.append(self.blank)

    def drop_col(self, col):
        for row in self._grid:
            row.remove(col)

    def drop_row(self, row):
        self._grid.remove(row)

    def __repr__(self):
        temp = AsciiTable(self._grid)
        temp.padding_left = 3
        temp.padding_right = 3
        temp.justify_columns = {col: 'center' for col in range(self.ncols)}
        return temp.table.__str__()
