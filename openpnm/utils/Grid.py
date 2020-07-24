from terminaltables import AsciiTable


class Grid():

    def __init__(self, rows=1, cols=1, blank='---'):
        super().__init__()
        self.blank = blank
        header = [blank for i in range(cols)]
        self._grid = AsciiTable([])
        [self._grid.table_data.append(header.copy()) for _ in range(rows)]
        self._grid.inner_heading_row_border = False
        self._grid.padding_left = 3
        self._grid.padding_right = 3
        self._grid.justify_columns = {col: 'center' for col in range(self.ncols)}

    def __getitem__(self, row):
        return self._grid.table_data[row]

    def __setitem__(self, row, cols):
        if len(cols) != self.ncols:
            raise Exception('Must write entire row')
        self._grid.table_data[row] = cols

    def _set_blank(self, blank):
        if hasattr(self, '_blank'):
            for i, row in enumerate(self._grid.table_data):
                for j, col in enumerate(self._grid.table_data[i]):
                    if self._grid.table_data[i][j] == self._blank:
                        self._grid.table_data[i][j] = blank
        self._blank = blank

    def _get_blank(self):
        if hasattr(self, '_blank'):
            return self._blank

    blank = property(fset=_set_blank, fget=_get_blank)

    @property
    def ncols(self):
        header = [self.blank for i in self._grid.table_data[0]]
        return len(header)

    @property
    def nrows(self):
        return len(self._grid.table_data)

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def size(self):
        return self.nrows * self.ncols

    def get_row(self, row):
        temp = self._grid.table_data[row]
        temp2 = Grid(1, self.ncols)
        for i, item in enumerate(temp):
            temp2[0][i] = item
        return temp2

    def get_col(self, col):
        temp = [self._grid.table_data[row][col] for row in range(self.nrows)]
        temp2 = Grid(1, 1)
        for row, item in enumerate(temp):
            if row == 0:
                temp2[0][0] = item
            else:
                temp2.add_row()
                temp2[row][0] = item
        return temp2

    def add_row(self, row=None):
        header = [self.blank for i in self._grid.table_data[0]]
        if row:
            self._grid.table_data.insert(row, header)
        else:
            self._grid.table_data.append(header)

    def add_col(self, col=None):
        for row in self._grid.table_data:
            if col:
                row.inert(col, self.blank)
            else:
                row.append(self.blank)

    def drop_col(self, col):
        for row in self._grid.table_data:
            row.remove(col)

    def drop_row(self, row):
        self._grid.table_data.remove(row)

    def __repr__(self):
        return self._grid.table.__str__()
